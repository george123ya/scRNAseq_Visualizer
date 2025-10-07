// violin.js

// **1. Shared State and Plotting Logic Scope**
// These variables are accessible by all Shiny message handlers and the binding.
let currentDatasets;
let plotElements = {};

// Utility function to convert RGB to HEX
function rgbToHex(rgb) {
  if (rgb.startsWith('#')) return rgb.toUpperCase();
  const match = rgb.match(/\d+/g);
  if (!match) return '#000000';
  const [r, g, b] = match.map(Number);
  return `#${((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1).toUpperCase()}`;
}

// Pickr-based color picker
async function createInlineColorPicker(initialColor, onColorChange, onClose, colorBoxContainer) {
  let Pickr;
  try {
    const pickrModule = await import('https://esm.sh/@simonwep/pickr@1.9.0');
    Pickr = pickrModule.default;
  } catch (error) {
    console.error('Failed to load Pickr:', error);
    return null;
  }

  if (!document.querySelector('link[href*="pickr"]')) {
    const link = document.createElement('link');
    link.rel = 'stylesheet';
    link.href = 'https://esm.sh/@simonwep/pickr@1.9.0/dist/themes/nano.min.css';
    document.head.appendChild(link);
  }

  const container = document.createElement('div');
  container.style.cssText = `
    width: 100%;
    height: 100%;
    border-radius: 2px;
  `;

  colorBoxContainer.innerHTML = '';
  colorBoxContainer.appendChild(container);

  const pickr = Pickr.create({
    el: container,
    theme: 'nano',
    default: initialColor,
    swatches: null,
    components: {
      preview: true,
      opacity: false,
      hue: true,
      interaction: {
        hex: true,
        rgba: false,
        input: true,
        save: true,
        clear: true
      }
    },
    defaultRepresentation: 'HEX',
    position: 'bottom-middle'
  });

  pickr.on('change', (color) => {
    const hex = color.toHEXA().toString().toUpperCase();
    onColorChange(hex, false);
  });

  pickr.on('save', (color) => {
    const hex = color ? color.toHEXA().toString().toUpperCase() : initialColor;
    onColorChange(hex, true);
    pickr.hide();
  });

  pickr.on('hide', () => {
    onClose();
  });

  pickr.on('clear', () => {
    onColorChange(initialColor, true);
    pickr.hide();
  });

  pickr.show();

  return container;
}


// **2. Core Plotting and Update Function**
function updatePlot() {
    if (!plotElements.ctx || !plotElements.svg || !currentDatasets) {
        console.warn("Plot elements or data not initialized.");
        return;
    }

    const { ctx, svg, width, height, margin, legendWidth, plotWidth, pointOffsets, showPoints, logScale, pointSize, legendPosition, groupColors } = plotElements;

    ctx.save();
    ctx.clearRect(0, 0, width, height);
    svg.selectAll('*').remove();

    // Re-create the clipping path
    svg.append('clipPath')
        .attr('id', 'violin-clip')
        .append('rect')
        .attr('x', margin.left)
        .attr('y', margin.top)
        .attr('width', width - margin.left - margin.right - legendWidth)
        .attr('height', height - margin.top - margin.bottom);

    // D3 Utility Functions
    function standardDeviation(arr) {
        if (!Array.isArray(arr) || arr.length < 2) return 0;
        const mean = d3.mean(arr);
        if (mean === undefined || mean === null) return 0;
        const variance = d3.sum(arr.map(x => (x - mean) ** 2)) / (arr.length - 1);
        return Math.sqrt(variance);
    }
    function iqr(arr) {
        const sorted = arr.slice().sort(d3.ascending);
        const q1 = d3.quantile(sorted, 0.25);
        const q3 = d3.quantile(sorted, 0.75);
        return q3 - q1;
    }
    function kernelDensityEstimator(kernel, X) {
        return function(V) {
            return X.map(x => [x, d3.mean(V, v => kernel(x - v))]);
        };
    }
    function kernelEpanechnikov(k) {
        return function(v) {
            return Math.abs(v /= k) <= 1 ? 0.75 * (1 - v*v) / k : 0;
        };
    }

    const allValues = currentDatasets.flatMap(d => d.values);
    if (allValues.length === 0) {
        console.warn("No valid values in datasets");
        return;
    }

    const dataExtent = d3.extent(allValues);
    const padding = (dataExtent[1] - dataExtent[0]) * 0.05;
    const y = logScale ?
        d3.scaleLog().domain([Math.max(d3.min(allValues, d => d > 0 ? d : 1), 1), d3.max(allValues)]).range([height - margin.bottom, margin.top]).nice() :
        d3.scaleLinear().domain([dataExtent[0] - padding, dataExtent[1] + padding]).range([height - margin.bottom, margin.top]);

    // Draw violins
    currentDatasets.forEach((dataset, i) => {
        if (!dataset.values || !Array.isArray(dataset.values) || dataset.values.length === 0) return;

        const n = dataset.values.length;
        const stdDev = standardDeviation(dataset.values);
        const dataIQR = iqr(dataset.values);
        const h = Math.max(0.9 * Math.min(stdDev, dataIQR / 1.34) * Math.pow(n, -1/5), 0.01);

        const dataExtentGroup = d3.extent(dataset.values);
        const kdePadding = (dataExtentGroup[1] - dataExtentGroup[0]) * 0.05;
        const kdeTicks = d3.range(dataExtentGroup[0] - kdePadding, dataExtentGroup[1] + kdePadding, (dataExtentGroup[1] - dataExtentGroup[0]) / 100);
        const kde = kernelDensityEstimator(kernelEpanechnikov(h), kdeTicks);
        let density = kde(dataset.values);

        const x = d3.scaleLinear().domain([0, d3.max(density, d => d[1]) || 1]).range([0, plotWidth * 0.4]);
        const xCenter = margin.left + plotWidth * (i + 0.5);

        svg.append('path')
            .datum(density)
            .attr('class', `violin-path-${i}`)
            .attr('fill', groupColors[i])
            .attr('fill-opacity', 0.6)
            .attr('stroke', '#333')
            .attr('stroke-width', 1)
            .attr('clip-path', 'url(#violin-clip)')
            .attr('d', d3.area()
                .x0(d => xCenter - x(d[1]))
                .x1(d => xCenter + x(d[1]))
                .y(d => y(d[0]))
                .curve(d3.curveBasis)
            );
    });

    // Draw points on canvas
    if (showPoints) {
        const offsets = pointOffsets;
        ctx.globalAlpha = 0.5;
        currentDatasets.forEach((dataset, i) => {
            const xCenter = margin.left + plotWidth * (i + 0.5);
            ctx.fillStyle = groupColors[i];
            dataset.values.forEach((v, j) => {
                ctx.beginPath();
                ctx.arc(
                    xCenter + offsets[i][j],
                    y(v),
                    pointSize, // Use pointSize directly
                    0,
                    2 * Math.PI
                );
                ctx.fill();
            });
        });
        ctx.globalAlpha = 1.0;
    }

    // Draw legend or labels based on legendPosition
    const legendX = width - margin.right - legendWidth + 5;
    const legendY = margin.top;
    const labelY = height - margin.bottom + 20;

    if (legendPosition === 'Legend') {
        const legend = svg.append('g')
            .attr('class', 'legend-group')
            .attr('transform', `translate(${legendX}, ${legendY})`);

        currentDatasets.forEach((dataset, i) => {
            if (!dataset.values || dataset.values.length === 0) return;
            
            const legendItem = legend.append('g')
                .attr('transform', `translate(0, ${i * 20})`);
                
            const colorBoxContainer = legendItem.append('foreignObject')
                .attr('width', 15)
                .attr('height', 15)
                .node();

            legendItem.append('rect')
                .attr('class', `legend-rect-${i}`)
                .attr('width', 15)
                .attr('height', 15)
                .attr('fill', groupColors[i])
                .attr('fill-opacity', 0.6)
                .attr('stroke', '#333')
                .attr('stroke-width', 1)
                .style('cursor', 'pointer')
                .on('click', async () => {
                    const initialColor = groupColors[i];
                    await createInlineColorPicker(
                        initialColor,
                        (hex, isFinal) => {
                            // Live change: only update the SVG path's color
                            svg.select(`.violin-path-${i}`).attr('fill', hex);
                            groupColors[i] = hex; 
                            plotElements.groupColors = groupColors;

                            if (isFinal) {
                                // Final change: re-render the points on the canvas
                                updatePlot();
                                Shiny.setInputValue('group_colors', groupColors, {priority: 'event'});
                            }
                        },
                        () => {
                            // On close, the SVG path is already updated. No further action needed.
                        },
                        colorBoxContainer
                    );
                });
                
            legendItem.append('text')
                .attr('x', 20)
                .attr('y', 12)
                .attr('fill', '#333')
                .attr('font-size', 12)
                .text(dataset.metadata?.name || `Group ${i + 1}`);
        });
    } else {
        currentDatasets.forEach((dataset, i) => {
            if (!dataset.values || dataset.values.length === 0) return;
            const xCenter = margin.left + plotWidth * (i + 0.5);
            svg.append('text')
                .attr('class', 'bottom-label')
                .attr('x', xCenter)
                .attr('y', labelY)
                .attr('text-anchor', 'middle')
                .attr('fill', '#333')
                .attr('font-size', 12)
                .text(dataset.metadata?.name || `Group ${i + 1}`);
        });
    }

    // Y-axis
    const yAxis = d3.axisLeft(y).ticks(5).tickFormat(logScale ? d3.format('.2f') : null);
    svg.select('.y-axis').remove();
    svg.append('g')
        .attr('class', 'y-axis')
        .attr('transform', `translate(${margin.left},0)`)
        .call(yAxis);

    // Mouseover handlers
    const hoverLine = svg.select('.hover-line').empty() ? svg.append('line') : svg.select('.hover-line');
    hoverLine
        .attr('class', 'hover-line')
        .attr('stroke', 'red')
        .attr('stroke-width', 1)
        .style('display', 'none');

    const hoverLabel = svg.select('.hover-label').empty() ? svg.append('text') : svg.select('.hover-label');
    hoverLabel
        .attr('class', 'hover-label')
        .attr('x', margin.left + (width - margin.left - margin.right - legendWidth) / 2)
        .attr('y', height - margin.bottom + 40)
        .attr('text-anchor', 'middle')
        .attr('fill', 'red')
        .style('font-size', '12px')
        .text('');

    const counters = currentDatasets.map((dataset, i) => {
        if (!dataset.values || !Array.isArray(dataset.values) || dataset.values.length === 0) return null;
        const xCenter = margin.left + plotWidth * (i + 0.5);
        return svg.select(`.counter-${i}`).empty() ? svg.append('text') : svg.select(`.counter-${i}`);
    }).filter(d => d !== null);

    counters.forEach((counter, i) => {
        counter
            .attr('class', `counter-${i}`)
            .attr('x', margin.left + plotWidth * (i + 0.5))
            .attr('y', margin.top - 10)
            .attr('text-anchor', 'middle')
            .attr('fill', '#333')
            .attr('font-size', 12)
            .style('display', 'none')
            .text('');
    });

    const allSorted = currentDatasets.map(d => d.values.slice().sort(d3.ascending));

    d3.select(plotElements.el)
        .on('mousemove', function(event) {
            const [mx, my] = d3.pointer(event);
            if (my < margin.top || my > height - margin.bottom) return;
            const yVal = y.invert(my);
            const totalCount = allSorted.reduce((sum, sorted) => sum + d3.bisectLeft(sorted, yVal), 0);
            
            hoverLine.attr('x1', margin.left).attr('x2', width - margin.right - legendWidth).attr('y1', my).attr('y2', my).style('display', null);
            hoverLabel.text(`Count ≤ ${yVal.toFixed(2)}`);
            allSorted.forEach((sorted, i) => {
                if (counters[i]) {
                    const count = d3.bisectLeft(sorted, yVal);
                    const pct = (count / sorted.length) * 100;

                    // Clear any previous tspans
                    counters[i].selectAll("tspan").remove();

                    // First line (count)
                    counters[i]
                        .append("tspan")
                        .attr("x", margin.left + plotWidth * (i + 0.5))
                        .attr("text-anchor", "middle")
                        .text(`n=${count}`);

                    // Second line (percentage)
                    counters[i]
                        .append("tspan")
                        .attr("x", margin.left + plotWidth * (i + 0.5))
                        .attr("dy", 12) // move down relative to previous line
                        .attr("text-anchor", "middle")
                        .text(`(${pct.toFixed(1)}%)`);

                    counters[i].style("display", null);
                }
            });
            Shiny.setInputValue('hover_count', { count: totalCount, value: yVal }, {priority: 'event'});
        })
        .on('mouseleave', function() {
            hoverLine.style('display', 'none');
            hoverLabel.text('');
            counters.forEach(counter => counter.style('display', 'none'));
            Shiny.setInputValue('hover_count', null, {priority: 'event'});
        });
    ctx.restore();
}

// **3. Shiny Output Binding**
var violinBinding = new Shiny.OutputBinding();
$.extend(violinBinding, {
    find: function(scope) {
        return $(scope).find('.violin_container');
    },

    renderValue: function(el, data) {
        if (!data || !data.datasets) return;
        
        // Clear the SVG when rendering the violin plot
        d3.select(el).select('#violin_overlay').selectAll('*').remove();
        
        const width = 800;
        const height = 400;
        const margin = { top: 40, right: 40, bottom: 40, left: 40 };
        const legendWidth = 150;
        const plotWidth = (width - margin.left - margin.right - legendWidth) / data.datasets.length;
        
        currentDatasets = data.datasets;

        plotElements = {
            el: el,
            ctx: d3.select(el).select('#violin_plot').node().getContext('2d'),
            svg: d3.select(el).select('#violin_overlay'),
            width: width,
            height: height,
            margin: margin,
            legendWidth: legendWidth,
            plotWidth: plotWidth,
            pointOffsets: currentDatasets.map(d => d.values.map(() => (Math.random() - 0.5) * (plotWidth * 0.2))),
            groupColors: data.datasets.map((d, i) => d.metadata?.color ? d.metadata.color : d3.schemeSet1[i % d3.schemeSet1.length]),
            showPoints: data.showPoints,
            logScale: data.logScale,
            pointSize: data.pointSize,
            legendPosition: data.legendPosition
        };

        updatePlot();
    }
});
Shiny.outputBindings.register(violinBinding, 'violinBinding');


// **4. Shiny Custom Message Handlers**
Shiny.addCustomMessageHandler('drawViolins', function(message) {
    const el = document.querySelector('.violin_container');
    if (el) {
        violinBinding.renderValue(el, message);
    }
});

Shiny.addCustomMessageHandler('updateViolinInputs', function(message) {
    if (!plotElements.el) return;
    
    // Update the state from the message
    plotElements.showPoints = message.showPoints;
    plotElements.logScale = message.logScale;
    plotElements.pointSize = message.pointSize;
    plotElements.legendPosition = message.legendPosition;
    
    // Re-render the plot with the new state
    updatePlot();
});

// **5. Shared State for Scatterplot** (Unchanged)
let scatterState = {
  scatterplot: null,
  normalizedPoints: [],
  originalPoints: [],
  xDomain: [0, 1],
  yDomain: [0, 1],
  currentView: null,
  drawOptions: {},
  xLabel: '',
  yLabel: '',
  logScale: false,
  xAxis: null,
  yAxis: null,
  xAxisG: null,
  yAxisG: null,
  viewSubscription: null
};

// **6. normalizePoints** (Add validation log)
function normalizePoints(points, xDomain, yDomain, logScale) {
  let xMin = xDomain[0], xMax = xDomain[1];
  let yMin = yDomain[0], yMax = yDomain[1];
  const epsilon = 1e-6;
  if (xMax - xMin < epsilon) { 
    xMax = xMin + epsilon; 
    console.warn('Zero x-range; added epsilon');
  }
  if (yMax - yMin < epsilon) { 
    yMax = yMin + epsilon; 
    console.warn('Zero y-range; added epsilon');
  }
  if (isNaN(xMin) || isNaN(xMax) || isNaN(yMin) || isNaN(yMax)) {
    console.error('Invalid domains in normalizePoints:', { xDomain, yDomain });
    return points.map(() => [-1, -1, 0, 0]);  // Fallback degenerate points
  }
  const xScale = logScale ? d3.scaleLog().domain([Math.max(xMin, 1e-6), xMax]).range([0, 1]) : d3.scaleLinear().domain([xMin, xMax]).range([0, 1]);
  const yScale = logScale ? d3.scaleLog().domain([Math.max(yMin, 1e-6), yMax]).range([0, 1]) : d3.scaleLinear().domain([yMin, yMax]).range([0, 1]);
  
  return points.map(p => [
    (logScale ? xScale(p[0]) : (p[0] - xMin) / (xMax - xMin)) * 2 - 1,
    (logScale ? yScale(p[1]) : (p[1] - yMin) / (yMax - yMin)) * 2 - 1,
    p[2],
    p[3]
  ]);
}

// **7. computeDrawOptions** (Unchanged)
function computeDrawOptions(colorBy, colorLevels) {
  const options = { colorBy: null, pointColor: ['#000000'] };
  if (colorBy !== "None" && colorLevels && colorLevels.length > 0) {
    options.colorBy = 'valueA';
    const colorScale = d3.scaleOrdinal(d3.schemeCategory10);
    options.pointColor = colorScale.domain(d3.range(colorLevels.length)).range();
  }
  return options;
}

// **UPDATED: Shared Axes Update Function** (Add NaN checks)
function updateAxesFromView(event, margin, width, height) {
  if (!scatterState.xDomain || !scatterState.yDomain || !event || !event.xScale || !event.yScale || 
      !scatterState.xAxis || !scatterState.yAxis || !scatterState.xAxisG || !scatterState.yAxisG) {
    console.warn('Skipping axes update: invalid inputs');
    return;
  }
  const xDom = event.xScale.domain();
  const yDom = event.yScale.domain();
  if (xDom.some(isNaN) || yDom.some(isNaN) || xDom[0] === undefined || yDom[0] === undefined) {
    console.warn('Skipping axes update: NaN/invalid event domains', { xDom, yDom });
    return;
  }
  const rangeX = scatterState.xDomain[1] - scatterState.xDomain[0];
  const rangeY = scatterState.yDomain[1] - scatterState.yDomain[0];
  const newXDomain = [
    scatterState.xDomain[0] + (xDom[0] + 1) / 2 * rangeX,
    scatterState.xDomain[0] + (xDom[1] + 1) / 2 * rangeX
  ];
  const newYDomain = [
    scatterState.yDomain[0] + (yDom[0] + 1) / 2 * rangeY,
    scatterState.yDomain[0] + (yDom[1] + 1) / 2 * rangeY
  ];
  if (newXDomain.some(isNaN) || newYDomain.some(isNaN)) {
    console.warn('Skipping axes update: NaN computed domains', { newXDomain, newYDomain });
    return;
  }
  const newXScale = d3.scaleLinear().domain(newXDomain).range([margin.left, width - margin.right]);
  const newYScale = d3.scaleLinear().domain(newYDomain).range([height - margin.bottom, margin.top]);
  scatterState.xAxis.scale(newXScale);
  scatterState.yAxis.scale(newYScale);
  scatterState.xAxisG.call(scatterState.xAxis);
  scatterState.yAxisG.call(scatterState.yAxis);
  console.log('Axes updated successfully');
}

// **UPDATED: RAF Updater** (Add validity checks; auto-stop after ~600ms transition)
function startAnimatedAxesUpdate(scatterplot, container, margin, width, height, xDomain, yDomain) {
  if (!scatterplot || !xDomain || !yDomain || xDomain.some(isNaN) || yDomain.some(isNaN)) {
    console.warn('Skipping animated axes: invalid params');
    return () => {};
  }

  let rafId;
  let latestEvent = null;
  const startTime = Date.now();
  const transitionDuration = 600;  // ms; match regl default

  const tempHandler = (event) => { 
    latestEvent = { xScale: event.xScale, yScale: event.yScale };  // FIXED: Keep scales, extract domain in updateAxes
  };
  const tempSub = scatterplot.subscribe('view', tempHandler);
  scatterState.viewSubscription = { handler: tempHandler, sub: tempSub };

  const updateAxes = () => {
    if (Date.now() - startTime > transitionDuration) {
      console.log('Animation loop stopped: transition complete');
      return;  // Auto-stop after duration
    }
    if (latestEvent && !latestEvent.xScale.domain().some(isNaN)) {
      updateAxesFromView(latestEvent, margin, width, height);
    }
    rafId = requestAnimationFrame(updateAxes);
  };
  updateAxes();

  return () => {
    cancelAnimationFrame(rafId);
    if (scatterState.viewSubscription) {
      const { handler, sub } = scatterState.viewSubscription;
      if (sub && typeof sub.unsubscribe === 'function') {
        sub.unsubscribe();
      } else if (scatterplot._pubSub) {
        scatterplot._pubSub.unsubscribe('view', handler);
      }
      scatterState.viewSubscription = null;
    }
  };
}

// **UPDATED: autoAdjustZoom** (Better validation; fallback to reset)
const autoAdjustZoom = (points, scatterplot, xDomain, yDomain) => {
  if (!points || points.length === 0) {
    console.warn('Invalid zoom: no points');
    return;
  }
  if (!xDomain || !yDomain || xDomain.some(d => d === undefined || isNaN(d)) || yDomain.some(d => d === undefined || isNaN(d))) {
    console.warn('Invalid zoom params; falling back to reset', { xDomain, yDomain });
    try {
      scatterplot.reset({ preventEvent: false });  // Triggers view event for sync
    } catch (e) {
      console.error('Reset error:', e);
    }
    return;
  }
  try {
    const rangeX = xDomain[1] - xDomain[0];
    const rangeY = yDomain[1] - yDomain[0];
    const padX = rangeX * 0.1;
    const padY = rangeY * 0.1;
    const normalizedX = -1 + (2 * padX / rangeX);
    const normalizedWidth = 2 - (2 * padX / rangeX);
    const normalizedY = -1 + (2 * padY / rangeY);
    const normalizedHeight = 2 - (2 * padY / rangeY);
    if (isNaN(normalizedX) || isNaN(normalizedY) || isNaN(normalizedWidth) || isNaN(normalizedHeight)) {
      console.warn('Invalid zoom bounds; falling back to reset');
      scatterplot.reset({ preventEvent: false });
      return;
    }
    scatterplot.zoomToArea({ x: normalizedX, y: normalizedY, width: normalizedWidth, height: normalizedHeight }, true);
    console.log('Zoom applied successfully');
  } catch (e) {
    console.error('Zoom error:', e);
    scatterplot.reset({ preventEvent: false });
  }
};

// **9. Updated scatterplotData Handler** (Add manual post-sync)
Shiny.addCustomMessageHandler('scatterplotData', function(message) {
  const { data, x_label, y_label, color_by, color_levels, point_size, logScale, metadata, forceRedraw } = message;
  const points = data || [];
  const container = document.getElementById('scatterplot');
  if (!container) return;
  const width = container.clientWidth;
  const height = container.clientHeight;
  const margin = { top: 20, right: 20, bottom: 50, left: 60 };

  let xDomain = points.length > 0 ? d3.extent(points, d => d[0]) : [0, 1];
  let yDomain = points.length > 0 ? d3.extent(points, d => d[1]) : [0, 1];
  if (xDomain[0] === undefined) xDomain = [0, 1];
  if (yDomain[0] === undefined) yDomain = [0, 1];
  if (xDomain.some(isNaN) || yDomain.some(isNaN)) {
    console.error('Invalid computed domains:', { xDomain, yDomain });
    xDomain = [0, 1]; yDomain = [0, 1];
  }
  const normalizedPoints = normalizePoints(points, xDomain, yDomain, logScale || false);

  const isFullRedraw = !scatterState.scatterplot || forceRedraw || scatterState.originalPoints.length !== points.length;

  try {
    if (isFullRedraw) {
      // Full recreation (unchanged, but ensure valid initial axes)
      container.innerHTML = '';
      const canvas = document.createElement('canvas');
      canvas.id = 'scatterplotCanvas';
      canvas.style.position = 'absolute';
      canvas.style.top = `${margin.top}px`;
      canvas.style.left = `${margin.left}px`;
      canvas.style.width = `${width - margin.left - margin.right}px`;
      canvas.style.height = `${height - margin.top - margin.bottom}px`;
      container.appendChild(canvas);

      const svg = d3.select(container).append('svg')
        .attr('width', width)
        .attr('height', height)
        .style('position', 'absolute')
        .style('top', 0)
        .style('left', 0)
        .style('pointer-events', 'none');

      const xAxis = d3.axisBottom().ticks(6);
      const xAxisG = svg.append('g').attr('transform', `translate(0, ${height - margin.bottom})`);
      // Initial call with full domain to show ticks
      const initXScale = d3.scaleLinear().domain(xDomain).range([margin.left, width - margin.right]);
      xAxis.scale(initXScale);
      xAxisG.call(xAxis);

      svg.append('text')
        .attr('class', 'x-label')
        .attr('x', margin.left + (width - margin.left - margin.right) / 2)
        .attr('y', height - 10)
        .attr('text-anchor', 'middle')
        .attr('fill', 'black')
        .text(x_label);

      const yAxis = d3.axisLeft().ticks(6);
      const yAxisG = svg.append('g').attr('transform', `translate(${margin.left},0)`);
      // Initial call with full domain
      const initYScale = d3.scaleLinear().domain(yDomain).range([height - margin.bottom, margin.top]);
      yAxis.scale(initYScale);
      yAxisG.call(yAxis);

      svg.append('text')
        .attr('class', 'y-label')
        .attr('transform', 'rotate(-90)')
        .attr('x', -(margin.top + (height - margin.top - margin.bottom) / 2))
        .attr('y', 15)
        .attr('text-anchor', 'middle')
        .attr('fill', 'black')
        .text(y_label);

      svg.selectAll('.domain').attr('stroke', 'black').attr('stroke-width', 1.5);

      let tooltip = document.getElementById('scatterplotTooltip');
      if (!tooltip) {
        tooltip = document.createElement('div');
        tooltip.id = 'scatterplotTooltip';
        tooltip.style.cssText = `position:absolute;background-color:rgba(0,0,0,0.8);color:white;padding:5px 10px;border-radius:4px;font-size:12px;pointer-events:none;z-index:1000;display:none;`;
        container.appendChild(tooltip);
      }

      import('https://esm.sh/regl-scatterplot@1.14.1').then(module => {
        const createScatterplot = module.default;
        const scatterplot = createScatterplot({
          canvas: canvas,
          width: width - margin.left - margin.right,
          height: height - margin.top - margin.bottom,
          xScale: d3.scaleLinear().domain([-1, 1]).range([0, width - margin.left - margin.right]),  // Note: regl uses normalized scales
          yScale: d3.scaleLinear().domain([-1, 1]).range([height - margin.top - margin.bottom, 0]),
          pointSize: point_size || 2
        });

        let drawOptions = computeDrawOptions(color_by, color_levels);
        scatterplot.set(drawOptions);
        scatterplot.draw(normalizedPoints);

        scatterState = { 
          ...scatterState, 
          scatterplot, 
          normalizedPoints, 
          originalPoints: points,
          xDomain, 
          yDomain, 
          drawOptions,
          xLabel: x_label,
          yLabel: y_label,
          logScale: logScale || false,
          xAxis,
          yAxis,
          xAxisG,
          yAxisG,
          viewSubscription: null
        };

        scatterplot.subscribe('pointOver', pointIndex => {
          const point = scatterState.originalPoints[pointIndex];
          const [px, py] = scatterplot.getScreenPosition(pointIndex);
          let tooltipContent = `X: ${point[0].toFixed(2)}<br>Y: ${point[1].toFixed(2)}`;
          if (metadata && metadata[pointIndex]) {
            tooltipContent += `<br>Metadata: ${metadata[pointIndex]}`;
          }
          tooltip.innerHTML = tooltipContent;
          tooltip.style.display = 'block';
          tooltip.style.left = `${px + margin.left + 10}px`;
          tooltip.style.top = `${py + margin.top + 10}px`;
        });

        scatterplot.subscribe('pointOut', () => {
          tooltip.style.display = 'none';
        });

        scatterplot.subscribe('view', event => {
          scatterState.currentView = { xScale: event.xScale, yScale: event.yScale };  // FIXED: Store scales
          updateAxesFromView(event, margin, width, height);
        });

        setTimeout(() => {
          if (!scatterState.currentView) {
            autoAdjustZoom(points, scatterplot, xDomain, yDomain);
          }
        }, 100);
      });
    } else {
      // Partial data update (add manual sync post-draw)
      console.log('Performing transitional data update...');
      const prevView = scatterState.currentView;
      const useTransition = points.length === scatterState.originalPoints.length;
      scatterState.originalPoints = points;
      scatterState.xDomain = xDomain;
      scatterState.yDomain = yDomain;
      const newNormalizedPoints = normalizePoints(points, xDomain, yDomain, logScale || false);
      scatterState.normalizedPoints = newNormalizedPoints;
      scatterState.logScale = logScale || false;

      const newDrawOptions = computeDrawOptions(color_by, color_levels);
      if (JSON.stringify(newDrawOptions) !== JSON.stringify(scatterState.drawOptions)) {
        scatterState.drawOptions = newDrawOptions;
        scatterState.scatterplot.set(newDrawOptions);
      }

      if (point_size !== undefined) {
        scatterState.scatterplot.set({ pointSize: point_size });
      }

      d3.select(container).select('.x-label').text(x_label);
      d3.select(container).select('.y-label').text(y_label);

      let targetXDomain, targetYDomain;
      if (prevView && useTransition) {
        const rangeX = xDomain[1] - xDomain[0];
        const rangeY = yDomain[1] - yDomain[0];
        targetXDomain = [
          xDomain[0] + (prevView.xScale.domain()[0] + 1) / 2 * rangeX,
          xDomain[0] + (prevView.xScale.domain()[1] + 1) / 2 * rangeX
        ];
        targetYDomain = [
          yDomain[0] + (prevView.yScale.domain()[0] + 1) / 2 * rangeY,
          yDomain[0] + (prevView.yScale.domain()[1] + 1) / 2 * rangeY
        ];
        const normX = -1 + 2 * ((targetXDomain[0] - xDomain[0]) / rangeX);
        const normW = 2 * ((targetXDomain[1] - targetXDomain[0]) / rangeX);
        const normY = -1 + 2 * ((targetYDomain[0] - yDomain[0]) / rangeY);
        const normH = 2 * ((targetYDomain[1] - targetYDomain[0]) / rangeY);
        scatterState.scatterplot.zoomToArea({ x: normX, y: normY, width: normW, height: normH }, false);
      }

      const drawPromise = scatterState.scatterplot.draw(newNormalizedPoints, { transition: useTransition });
      if (useTransition) {
        const stopAnimation = startAnimatedAxesUpdate(scatterState.scatterplot, container, margin, width, height, xDomain, yDomain);
        drawPromise.then(() => {
          stopAnimation();
          if (scatterState.viewSubscription) {
            const { handler, sub } = scatterState.viewSubscription;
            if (sub && typeof sub.unsubscribe === 'function') {
              sub.unsubscribe();
            } else if (scatterState.scatterplot._pubSub) {
              scatterState.scatterplot._pubSub.unsubscribe('view', handler);
            }
            scatterState.viewSubscription = null;
          }
          scatterState.scatterplot.refresh();
          // NEW: Manual post-sync if no recent view
          setTimeout(() => {
            if (scatterState.currentView && !scatterState.currentView.xScale.domain().some(isNaN)) {
              updateAxesFromView(scatterState.currentView, margin, width, height);
            } else if (targetXDomain && targetYDomain) {
              autoAdjustZoom(points, scatterState.scatterplot, targetXDomain, targetYDomain);
            } else {
              autoAdjustZoom(points, scatterState.scatterplot, xDomain, yDomain);
            }
          }, 50);
        });
      } else {
        drawPromise.then(() => {
          scatterState.scatterplot.refresh();
          // NEW: Manual sync
          setTimeout(() => {
            if (!prevView) {
              autoAdjustZoom(points, scatterState.scatterplot, xDomain, yDomain);
            } else if (targetXDomain && targetYDomain) {
              autoAdjustZoom(points, scatterState.scatterplot, targetXDomain, targetYDomain);
            }
          }, 50);
        });
      }
    }
  } catch (e) {
    console.error('Scatterplot update error:', e);
  }
});

// **10. Updated updateScatterplotInputs Handler** (Add manual sync for log/color)
Shiny.addCustomMessageHandler('updateScatterplotInputs', function(message) {
  if (!scatterState.scatterplot) return;

  const { point_size, logScale, color_by, color_levels, x_label, y_label } = message;
  const container = document.getElementById('scatterplot');
  if (!container) return;
  const width = container.clientWidth;
  const height = container.clientHeight;
  const margin = { top: 20, right: 20, bottom: 50, left: 60 };
  let needsRedraw = false;
  let needsRefresh = false;

  if (point_size !== undefined) {
    scatterState.scatterplot.set({ pointSize: point_size });
    needsRefresh = true;
  }

  if (color_by !== undefined || color_levels !== undefined) {
    const newDrawOptions = computeDrawOptions(color_by, color_levels);
    if (JSON.stringify(newDrawOptions) !== JSON.stringify(scatterState.drawOptions)) {
      scatterState.drawOptions = newDrawOptions;
      scatterState.scatterplot.set(newDrawOptions);
      const prevView = scatterState.currentView;
      scatterState.scatterplot.draw(scatterState.normalizedPoints, { transition: true, transitionDuration: 200 }).then(() => {
        const stopAnim = startAnimatedAxesUpdate(scatterState.scatterplot, container, margin, width, height, scatterState.xDomain, scatterState.yDomain);
        setTimeout(() => {
          stopAnim();
          if (scatterState.viewSubscription) {
            const { handler, sub } = scatterState.viewSubscription;
            if (sub && typeof sub.unsubscribe === 'function') {
              sub.unsubscribe();
            } else if (scatterState.scatterplot._pubSub) {
              scatterState.scatterplot._pubSub.unsubscribe('view', handler);
            }
            scatterState.viewSubscription = null;
          }
          scatterState.scatterplot.refresh();
          // NEW: Manual sync for color transition
          setTimeout(() => {
            if (prevView && !prevView.xScale.domain().some(isNaN)) {
              updateAxesFromView(prevView, margin, width, height);
            } else {
              autoAdjustZoom(scatterState.originalPoints, scatterState.scatterplot, scatterState.xDomain, scatterState.yDomain);
            }
          }, 50);
        }, 200);
      });
      needsRefresh = true;
    }
  }

  if (x_label !== undefined) d3.select(container).select('.x-label')?.text(x_label);
  if (y_label !== undefined) d3.select(container).select('.y-label')?.text(y_label);

  if (logScale !== undefined && logScale !== scatterState.logScale) {
    const newNormalizedPoints = normalizePoints(scatterState.originalPoints, scatterState.xDomain, scatterState.yDomain, logScale);
    scatterState.normalizedPoints = newNormalizedPoints;
    scatterState.logScale = logScale;
    needsRedraw = true;
    const prevView = scatterState.currentView;
    const drawPromise = scatterState.scatterplot.draw(newNormalizedPoints, { transition: true });
    const stopAnim = startAnimatedAxesUpdate(scatterState.scatterplot, container, margin, width, height, scatterState.xDomain, scatterState.yDomain);
    drawPromise.then(() => {
      stopAnim();
      if (scatterState.viewSubscription) {
        const { handler, sub } = scatterState.viewSubscription;
        if (sub && typeof sub.unsubscribe === 'function') {
          sub.unsubscribe();
        } else if (scatterState.scatterplot._pubSub) {
          scatterState.scatterplot._pubSub.unsubscribe('view', handler);
        }
        scatterState.viewSubscription = null;
      }
      scatterState.scatterplot.refresh();
      // NEW: Manual sync for log toggle
      setTimeout(() => {
        if (prevView && !prevView.xScale.domain().some(isNaN)) {
          updateAxesFromView(prevView, margin, width, height);
        } else {
          autoAdjustZoom(scatterState.originalPoints, scatterState.scatterplot, scatterState.xDomain, scatterState.yDomain);
        }
      }, 50);
    });
  }

  if (needsRefresh && !needsRedraw) {
    scatterState.scatterplot.refresh();
    // Quick sync for size-only changes
    setTimeout(() => updateAxesFromView(scatterState.currentView || { xScale: d3.scaleLinear().domain([-1,1]), yScale: d3.scaleLinear().domain([-1,1]) }, margin, width, height), 0);
  }
});