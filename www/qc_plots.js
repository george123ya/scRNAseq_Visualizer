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
                    counters[i].text(`n=${count}`).style('display', null);
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

// This function needs to be outside the message handler to be accessible
const autoAdjustZoom = (points, scatterplot, xDomain, yDomain) => {
    if (!points || points.length === 0) {
      console.warn('No points to zoom to.');
      return;
    }
    if (!xDomain[0] || !xDomain[1] || !yDomain[0] || !yDomain[1]) {
      console.warn('Invalid domain values:', xDomain, yDomain);
      return;
    }
    const rangeX = xDomain[1] - xDomain[0] || 2;
    const rangeY = yDomain[1] - yDomain[0] || 2;
    const padX = rangeX * 0.1;
    const padY = rangeY * 0.1;
    const normalizedX = ((xDomain[0] - padX) - xDomain[0]) / (xDomain[1] - xDomain[0]) * 2 - 1;
    const normalizedWidth = (rangeX + 2 * padX) / (xDomain[1] - xDomain[0]) * 2;
    const normalizedY = ((yDomain[0] - padY) - yDomain[0]) / (yDomain[1] - yDomain[0]) * 2 - 1;
    const normalizedHeight = (rangeY + 2 * padY) / (yDomain[1] - yDomain[0]) * 2;
    if (isNaN(normalizedX) || isNaN(normalizedY) || isNaN(normalizedWidth) || isNaN(normalizedHeight)) {
      console.warn('Invalid zoom bounds:', { normalizedX, normalizedY, normalizedWidth, normalizedHeight });
      return;
    }
    scatterplot.zoomToArea({ x: normalizedX, y: normalizedY, width: normalizedWidth, height: normalizedHeight }, true);
};


Shiny.addCustomMessageHandler('scatterplotData', function(message) {
  // Clear the container before re-drawing
  const container = document.getElementById('scatterplot');
  container.innerHTML = '';

  const { data, x_label, y_label, color_by, color_levels, point_size, metadata } = message; // Add metadata to destructuring
  const points = data;
    
  const width = container.clientWidth;
  const height = container.clientHeight;
  const margin = { top: 20, right: 20, bottom: 50, left: 60 };
  
  const xDomain = d3.extent(points, d => d[0]);
  const yDomain = d3.extent(points, d => d[1]);
  
  const normalizedPoints = points.map(p => [
      ((p[0] - xDomain[0]) / (xDomain[1] - xDomain[0])) * 2 - 1,
      ((p[1] - yDomain[0]) / (yDomain[1] - yDomain[0])) * 2 - 1,
      p[2],
      p[3]
  ]);

  const canvas = document.createElement('canvas');
  canvas.id = 'scatterplotCanvas';
  canvas.style.position = 'absolute';
  canvas.style.top = margin.top + 'px';
  canvas.style.left = margin.left + 'px';
  canvas.style.width = (width - margin.left - margin.right) + 'px';
  canvas.style.height = (height - margin.top - margin.bottom) + 'px';
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
  svg.append('text')
    .attr('x', margin.left + (width - margin.left - margin.right) / 2)
    .attr('y', height - 10)
    .attr('text-anchor', 'middle')
    .attr('fill', 'black')
    .text(x_label);
  
  const yAxis = d3.axisLeft().ticks(6);
  const yAxisG = svg.append('g').attr('transform', `translate(${margin.left},0)`);
  svg.append('text')
    .attr('transform', 'rotate(-90)')
    .attr('x', -(margin.top + (height - margin.top - margin.bottom) / 2))
    .attr('y', 15)
    .attr('text-anchor', 'middle')
    .attr('fill', 'black')
    .text(y_label);
  
  svg.selectAll('.domain').attr('stroke', 'black').attr('stroke-width', 1.5);
  
  const tooltip = document.createElement('div');
  tooltip.id = 'scatterplotTooltip';
  tooltip.style.cssText = `position:absolute;background-color:rgba(0,0,0,0.8);color:white;padding:5px 10px;border-radius:4px;font-size:12px;pointer-events:none;z-index:1000;display:none;`;
  container.appendChild(tooltip);

  import('https://esm.sh/regl-scatterplot@1.14.1').then(module => {
    const createScatterplot = module.default;
    const scatterplot = createScatterplot({
        canvas: canvas,
        width: width - margin.left - margin.right,
        height: height - margin.top - margin.bottom,
        xScale: d3.scaleLinear().domain([-1, 1]).range([0, width - margin.left - margin.right]),
        yScale: d3.scaleLinear().domain([-1, 1]).range([height - margin.top - margin.bottom, 0]),
            pointSize: message.point_size // Direct use of the point_size value
    });

    const drawOptions = {
        colorBy: null,
        sizeBy: null
    };
    
    if (message.color_by !== "None" && message.color_levels && message.color_levels.length > 0) {
        drawOptions.colorBy = 'valueA';
        const colorScale = d3.scaleOrdinal(d3.schemeCategory10);
        const colorMap = colorScale.domain(d3.range(message.color_levels.length)).range();
        drawOptions.pointColor = colorMap;
    } else {
        drawOptions.pointColor = ['#000000'];
    }

    console.log(drawOptions);

    // We can now remove the sizeBy and pointSizeMap since pointSize is set directly
    // drawOptions.sizeBy = 'valueB'; // This is no longer needed
    // const pointSizeMap = d3.scaleLinear() ...
    // drawOptions.pointSize = [pointSizeMap(message.point_size)]; // Not needed

    scatterplot.set(drawOptions);
    scatterplot.draw(normalizedPoints);
      
     // Subscribe to events and update the tooltip
    scatterplot.subscribe('pointOver', pointIndex => {
        const point = data[pointIndex];
        const [px, py] = scatterplot.getScreenPosition(pointIndex);
        
        // Update tooltip content
        let tooltipContent = `X: ${point[0].toFixed(2)}<br>Y: ${point[1].toFixed(2)}`;
        
        // If metadata exists, add it to the tooltip
        if (metadata && metadata[pointIndex]) {
            tooltipContent += `<br>Metadata: ${metadata[pointIndex]}`;
        }
        
        tooltip.innerHTML = tooltipContent;
        tooltip.style.display = 'block';
        tooltip.style.left = (px + margin.left + 10) + 'px';
        tooltip.style.top = (py + margin.top + 10) + 'px';
    });

    scatterplot.subscribe('pointOut', () => {
    tooltip.style.display = 'none';
    });

    scatterplot.subscribe('view', event => {
    const newXDomain = [
        xDomain[0] + (event.xScale.domain()[0] + 1) / 2 * (xDomain[1] - xDomain[0]),
        xDomain[0] + (event.xScale.domain()[1] + 1) / 2 * (xDomain[1] - xDomain[0])
    ];
    const newYDomain = [
        yDomain[0] + (event.yScale.domain()[0] + 1) / 2 * (yDomain[1] - yDomain[0]),
        yDomain[0] + (event.yScale.domain()[1] + 1) / 2 * (yDomain[1] - yDomain[0])
    ];

    const newXScale = d3.scaleLinear().domain(newXDomain).range([margin.left, width - margin.right]);
    const newYScale = d3.scaleLinear().domain(newYDomain).range([height - margin.bottom, margin.top]);

    xAxis.scale(newXScale);
    yAxis.scale(newYScale);

    xAxisG.call(xAxis);
    yAxisG.call(yAxis);
    });

    setTimeout(() => {
    autoAdjustZoom(points, scatterplot, xDomain, yDomain);
    }, 100);
  });
});