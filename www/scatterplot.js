// Global variables
let scatterplot = null;
let genePlots = {}; // Store gene expression plots
let points = [];
let filteredPoints = [];
let annotations = {};
let currentColorBy = null;
let clusters = [];
let clusterColors = [];
let geneExprRange = [];
let useTransition = false;
let pointsXY = [];
let pointsMAGIC = [];
let isSyncMode = false;
let syncing = false;
let renderer = null; // Shared renderer for performance
let gene_number = 0; // Number of gene expression columns
let activeGeneExpressions = new Set(); // Track active gene expressions
let syncCooldownTimeout = null;

let geneExpressionCache = new Map(); // Cache for gene expression data
let geneDataMatrix = null; // Store the full gene expression matrix
let currentGeneNames = []; // Track current gene names
let geneDataInfo = null; // Store gene data dimensions and info


// Updated global variables to handle MAGIC data
let geneExpressionOriginal = new Map(); // Original gene expression data
let geneExpressionMAGIC = new Map();    // MAGIC-imputed gene expression data
let isMAGICActive = false;

let initScatterplotPromise = null; // Promise to track scatterplot initialization

let mainAnnotation = null;

// Viridis color palette for gene expression
const viridisColors = [
  '#440154', '#482777', '#3f4a8a', '#31678e', '#26838f', '#1f9d8a', 
  '#6cce5a', '#b6de2b', '#fee825', '#fcce25'
];

function createSyncButton() {
  if (document.getElementById('syncButton')) return;

  const buttonContainer = document.createElement('div');
  buttonContainer.id = 'syncButtonContainer';
  buttonContainer.style.cssText = `
    position: absolute;
    top: 10px;
    right: 10px;
    z-index: 1000;
  `;

  const syncButton = document.createElement('button');
  syncButton.id = 'syncButton';
  syncButton.textContent = 'Disable Sync Mode';
  syncButton.style.cssText = `
    padding: 8px 16px;
    background: #dc3545;
    color: white;
    border: none;
    border-radius: 4px;
    cursor: pointer;
    font-size: 14px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.2);
  `;
  syncButton.onclick = toggleSyncMode;
  buttonContainer.appendChild(syncButton);

  const plotGrid = document.getElementById('plotGrid');
  plotGrid.style.position = 'relative';
  plotGrid.appendChild(buttonContainer);

  // Explicitly start sync mode now
  setTimeout(() => {
    toggleSyncMode();  // will set text to Disable + red, and call setupSynchronization()
  }, 0);
}


function removeSyncButton() {
  const buttonContainer = document.getElementById('syncButtonContainer');
  if (buttonContainer) {
    buttonContainer.remove();
  }
}

function toggleSyncMode() {
  const button = document.getElementById('syncButton');
  console.log(`Toggling sync mode: ${isSyncMode ? 'Disabling' : 'Enabling'}`);

  if (!isSyncMode) {
    isSyncMode = true;
    button.textContent = 'Disable Sync Mode';
    button.style.background = '#dc3545';
    setupSynchronization();
  } else {
    isSyncMode = false;
    button.textContent = 'Enable Sync Mode';
    button.style.background = '#007bff';
    removeSynchronization();
  }
}

function createPlotGrid() {
  let plotGrid = document.getElementById('plotGrid');
  if (!plotGrid) {
    // Find the main canvas and its parent container
    const mainCanvas = document.getElementById('scatterplot_canvas');
    const parent = mainCanvas.parentElement;

    // Capture parent's explicit dimensions before clearing it
    const parentRect = parent.getBoundingClientRect();

    // Create grid container
    plotGrid = document.createElement('div');
    plotGrid.id = 'plotGrid';
    plotGrid.style.cssText = `
      width: ${parentRect.width}px;   /* lock width */
      height: ${parentRect.height}px; /* lock height */
      display: grid;
      gap: 5px;
      padding: 5px;
      box-sizing: border-box;
      overflow: hidden;
    `;

    // Replace parent's content with the grid
    parent.innerHTML = '';
    parent.appendChild(plotGrid);
    plotGrid.appendChild(mainCanvas);
  }
  return plotGrid;
}

function updatePlotLayout() {
  const plotGrid = createPlotGrid();
  const totalPlots = 1 + activeGeneExpressions.size;

  // Force max 2 rows
  const maxCols = 2;
  const cols = Math.min(maxCols, totalPlots);
  const rows = totalPlots <= maxCols ? 1 : 2;

  plotGrid.style.gridTemplateColumns = `repeat(${cols}, 1fr)`;
  plotGrid.style.gridTemplateRows = `repeat(${rows}, 1fr)`;

  // --- Compute cell sizes (full space since legends are overlays) ---
  const gridRect = plotGrid.getBoundingClientRect();
  const cellWidth = gridRect.width / cols;
  const cellHeight = gridRect.height / rows;
  const dpr = window.devicePixelRatio || 1;

  const allCanvases = [
    document.getElementById('scatterplot_canvas'),
    ...Object.keys(genePlots).map(id => document.getElementById(`gene_canvas_${id}`))
  ].filter(Boolean);

  allCanvases.forEach(canvas => {
    canvas.style.width = '100%';
    canvas.style.height = '100%';
    canvas.width = cellWidth * dpr;
    canvas.height = cellHeight * dpr;
  });

  // Resize plots
  setTimeout(() => {
    if (scatterplot) scatterplot.set({ width: cellWidth, height: cellHeight });
    Object.values(genePlots).forEach(plot => plot.set({ width: cellWidth, height: cellHeight }));
  }, 50);
}

function createPlotTitle(canvas, title, plotId) {
  // Ensure each canvas is wrapped in its own relative container
  let container = canvas.parentElement;
  if (!container || !container.classList.contains('canvas-container')) {
    container = document.createElement('div');
    container.className = 'canvas-container';
    container.style.cssText = `
      position: relative;
      width: 100%;
      height: 100%;
      display: flex;
      align-items: stretch;
    `;
    canvas.parentElement.insertBefore(container, canvas);
    container.appendChild(canvas);
  }

  // Remove any old title for this canvas
  const oldTitle = container.querySelector('.plot-title');
  if (oldTitle) oldTitle.remove();

  // Create title overlay
  const titleOverlay = document.createElement('div');
  titleOverlay.className = 'plot-title';
  titleOverlay.style.cssText = `
    position: absolute;
    top: 5px;
    left: 5px;
    background: rgba(255, 255, 255, 0.9);
    padding: 4px 8px;
    border-radius: 4px;
    font-size: 12px;
    font-weight: bold;
    color: #333;
    pointer-events: auto;
    z-index: 10;
    cursor: pointer;
  `;
  titleOverlay.textContent = title || `Plot ${plotId}`;
  titleOverlay.title = 'Click to download plot as PNG or SVG';

  // Create dropdown for format selection
  const formatDropdown = document.createElement('div');
  formatDropdown.className = 'format-dropdown';
  formatDropdown.style.cssText = `
    position: absolute;
    top: 100%;
    left: 0;
    background: rgba(255, 255, 255, 0.95);
    border: 1px solid #ddd;
    border-radius: 4px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.15);
    z-index: 11;
    display: none;
    padding: 5px;
  `;

  const pngOption = document.createElement('div');
  pngOption.textContent = 'Download as PNG';
  pngOption.style.cssText = `
    padding: 5px 10px;
    cursor: pointer;
    font-size: 12px;
    color: #333;
  `;
  pngOption.addEventListener('mouseover', () => { pngOption.style.background = '#f0f0f0'; });
  pngOption.addEventListener('mouseout', () => { pngOption.style.background = 'none'; });

  const svgOption = document.createElement('div');
  svgOption.textContent = 'Download as SVG';
  svgOption.style.cssText = `
    padding: 5px 10px;
    cursor: pointer;
    font-size: 12px;
    color: #333;
  `;
  svgOption.addEventListener('mouseover', () => { svgOption.style.background = '#f0f0f0'; });
  svgOption.addEventListener('mouseout', () => { svgOption.style.background = 'none'; });

  formatDropdown.appendChild(pngOption);
  formatDropdown.appendChild(svgOption);
  titleOverlay.appendChild(formatDropdown);

  // Toggle dropdown on click
  titleOverlay.addEventListener('click', (e) => {
    e.stopPropagation();
    formatDropdown.style.display = formatDropdown.style.display === 'block' ? 'none' : 'block';
  });

  // Handle format selection
  const handleSelection = (format) => {
    const isMainPlot = plotId === 'main';
    const legendKey = isMainPlot ? 'main_categorical' : `${plotId}_gene`;
    const legendData = currentLegendData.get(legendKey);
    console.log('Selected format:', format, 'Plot ID:', plotId, 'Legend key:', legendKey, 'Legend data:', legendData);
    downloadPlot(canvas, plotId, legendData, format, title);
    formatDropdown.style.display = 'none';
  };

  pngOption.addEventListener('click', () => handleSelection('png'));
  svgOption.addEventListener('click', () => handleSelection('svg'));

  // Close dropdown when clicking outside
  document.addEventListener('click', (e) => {
    if (!titleOverlay.contains(e.target)) {
      formatDropdown.style.display = 'none';
    }
  });

  container.appendChild(titleOverlay);
}

async function downloadPlot(canvas, plotId, legendData, format, titleText) {
  console.log('downloadPlot called with:', { plotId, format, legendData, titleText });
  const plotContainer = canvas.parentElement;
  const isMainPlot = plotId === 'main';
  const fileName = `${isMainPlot ? 'main_plot' : `gene_plot_${plotId}`}.${format}`;

  // Helper function to escape SVG text
  const escapeSvgText = (text) => {
    if (!text) return '';
    return text.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;').replace(/'/g, '&apos;');
  };

  if (format === 'png') {
    const tempCanvas = document.createElement('canvas');
    const ctx = tempCanvas.getContext('2d');

    const canvasWidth = canvas.width;
    const canvasHeight = canvas.height;

    const titleHeight = titleText ? 30 : 0;
    let legendWidth = 0;
    let legendHeight = 0;
    const isCategorical = legendData && Array.isArray(legendData.names) && Array.isArray(legendData.colors) && legendData.names.length === legendData.colors.length;

    console.log('isCategorical:', isCategorical, 'legendData:', legendData);

    if (isCategorical) {
      const itemHeight = 20;
      const legendTitleHeight = 30;
      legendHeight = legendData.names.length * itemHeight + legendTitleHeight;
      legendWidth = 200;
    } else if (legendData && typeof legendData.geneName === 'string') {
      legendHeight = 120;
      legendWidth = 150;
    } else {
      console.warn('No valid legend data for plot:', plotId, 'Legend data:', legendData);
    }

    tempCanvas.width = canvasWidth + (legendWidth ? legendWidth + 20 : 0);
    tempCanvas.height = Math.max(canvasHeight, legendHeight || 0) + titleHeight;

    ctx.fillStyle = '#FFFFFF';
    ctx.fillRect(0, 0, tempCanvas.width, tempCanvas.height);

    if (titleText) {
      ctx.fillStyle = '#333';
      ctx.font = 'bold 12px sans-serif';
      ctx.textAlign = 'left';
      ctx.fillText(titleText, 10, 20);
    }

    try {
      ctx.drawImage(canvas, 0, titleHeight, canvasWidth, canvasHeight);
    } catch (error) {
      console.error('Failed to draw canvas image:', error);
      tempCanvas.remove();
      return;
    }

    if (legendData) {
      if (isCategorical) {
        ctx.fillStyle = 'rgba(255, 255, 255, 0.95)';
        ctx.fillRect(canvasWidth + 10, titleHeight, legendWidth, legendHeight);
        ctx.strokeStyle = '#ddd';
        ctx.strokeRect(canvasWidth + 10, titleHeight, legendWidth, legendHeight);

        ctx.fillStyle = '#333';
        ctx.font = 'bold 12px sans-serif';
        ctx.textAlign = 'center';
        const legendTitle = legendData.annotationName ? legendData.annotationName.replace(/_/g, ' ') : 'Legend';
        ctx.fillText(legendTitle, canvasWidth + 10 + legendWidth / 2, titleHeight + 20);

        ctx.font = '10px sans-serif';
        ctx.textAlign = 'left';
        legendData.names.forEach((name, i) => {
          const y = titleHeight + 40 + i * 20;
          const color = legendData.colors[i] || '#000000';
          if (!/^#[0-9A-F]{6}$/i.test(color)) {
            console.warn(`Invalid legend color at index ${i}: ${color}`);
            ctx.fillStyle = '#000000';
          } else {
            ctx.fillStyle = color;
          }
          ctx.fillRect(canvasWidth + 20, y - 8, 12, 12);
          ctx.strokeStyle = 'rgba(0,0,0,0.2)';
          ctx.strokeRect(canvasWidth + 20, y - 8, 12, 12);
          ctx.fillStyle = '#333';
          ctx.fillText(name || `Item ${i}`, canvasWidth + 38, y);
        });
      } else if (typeof legendData.geneName === 'string') {
        ctx.fillStyle = 'rgba(255, 255, 255, 0.95)';
        ctx.fillRect(canvasWidth + 10, titleHeight, legendWidth, legendHeight);
        ctx.strokeStyle = '#ddd';
        ctx.strokeRect(canvasWidth + 10, titleHeight, legendWidth, legendHeight);

        ctx.fillStyle = '#333';
        ctx.font = 'bold 12px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText(legendData.geneName, canvasWidth + 10 + legendWidth / 2, titleHeight + 20);

        const gradientHeight = 80;
        const gradient = ctx.createLinearGradient(0, titleHeight + gradientHeight + 30, 0, titleHeight + 30);
        gradient.addColorStop(0, viridisColors[0] || '#000000');
        gradient.addColorStop(0.5, viridisColors[Math.floor(viridisColors.length / 2)] || '#888888');
        gradient.addColorStop(1, viridisColors[viridisColors.length - 1] || '#FFFFFF');
        ctx.fillStyle = gradient;
        ctx.fillRect(canvasWidth + 20, titleHeight + 30, 20, gradientHeight);
        ctx.strokeStyle = '#ccc';
        ctx.strokeRect(canvasWidth + 20, titleHeight + 30, 20, gradientHeight);

        ctx.fillStyle = '#333';
        ctx.font = '9px sans-serif';
        ctx.textAlign = 'left';
        ctx.fillText(legendData.maxVal || 'Max', canvasWidth + 50, titleHeight + 35);
        ctx.fillText(legendData.midVal || 'Mid', canvasWidth + 50, titleHeight + 30 + gradientHeight / 2);
        ctx.fillText(legendData.minVal || 'Min', canvasWidth + 50, titleHeight + 30 + gradientHeight - 5);

        ctx.fillStyle = '#666';
        ctx.font = '9px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText(`Data range: ${legendData.realMin || '0'} – ${legendData.realMax || '0'}`, canvasWidth + 10 + legendWidth / 2, titleHeight + gradientHeight + 40);
      }
    } else {
      console.warn('No legend data provided for PNG download:', plotId);
    }

    try {
      const link = document.createElement('a');
      link.download = fileName;
      link.href = tempCanvas.toDataURL('image/png');
      link.click();
      link.remove();
      tempCanvas.remove();
    } catch (error) {
      console.error('Failed to download PNG:', error);
      tempCanvas.remove();
    }
  } else if (format === 'svg') {
    let svgContent = '<svg xmlns="http://www.w3.org/2000/svg" ';

    const canvasWidth = canvas.width;
    const canvasHeight = canvas.height;

    const titleHeight = titleText ? 30 : 0;
    let legendWidth = 0;
    let legendHeight = 0;
    const isCategorical = legendData && Array.isArray(legendData.names) && Array.isArray(legendData.colors) && legendData.names.length === legendData.colors.length;

    console.log('isCategorical:', isCategorical, 'legendData:', legendData);

    if (isCategorical) {
      const itemHeight = 20;
      const legendTitleHeight = 30;
      legendHeight = legendData.names.length * itemHeight + legendTitleHeight;
      legendWidth = 200;
    } else if (legendData && typeof legendData.geneName === 'string') {
      legendHeight = 120;
      legendWidth = 150;
    }

    const svgWidth = canvasWidth + (legendWidth ? legendWidth + 20 : 0);
    const svgHeight = Math.max(canvasHeight, legendHeight || 0) + titleHeight;
    svgContent += `width="${svgWidth}" height="${svgHeight}">`;

    if (titleText) {
      svgContent += `
        <text x="10" y="20" font-family="sans-serif" font-size="12px" font-weight="bold" fill="#333">
          ${escapeSvgText(titleText)}
        </text>
      `;
    }

    const plot = isMainPlot ? scatterplot : genePlots[plotId];
    console.log('Plot object:', plot, 'genePlots:', genePlots, 'plotId:', plotId);
    if (!plot) {
      console.error('Plot not found for ID:', plotId);
      console.log('Available genePlots keys:', Object.keys(genePlots));
      alert('Error: Plot not found. Using rasterized image in SVG.');
      const imageData = canvas.toDataURL('image/png');
      svgContent += `<image x="0" y="${titleHeight}" width="${canvasWidth}" height="${canvasHeight}" href="${imageData}" />`;
    } else {
      const pointsToRender = isMainPlot ? filteredPoints : filteredPoints;
      let colors = [];
      if (isCategorical) {
        colors = legendData.colors;
      } else if (legendData && typeof legendData.geneName === 'string') {
        colors = viridisColors;
      } else {
        colors = ['#000000'];
      }

      console.log('Points to render length:', pointsToRender.length, 'Colors length:', colors.length);
      console.log('Points sample:', pointsToRender.slice(0, 5));

      let view;
      try {
        view = plot.get('cameraView') || { x: 0, y: 0, scale: 1 };
        console.log('Camera view:', view);
      } catch (error) {
        console.warn('Could not get cameraView, using default scaling:', error);
        view = { x: 0, y: 0, scale: 1 };
      }

      svgContent += `<g transform="translate(0,${titleHeight})">`;
      if (pointsToRender && pointsToRender.length > 0) {
        const maxPoints = 5000;
        if (pointsToRender.length > maxPoints) {
          console.warn(`Dataset too large (${pointsToRender.length} points). Using rasterized image for SVG.`);
          alert(`Dataset too large (${pointsToRender.length} points). Rendering as rasterized image in SVG.`);
          const imageData = canvas.toDataURL('image/png');
          svgContent += `<image x="0" y="0" width="${canvasWidth}" height="${canvasHeight}" href="${imageData}" />`;
        } else {
          const step = Math.ceil(pointsToRender.length / maxPoints);
          console.log(`Rendering ${Math.min(pointsToRender.length, maxPoints)} points for SVG (total: ${pointsToRender.length}, step: ${step})`);

          try {
            pointsToRender.forEach((point, i) => {
              if (i % step !== 0) return;
              const x = (point[0] * view.scale + view.x + 1) * canvasWidth / 2;
              const y = canvasHeight - ((point[1] * view.scale + view.y + 1) * canvasHeight / 2);
              if (isNaN(x) || isNaN(y) || !isFinite(x) || !isFinite(y)) {
                console.warn(`Invalid point coordinates at index ${i}: [${point[0]}, ${point[1]}]`);
                return;
              }
              let color;
              if (isCategorical) {
                color = colors[i % colors.length] || '#000000';
                if (!/^#[0-9A-F]{6}$/i.test(color)) {
                  console.warn(`Invalid color at index ${i}: ${color}`);
                  color = '#000000';
                }
              } else if (legendData && typeof legendData.geneName === 'string' && point[2] !== undefined) {
                const normalizedValue = (point[2] - (legendData.realMin || 0)) / ((legendData.realMax || 1) - (legendData.realMin || 0));
                const colorIndex = Math.min(colors.length - 1, Math.max(0, Math.floor(normalizedValue * colors.length)));
                color = colors[colorIndex] || '#000000';
                if (!/^#[0-9A-F]{6}$/i.test(color)) {
                  console.warn(`Invalid gene color at index ${i}: ${color}`);
                  color = '#000000';
                }
              } else {
                color = colors[0] || '#000000';
              }
              const pointSize = plot.get('pointSize') || 3;
              svgContent += `<circle cx="${x}" cy="${y}" r="${pointSize}" fill="${color}" opacity="0.8" />`;
            });
          } catch (error) {
            console.error('Error rendering SVG points:', error);
            alert('Error rendering SVG points. Using rasterized image.');
            svgContent += `<image x="0" y="0" width="${canvasWidth}" height="${canvasHeight}" href="${canvas.toDataURL('image/png')}" />`;
          }
        }
      } else {
        console.warn('No points available to render in SVG for plot:', plotId);
        svgContent += `<text x="10" y="20" font-family="sans-serif" font-size="12px" fill="#333">No points available</text>`;
      }
      svgContent += `</g>`;
    }

    if (legendData) {
      if (isCategorical) {
        const itemHeight = 20;
        svgContent += `
          <g transform="translate(${canvasWidth + 10},${titleHeight})">
            <rect x="0" y="0" width="${legendWidth}" height="${legendHeight}" fill="rgba(255,255,255,0.95)" stroke="#ddd" />
            <text x="${legendWidth / 2}" y="20" font-family="sans-serif" font-size="12px" font-weight="bold" text-anchor="middle">
              ${escapeSvgText(legendData.annotationName ? legendData.annotationName.replace(/_/g, ' ') : 'Legend')}
            </text>
        `;
        legendData.names.forEach((name, i) => {
          const y = 40 + i * itemHeight;
          const color = legendData.colors[i] || '#000000';
          if (!/^#[0-9A-F]{6}$/i.test(color)) {
            console.warn(`Invalid legend color at index ${i}: ${color}`);
            svgContent += `
              <rect x="10" y="${y - 8}" width="12" height="12" fill="#000000" stroke="rgba(0,0,0,0.2)" />
              <text x="30" y="${y}" font-family="sans-serif" font-size="10px" fill="#333">${escapeSvgText(name || `Item ${i}`)}</text>
            `;
          } else {
            svgContent += `
              <rect x="10" y="${y - 8}" width="12" height="12" fill="${color}" stroke="rgba(0,0,0,0.2)" />
              <text x="30" y="${y}" font-family="sans-serif" font-size="10px" fill="#333">${escapeSvgText(name || `Item ${i}`)}</text>
            `;
          }
        });
        svgContent += `</g>`;
      } else if (typeof legendData.geneName === 'string') {
        const gradientId = `gradient_${plotId}`;
        svgContent += `
          <defs>
            <linearGradient id="${gradientId}" x1="0" y1="1" x2="0" y2="0">
              <stop offset="0%" stop-color="${viridisColors[0] || '#000000'}" />
              <stop offset="50%" stop-color="${viridisColors[Math.floor(viridisColors.length / 2)] || '#888888'}" />
              <stop offset="100%" stop-color="${viridisColors[viridisColors.length - 1] || '#FFFFFF'}" />
            </linearGradient>
          </defs>
          <g transform="translate(${canvasWidth + 10},${titleHeight})">
            <rect x="0" y="0" width="${legendWidth}" height="${legendHeight}" fill="rgba(255,255,255,0.95)" stroke="#ddd" />
            <text x="${legendWidth / 2}" y="20" font-family="sans-serif" font-size="12px" font-weight="bold" text-anchor="middle">
              ${escapeSvgText(legendData.geneName)}
            </text>
            <rect x="20" y="30" width="20" height="80" fill="url(#${gradientId})" stroke="#ccc" />
            <text x="50" y="35" font-family="sans-serif" font-size="9px" fill="#333">${escapeSvgText(legendData.maxVal || 'Max')}</text>
            <text x="50" y="${30 + 80 / 2}" font-family="sans-serif" font-size="9px" fill="#333">${escapeSvgText(legendData.midVal || 'Mid')}</text>
            <text x="50" y="${30 + 80 - 5}" font-family="sans-serif" font-size="9px" fill="#333">${escapeSvgText(legendData.minVal || 'Min')}</text>
            <text x="${legendWidth / 2}" y="${30 + 80 + 15}" font-family="sans-serif" font-size="9px" fill="#666" text-anchor="middle">
              Data range: ${escapeSvgText(`${legendData.realMin || '0'} – ${legendData.realMax || '0'}`)}
            </text>
          </g>
        `;
      }
    } else {
      console.warn('No legend data provided for SVG download:', plotId);
    }

    svgContent += `</svg>`;

    try {
      console.log('SVG content length:', svgContent.length);
      const blob = new Blob([svgContent], { type: 'image/svg+xml' });
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.download = fileName;
      link.href = url;
      link.click();
      URL.revokeObjectURL(url);
      link.remove();
    } catch (error) {
      console.error('Failed to download SVG:', error);
      alert('Error downloading SVG: Invalid format or too large. Try PNG format.');
    }
  }
}

let lastClicked = null;

// Add these global variables at the top of your file
let isUpdatingVisibility = false;
let currentLegendData = new Map(); // Track current legend data to avoid unnecessary recreations

// Add these helper functions first (if not already added)
let customColors = new Map(); // Add to globals if not already there

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
  // Dynamically import Pickr
  let Pickr;
  try {
    const pickrModule = await import('https://esm.sh/@simonwep/pickr@1.9.0');
    Pickr = pickrModule.default;
  } catch (error) {
    console.error('Failed to load Pickr:', error);
    return null; // Fallback to prevent breaking the app
  }

  // Load Pickr CSS dynamically (only once)
  if (!document.querySelector('link[href*="pickr"]')) {
    const link = document.createElement('link');
    link.rel = 'stylesheet';
    link.href = 'https://esm.sh/@simonwep/pickr@1.9.0/dist/themes/nano.min.css';
    document.head.appendChild(link);
  }

  // Create container for the color picker
  const container = document.createElement('div');
  container.style.cssText = `
    width: 100%;
    height: 100%;
    border-radius: 2px;
  `;

  // Append container to colorBoxContainer (already in DOM) before Pickr initialization
  colorBoxContainer.innerHTML = '';
  colorBoxContainer.appendChild(container);

  // Initialize Pickr
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
    position: 'bottom-middle' // Ensure popup stays in view
  });

  // Live color changes
  pickr.on('change', (color) => {
    const hex = color.toHEXA().toString().toUpperCase();
    onColorChange(hex, false);
  });

  // Final color selection
  pickr.on('save', (color) => {
    const hex = color ? color.toHEXA().toString().toUpperCase() : initialColor;
    onColorChange(hex, true);
    pickr.hide();
  });

  // Picker close
  pickr.on('hide', () => {
    onClose();
  });

  // Open the picker immediately
  pickr.show();

  return container;
}

// Rewritten createPlotLegend function
function createPlotLegend(plotId, legendData, type = 'categorical') {
  // Generate unique legend key
  const legendKey = `${plotId}_${type}`;
  const existingData = currentLegendData.get(legendKey);

  // Handle visibility updates for existing categorical legend
  if (isUpdatingVisibility && type === 'categorical' && existingData) {
    const container = document.getElementById(plotId === 'main' ? 'scatterplot_canvas' : `gene_canvas_${plotId}`);
    if (!container) return;

    const plotContainer = container.parentElement;
    if (!plotContainer) return;

    const legend = plotContainer.querySelector('.plot-legend');
    if (legend) {
      const itemsContainer = legend.children[1];
      if (itemsContainer && itemsContainer.classList) {
        itemsContainer.querySelectorAll('.legend-item').forEach((el, idx) => {
          el.style.opacity = legendData.visible.has(idx) ? '1' : '0.4';
        });
        return;
      }
    }
  }

  // Store legend data
  currentLegendData.set(legendKey, legendData);

  // Get plot container
  const container = document.getElementById(plotId === 'main' ? 'scatterplot_canvas' : `gene_canvas_${plotId}`);
  if (!container) return;

  const plotContainer = container.parentElement;
  if (!plotContainer) return;

  // Remove existing legend
  const existingLegend = plotContainer.querySelector('.plot-legend');
  if (existingLegend) existingLegend.remove();

  // Ensure plot container is positioned
  plotContainer.style.position = 'relative';

  // Create legend container
  const legendContainer = document.createElement('div');
  legendContainer.className = 'plot-legend';
  const legendId = `legend_${plotId}_${type}`;
  legendContainer.style.cssText = `
    position: absolute;
    top: 8px;
    right: 8px;
    width: 150px;
    background: rgba(255, 255, 255, 0.95);
    border: 1px solid #ddd;
    font-size: 11px;
    padding: 0;
    z-index: 10;
    border-radius: 4px;
    pointer-events: auto;
    box-shadow: 0 2px 8px rgba(0,0,0,0.15);
    min-width: 120px;
    min-height: 100px;
  `;

  if (type === 'gene') {
    // Gene legend (gradient-based)
    const { minVal, maxVal, midVal, realMin, realMax, geneName } = legendData;
    const contentContainer = document.createElement('div');
    contentContainer.style.cssText = 'padding: 6px;';
    contentContainer.innerHTML = `
      <div style="font-weight: bold; margin-bottom: 6px; text-align: center; margin-top: 4px;">
        ${geneName}
      </div>
      <div style="display: flex; justify-content: center; align-items: center;">
        <div style="
          width: 20px;
          height: 80px;
          background: linear-gradient(to top, 
            ${viridisColors[0]}, 
            ${viridisColors[Math.floor(viridisColors.length/2)]}, 
            ${viridisColors[viridisColors.length-1]}
          );
          border: 1px solid #ccc;
          margin-right: 6px;">
        </div>
        <div style="
          display: flex;
          flex-direction: column;
          justify-content: space-between;
          height: 80px;
          font-size: 9px;
          color: #333;
          text-align: left;">
          <div>${maxVal}</div>
          <div>${midVal}</div>
          <div>${minVal}</div>
        </div>
      </div>
      <div style="
        font-size: 9px;
        text-align: center;
        margin-top: 4px;
        color: #666;">
        Data range: ${realMin} – ${realMax}
      </div>
    `;
    legendContainer.appendChild(contentContainer);
  } else {
    // Categorical legend
    const { names, colors, visible, annotationName } = legendData;

    // Title
    const titleDiv = document.createElement('div');
    titleDiv.style.cssText = `
      font-weight: bold;
      margin: 6px;
      text-align: center;
      border-bottom: 1px solid #eee;
      padding-bottom: 4px;
    `;
    titleDiv.textContent = annotationName.replace(/_/g, ' ');
    legendContainer.appendChild(titleDiv);

    // Items container
    const itemsContainer = document.createElement('div');
    itemsContainer.style.cssText = `
      max-height: 150px;
      overflow-y: auto;
      margin: 0 6px 6px 6px;
    `;

    names.forEach((name, i) => {
      const item = document.createElement('div');
      item.className = 'legend-item';
      item.style.cssText = `
        display: flex;
        align-items: center;
        padding: 3px;
        border-radius: 3px;
        user-select: none;
        opacity: ${visible.has(i) ? '1' : '0.4'};
        transition: opacity 0.2s;
      `;

      // Color box container
      const colorBoxContainer = document.createElement('div');
      colorBoxContainer.style.cssText = `
        flex: 0 0 12px;
        height: 12px;
        margin-right: 6px;
        position: relative;
        border-radius: 2px;
        overflow: hidden;
      `;

      // Color box
      const colorBox = document.createElement('div');
      const customColor = customColors.get(`${annotationName}_${i}`) || colors[i];
      colorBox.style.cssText = `
        width: 100%;
        height: 100%;
        background: ${customColor};
        border-radius: 2px;
        cursor: pointer;
        border: 1px solid rgba(0,0,0,0.2);
        transition: transform 0.2s;
        box-sizing: border-box;
      `;

      // Hover effects
      colorBox.addEventListener('mouseenter', () => {
        if (!colorBox.dataset.editing) {
          colorBox.style.transform = 'scale(1.1)';
          colorBox.style.boxShadow = '0 2px 4px rgba(0,0,0,0.3)';
        }
      });
      colorBox.addEventListener('mouseleave', () => {
        if (!colorBox.dataset.editing) {
          colorBox.style.transform = 'scale(1)';
          colorBox.style.boxShadow = 'none';
        }
      });

      // Color picker handler
      let originalColor = rgbToHex(customColor);
      colorBox.addEventListener('click', async (e) => {
        e.preventDefault();
        e.stopPropagation();

        if (colorBox.dataset.editing === 'true') return;
        colorBox.dataset.editing = 'true';
        originalColor = rgbToHex(colorBox.style.backgroundColor || customColor);

        try {
          const colorPicker = await createInlineColorPicker(
            originalColor,
            (newColor, isFinal) => {
              colorBox.style.backgroundColor = newColor;

              if (isFinal) {
                customColors.set(`${annotationName}_${i}`, newColor);
                if (mainAnnotation && mainAnnotation.colors) {
                  mainAnnotation.colors[i] = newColor;
                }
                Shiny.setInputValue('colorChange', {
                  annotation: annotationName,
                  categoryIndex: i,
                  categoryName: name,
                  newColor: newColor,
                  allColors: mainAnnotation ? mainAnnotation.colors : colors,
                  timestamp: Date.now()
                }, { priority: 'event' });
              } else {
                if (mainAnnotation && mainAnnotation.colors) {
                  mainAnnotation.colors[i] = newColor;
                }
              }

              if (scatterplot && mainAnnotation && mainAnnotation.colorBy === annotationName) {
                scatterplot.set({ pointColor: mainAnnotation.colors });
                scatterplot.draw(filteredPoints, { transition: false });
              }
            },
            () => {
              colorBox.dataset.editing = 'false';
              colorBoxContainer.innerHTML = '';
              colorBoxContainer.appendChild(colorBox);
              colorBox.style.transform = 'scale(1)';
              colorBox.style.boxShadow = 'none';
            },
            colorBoxContainer // Pass colorBoxContainer to append directly
          );

          if (!colorPicker) {
            colorBox.dataset.editing = 'false';
          }
        } catch (error) {
          console.error('Failed to create color picker:', error);
          colorBox.dataset.editing = 'false';
        }
      });

      // Text span for visibility toggle
      const textSpan = document.createElement('span');
      textSpan.style.cssText = `
        font-size: 10px;
        color: #333;
        flex: 1;
        word-break: break-word;
        overflow-wrap: anywhere;
        white-space: normal;
        cursor: pointer;
      `;
      textSpan.textContent = name;

      // Visibility toggle handler
      textSpan.addEventListener('click', (e) => {
        e.preventDefault();
        e.stopPropagation();

        const scrollTop = itemsContainer.scrollTop;
        const currentlyVisible = Array.from(visible);
        isUpdatingVisibility = true;

        // Handle visibility toggling (unchanged)
        if (e.shiftKey && lastClicked !== null) {
          const start = Math.min(lastClicked, i);
          const end = Math.max(lastClicked, i);
          if (e.ctrlKey || e.metaKey) {
            for (let idx = start; idx <= end; idx++) {
              visible.add(idx);
            }
          } else {
            visible.clear();
            for (let idx = start; idx <= end; idx++) {
              visible.add(idx);
            }
          }
        } else if (e.ctrlKey || e.metaKey) {
          if (visible.has(i)) {
            visible.delete(i);
          } else {
            visible.add(i);
          }
          if (visible.size === 0) {
            names.forEach((_, idx) => visible.add(idx));
          }
        } else {
          if (currentlyVisible.length === names.length) {
            visible.clear();
            visible.add(i);
          } else if (currentlyVisible.length === 1 && currentlyVisible[0] === i) {
            visible.clear();
            names.forEach((_, idx) => visible.add(idx));
          } else {
            visible.clear();
            visible.add(i);
          }
        }

        lastClicked = i;
        itemsContainer.querySelectorAll('.legend-item').forEach((el, idx) => {
          el.style.opacity = visible.has(idx) ? '1' : '0.4';
        });
        itemsContainer.scrollTop = scrollTop;

        // Send visible category indices to Shiny
        Shiny.setInputValue('visibleCategories', {
          plotId: plotId,
          annotationName: annotationName,
          visibleIndices: Array.from(visible),
          timestamp: Date.now()
        }, { priority: 'event' });

        redrawAllPlots(mainAnnotation).finally(() => {
          isUpdatingVisibility = false;
          setTimeout(() => {
            itemsContainer.scrollTop = scrollTop;
          }, 0);
        });
      });

      // Assemble legend item
      colorBoxContainer.appendChild(colorBox);
      item.appendChild(colorBoxContainer);
      item.appendChild(textSpan);
      itemsContainer.appendChild(item);
    });

    legendContainer.appendChild(itemsContainer);

    // Send initial visible category indices to Shiny
    Shiny.setInputValue('visibleCategories', {
      plotId: plotId,
      annotationName: annotationName,
      visibleIndices: Array.from(visible),
      timestamp: Date.now()
    }, { priority: 'event' });
  }

  // Append and make draggable/resizable
  plotContainer.appendChild(legendContainer);
  makeLegendDraggableAndResizable(legendContainer, legendId);
}

// Add these global variables for legend positioning
let legendPositions = new Map(); // Store positions for each legend
let legendSizes = new Map(); // Store sizes for each legend

// Function to update legend content layout when resized
// Function to update legend content layout when resized
function updateLegendContentLayout(legendContainer, width, height, type) {
  const dragHandleHeight = 20;
  const padding = 12;
  const availableHeight = height - dragHandleHeight - padding;
  const availableWidth = width - padding;

  if (type === 'gene') {
    // Update gene legend layout
    const contentContainer = legendContainer.querySelector('div:not(.legend-drag-handle):not(.legend-resize-handle)');
    if (contentContainer) {
      // Scale the colorbar and text based on available space
      const colorbarContainer = contentContainer.querySelector('div:nth-child(2)');
      if (colorbarContainer) {
        const colorbar = colorbarContainer.querySelector('div:first-child');
        const textContainer = colorbarContainer.querySelector('div:last-child');
        
        if (colorbar && textContainer) {
          // Scale colorbar height based on available space
          const maxColorbarHeight = Math.min(120, availableHeight - 60); // Leave space for title and range text
          const colorbarHeight = Math.max(60, maxColorbarHeight);
          
          colorbar.style.height = colorbarHeight + 'px';
          textContainer.style.height = colorbarHeight + 'px';
          
          // Adjust font sizes for smaller widths
          if (availableWidth < 130) {
            contentContainer.style.fontSize = '10px';
            const rangeText = contentContainer.querySelector('div:last-child');
            if (rangeText) rangeText.style.fontSize = '8px';
          } else {
            contentContainer.style.fontSize = '11px';
            const rangeText = contentContainer.querySelector('div:last-child');
            if (rangeText) rangeText.style.fontSize = '9px';
          }
        }
      }
    }
  } else {
    // Update categorical legend layout
    // Find the items container more reliably
    let itemsContainer = null;
    
    // Method 1: Look for container with legend-item children
    const containers = legendContainer.querySelectorAll('div');
    for (let container of containers) {
      if (container.classList.contains('legend-drag-handle') || 
          container.classList.contains('legend-resize-handle')) {
        continue;
      }
      
      // Check if this container has legend-item children
      if (container.querySelectorAll('.legend-item').length > 0) {
        itemsContainer = container;
        break;
      }
    }
    
    // Method 2: If not found, look for container with overflow-y auto or scroll
    if (!itemsContainer) {
      for (let container of containers) {
        if (container.classList.contains('legend-drag-handle') || 
            container.classList.contains('legend-resize-handle')) {
          continue;
        }
        
        const computedStyle = window.getComputedStyle(container);
        if (computedStyle.overflowY === 'auto' || computedStyle.overflowY === 'scroll') {
          itemsContainer = container;
          break;
        }
      }
    }
    
    console.log('Resizing legend - Available height:', availableHeight);
    console.log('Items container found:', !!itemsContainer);
    
    if (itemsContainer) {
      // Calculate the height for items container
      // Account for title div (if present)
      const titleDiv = legendContainer.querySelector('div:not(.legend-drag-handle):not(.legend-resize-handle)');
      const titleHeight = (titleDiv && titleDiv !== itemsContainer) ? titleDiv.offsetHeight + 10 : 20; // 10px for margins
      
      const newHeight = Math.max(50, availableHeight - titleHeight);
      console.log('Setting items container height to:', newHeight);
      
      // Apply the new height and ensure scrolling
      itemsContainer.style.maxHeight = newHeight + 'px';
      itemsContainer.style.height = newHeight + 'px';
      itemsContainer.style.overflowY = 'auto';
      itemsContainer.style.boxSizing = 'border-box';
      
      // Force a layout recalculation
      itemsContainer.offsetHeight;
      
      // Adjust item layout for different widths
      const items = itemsContainer.querySelectorAll('.legend-item');
      items.forEach(item => {
        const colorBox = item.querySelector('div:first-child');
        const textSpan = item.querySelector('span');
        
        if (availableWidth < 120) {
          // Compact layout for very narrow legends
          item.style.padding = '2px';
          item.style.fontSize = '9px';
          if (colorBox) {
            colorBox.style.width = '10px';
            colorBox.style.height = '10px';
            colorBox.style.marginRight = '4px';
          }
          if (textSpan) {
            textSpan.style.fontSize = '9px';
            textSpan.style.lineHeight = '1.1';
          }
        } else if (availableWidth < 140) {
          // Medium compact layout
          item.style.padding = '2px';
          item.style.fontSize = '9px';
          if (colorBox) {
            colorBox.style.width = '11px';
            colorBox.style.height = '11px';
            colorBox.style.marginRight = '5px';
          }
          if (textSpan) {
            textSpan.style.fontSize = '9px';
            textSpan.style.lineHeight = '1.2';
          }
        } else {
          // Normal layout
          item.style.padding = '3px';
          item.style.fontSize = '10px';
          if (colorBox) {
            colorBox.style.width = '12px';
            colorBox.style.height = '12px';
            colorBox.style.marginRight = '6px';
          }
          if (textSpan) {
            textSpan.style.fontSize = '10px';
            textSpan.style.lineHeight = 'normal';
          }
        }
      });
      
      console.log('Final items container style:', {
        maxHeight: itemsContainer.style.maxHeight,
        height: itemsContainer.style.height,
        overflowY: itemsContainer.style.overflowY
      });
    } else {
      console.warn('Could not find items container in legend');
    }
  }
}

// Enhanced draggable and resizable functionality
function makeLegendDraggableAndResizable(legendContainer, legendId) {
  // Get the legend type from the legendId
  const type = legendId.includes('_gene') ? 'gene' : 'categorical';
  
  // Get or set initial position and size
  const savedPosition = legendPositions.get(legendId) || { x: null, y: null };
  const savedSize = legendSizes.get(legendId) || { width: 150, height: null };

  // Apply saved size
  legendContainer.style.width = savedSize.width + 'px';
  if (savedSize.height) {
    legendContainer.style.height = savedSize.height + 'px';
  }

  // Apply saved position or default
  if (savedPosition.x !== null && savedPosition.y !== null) {
    legendContainer.style.right = 'auto';
    legendContainer.style.left = savedPosition.x + 'px';
    legendContainer.style.top = savedPosition.y + 'px';
  }

  // Create drag handle
  const dragHandle = document.createElement('div');
  dragHandle.className = 'legend-drag-handle';
  dragHandle.style.cssText = `
    width: 100%;
    height: 20px;
    background: rgba(0, 0, 0, 0.1);
    cursor: move;
    border-radius: 4px 4px 0 0;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 10px;
    color: #666;
    user-select: none;
    border-bottom: 1px solid #ddd;
    flex-shrink: 0;
  `;
  dragHandle.innerHTML = '⋮⋮⋮';
  dragHandle.title = 'Drag to move';

  // Insert drag handle at the top
  legendContainer.insertBefore(dragHandle, legendContainer.firstChild);

  // Create resize handle
  const resizeHandle = document.createElement('div');
  resizeHandle.className = 'legend-resize-handle';
  resizeHandle.style.cssText = `
    position: absolute;
    bottom: 0;
    right: 0;
    width: 12px;
    height: 12px;
    background: rgba(0, 0, 0, 0.2);
    cursor: nw-resize;
    border-radius: 0 0 4px 0;
    z-index: 10;
    flex-shrink: 0;
  `;
  resizeHandle.innerHTML = '◢';
  resizeHandle.style.fontSize = '8px';
  resizeHandle.style.color = '#999';
  resizeHandle.style.lineHeight = '12px';
  resizeHandle.style.textAlign = 'center';
  resizeHandle.title = 'Drag to resize';
  legendContainer.appendChild(resizeHandle);

  // Make legend container relative for absolute positioning of children
  legendContainer.style.position = 'absolute';
  legendContainer.style.display = 'flex';
  legendContainer.style.flexDirection = 'column';

  // Dragging functionality
  let isDragging = false;
  let dragStartX, dragStartY, legendStartX, legendStartY;

  dragHandle.addEventListener('mousedown', (e) => {
    isDragging = true;
    dragStartX = e.clientX;
    dragStartY = e.clientY;
    
    const rect = legendContainer.getBoundingClientRect();
    const containerRect = legendContainer.offsetParent.getBoundingClientRect();
    legendStartX = rect.left - containerRect.left;
    legendStartY = rect.top - containerRect.top;
    
    e.preventDefault();
    dragHandle.style.cursor = 'grabbing';
  });

  document.addEventListener('mousemove', (e) => {
    if (!isDragging) return;
    
    const deltaX = e.clientX - dragStartX;
    const deltaY = e.clientY - dragStartY;
    
    const newX = legendStartX + deltaX;
    const newY = legendStartY + deltaY;
    
    // Keep within bounds
    const containerRect = legendContainer.offsetParent.getBoundingClientRect();
    const legendRect = legendContainer.getBoundingClientRect();
    
    const maxX = containerRect.width - legendRect.width;
    const maxY = containerRect.height - legendRect.height;
    
    const clampedX = Math.max(0, Math.min(newX, maxX));
    const clampedY = Math.max(0, Math.min(newY, maxY));
    
    legendContainer.style.left = clampedX + 'px';
    legendContainer.style.top = clampedY + 'px';
    legendContainer.style.right = 'auto';
    
    // Save position
    legendPositions.set(legendId, { x: clampedX, y: clampedY });
  });

  document.addEventListener('mouseup', () => {
    if (isDragging) {
      isDragging = false;
      dragHandle.style.cursor = 'move';
    }
  });

  // Resizing functionality
  let isResizing = false;
  let resizeStartX, resizeStartY, legendStartWidth, legendStartHeight;

  resizeHandle.addEventListener('mousedown', (e) => {
    isResizing = true;
    resizeStartX = e.clientX;
    resizeStartY = e.clientY;
    legendStartWidth = legendContainer.offsetWidth;
    legendStartHeight = legendContainer.offsetHeight;
    
    e.preventDefault();
    e.stopPropagation();
    resizeHandle.style.cursor = 'nw-resize';
  });

  document.addEventListener('mousemove', (e) => {
    if (!isResizing) return;
    
    const deltaX = e.clientX - resizeStartX;
    const deltaY = e.clientY - resizeStartY;
    
    const newWidth = Math.max(120, legendStartWidth + deltaX); // Min width 120px
    const newHeight = Math.max(100, legendStartHeight + deltaY); // Min height 100px
    
    legendContainer.style.width = newWidth + 'px';
    legendContainer.style.height = newHeight + 'px';
    
    // Update content layout based on legend type
    updateLegendContentLayout(legendContainer, newWidth, newHeight, type);
    
    // Save size
    legendSizes.set(legendId, { width: newWidth, height: newHeight });
  });

  document.addEventListener('mouseup', () => {
    if (isResizing) {
      isResizing = false;
      resizeHandle.style.cursor = 'nw-resize';
    }
  });

  // Add double-click to reset position and size
  dragHandle.addEventListener('dblclick', () => {
    legendContainer.style.width = '150px';
    legendContainer.style.height = 'auto';
    legendContainer.style.right = '8px';
    legendContainer.style.left = 'auto';
    legendContainer.style.top = '8px';
    
    // Reset content layout to original
    updateLegendContentLayout(legendContainer, 150, null, type);
    
    // Clear saved position and size
    legendPositions.delete(legendId);
    legendSizes.delete(legendId);
  });
}

async function createGenePlot(geneId) {
  console.log(`Creating gene plot for ${geneId}`);
  if (genePlots[geneId]) {
    console.log(`Gene plot for ${geneId} already exists - skipping creation`);
    return; // Already exists
  }

  // Create canvas
  const canvas = document.createElement('canvas');
  canvas.id = `gene_canvas_${geneId}`;
  canvas.style.cssText = `
    width: auto;
    height: auto;
    border-radius: 4px;
    border: 1px solid #ddd;
    box-sizing: border-box;
  `;

  // Add canvas to plot grid
  const plotGrid = document.getElementById('plotGrid');
  plotGrid.appendChild(canvas);
  
  // Add title overlay - use the actual gene name
  createPlotTitle(canvas, geneId, geneId);
  
  // Import regl-scatterplot if not already imported
  const module = await import('https://esm.sh/regl-scatterplot@1.14.1');
  const createScatterplot = module.default;
  const { createRenderer } = module;

  // Use shared renderer or create new one
  if (!renderer) {
    renderer = createRenderer();
  }

  // Get optimal settings based on point count
  const pointCount = points.length;
  let pointSize = 3;
  let opacity = 0.8;
  let performanceMode = false;

  // if (pointCount >= 1000000) {
  //   pointSize = 2;
  //   opacity = 0.6;
  //   performanceMode = true;
  // } else if (pointCount >= 500000) {
  //   pointSize = 1;
  //   opacity = 0.4;
  // }

  if (pointCount < 300) {
    pointSize = 100;
  } else if (pointCount < 50000) {
    pointSize = 30;
  } else {
    pointSize = 3;
  }

  const plot = createScatterplot({
    renderer,
    canvas,
    width: canvas.clientWidth || 400,
    height: canvas.clientHeight || 400,
    pointSize,
    opacity,
    performanceMode,
    lassoOnLongPress: true,
    lassoType: 'freeform',
  });

  // Subscribe to events
  plot.subscribe('select', ({ points: selectedIndices }) => {
    if (!syncing) {
      const eps = 1e-8;
      const pointToIndex = new Map();

      // Build the map once: O(M)
      points.forEach((p, index) => {
        // Round to ~10 decimal places to match eps tolerance and avoid float quirks
        const key = `${p[0].toFixed(10)},${p[1].toFixed(10)}`;
        if (!pointToIndex.has(key)) {
          pointToIndex.set(key, index);
        }
      });

      // Now map with fast lookups: O(N)
      const originalIndices = selectedIndices
        .map(i => {
          const t = filteredPoints[i];
          const key = `${t[0].toFixed(10)},${t[1].toFixed(10)}`;
          return pointToIndex.get(key) ?? -1;
        })
        .filter(i => i !== -1);

      Shiny.setInputValue('selectedPoints', originalIndices);
      
      // Sync selection to other plots if in sync mode
      if (isSyncMode) {
        syncSelectionToAllPlots(selectedIndices, geneId);
      }
    }
  });
  
  plot.subscribe('deselect', () => {
    if (!syncing) {
      Shiny.setInputValue('selectedPoints', []);
      
      // Sync deselection to other plots if in sync mode
      if (isSyncMode) {
        syncDeselectionToAllPlots(geneId);
      }
    }
  });

  // Store the plot BEFORE trying to inherit view
  genePlots[geneId] = plot;
  console.log(`Gene plot ${geneId} successfully created and stored`);
  
  // **Inherit the main plot's current view ONLY during creation**
  // Wait a frame to ensure plot is fully initialized
  await new Promise(resolve => requestAnimationFrame(resolve));
  
  if (scatterplot && !scatterplot._destroyed) {
    try {
      const mainCameraView = scatterplot.get('cameraView');
      
      if (mainCameraView) {
        console.log(`Inheriting main plot view for ${geneId} during creation:`, mainCameraView);
        
        // Apply the view to the newly created plot without triggering events
        plot.set({ cameraView: mainCameraView }, { preventEvent: true });
        
        // Also sync ALL existing gene plots to maintain consistency
        Object.keys(genePlots).forEach(existingGeneId => {
          if (existingGeneId !== geneId && genePlots[existingGeneId] && !genePlots[existingGeneId]._destroyed) {
            try {
              console.log(`Also syncing existing plot ${existingGeneId} to main view`);
              genePlots[existingGeneId].set({ cameraView: mainCameraView }, { preventEvent: true });
            } catch (error) {
              console.warn(`Failed to sync existing plot ${existingGeneId}:`, error);
            }
          }
        });

        let keepSync = false;

        if (isSyncMode) {
          keepSync = true;
        }

        removeSynchronization();

        if (keepSync) {
          setTimeout(() => {
            toggleSyncMode();  // will set text to Disable + red, and call setupSynchronization()
          }, 0);
        }
      }
    } catch (error) {
      console.warn(`Failed to inherit main plot view for ${geneId}:`, error);
    }
  }

  // Update layout after creating plot
  updatePlotLayout();
}



function removeGenePlot(geneId) {
  if (!genePlots[geneId]) return;

  console.log(`Removing gene plot for: ${geneId}`);

  const plot = genePlots[geneId];
  
  // Remove sync handler for this specific plot
  if (syncHandlers.has(geneId)) {
    try {
      plot.unsubscribe('view', syncHandlers.get(geneId));
      syncHandlers.delete(geneId);
    } catch (error) {
      console.warn(`Error removing sync handler for ${geneId}:`, error);
    }
  }
  
  // Also unsubscribe from other events
  try {
    plot.unsubscribe('select');
    plot.unsubscribe('deselect');
  } catch (error) {
    console.warn('Error unsubscribing from plot events:', error);
  }

  // Destroy WebGL plot
  plot.destroy();
  delete genePlots[geneId];

  // Remove canvas + title safely
  const canvas = document.getElementById(`gene_canvas_${geneId}`);
  if (canvas) {
    const container = canvas.parentElement;

    // Remove title
    const titleOverlay = container.querySelector('.plot-title');
    if (titleOverlay) titleOverlay.remove();

    // Remove legend (sibling in same container)
    const legend = container.querySelector('.plot-legend');
    if (legend) legend.remove();

    // Remove canvas
    canvas.remove();

    // If container empty, remove it
    if (container && container.classList.contains('canvas-container') && container.children.length === 0) {
      container.remove();
    }
  }

  // Re-setup synchronization with remaining plots if sync is active
  if (isSyncMode && activeGeneExpressions.size > 0) {
    setupSynchronization();
  }

  // Ensure layout is updated AFTER DOM removal
  updatePlotLayout();

  // Extra cleanup: if no active genes left
  if (activeGeneExpressions.size === 0) {
    console.log("No active genes remaining. Clearing legends & sync.");
    removeSyncButton();
    isSyncMode = false;
    clearAllSyncHandlers();
  }
}

// Store sync handlers to manage them properly
let syncHandlers = new Map();

function clearAllSyncHandlers() {
  console.log("Clearing all sync handlers...");
  
  // Clear main plot handler
  if (scatterplot && syncHandlers.has('main')) {
    try {
      scatterplot.unsubscribe('view', syncHandlers.get('main'));
      syncHandlers.delete('main');
    } catch (error) {
      console.warn('Error clearing main plot sync handler:', error);
    }
  }
  
  // Clear gene plot handlers
  Object.keys(genePlots).forEach(geneId => {
    const plot = genePlots[geneId];
    if (plot && syncHandlers.has(geneId)) {
      try {
        plot.unsubscribe('view', syncHandlers.get(geneId));
        syncHandlers.delete(geneId);
      } catch (error) {
        console.warn(`Error clearing gene plot sync handler for ${geneId}:`, error);
      }
    }
  });
}

function setupSynchronization() {
  console.log("Setting up synchronization...");
  
  if (!isSyncMode) return;
  
  // Clear any existing handlers first
  clearAllSyncHandlers();

  // Get only valid, non-destroyed plots
  const allPlots = [scatterplot, ...Object.values(genePlots)].filter(plot => {
    return plot && !plot._destroyed && typeof plot.subscribe === 'function';
  });

  console.log(`Syncing ${allPlots.length} plots`);

  const createSyncHandler = (sourceId) => {
    return () => {
      if (!isSyncMode || syncing) return;
      syncing = true;
      
      try {
        const sourcePlot = sourceId === 'main' ? scatterplot : genePlots[sourceId];
        if (!sourcePlot || sourcePlot._destroyed) {
          syncing = false;
          return;
        }
        
        const cameraView = sourcePlot.get('cameraView');

        // Get all valid target plots
        const validTargets = allPlots.filter(target => {
          return target && target !== sourcePlot && !target._destroyed && typeof target.set === 'function';
        });

        // Direct, synchronous updates for best performance
        validTargets.forEach(target => {
          try {
            // console.log(`Syncing camera view to target: ${target.id}`);
            target.set({ cameraView }, { preventEvent: true });
          } catch (error) {
            console.warn('Failed to sync camera to target:', error);
          }
        });
      } catch (error) {
        console.warn('Error in sync handler:', error);
      }

      // Reset sync flag immediately
      syncing = false;
    };
  };

  // Set up sync handlers for each plot
  if (scatterplot && !scatterplot._destroyed) {
    const handler = createSyncHandler('main');
    syncHandlers.set('main', handler);
    scatterplot.subscribe('view', handler);
  }

  Object.keys(genePlots).forEach(geneId => {
    const plot = genePlots[geneId];
    if (plot && !plot._destroyed) {
      const handler = createSyncHandler(geneId);
      syncHandlers.set(geneId, handler);
      plot.subscribe('view', handler);
    }
  });
}

function removeSynchronization() {
  console.log("Removing synchronization...");

  isSyncMode = false;
  syncing = false;

  // Use the centralized handler clearing function
  clearAllSyncHandlers();

  // Reset camera views independently for all valid plots
  const allPlots = [scatterplot, ...Object.values(genePlots)].filter(plot => {
    return plot && !plot._destroyed && typeof plot.get === 'function' && typeof plot.set === 'function';
  });
  
  allPlots.forEach(plot => {
    try {
      const view = plot.get('cameraView');
      // Re-apply view but break any shared internal references
      plot.set({ cameraView: { ...view } }, { preventEvent: true });
    } catch (error) {
      console.warn('Failed to reset camera view for plot:', error);
    }
  });
}

function syncSelectionToAllPlots(selectedIndices, sourceId) {
  if (!isSyncMode || syncing) return;
  
  syncing = true;
  
  // Sync to main plot
  if (sourceId !== 'main' && scatterplot && !scatterplot._destroyed) {
    try {
      scatterplot.select(selectedIndices, { preventEvent: true });
    } catch (error) {
      console.warn('Failed to sync selection to main plot:', error);
    }
  }
  
  // Sync to gene plots
  Object.keys(genePlots).forEach(geneId => {
    if (geneId !== sourceId && genePlots[geneId] && !genePlots[geneId]._destroyed) {
      try {
        genePlots[geneId].select(selectedIndices, { preventEvent: true });
      } catch (error) {
        console.warn(`Failed to sync selection to gene plot ${geneId}:`, error);
      }
    }
  });
  
  Promise.resolve().then(() => { syncing = false; });
}

function syncDeselectionToAllPlots(sourceId) {
  if (!isSyncMode || syncing) return;
  
  syncing = true;
  
  // Sync to main plot
  if (sourceId !== 'main' && scatterplot && !scatterplot._destroyed) {
    try {
      scatterplot.deselect({ preventEvent: true });
    } catch (error) {
      console.warn('Failed to sync deselection to main plot:', error);
    }
  }
  
  // Sync to gene plots
  Object.keys(genePlots).forEach(geneId => {
    if (geneId !== sourceId && genePlots[geneId] && !genePlots[geneId]._destroyed) {
      try {
        genePlots[geneId].deselect({ preventEvent: true });
      } catch (error) {
        console.warn(`Failed to sync deselection to gene plot ${geneId}:`, error);
      }
    }
  });
  
  Promise.resolve().then(() => { syncing = false; });
}

async function initScatterplot() {
  // Wait until DOM is fully loaded
  if (document.readyState !== 'complete') {
    await new Promise(resolve => window.addEventListener('load', resolve));
  }

  const canvas = document.getElementById('scatterplot_canvas');
  if (!canvas) {
    console.error("Canvas element 'scatterplot_canvas' not found!");
    return;
  }

  // Style the main canvas
  canvas.style.cssText = `
    width: 100%;
    height: 100%;
    border-radius: 4px;
    border: 1px solid #ddd;
    box-sizing: border-box;
  `;

  // Create the plot grid and add title to main plot
  createPlotGrid();
  createPlotTitle(canvas, 'Main View', 'main');

  // Import regl-scatterplot
  const module = await import('https://esm.sh/regl-scatterplot@1.14.1');
  const createScatterplot = module.default;
  const { createRenderer } = module;

  // Create shared renderer
  renderer = createRenderer();

  scatterplot = createScatterplot({
    renderer,
    canvas,
    width: canvas.clientWidth || 800,
    height: canvas.clientHeight || 400,
    // pointSize: 2.5,
    pointScaleMode: 'asinh',
    opacity: 0.8,
    lassoOnLongPress: true,
    lassoType: 'freeform',
  });

  subscribeMainPlotEvents();
  updatePlotLayout();

  initScatterplotPromise = Promise.resolve();
}

function subscribeMainPlotEvents() {
  // Modified plot.subscribe for selection
  scatterplot.subscribe('select', ({ points: selectedIndices }) => {
    if (!syncing) {
      const eps = 1e-8;
      const pointToIndex = new Map();

      // Build the map once: O(M)
      points.forEach((p, index) => {
        // Round to ~10 decimal places to match eps tolerance and avoid float quirks
        const key = `${p[0].toFixed(10)},${p[1].toFixed(10)}`;
        if (!pointToIndex.has(key)) {
          pointToIndex.set(key, index);
        }
      });

      // Now map with fast lookups: O(N)
      const originalIndices = selectedIndices
        .map(i => {
          const t = filteredPoints[i];
          const key = `${t[0].toFixed(10)},${t[1].toFixed(10)}`;
          return pointToIndex.get(key) ?? -1;
        })
        .filter(i => i !== -1);

      console.log('Original indices of selected points:', originalIndices);
      console.log('Selected indices in filteredPoints:', selectedIndices);

      
      // Get the visible categories for the current plot
      let plotId = 'main';
      const legendKey = `${plotId}_categorical`;
      const legendData = currentLegendData.get(legendKey);
      const visibleIndices = legendData ? Array.from(legendData.visible) : [];

      // Send both selected points and visible categories to Shiny
      Shiny.setInputValue('selectedPoints', {
        plotId: plotId,
        selectedIndices: originalIndices,
        originalIndices: originalIndices,
        visibleCategories: visibleIndices,
        annotationName: legendData ? legendData.annotationName : null,
        timestamp: Date.now()
      }, { priority: 'event' });
      
      console.log('Selected points:', selectedIndices);

      // Sync selection to other plots if in sync mode
      if (isSyncMode) {
        syncSelectionToAllPlots(selectedIndices, geneId);
      }
    }
  });

  scatterplot.subscribe('deselect', () => {
    if (!syncing) {
      // Get the visible categories for the current plot
      let plotId = 'main';
      const legendKey = `${plotId}_categorical`;
      const legendData = currentLegendData.get(legendKey);
      const visibleIndices = legendData ? Array.from(legendData.visible) : [];

      // Send deselection info with visible categories
      Shiny.setInputValue('selectedPoints', {
        plotId: plotId,
        selectedIndices: [],
        originalIndices: [],
        visibleCategories: visibleIndices,
        annotationName: legendData ? legendData.annotationName : null,
        timestamp: Date.now()
      }, { priority: 'event' });

      // Sync deselection to other plots if in sync mode
      if (isSyncMode) {
        syncDeselectionToAllPlots(geneId);
      }
    }
  });
}

function filterPoints() {

  if (mainAnnotation == null) return points;

  if (mainAnnotation.var_type != 'categorical') {

    let minVal = Infinity;
    let maxVal = -Infinity;

    const pointsColored = points.map((p, i) => {
      const newP = [...p];
      const rawVal = mainAnnotation.annotation[i] - 1;

      // track min and max while building new array
      if (rawVal < minVal) minVal = rawVal;
      if (rawVal > maxVal) maxVal = rawVal;

      newP[2] = rawVal; // temporarily store raw
      return newP;
    });

    const range = maxVal - minVal || 1;

    // second pass to normalize
    pointsColored.forEach(p => {
      p[2] = (p[2] - minVal) / range;
    });

    return pointsColored;

  } else {

    const pointsColored = points.map((p, i) => {
      const newP = [...p];
      newP[2] = mainAnnotation.annotation[i] - 1; // match index
      return newP;
    });

    console.log(mainAnnotation);
  
  
    return pointsColored.filter(p => {
      const valueIndex = 2 + gene_number;
      const value = p[valueIndex];
      return mainAnnotation.visible.has(value);
    });

  }
  
}



async function redrawAllPlots(mainAnnotation=null) {
  if (points.length === 0) return;
  
  // Always update filteredPoints before redrawing
  filteredPoints = filterPoints();
  
  // Start all plot redraws simultaneously
  const allDrawPromises = [];
  
  // Add main plot redraw
  allDrawPromises.push(redrawMainPlot(mainAnnotation));
  
  // Add gene plot redraws
  Object.keys(genePlots).forEach(geneId => {
    allDrawPromises.push(Promise.resolve(redrawGenePlot(geneId)));
  });
  
  // Wait for all plots to finish drawing simultaneously
  await Promise.all(allDrawPromises);
}

// Additional debugging function to verify coordinate consistency
function verifyCoordinateConsistency() {
  console.log("=== Coordinate Consistency Check ===");
  console.log(`Points array length: ${points.length}`);
  console.log(`FilteredPoints length: ${filteredPoints.length}`);
  console.log(`PointsXY length: ${pointsXY.length}`);
  console.log(`PointsMAGIC length: ${pointsMAGIC.length}`);
  console.log(`MAGIC active: ${isMAGICActive}`);
  
  if (points.length > 0) {
    const firstPoint = points[0];
    const firstXY = pointsXY[0];
    const firstMAGIC = pointsMAGIC[0];
    
    console.log(`First point coordinates: [${firstPoint[0]}, ${firstPoint[1]}]`);
    console.log(`First XY coordinates: [${firstXY[0]}, ${firstXY[1]}]`);
    console.log(`First MAGIC coordinates: [${firstMAGIC[0]}, ${firstMAGIC[1]}]`);
    
    if (isMAGICActive) {
      const coordsMatch = Math.abs(firstPoint[0] - firstMAGIC[0]) < 0.001 && 
                         Math.abs(firstPoint[1] - firstMAGIC[1]) < 0.001;
      console.log(`Coordinates match MAGIC: ${coordsMatch}`);
    } else {
      const coordsMatch = Math.abs(firstPoint[0] - firstXY[0]) < 0.001 && 
                         Math.abs(firstPoint[1] - firstXY[1]) < 0.001;
      console.log(`Coordinates match XY: ${coordsMatch}`);
    }
  }
  console.log("================================");
}


async function redrawMainPlot(annotation=null) {
  if (!scatterplot || filteredPoints.length === 0) return;

  if (annotation !== null) {
    // ... your existing color and pointSize setup code ...
    
    console.log(annotation);
    
    if (annotation.var_type != 'categorical') {
      colors = viridisColors
    } else {
      colors = annotation.colors;
    }

    if (filteredPoints.length < 300) {
      pointSize = 100;
    } else if (filteredPoints.length < 100000) {
      pointSize = 30;
    } else {
      pointSize = 3;
    }
    
    console.log(points[10]);
    console.log(filteredPoints[10]);
    
    const config = { 
      colorBy: 'category',
      pointColor: colors,
      sizeBy: 'category',
      pointSize: Array(filteredPoints.length).fill(pointSize),
    };
    
    scatterplot.set(config);
    await scatterplot.draw(filteredPoints, { transition: useTransition });
    
    // Your existing bounds/resizing code...
    if (initScatterplotPromise) {
      console.log("Resizing plot");
      const xs = filteredPoints.map(p => p[0]);
      const ys = filteredPoints.map(p => p[1]);
      const minX = Math.min(...xs);
      const maxX = Math.max(...xs);
      const minY = Math.min(...ys);
      const maxY = Math.max(...ys);
      const padX = (maxX - minX) * 0.1 || 1;
      const padY = (maxY - minY) * 0.1 || 1;
      
      scatterplot.zoomToArea({
        x: minX - padX, y: minY - padY,
        width: maxX + padX - (minX - padX),
        height: maxY + padY - (minY - padY)
      });
      
      initScatterplotPromise = null;
    }

    // Only create/update legend if not in visibility update mode
    if (annotation.var_type == 'categorical') {
      createPlotLegend('main', {
        names: annotation.names,
        colors: annotation.colors,
        visible: annotation.visible,
        annotationName: annotation.colorBy
      }, 'categorical');
    } else {
      const values = Array.from(annotation.annotation);
      const minVal = Math.min(...values);
      const maxVal = Math.max(...values);
      const meanVal = values.reduce((a, b) => a + b, 0) / values.length;

      createPlotLegend('main', {
        minVal: minVal.toFixed(2),
        maxVal: maxVal.toFixed(2),
        midVal: meanVal.toFixed(2),
        realMin: minVal.toFixed(2),
        realMax: maxVal.toFixed(2),
        geneName: annotation.colorBy
      }, 'gene');
    }

  } else {
    // Your existing else block code...
    if (points.length < 300) {
      pointSize = 100;
    } else {
      pointSize = 30;
    }
    
    const config = { 
      pointColor: '#000000',
    };
    
    scatterplot.set(config);
    await scatterplot.draw(points, { transition: useTransition });
    
    if (initScatterplotPromise) {
      const xs = points.map(p => p[0]);
      const ys = points.map(p => p[1]);
      const minX = xs.reduce((a, b) => Math.min(a, b), Infinity);
      const maxX = xs.reduce((a, b) => Math.max(a, b), -Infinity);
      const minY = ys.reduce((a, b) => Math.min(a, b), Infinity);
      const maxY = ys.reduce((a, b) => Math.max(a, b), -Infinity);
      const padX = (maxX - minX) * 0.1 || 1;
      const padY = (maxY - minY) * 0.1 || 1;
      
      scatterplot.zoomToArea({
        x: minX - padX, y: minY - padY,
        width: maxX + padX - (minX - padX),
        height: maxY + padY - (minY - padY)
      });

      initScatterplotPromise = null;
    }
  }
}

// Add a function to calculate percentiles
function calculatePercentile(values, percentile) {
  if (!values || values.length === 0) return 0;
  
  // Filter out any invalid values and sort
  const validValues = values.filter(v => !isNaN(v) && isFinite(v)).sort((a, b) => a - b);
  if (validValues.length === 0) return 0;
  
  const index = (percentile / 100) * (validValues.length - 1);
  
  if (index === Math.floor(index)) {
    return validValues[index];
  } else {
    const lower = validValues[Math.floor(index)];
    const upper = validValues[Math.ceil(index)];
    return lower + (upper - lower) * (index - Math.floor(index));
  }
}

// Update the processGeneExpressionData function to calculate percentiles
function processGeneExpressionData(expressionData) {
  console.log("Raw expression data received:", expressionData);
  
  // Check if we have valid data
  if (!expressionData || !expressionData.data || !expressionData.genes || 
      !Array.isArray(expressionData.genes) || expressionData.genes.length === 0) {
    console.log("No valid gene expression data provided");
    geneDataMatrix = null;
    currentGeneNames = [];
    geneDataInfo = null;
    return;
  }

  try {
    console.log(`Processing gene expression data for ${expressionData.genes.length} genes`);
    console.log("Gene names:", expressionData.genes);
    
    // Check if we have actual expression data
    if (!expressionData.data || expressionData.data === "" || 
        !expressionData.nrows || !expressionData.ncols) {
      console.log("No expression matrix data available");
      geneDataMatrix = null;
      currentGeneNames = [];
      geneDataInfo = null;
      return;
    }
    
    // Decode binary data for original expression
    const binaryStr = atob(expressionData.data);
    const len = binaryStr.length;
    const buffer = new ArrayBuffer(len);
    const view = new Uint8Array(buffer);
    for (let i = 0; i < len; i++) {
      view[i] = binaryStr.charCodeAt(i);
    }
    const floatArray = new Float64Array(buffer);
    
    // Process MAGIC data if available
    let magicFloatArray = null;
    if (expressionData.magic_data && expressionData.magic_data !== "") {
      const magicBinaryStr = atob(expressionData.magic_data);
      const magicLen = magicBinaryStr.length;
      const magicBuffer = new ArrayBuffer(magicLen);
      const magicView = new Uint8Array(magicBuffer);
      for (let i = 0; i < magicLen; i++) {
        magicView[i] = magicBinaryStr.charCodeAt(i);
      }
      magicFloatArray = new Float64Array(magicBuffer);
      console.log("MAGIC expression data processed");
    }
    
    // Store gene data info
    geneDataInfo = {
      nrows: expressionData.nrows,
      ncols: expressionData.ncols,
      genes: expressionData.genes,
      ranges: expressionData.ranges || {},
      magic_ranges: expressionData.magic_ranges || {}
    };
    
    console.log("Gene data info:", geneDataInfo);
    
    // Reshape into matrix (genes as columns, cells as rows) - Original data
    geneDataMatrix = [];
    for (let row = 0; row < geneDataInfo.nrows; row++) {
      geneDataMatrix[row] = [];
      for (let col = 0; col < geneDataInfo.ncols; col++) {
        const index = col * geneDataInfo.nrows + row; // Column-major order
        geneDataMatrix[row][col] = floatArray[index];
      }
    }

    // Reshape MAGIC matrix if available
    let magicDataMatrix = null;
    if (magicFloatArray) {
      magicDataMatrix = [];
      for (let row = 0; row < geneDataInfo.nrows; row++) {
        magicDataMatrix[row] = [];
        for (let col = 0; col < geneDataInfo.ncols; col++) {
          const index = col * geneDataInfo.nrows + row;
          magicDataMatrix[row][col] = magicFloatArray[index];
        }
      }
    }
    
    currentGeneNames = expressionData.genes;
    
    // Cache both original and MAGIC data for each gene
    geneExpressionOriginal.clear();
    geneExpressionMAGIC.clear();
    geneExpressionCache.clear();
    
    expressionData.genes.forEach((geneName, geneIndex) => {
      const geneValues = new Array(geneDataInfo.nrows);
      let min = Infinity, max = -Infinity, sum = 0;

      for (let row = 0; row < geneDataInfo.nrows; row++) {
        const val = geneDataMatrix[row][geneIndex];
        geneValues[row] = val;
        if (val < min) min = val;
        if (val > max) max = val;
        sum += val;
      }

      const mean = sum / geneValues.length;
      const p90 = calculatePercentile(geneValues, 90);
      const p95 = calculatePercentile(geneValues, 95);
      const p99 = calculatePercentile(geneValues, 99);

      const originalGeneData = {
        values: geneValues,
        range: expressionData.ranges[geneName] || { min, max, mean },
        percentiles: { p90, p95, p99, min, max, mean }
      };
      geneExpressionOriginal.set(geneName, originalGeneData);

      // MAGIC data if available
      if (magicDataMatrix) {
        const magicGeneValues = new Array(geneDataInfo.nrows);
        let mMin = Infinity, mMax = -Infinity, mSum = 0;

        for (let row = 0; row < geneDataInfo.nrows; row++) {
          const val = magicDataMatrix[row][geneIndex];
          magicGeneValues[row] = val;
          if (val < mMin) mMin = val;
          if (val > mMax) mMax = val;
          mSum += val;
        }

        const mMean = mSum / magicGeneValues.length;
        const magicP90 = calculatePercentile(magicGeneValues, 90);
        const magicP95 = calculatePercentile(magicGeneValues, 95);
        const magicP99 = calculatePercentile(magicGeneValues, 99);

        const magicGeneData = {
          values: magicGeneValues,
          range: expressionData.magic_ranges[geneName] || { min: mMin, max: mMax, mean: mMean },
          percentiles: { p90: magicP90, p95: magicP95, p99: magicP99, min: mMin, max: mMax, mean: mMean }
        };
        geneExpressionMAGIC.set(geneName, magicGeneData);
      }
    });

    
    // Set the active cache based on current MAGIC state
    updateActiveGeneExpressionCache();
    
    console.log(`Successfully processed gene expression data: ${geneDataInfo.nrows} cells x ${geneDataInfo.ncols} genes`);
    console.log("Available genes:", currentGeneNames);
    console.log("MAGIC data available:", magicDataMatrix !== null);
    
  } catch (error) {
    console.error('Error processing gene expression data:', error);
    console.error("Error details:", error.stack);
    geneDataMatrix = null;
    currentGeneNames = [];
    geneDataInfo = null;
  }
}

// function updateActiveGeneExpressionCache() {
//   geneExpressionCache.clear();
  
//   // Original or MAGIC depending on isMAGICActive
//   const sourceCache = isMAGICActive ? geneExpressionMAGIC : geneExpressionOriginal;

//   sourceCache.forEach((geneData, geneName) => {
//     geneExpressionCache.set(geneName, geneData);
//   });

//   console.log(`Active cache updated with ${isMAGICActive ? "MAGIC" : "original"} expression data`);
// }

// Add a global setting for vmax mode
let vmaxMode = 'p90'; // Options: 'max', 'p90', 'p95', 'p99'

// Function to set vmax mode
function setVmaxMode(mode) {
  vmaxMode = mode;
  console.log(`Setting vmax mode to: ${mode}`);
  // Redraw all gene plots with new scaling
  Object.keys(genePlots).forEach(geneId => {
    redrawGenePlot(geneId);
  });
}

function getGeneExpressionValue(geneId, cellIndex) {
  const geneData = geneExpressionCache.get(geneId);
  if (!geneData || cellIndex >= geneData.values.length) {
    return 0;
  }
  return geneData.values[cellIndex];
}

// Add a function to get available genes (useful for debugging)
function getAvailableGenes() {
  return Array.from(geneExpressionCache.keys());
}

// Add performance monitoring
function logPerformanceStats() {
  if (geneDataInfo) {
    console.log("Gene Expression Data Stats:");
    console.log(`- Cells: ${geneDataInfo.nrows}`);
    console.log(`- Genes: ${geneDataInfo.ncols}`);
    console.log(`- Memory usage: ~${(geneDataInfo.nrows * geneDataInfo.ncols * 4 / 1024 / 1024).toFixed(2)} MB`);
    console.log(`- Cached genes: ${geneExpressionCache.size}`);
    console.log(`- Active plots: ${activeGeneExpressions.size}`);
  }
}


function getColorLimits(geneData, vmaxMode) {
  // pick vmax from percentiles
  let vmax;
  switch (vmaxMode) {
    case 'p90': vmax = geneData.percentiles.p90; break;
    case 'p95': vmax = geneData.percentiles.p95; break;
    case 'p99': vmax = geneData.percentiles.p99; break;
    case 'max': 
    default:    vmax = geneData.percentiles.max; break;
  }

  // vmin = min(0, actual min)
  let vmin = Math.min(0, geneData.percentiles.min);

  // If vmax and vmin are almost identical (like all zeros)
  if (Math.abs(vmax - vmin) < 1e-9) {
    vmin = -0.1;
    vmax =  0.1;
  }

  return { vmin, vmax };
}


// Updated redrawGenePlot function with percentile scaling
function redrawGenePlot(geneId) {
  const plot = genePlots[geneId];
  if (!plot) return;

  const geneData = geneExpressionCache.get(geneId);
  if (!geneData) {
    console.warn(`No expression data found for gene: ${geneId}`);
    return;
  }

  // 1. Determine vmin/vmax using helper
  const { vmin, vmax } = getColorLimits(geneData, vmaxMode);
  const range = vmax - vmin;

  // 2. Create gene points - handle categorical filtering properly
  const genePoints = [];

  if (mainAnnotation && mainAnnotation.var_type === 'categorical') {
    // For categorical filtering, build points from visible categories only
    for (let originalIndex = 0; originalIndex < points.length; originalIndex++) {
      const categoryIndex = mainAnnotation.annotation[originalIndex] - 1;
      if (mainAnnotation.visible.has(categoryIndex)) {
        const rawVal = geneData.values[originalIndex];
        const normalized = Math.min(Math.max((rawVal - vmin) / range, 0.0), 1.001);
        genePoints.push([
          points[originalIndex][0], 
          points[originalIndex][1], 
          normalized
        ]);
      }
    }
  } else {
    // For non-categorical or no filtering, use all points
    for (let i = 0; i < points.length; i++) {
      const rawVal = geneData.values[i];
      const normalized = Math.min(Math.max((rawVal - vmin) / range, 0.0), 1.001);
      genePoints.push([points[i][0], points[i][1], normalized]);
    }
  }

  if (genePoints.length === 0) {
    console.warn(`No valid points to plot for gene: ${geneId}`);
    return;
  }

  // 3. Set plot config
  const config = {
    colorBy: 'valueA',
    pointColor: viridisColors,
  };

  if (!useTransition) plot.clear();
  plot.set(config);
  plot.draw(genePoints, { transition: useTransition });

  // 4. Update legend
  createPlotLegend(geneId, {
    minVal: vmin.toFixed(2),
    maxVal: vmax.toFixed(2),
    midVal: ((vmin + vmax) / 2).toFixed(2),
    realMin: geneData.percentiles.min.toFixed(2),
    realMax: geneData.percentiles.max.toFixed(2),
    geneName: geneId
  }, 'gene');
}

// Add controls for vmax mode (you can integrate this into your UI)
function createVmaxControls() {
  // Check if controls already exist
  if (document.getElementById('vmaxControls')) return;

  const controlsContainer = document.createElement('div');
  controlsContainer.id = 'vmaxControls';
  controlsContainer.style.cssText = `
    position: absolute;
    top: 50px;
    right: 10px;
    z-index: 1000;
    background: rgba(255, 255, 255, 0.9);
    padding: 8px;
    border-radius: 4px;
    border: 1px solid #ddd;
    font-size: 12px;
  `;

  controlsContainer.innerHTML = `
    <div style="margin-bottom: 4px; font-weight: bold;">Color Scale Max:</div>
    <select id="vmaxSelect" style="width: 100%; font-size: 12px;">
      <option value="p90" selected>90th Percentile</option>
      <option value="p95">95th Percentile</option>
      <option value="p99">99th Percentile</option>
      <option value="max">Maximum</option>
    </select>
  `;

  const select = controlsContainer.querySelector('#vmaxSelect');
  select.onchange = (e) => {
    setVmaxMode(e.target.value);
  };

  // Add to the plot grid container
  const plotGrid = document.getElementById('plotGrid');
  if (plotGrid) {
    plotGrid.style.position = 'relative';
    plotGrid.appendChild(controlsContainer);
  }
}

function setSpinner(visible) {
  const overlay = document.getElementById('loadingOverlay');
  if (overlay) {
    overlay.style.display = visible ? 'flex' : 'none';
  }
}

// Shiny message handlers
Shiny.addCustomMessageHandler('showSpinner', function(show) {
  setSpinner(show);
});

// Replace your current updateData handler with this version
Shiny.addCustomMessageHandler('updateData', async function(message) {
  // Capture timing immediately
  const jsReceiveTime = Date.now();
  const jsProcessStart = performance.now();
  
  try {
    const { 
      base64, 
      annotationData, 
      numCols, 
      clusters: cl, 
      colors, 
      metacellColors, 
      geneExprRanges: ge,
      sendTime,  // R timestamp when message was sent
      timingId   // Unique ID for this timing measurement
    } = message;

    // Calculate transfer time if sendTime is provided
    let transferTime = null;
    if (sendTime) {
      transferTime = jsReceiveTime - sendTime;
      console.log('⏱ R→JS transfer time:', transferTime.toFixed(2), 'ms');
    }

    // === YOUR EXISTING CODE STARTS HERE ===
    geneExprRange = ge;
    if (geneExprRange && Object.keys(geneExprRange).length > 0) {
      console.log("Gene expression range:", geneExprRange);
      gene_number = Object.keys(geneExprRange).length;
    } else {
      console.warn("No gene expression range provided, defaulting to 0 genes.");
      gene_number = 0;
    }

    clusters = cl;
    
    if (!base64 || !numCols || numCols <= 0) {
      console.error('Invalid message data:', message);
      setSpinner(false);
      return;
    }
    
    // Parse binary data
    const binaryStr = atob(base64);
    const len = binaryStr.length;
    const buffer = new ArrayBuffer(len);
    const view = new Uint8Array(buffer);
    for (let i = 0; i < len; i++) {
      view[i] = binaryStr.charCodeAt(i);
    }
    const floatArray = new Float32Array(buffer);
    const pointCount = Math.floor(floatArray.length / numCols);
    
    if (pointCount <= 0 || pointCount > 3000000) {
      console.error('Invalid point count:', pointCount);
      setSpinner(false);
      return;
    }
    
    points = [];
    pointsXY = [];
    pointsMAGIC = [];

    for (let i = 0; i < pointCount; i++) {
      const baseIdx = i * numCols;
      const x = floatArray[baseIdx];
      const y = floatArray[baseIdx + 1];
      const magicX = floatArray[baseIdx + 2];
      const magicY = floatArray[baseIdx + 3];
      
      const rest = [];
      for (let j = 4; j < numCols; j++) {
        rest.push(floatArray[baseIdx + j]);
      }

      points.push([x, y, ...rest]);
      pointsXY.push([x, y]);
      pointsMAGIC.push([x, y]);
    }

    if ('seacell' in annotationData) {
      clusterColors = [...colors, ...metacellColors];
    }

    mainAnnotation = null;

    // Set up annotations
    annotations = {};
    if (annotationData) {
      Object.keys(annotationData).forEach(key => {
        const names = annotationData[key].names;
        const colors = annotationData[key].colors;

        if (Array.isArray(names)) {
          annotations[key] = {
            names,
            colors,
            visible: new Set(names.map((_, i) => i))
          };
        } else {
          console.warn(`annotationData[${key}].names is not an array:`, names);
          annotations[key] = { names: [], colors: [], visible: new Set() };
        }
      });
    }

    console.log("Annotations set:", annotations);

    // Update all plot settings based on point count
    const pointCount2 = points.length;
    let pointSize = 2.5;
    let opacity = 0.8;
    let performanceMode = false;

    if (pointCount2 >= 1000000) {
      pointSize = 2;
      opacity = 0.6;
      performanceMode = true;
    } else if (pointCount2 >= 500000) {
      pointSize = 1;
      opacity = 0.4;
    }

    if (scatterplot) {
      scatterplot.set({ pointSize, opacity, performanceMode });
    }

    Object.values(genePlots).forEach(plot => {
      plot.set({ pointSize, opacity, performanceMode });
    });

    initScatterplotPromise = Promise.resolve();

    // Time the redraw operation specifically
    const redrawStart = performance.now();
    await redrawAllPlots(mainAnnotation);
    const redrawEnd = performance.now();
    
    setSpinner(false);
    // === YOUR EXISTING CODE ENDS HERE ===
    
    // Calculate total processing time
    const jsProcessEnd = performance.now();
    const jsProcessingTime = jsProcessEnd - jsProcessStart;
    const redrawTime = redrawEnd - redrawStart;
    
    console.log('⏱ Total JS processing:', jsProcessingTime.toFixed(2), 'ms');
    console.log('⏱ redrawAllPlots took:', redrawTime.toFixed(2), 'ms');
    
    // Send timing back to R in the format it expects
    const completionTime = Date.now();
    
    Shiny.setInputValue("updateData_done", {
      timingId: timingId,
      jsReceiveTime: jsReceiveTime,
      jsCompletionTime: completionTime,
      transferTime: transferTime,
      jsProcessingTime: jsProcessingTime,
      redrawTime: redrawTime,
      totalEndToEndTime: sendTime ? (completionTime - sendTime) : null,
      dataStats: {
        pointCount: pointCount,
        numCols: numCols,
        base64Length: base64?.length || 0
      },
      timestamp: completionTime
    }, {priority: 'event'});
    
    console.log('📤 Sent timing data back to R');
    
  } catch (error) {
    console.error('❌ Error in updateData:', error);
    setSpinner(false);
    
    // Send error info back to R
    Shiny.setInputValue("updateData_done", {
      timingId: message.timingId || "unknown",
      error: error.message,
      timestamp: Date.now()
    }, {priority: 'event'});
  }
});

Shiny.addCustomMessageHandler('colorByChange', async function(message) {
  // Decode binary annotation data from base64
  if (message.annotation_raw) {
    try {
      // Decode base64 to binary string
      const binaryString = atob(message.annotation_raw);
      
      // Create ArrayBuffer from binary string
      const buffer = new ArrayBuffer(binaryString.length);
      const bytes = new Uint8Array(buffer);
      
      for (let i = 0; i < binaryString.length; i++) {
        bytes[i] = binaryString.charCodeAt(i);
      }
      
      // Interpret as Float32 array (matches R's writeBin with size=4)
      message.annotation = new Float32Array(buffer);
      
      // Validate length if provided
      if (message.annotation_length && message.annotation.length !== message.annotation_length) {
        console.warn('Binary data length mismatch:', message.annotation.length, 'vs expected:', message.annotation_length);
      }
      
      console.log('Binary transfer successful:', message.annotation.length, 'values');
      
    } catch (error) {
      console.error('Failed to decode binary annotation data:', error);
      // Fallback: could request non-binary version or show error
      return;
    }
  } else {
    console.warn('No binary annotation data received');
    return;
  }

  // Visible set = all indices by default
  message.visible = new Set(message.names.map((_, i) => i));

  // Remove sync button for non-gene expression modes
  if (!activeGeneExpressions.size) {
    removeSyncButton();
    if (isSyncMode) {
      toggleSyncMode();
    }
  }

  mainAnnotation = message;
  await redrawAllPlots(mainAnnotation);
});

// Shiny.addCustomMessageHandler('mainPlotAnnotationChange', function(message) {
//   // from base64 to array
//   const annotation = JSON.parse(atob(message.annotation));
//   currentAnnotation = annotation;
//   redrawMainPlot(currentAnnotation);
// });

// Global variable to track the active layer
let activeLayer = 'X'; // Default layer

// New handler for layer changes
Shiny.addCustomMessageHandler('updateLayer', function(message) {
  console.log(`Updating layer to: ${message.layer}`);
  
  activeLayer = message.layer;
  
  // Update coordinates based on the selected layer
  let coordSource;
  if (activeLayer === 'MAGIC') {
    coordSource = pointsMAGIC;
    console.log('Using MAGIC coordinates');
  } else {
    coordSource = pointsXY; // Default to X or other layer coordinates
    console.log(`Using coordinates for layer: ${activeLayer}`);
  }
  
  // Update points array with new coordinates
  for (let i = 0; i < points.length; i++) {
    points[i][0] = coordSource[i][0]; // Replace x
    points[i][1] = coordSource[i][1]; // Replace y
  }
  
  // Update filteredPoints
  filteredPoints = filterPoints();
  
  // Update gene expression cache based on layer
  updateActiveGeneExpressionCache();
  
  // Check if expression data is available for the layer
  if (activeLayer === 'MAGIC' && geneExpressionMAGIC.size === 0) {
    console.warn(`No expression data available for layer: ${activeLayer}`);
  }
  
  // Enable transition for smooth animation
  useTransition = true;
  
  // Redraw all plots with new coordinates and expression data
  redrawAllPlots(mainAnnotation);
  
  // Reset transition after animation
  setTimeout(() => {
    useTransition = false;
    console.log(`Layer transition to ${activeLayer} completed`);
  }, 1500);
});

// Modified geneSearchChange handler
Shiny.addCustomMessageHandler('geneSearchChange', function(message) {
  console.log("Gene search change message:", message);

  // Normalize gene names to always be an array
  const geneNames = message.genes 
    ? (Array.isArray(message.genes) ? message.genes : [message.genes])
    : [];

  // Normalize expression data and its genes
  const expressionData = message.expression_data || {};
  expressionData.genes = expressionData.genes 
    ? (Array.isArray(expressionData.genes) ? expressionData.genes : [expressionData.genes])
    : [];

  console.log("Processed gene names:", geneNames);
  console.log("Expression data:", expressionData);

  // Process both original and MAGIC gene expression data
  processGeneExpressionData(expressionData);

  console.log("Active genes before update:", Array.from(activeGeneExpressions));

  // Filter out empty strings and ensure valid names
  const newGenes = new Set(
    geneNames.filter(name => name && typeof name === 'string' && name.trim() !== '')
  );

  console.log("New genes to process:", Array.from(newGenes));

  // Remove deselected genes
  for (const gene of Array.from(activeGeneExpressions)) {
    if (!newGenes.has(gene)) {
      console.log(`Removing gene: ${gene}`);
      activeGeneExpressions.delete(gene);
      removeGenePlot(gene);
    }
  }

  // Add new genes (with duplicate protection)
  const genesToAdd = [];
  for (const gene of newGenes) {
    if (!activeGeneExpressions.has(gene) && !genePlots[gene]) {
      // Check if we have data in the active layer's cache
      const hasData = geneExpressionCache.has(gene);
      
      if (hasData) {
        console.log(`Adding gene: ${gene} for layer: ${activeLayer}`);
      } else {
        console.warn(`No expression data for gene: ${gene} in layer: ${activeLayer}`);
      }
      activeGeneExpressions.add(gene);
      genesToAdd.push(gene);
    } else if (genePlots[gene]) {
      console.log(`Gene plot for ${gene} already exists - skipping`);
    }
  }

  // Create plots for new genes
  console.log(`About to create ${genesToAdd.length} new gene plots:`, genesToAdd);
  const createPromises = genesToAdd.map(gene => createGenePlot(gene));

  Promise.all(createPromises).then(() => {
    console.log(`Successfully created ${genesToAdd.length} gene plots`);
    
    // Create sync button only when we have gene plots available
    if (activeGeneExpressions.size >= 1 && !document.getElementById('syncButton')) {
      createSyncButton();
    }

    // Draw the gene plots with data from the active layer
    genesToAdd.forEach(gene => {
      console.log(`Drawing data for gene plot: ${gene} in layer: ${activeLayer}`);
      redrawGenePlot(gene);
    });

    console.log(`Current sync mode state: ${isSyncMode}`);
    console.log("View inheritance complete - plots are now independent unless sync is manually enabled");

    updatePlotLayout();
  });

  // Cleanup if no genes are active
  if (activeGeneExpressions.size === 0) {
    console.log("No active genes left. Clearing legends and sync button.");
    removeSyncButton();
    isSyncMode = false;
    clearAllSyncHandlers();
  }
});

// Updated active gene expression cache based on layer
function updateActiveGeneExpressionCache() {
  geneExpressionCache.clear();
  
  // Select source cache based on activeLayer
  const sourceCache = activeLayer === 'MAGIC' ? geneExpressionMAGIC : geneExpressionOriginal;

  sourceCache.forEach((geneData, geneName) => {
    geneExpressionCache.set(geneName, geneData);
  });

  console.log(`Active cache updated with ${activeLayer} expression data`);
}

// **NEW: Helper function to manually sync all gene plots to main plot view**
function syncAllGenePlotViewsToMain() {
  if (!scatterplot || scatterplot._destroyed) {
    console.warn('Main plot not available for view sync');
    return;
  }

  try {
    const mainCameraView = scatterplot.get('cameraView');
    if (!mainCameraView) {
      console.warn('Main plot camera view not available');
      return;
    }

    console.log('Syncing all gene plot views to main plot');
    
    Object.keys(genePlots).forEach(geneId => {
      const plot = genePlots[geneId];
      if (plot && !plot._destroyed) {
        try {
          plot.set({ cameraView: mainCameraView }, { preventEvent: true });
          console.log(`Synced view for ${geneId}`);
        } catch (error) {
          console.warn(`Failed to sync view for ${geneId}:`, error);
        }
      }
    });
  } catch (error) {
    console.error('Error syncing gene plot views to main:', error);
  }
}


function createMAGICIndicator() {
  if (document.getElementById('magicIndicator')) return;

  const indicator = document.createElement('div');
  indicator.id = 'magicIndicator';
  indicator.style.cssText = `
    position: absolute;
    top: 10px;
    left: 10px;
    z-index: 1000;
    padding: 4px 8px;
    background: ${isMAGICActive ? '#28a745' : '#6c757d'};
    color: white;
    border-radius: 4px;
    font-size: 12px;
    font-weight: bold;
  `;
  indicator.textContent = isMAGICActive ? 'MAGIC ON' : 'MAGIC OFF';

  const plotGrid = document.getElementById('plotGrid');
  if (plotGrid) {
    plotGrid.appendChild(indicator);
  }
}

// Update indicator when MAGIC state changes
function updateMAGICIndicator() {
  const indicator = document.getElementById('magicIndicator');
  if (indicator) {
    indicator.style.background = isMAGICActive ? '#28a745' : '#6c757d';
    indicator.textContent = isMAGICActive ? 'MAGIC ON' : 'MAGIC OFF';
  }
}

// Enhanced toggle MAGIC function
function toggleMAGIC(isActivated) {
  console.log(`Toggling MAGIC: ${isActivated ? 'ON' : 'OFF'}`);
  
  isMAGICActive = isActivated;
  
  // Update coordinates in the points array
  const coordSource = isActivated ? pointsMAGIC : pointsXY;
  for (let i = 0; i < points.length; i++) {
    points[i][0] = coordSource[i][0];  // Replace x
    points[i][1] = coordSource[i][1];  // Replace y
  }

  // Update filteredPoints after coordinate changes
  // This ensures all plots work with the same updated coordinate set
  filteredPoints = filterPoints();

  // Update gene expression cache to use appropriate data
  updateActiveGeneExpressionCache();
  
  // Check if we have MAGIC expression data when activating
  if (isActivated && geneExpressionMAGIC.size === 0) {
    console.warn("MAGIC coordinates activated but no MAGIC expression data available");
  }
  
  // Update indicator
  updateMAGICIndicator();
  
  // Enable transition for smooth animation
  useTransition = true;
  
  // Redraw all plots with new coordinates and expression data
  redrawAllPlots(mainAnnotation);

  // Extend timeout to ensure transitions complete
  setTimeout(() => {
    useTransition = false;
    console.log("MAGIC transition completed");
  }, 1500); // Increased from 1000ms to 1500ms
}


// Update the Shiny message handler
Shiny.addCustomMessageHandler('toggleMAGIC', function(isActivated) {
  toggleMAGIC(isActivated);
});

// Handle window resize
window.addEventListener('resize', () => {
  if (scatterplot) scatterplot.refresh();
  Object.values(genePlots).forEach(plot => plot.refresh());
});

window.geneDebug = {
  getAvailableGenes,
  logPerformanceStats,
  getGeneExpressionValue,
  geneExpressionCache: () => geneExpressionCache,
  geneDataInfo: () => geneDataInfo
};

window.geneDebug = {
  ...window.geneDebug,
  setVmaxMode,
  createVmaxControls,
  calculatePercentile
};

window.geneDebug = {
  ...window.geneDebug,
  geneExpressionOriginal: () => geneExpressionOriginal,
  geneExpressionMAGIC: () => geneExpressionMAGIC,
  isMAGICActive: () => isMAGICActive,
  createMAGICIndicator,
  toggleMAGIC
};

window.geneDebug = {
  ...window.geneDebug,
  verifyCoordinateConsistency
};

window.geneDebug = {
  ...window.geneDebug,
  syncAllGenePlotViewsToMain,
  getCurrentMainView: () => {
    if (scatterplot && !scatterplot._destroyed) {
      return scatterplot.get('cameraView');
    }
    return null;
  }
};

// Initialize
initScatterplot();


// www/capture_canvas.js
Shiny.addCustomMessageHandler('captureCanvases', function(message) {
  console.log("Received captureCanvases message:", message);
  console.log("Starting canvas capture...");
  
  // Add small delay to ensure all DOM elements are ready
  setTimeout(() => {
    const canvases = document.querySelectorAll('canvas[id^="scatterplot_canvas"], canvas[id^="gene_canvas_"]');
    console.log("Found", canvases.length, "canvases");

    const canvasData = {};

    // Ensure html2canvas is loaded
    if (typeof html2canvas === 'undefined') {
      console.error("html2canvas is not loaded. Please include the library.");
      Shiny.setInputValue('canvas_data', canvasData, { priority: 'event' });
      return;
    }

    if (canvases.length === 0) {
      console.warn("No canvases found to capture");
      Shiny.setInputValue('canvas_data', canvasData, { priority: 'event' });
      return;
    }

    // Process canvases sequentially to avoid rendering issues
    const capturePromises = Array.from(canvases).map((canvas, index) => {
      return new Promise(resolve => {
        // Add small delay between captures to prevent conflicts
        setTimeout(() => {
          console.log("Capturing canvas:", canvas.id, `(${index + 1}/${canvases.length})`);
          
          // Get the plot container (parent of canvas, which includes the legend)
          const plotContainer = canvas.parentElement;
          if (!plotContainer) {
            console.error("No plot container found for canvas:", canvas.id);
            resolve();
            return;
          }

          // Ensure the canvas is visible and has content
          if (canvas.width === 0 || canvas.height === 0) {
            console.warn("Canvas has zero dimensions:", canvas.id);
            resolve();
            return;
          }

          // Use html2canvas to capture the entire plot container (canvas + legend)
          html2canvas(plotContainer, {
            backgroundColor: '#FFFFFF', // Ensure white background
            scale: 2, // Increase resolution for better quality
            logging: false, // Disable for production
            useCORS: true,
            allowTaint: true,
            width: plotContainer.offsetWidth,
            height: plotContainer.offsetHeight,
            scrollX: 0,
            scrollY: 0
          }).then(offscreenCanvas => {
            // Verify the captured canvas has content
            if (offscreenCanvas.width === 0 || offscreenCanvas.height === 0) {
              console.error("Captured canvas has zero dimensions:", canvas.id);
              resolve();
              return;
            }
            
            // Export as JPEG
            const dataURL = offscreenCanvas.toDataURL('image/jpeg', 0.8);
            canvasData[canvas.id] = dataURL;
            console.log("Captured canvas with legend:", canvas.id, 
                       `(${offscreenCanvas.width}x${offscreenCanvas.height})`);
            resolve();
          }).catch(error => {
            console.error("Error capturing canvas:", canvas.id, error);
            // Fallback to capturing only the canvas
            try {
              const offscreen = document.createElement('canvas');
              offscreen.width = canvas.width;
              offscreen.height = canvas.height;
              const offctx = offscreen.getContext('2d');
              offctx.fillStyle = '#FFFFFF';
              offctx.fillRect(0, 0, offscreen.width, offscreen.height);
              offctx.drawImage(canvas, 0, 0);
              const dataURL = offscreen.toDataURL('image/jpeg', 0.8);
              canvasData[canvas.id] = dataURL;
              console.log("Fallback: Captured canvas without legend:", canvas.id);
            } catch (fallbackError) {
              console.error("Fallback capture also failed:", canvas.id, fallbackError);
            }
            resolve();
          });
        }, index * 100); // Stagger captures by 100ms each
      });
    });

    // Wait for all captures to complete
    Promise.all(capturePromises).then(() => {
      console.log("Sending canvas data to Shiny:", Object.keys(canvasData));
      console.log("Total canvases captured:", Object.keys(canvasData).length);
      
      // Always send with unique timestamp to force Shiny to detect changes
      const uniqueData = {
        ...canvasData,
        timestamp: Date.now(),
        capture_id: Math.random().toString(36).substr(2, 9) // Additional uniqueness
      };
      
      Shiny.setInputValue('canvas_data', uniqueData, { priority: 'event' });
    }).catch(error => {
      console.error("Error in Promise.all:", error);
      // Still send data even if there are errors, with unique identifier
      const uniqueData = {
        ...canvasData,
        timestamp: Date.now(),
        capture_id: Math.random().toString(36).substr(2, 9),
        error: true
      };
      Shiny.setInputValue('canvas_data', uniqueData, { priority: 'event' });
    });
    
  }, 100); // Initial delay of 100ms
});