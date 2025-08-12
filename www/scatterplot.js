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

// Viridis color palette for gene expression
const viridisColors = [
  '#440154', '#482777', '#3f4a8a', '#31678e', '#26838f', '#1f9d8a', 
  '#6cce5a', '#b6de2b', '#fee825', '#fcce25'
];

function createSyncButton() {
  // Check if button already exists
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
  syncButton.textContent = 'Enable Sync Mode';
  syncButton.style.cssText = `
    padding: 8px 16px;
    background: #007bff;
    color: white;
    border: none;
    border-radius: 4px;
    cursor: pointer;
    font-size: 14px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.2);
  `;

  syncButton.onclick = toggleSyncMode;
  buttonContainer.appendChild(syncButton);
  
  // Add to the plot grid container
  const plotGrid = document.getElementById('plotGrid');
  plotGrid.style.position = 'relative';
  plotGrid.appendChild(buttonContainer);
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

function createPlotTitle(canvas, title) {
  // Ensure each canvas is wrapped in its own relative container
  let container = canvas.parentElement;
  if (!container || !container.classList.contains('canvas-container')) {
    // Create wrapper div
    container = document.createElement('div');
    container.className = 'canvas-container';
    container.style.cssText = `
      position: relative;
      width: 100%;
      height: 100%;
      display: flex;
      align-items: stretch;
    `;

    // Insert container before canvas and move canvas inside it
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
    pointer-events: none;
    z-index: 10;
  `;
  titleOverlay.textContent = title;

  // Append to this canvas container
  container.appendChild(titleOverlay);
}

function createPlotLegend(plotId, legendData, type = 'categorical') {
  const container = document.getElementById(plotId === 'main' ? 'scatterplot_canvas' : `gene_canvas_${plotId}`);
  if (!container) return;

  const plotContainer = container.parentElement;
  if (!plotContainer) return;

  // Remove existing legend
  const existingLegend = plotContainer.querySelector('.plot-legend');
  if (existingLegend) existingLegend.remove();

  // Ensure plotContainer is relative for absolute legend placement
  plotContainer.style.position = 'relative';

  // Create overlay legend
  const legendContainer = document.createElement('div');
  legendContainer.className = 'plot-legend';
  legendContainer.style.cssText = `
    position: absolute;
    top: 8px;
    right: 8px; /* Change to 'left: 8px;' for top-left placement */
    width: 120px;
    background: rgba(255, 255, 255, 0.9);
    border: 1px solid #ddd;
    font-size: 11px;
    padding: 6px;
    z-index: 10; /* ensure it floats above plot */
    border-radius: 4px;
    pointer-events: auto; /* allow interaction */
  `;

  if (type === 'gene') {
    const { minVal, maxVal, midVal, geneName } = legendData;
    legendContainer.innerHTML = `
      <div style="font-weight: bold; margin-bottom: 6px; text-align: center;">${geneName}</div>
      <div style="display: flex; flex-direction: column; align-items: center;">
        <div style="font-size: 10px;">High</div>
        <div style="width: 20px; height: 80px;
          background: linear-gradient(to top, ${viridisColors[0]}, ${viridisColors[Math.floor(viridisColors.length/2)]}, ${viridisColors[viridisColors.length-1]});
          border: 1px solid #ccc; margin: 4px 0;">
        </div>
        <div style="font-size: 10px;">Low</div>
      </div>
      <div style="font-size: 9px; text-align: center; margin-top: 4px; color: #666;">
        Max: ${maxVal}<br>Mid: ${midVal}<br>Min: ${minVal}
      </div>
    `;
  } else {
    const { names, colors, visible, annotationName } = legendData;

    const titleDiv = document.createElement('div');
    titleDiv.style.cssText = `
      font-weight: bold;
      margin-bottom: 6px;
      text-align: center;
      border-bottom: 1px solid #eee;
      padding-bottom: 4px;
    `;
    titleDiv.textContent = annotationName.replace(/_/g, ' ');
    legendContainer.appendChild(titleDiv);

    const itemsContainer = document.createElement('div');
    itemsContainer.style.cssText = `
      max-height: 150px;
      overflow-y: auto;
    `;

    names.forEach((name, i) => {
      const item = document.createElement('div');
      item.className = 'legend-item';
      item.style.cssText = `
        display: flex;
        align-items: center;
        padding: 3px;
        cursor: pointer;
        border-radius: 3px;
        opacity: ${visible.has(i) ? '1' : '0.4'};
        transition: opacity 0.2s;
      `;
      item.innerHTML = `
        <div style="width: 12px; height: 12px; background: ${colors[i]};
          margin-right: 6px; border-radius: 2px;"></div>
        <span style="font-size: 10px; color: #333;">${name}</span>
      `;

      item.onclick = () => {
        if (visible.has(i)) {
          visible.delete(i);
          item.style.opacity = '0.4';
        } else {
          visible.add(i);
          item.style.opacity = '1';
        }
        redrawAllPlots();
      };

      itemsContainer.appendChild(item);
    });

    legendContainer.appendChild(itemsContainer);
  }

  // Add legend as overlay
  plotContainer.appendChild(legendContainer);
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
  createPlotTitle(canvas, geneId);
  
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
  let pointSize = 2.5;
  let opacity = 0.8;
  let performanceMode = false;

  if (pointCount >= 1000000) {
    pointSize = 2;
    opacity = 0.6;
    performanceMode = true;
  } else if (pointCount >= 500000) {
    pointSize = 1;
    opacity = 0.4;
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
      const originalIndices = selectedIndices.map(i => points.indexOf(filteredPoints[i])).filter(i => i !== -1);
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

  genePlots[geneId] = plot;
  console.log(`Gene plot ${geneId} successfully created and stored`);
  
  // **Inherit the main plot's current view ONLY during creation**
  // This is completely separate from sync system
  if (scatterplot && !scatterplot._destroyed) {
    try {

      const mainCameraView = scatterplot.get('cameraView');
      if (mainCameraView) {
        console.log(`Inheriting main plot view for ${geneId} during creation (one-time only)`);
        // Apply the view without triggering events
        plot.set({ cameraView: mainCameraView }, { preventEvent: true });

        console.log("awi");
        removeSynchronization();
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
  createPlotTitle(canvas, 'Main View');

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
    pointSize: 2.5,
    opacity: 0.8,
    lassoOnLongPress: true,
    lassoType: 'freeform',
  });

  subscribeMainPlotEvents();
  updatePlotLayout();

  initScatterplotPromise = Promise.resolve();
}

function subscribeMainPlotEvents() {
  scatterplot.subscribe('select', ({ points: selectedIndices }) => {
    if (!syncing) {
      const originalIndices = selectedIndices.map(i => points.indexOf(filteredPoints[i])).filter(i => i !== -1);
      Shiny.setInputValue('selectedPoints', originalIndices);
      
      // Sync selection to other plots if in sync mode
      if (isSyncMode) {
        syncSelectionToAllPlots(selectedIndices, 'main');
      }
    }
  });
  
  scatterplot.subscribe('deselect', () => {
    if (!syncing) {
      Shiny.setInputValue('selectedPoints', []);
      
      // Sync deselection to other plots if in sync mode
      if (isSyncMode) {
        syncDeselectionToAllPlots('main');
      }
    }
  });
}

function filterPoints() {
  if (Object.keys(annotations).length === 0) return points;

  // console.log(gene_number, annotations);
  // console.log(points.length, points[0]);
  
  return points.filter(p => {
    return !Object.keys(annotations).some((key, idx) => {
      const valueIndex = 2 + gene_number + idx;
      const value = p[valueIndex] - 1;
      // console.log(`Checking annotation ${key} with value ${value} (index ${valueIndex})`);
      if (!annotations[key].visible.has(value)) {
        return true;
      }
      return false;
    });
  });
}

function redrawAllPlots() {
  if (points.length === 0) return;
  
  // Always update filteredPoints before redrawing
  // This ensures all plots use the most current coordinates
  filteredPoints = filterPoints();
  
  // Redraw main plot
  redrawMainPlot();
  
  // Redraw gene plots
  Object.keys(genePlots).forEach(geneId => {
    redrawGenePlot(geneId);
  });
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


async function redrawMainPlot() {
  if (!scatterplot || filteredPoints.length === 0) return;

  if (currentColorBy === 'seacell') {
    const annotationKeys = Object.keys(annotations);
    const annotationIndex = annotationKeys.indexOf('seacell');
    const clusterIndex = annotationKeys.indexOf('cluster');

    if (annotationIndex === -1) return;

    const colorIdx = 2 + gene_number + annotationIndex;
    const clusterIdx = 2 + gene_number + clusterIndex;

    const regularCells = filteredPoints.filter(p => p[colorIdx] === 1);
    const metacells = filteredPoints.filter(p => p[colorIdx] === 2);

    let allPoints = [];
    regularCells.forEach(p => {
      allPoints.push([p[0], p[1], p[clusterIdx] - 1, 0]);
    });
    metacells.forEach(p => {
      allPoints.push([p[0], p[1], (p[clusterIdx] - 1) + clusters.length, 1]);
    });

    const config = { 
      colorBy: 'category', 
      sizeBy: 'value',
      pointSize: [5, 10],
      pointColor: clusterColors,
      opacity: 0.8
    };

    scatterplot.set(config);
    scatterplot.draw(allPoints);

    // Create legend for seacell mode (showing clusters)
    if (annotations['cluster']) {
      createPlotLegend('main', {
        names: annotations['cluster'].names,
        colors: annotations['cluster'].colors,
        visible: annotations['cluster'].visible,
        annotationName: 'Clusters'
      }, 'categorical');
    }
    
  } else if (currentColorBy && annotations[currentColorBy]) {
    const annotationKeys = Object.keys(annotations);
    const annotationIndex = annotationKeys.indexOf(currentColorBy);
    const colorIdx = 2 + gene_number + annotationIndex;
    let colors = annotations[currentColorBy].colors;

    const pointsColored = filteredPoints.map(p => {
      const newP = [...p];
      newP[2] = p[colorIdx] - 1;
      return newP;
    });

    if (pointsColored.length < 300) {
      pointSize = 100;
    } else {
      pointSize = 30;
    }

    
    const config = { 
      colorBy: 'category',
      pointColor: colors,
      sizeBy: 'category',
      pointSize: Array(pointsColored.length).fill(pointSize),
    };
    
    scatterplot.set(config);
    scatterplot.draw(pointsColored, { transition: useTransition });
    
    // if initScatterplot has been initialized, set x and y boundaries
    if (initScatterplotPromise) {

      // Compute bounds for x and y
      const xs = filteredPoints.map(p => p[0]);
      const ys = filteredPoints.map(p => p[1]);
      const minX = Math.min(...xs);
      const maxX = Math.max(...xs);
      const minY = Math.min(...ys);
      const maxY = Math.max(...ys);
      const padX = (maxX - minX) * 0.1 || 1;
      const padY = (maxY - minY) * 0.1 || 1;
      // let indices = Array.from({ length: pointsColored.length }, (_, i) => i);
      // console.log(indices);
      // scatterplot.zoomToPoints(indices);
      scatterplot.zoomToArea(
        { x: minX - padX, y: minY - padY,
          width: maxX + padX - (minX - padX),
          height: maxY + padY - (minY - padY)
        }
      )

      // set initScatterplotPromise to null
      initScatterplotPromise = null;
    }

    // Create legend for current annotation
    createPlotLegend('main', {
      names: annotations[currentColorBy].names,
      colors: annotations[currentColorBy].colors,
      visible: annotations[currentColorBy].visible,
      annotationName: currentColorBy
    }, 'categorical');
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
      // Original data
      const geneValues = [];
      for (let row = 0; row < geneDataInfo.nrows; row++) {
        geneValues[row] = geneDataMatrix[row][geneIndex];
      }
      
      const min = Math.min(...geneValues);
      const max = Math.max(...geneValues);
      const mean = geneValues.reduce((a, b) => a + b, 0) / geneValues.length;
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
        const magicGeneValues = [];
        for (let row = 0; row < geneDataInfo.nrows; row++) {
          magicGeneValues[row] = magicDataMatrix[row][geneIndex];
        }
        
        const magicMin = Math.min(...magicGeneValues);
        const magicMax = Math.max(...magicGeneValues);
        const magicMean = magicGeneValues.reduce((a, b) => a + b, 0) / magicGeneValues.length;
        const magicP90 = calculatePercentile(magicGeneValues, 90);
        const magicP95 = calculatePercentile(magicGeneValues, 95);
        const magicP99 = calculatePercentile(magicGeneValues, 99);
        
        const magicGeneData = {
          values: magicGeneValues,
          range: expressionData.magic_ranges[geneName] || { min: magicMin, max: magicMax, mean: magicMean },
          percentiles: { p90: magicP90, p95: magicP95, p99: magicP99, min: magicMin, max: magicMax, mean: magicMean }
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

function updateActiveGeneExpressionCache() {
  geneExpressionCache.clear();
  
  // Always use original expression data, regardless of MAGIC state
  // MAGIC only affects coordinates, not gene expression values
  const sourceCache = geneExpressionOriginal;
  
  sourceCache.forEach((geneData, geneName) => {
    geneExpressionCache.set(geneName, geneData);
  });
  
  console.log(`Active cache updated with original expression data (MAGIC affects coordinates only)`);
}

// Add a global setting for vmax mode
let vmaxMode = 'max'; // Options: 'max', 'p90', 'p95', 'p99'

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

// Updated redrawGenePlot function with percentile scaling
function redrawGenePlot(geneId) {
  const plot = genePlots[geneId];
  if (!plot || filteredPoints.length === 0) return;
  
  // Check if we have cached gene expression data for this gene
  const geneData = geneExpressionCache.get(geneId);
  if (!geneData) {
    console.warn(`No expression data found for gene: ${geneId}`);
    return;
  }
  
  console.log(`Redrawing plot for ${geneId} with ${geneData.values.length} expression values`);
  
  // Determine vmax based on mode first
  let vmax;
  let vmaxLabel;
  
  switch(vmaxMode) {
    case 'p90':
      vmax = geneData.percentiles.p90;
      vmaxLabel = `P90: ${vmax.toFixed(2)}`;
      break;
    case 'p95':
      vmax = geneData.percentiles.p95;
      vmaxLabel = `P95: ${vmax.toFixed(2)}`;
      break;
    case 'p99':
      vmax = geneData.percentiles.p99;
      vmaxLabel = `P99: ${vmax.toFixed(2)}`;
      break;
    case 'max':
    default:
      vmax = geneData.percentiles.max;
      vmaxLabel = `Max: ${vmax.toFixed(2)}`;
      break;
  }
  
  // Create points with gene expression values - normalize to [0,1] range
  const genePoints = [];
  
  for (let i = 0; i < filteredPoints.length; i++) {
    const point = filteredPoints[i];
    const originalIndex = points.indexOf(point);
    
    if (originalIndex !== -1 && originalIndex < geneData.values.length) {
      let expressionValue = geneData.values[originalIndex];
      
      // Normalize expression value to [0,1] based on vmax
      const normalizedValue = Math.min(expressionValue / vmax, 1.0);
      
      // Use coordinates from the points array (which gets updated by toggleMAGIC)
      // This ensures gene plots use the same coordinates as the main plot
      const x = point[0];  // This now contains the updated coordinates
      const y = point[1];  // This now contains the updated coordinates
      
      genePoints.push([x, y, normalizedValue]);
    }
  }
  
  if (genePoints.length === 0) {
    console.warn(`No valid points to plot for gene: ${geneId}`);
    return;
  }

  // Set up the plot configuration
  const config = {
    colorBy: 'valueA',
    pointColor: viridisColors,
    pointSize: 50,
    colorRange: [0, 1]
  };

  // Always use transition when useTransition is true
  // Don't clear the plot when transitioning to allow smooth animation
  if (!useTransition) {
    plot.clear();
  }
  
  plot.set(config);
  plot.draw(genePoints, { transition: useTransition });

  // Create gene expression legend with percentile info
  createPlotLegend(geneId, {
    minVal: geneData.percentiles.min.toFixed(2),
    maxVal: vmaxLabel,
    midVal: geneData.percentiles.mean.toFixed(2),
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

Shiny.addCustomMessageHandler('updateData', function(message) {
  try {
    const { base64, annotationData, numCols, clusters: cl, colors, metacellColors, geneExprRanges: ge } = message;

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
      pointsMAGIC.push([magicX, magicY]);
    }

    if ('seacell' in annotationData) {
      clusterColors = [...colors, ...metacellColors];
    }

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

    redrawAllPlots();
    setSpinner(false);
    
  } catch (error) {
    console.error('Error in updateData:', error);
    setSpinner(false);
  }
});

Shiny.addCustomMessageHandler('colorByChange', function(colorBy) {
  //  check if seacell is in annotations
  if (colorBy && !annotations[colorBy] && 'seacell' in annotations) {
    if (colorBy === 'seacell') {
      annotations['seacell'].visible.add(1);
    } else {
      annotations['seacell'].visible.delete(1);
    }
  }

  currentColorBy = colorBy;
  
  // Remove sync button for non-gene expression modes
  if (!activeGeneExpressions.size) {
    removeSyncButton();
    
    // If currently in sync mode, disable it
    if (isSyncMode) {
      toggleSyncMode();
    }
  }

  redrawAllPlots();
});

// Modify the existing geneSearchChange handler
Shiny.addCustomMessageHandler('geneSearchChange', function(message) {
  console.log("Gene search change message:", message);

  // --- Normalize gene names to always be an array ---
  const geneNames = message.genes 
    ? (Array.isArray(message.genes) ? message.genes : [message.genes])
    : [];

  // --- Normalize expression data and its genes ---
  const expressionData = message.expression_data || {};
  expressionData.genes = expressionData.genes 
    ? (Array.isArray(expressionData.genes) ? expressionData.genes : [expressionData.genes])
    : [];

  console.log("Processed gene names:", geneNames);
  console.log("Expression data:", expressionData);

  // --- Process both original and MAGIC gene expression data ---
  processGeneExpressionData(expressionData);

  console.log("Active genes before update:", Array.from(activeGeneExpressions));

  // Filter out empty strings and ensure valid names
  const newGenes = new Set(
    geneNames.filter(name => name && typeof name === 'string' && name.trim() !== '')
  );

  console.log("New genes to process:", Array.from(newGenes));

  // --- Remove deselected genes ---
  for (const gene of Array.from(activeGeneExpressions)) {
    if (!newGenes.has(gene)) {
      console.log(`Removing gene: ${gene}`);
      activeGeneExpressions.delete(gene);
      removeGenePlot(gene);
    }
  }

  // --- Add new genes (with duplicate protection) ---
  const genesToAdd = [];
  for (const gene of newGenes) {
    if (!activeGeneExpressions.has(gene) && !genePlots[gene]) { // Extra check to prevent duplicates
      // Check if we have data in either original or MAGIC cache
      const hasOriginal = geneExpressionOriginal.has(gene);
      const hasMAGIC = geneExpressionMAGIC.has(gene);
      
      if (hasOriginal || hasMAGIC) {
        console.log(`Adding gene: ${gene} (Original: ${hasOriginal}, MAGIC: ${hasMAGIC})`);
      } else {
        console.warn(`No expression data for gene: ${gene}`);
      }
      activeGeneExpressions.add(gene);
      genesToAdd.push(gene);
    } else if (genePlots[gene]) {
      console.log(`Gene plot for ${gene} already exists - skipping`);
    }
  }

  // --- Create plots for new genes (view inheritance happens inside createGenePlot) ---
  console.log(`About to create ${genesToAdd.length} new gene plots:`, genesToAdd);
  const createPromises = genesToAdd.map(gene => createGenePlot(gene));

  Promise.all(createPromises).then(() => {
    console.log(`Successfully created ${genesToAdd.length} gene plots`);
    
    // Create sync button only when we have gene plots available
    if (activeGeneExpressions.size >= 1 && !document.getElementById('syncButton')) {
      createSyncButton();
    }

    // Draw the gene plots with data
    genesToAdd.forEach(gene => {
      console.log(`Drawing data for gene plot: ${gene}`);
      redrawGenePlot(gene);
    });

    // COMPLETELY AVOID sync setup during creation
    // Sync will only happen when user manually toggles the sync button
    console.log(`Current sync mode state: ${isSyncMode}`);
    console.log("View inheritance complete - plots are now independent unless sync is manually enabled");

    updatePlotLayout();
  });

  // --- Cleanup if no genes are active ---
  if (activeGeneExpressions.size === 0) {
    console.log("No active genes left. Clearing legends and sync button.");
    removeSyncButton();
    isSyncMode = false;
    clearAllSyncHandlers();
  }
});

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
  redrawAllPlots();

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