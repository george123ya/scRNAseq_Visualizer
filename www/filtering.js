// filter.js - FIXED: Proper filter removal and clear functionality
// Place in www/filter.js

(function() {
  var numericVars = [];
  var filterRows = [];
  var nextRowId = 1;
  var varData = {};
  var currentFilters = [];

  // 🆕 Clear category filters
  Shiny.addCustomMessageHandler("clearCategoryFilters", function(data) {
    console.log("🗑️ Clearing category filters");
    
    if (window.categoryFilterManager) {
      window.categoryFilterManager.clearAll();
    }
  });

  // Handle restoreState with optional filtering
  Shiny.addCustomMessageHandler("restoreState", function(data) {
    console.log("🔄 Restoring state:", data);

    // 🆕 HANDLE ORIGINAL VERSION - COMPLETE RESET
    if (data.isOriginal) {
      console.log("   🌍 ORIGINAL VERSION - Full reset");
      
      // Clear all filters
      if (window.scatterplotManager) {
        window.scatterplotManager.clearFilter();
      }
      
      // Clear category filters
      if (window.categoryFilterManager) {
        window.categoryFilterManager.clearAll();
      }
      
      // Clear lasso selection
      if (window.lassoManager) {
        window.lassoManager.clearSelection();
      }
      
      // Reset plots to full data
      redrawAllPlots();
      
      return;
    }
    
    if (data.filterToSelection && data.cellMask) {
      // FILTER MODE: Show only selected cells
      const cellMask = data.cellMask;
      const selectedIndices = [];
      
      for (let i = 0; i < cellMask.length; i++) {
        if (cellMask[i]) {
          selectedIndices.push(i);
        }
      }
      
      console.log("   🔍 FILTERING to", selectedIndices.length, "cells");
      
      // Update all plots to show only these cells
      if (window.scatterplotManager) {
        window.scatterplotManager.filterToIndices(selectedIndices);
      }
      
    } else if (data.cellMask) {
      // HIGHLIGHT MODE: Show all cells, highlight filtered
      console.log("   ✨ HIGHLIGHTING filtered cells");
      
      // Restore category visibility
      if (data.visibleCategories && data.annotationName) {
        restoreCategoryVisibility(data.annotationName, data.visibleCategories);
      }
      
      // Restore lasso selection if any
      if (data.selectedPoints) {
        restoreLassoSelection(data.selectedPoints);
      }
      
    } else {
      // NORMAL MODE: Show everything
      console.log("   🌍 Showing all cells");
      
      if (window.scatterplotManager) {
        window.scatterplotManager.clearFilter();
      }
    }
  });

  // Restore numeric filters in UI
  Shiny.addCustomMessageHandler("restoreNumericFilters", function(data) {
    console.log("🔢 Restoring numeric filters:", data.filters);
    
    // Rebuild filter UI (you'll need to implement this based on your filter.js)
    if (window.filterManager) {
      window.filterManager.restoreFilters(data.filters);
    }
  });

  Shiny.addCustomMessageHandler("updateNumericVars", function(message) {
    numericVars = message;
    initFilterMenu();
  });

  Shiny.addCustomMessageHandler("updateVarData", function(message) {
    varData[message.var] = message.data;
    updateHistogramsForVar(message.var);
  });

  function initFilterMenu() {
    var container = document.getElementById('filterControls');
    if (!container) return;
    container.innerHTML = `
      <style>
        body {
          font-family: 'Roboto', sans-serif;
        }
        #toggleMenu {
          position: fixed;
          top: 16px;
          right: 16px;
          z-index: 1001;
          background: linear-gradient(135deg, #667eea, #764ba2);
          border: none;
          color: white;
          font-size: 24px;
          font-weight: 500;
          width: 48px;
          height: 48px;
          border-radius: 50%;
          display: flex;
          align-items: center;
          justify-content: center;
          cursor: pointer;
          box-shadow: 0 2px 8px rgba(0, 0, 0, 0.15);
          transition: transform 0.2s ease, box-shadow 0.2s ease;
        }
        #toggleMenu:hover {
          transform: translateY(-2px);
          box-shadow: 0 4px 12px rgba(0, 0, 0, 0.2);
        }
        #toggleMenu.active::after { content: '▼'; }
        #toggleMenu::after { content: '▲'; }
        #filterMenu {
          display: none;
          position: fixed;
          top: 80px;
          right: 16px;
          background: #ffffff;
          border-radius: 12px;
          padding: 20px;
          box-shadow: 0 4px 16px rgba(0, 0, 0, 0.1);
          z-index: 1000;
          width: 480px;
          max-width: 95vw;
          opacity: 0;
          transform: translateY(-10px);
          transition: opacity 0.3s ease, transform 0.3s ease;
          border: 1px solid #e0e0e0;
        }
        #filterMenu.open {
          display: block;
          opacity: 1;
          transform: translateY(0);
        }
        #filtersContainer {
          display: grid;
          gap: 12px;
          max-height: 400px;
          overflow-y: auto;
          padding-right: 8px;
        }
        .filterRow {
          display: grid;
          grid-template-columns: 120px 80px 120px 80px 40px;
          align-items: center;
          gap: 8px;
          background: #f9f9f9;
          padding: 8px;
          border-radius: 8px;
          transition: background 0.2s ease;
        }
        .filterRow:hover {
          background: #f0f2f5;
        }
        .varSelect, .threshInput {
          border: 1px solid #d0d0d0;
          border-radius: 6px;
          padding: 6px;
          font-size: 14px;
          background: white;
          color: #333;
          outline: none;
          transition: border-color 0.2s ease;
        }
        .varSelect:focus, .threshInput:focus {
          border-color: #B0BEC5;
        }
        .varSelect {
          width: 100%;
        }
        .threshInput {
          width: 100%;
        }
        .threshInput::placeholder {
          color: #888;
        }
        .histSvg {
          width: 120px;
          height: 100px;
        }
        .removeRow {
          width: 40px;
          height: 40px;
          border: none;
          border-radius: 6px;
          background: #e74c3c;
          color: white;
          font-size: 16px;
          cursor: pointer;
          display: flex;
          align-items: center;
          justify-content: center;
          transition: transform 0.2s ease, box-shadow 0.2s ease;
        }
        .removeRow:hover {
          transform: translateY(-1px);
          box-shadow: 0 2px 8px rgba(0, 0, 0, 0.15);
        }
        #buttonContainer {
          display: flex;
          justify-content: flex-end;
          gap: 8px;
          margin-top: 16px;
          padding-top: 8px;
          border-top: 1px solid #e0e0e0;
        }
        .btn-primary {
          padding: 8px 16px;
          border: none;
          border-radius: 6px;
          background: linear-gradient(135deg, #667eea, #764ba2);
          color: white;
          font-size: 14px;
          font-weight: 500;
          cursor: pointer;
          transition: transform 0.2s ease, box-shadow 0.2s ease;
        }
        .btn-primary:hover {
          transform: translateY(-1px);
          box-shadow: 0 2px 8px rgba(0, 0, 0, 0.15);
        }
        #output {
          margin-top: 20px;
          font-size: 14px;
          color: #333;
        }
        .brush .selection {
          fill: rgba(207, 216, 220, 0.2);
          stroke: #B0BEC5;
          stroke-width: 1;
        }
        .brush .handle {
          fill: #B0BEC5;
          cursor: ew-resize;
        }
      </style>
      <button id="toggleMenu"></button>
      <div id="filterMenu">
        <div id="filtersContainer"></div>
        <div id="buttonContainer">
          <button id="addFilter" class="btn-primary">+ Add Filter</button>
          <button id="clearFilters" class="btn-primary" style="background: #6c757d;">Clear</button>
          <button id="applyFilters" class="btn-primary">Apply</button>
        </div>
      </div>
    `;
    addFilterRow();
    document.getElementById('toggleMenu').onclick = toggleMenu;
    document.getElementById('addFilter').onclick = addFilterRow;
    document.getElementById('applyFilters').onclick = applyFilters;
    
    // Clear button: reset UI and send empty array
    var clearButton = document.getElementById('clearFilters');
    if (clearButton) {
      clearButton.onclick = function() {
        console.log('\n🗑️ ==================== CLEAR FILTERS ====================');
        console.log('🗑️ Clear button clicked');
        
        // Clear all filter rows and reset UI
        var container = document.getElementById('filtersContainer');
        container.innerHTML = '';
        filterRows = [];
        nextRowId = 1;
        currentFilters = [];
        
        // Add one empty row back
        addFilterRow();
        
        // Send EMPTY array to Shiny (this will clear all numeric filters)
        console.log('📤 Sending EMPTY filter array [] to Shiny');
        Shiny.setInputValue('applyFilters', JSON.stringify([]), {priority: 'event'});
        
        console.log('✅ All filters cleared');
        console.log('🗑️ =====================================================\n');
      };
    }
  }

  function toggleMenu() {
    var menu = document.getElementById('filterMenu');
    var toggle = document.getElementById('toggleMenu');
    if (menu.classList.contains('open')) {
      menu.classList.remove('open');
      toggle.classList.remove('active');
    } else {
      menu.classList.add('open');
      toggle.classList.add('active');
    }
  }

  function addFilterRow() {
    var container = document.getElementById('filtersContainer');
    var rowId = nextRowId++;
    var row = document.createElement('div');
    row.className = 'filterRow';
    row.id = 'filterRow_' + rowId;
    row.innerHTML = `
      <select id="varSelect${rowId}" class="varSelect">
        <option value="">Select variable</option>
      </select>
      <input type="number" id="minThresh${rowId}" class="threshInput" placeholder="Min" step="any">
      <svg id="hist${rowId}" class="histSvg"></svg>
      <input type="number" id="maxThresh${rowId}" class="threshInput" placeholder="Max" step="any">
      <button class="removeRow" type="button" data-row-id="${rowId}">−</button>
    `;
    container.appendChild(row);
    
    var select = document.getElementById(`varSelect${rowId}`);
    numericVars.forEach(function(varName) {
      var opt = document.createElement('option');
      opt.value = varName;
      opt.text = varName;
      select.add(opt);
    });
    select.onchange = function() {
      var col = select.value;
      if (col && !varData[col]) {
        console.log('📥 Requesting data for column:', col);
        Shiny.setInputValue('requestVarData', {var: col, rowId: rowId}, {priority: 'event'});
      } else if (col && varData[col]) {
        console.log('💾 Using cached data for column:', col);
        updateHistogram(rowId);
      }
    };
    
    var minInput = document.getElementById(`minThresh${rowId}`);
    var maxInput = document.getElementById(`maxThresh${rowId}`);
    minInput.onchange = function() { syncBrushWithInputs(rowId); };
    maxInput.onchange = function() { syncBrushWithInputs(rowId); };
    
    // ⭐ FIX: Attach remove handler to button
    var removeBtn = row.querySelector('.removeRow');
    removeBtn.onclick = function(e) {
      e.preventDefault();
      removeRow(rowId);
    };
    
    filterRows.push(rowId);
    console.log('✅ Added filter row:', rowId, 'Total rows:', filterRows.length);
  }

  // ⭐ FIX: Proper removeRow function
  function removeRow(rowId) {
    console.log('🗑️ Removing row:', rowId);
    
    // Remove from DOM
    var rowElement = document.getElementById('filterRow_' + rowId);
    if (rowElement) {
      rowElement.remove();
      console.log('✅ Removed from DOM');
    }
    
    // Remove from tracking array
    var oldLength = filterRows.length;
    filterRows = filterRows.filter(id => id !== rowId);
    console.log('✅ Removed from tracking:', oldLength, '→', filterRows.length, 'rows');
    
    // If no rows left, add one empty row
    if (filterRows.length === 0) {
      console.log('ℹ️ No rows left - adding empty row');
      addFilterRow();
    }
  }

  function updateHistogramsForVar(varName) {
    filterRows.forEach(function(rowId) {
      var select = document.getElementById(`varSelect${rowId}`);
      if (select && select.value === varName) {
        updateHistogram(rowId);
      }
    });
  }

  function updateHistogram(rowId) {
    var select = document.getElementById(`varSelect${rowId}`);
    var col = select.value;
    var svg = d3.select(`#hist${rowId}`);
    svg.selectAll('*').remove();
    if (!col || !varData[col]) {
      svg.style('display', 'none');
      return;
    }
    svg.style('display', 'block');
    var values = varData[col].map(v => parseFloat(v)).filter(v => !isNaN(v));
    if (values.length === 0) {
      svg.style('display', 'none');
      return;
    }
    var width = 120, height = 100, margin = { top: 5, right: 5, bottom: 20, left: 5 };
    svg.attr('width', width).attr('height', height + margin.bottom);
    var minV = Math.min(...values);
    var maxV = Math.max(...values);
    var bins = 20;
    var binSize = (maxV - minV) / bins || 1;
    var histogram = d3.histogram()
      .value(d => d)
      .domain([minV, maxV])
      .thresholds(d3.range(minV, maxV + binSize, binSize));
    var binsData = histogram(values);
    var x = d3.scaleLinear()
      .domain([minV, maxV])
      .range([0, width - margin.left - margin.right]);
    var y = d3.scaleLinear()
      .domain([0, d3.max(binsData, d => d.length)])
      .range([height, 0]);
    var container = svg.append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);
    container.selectAll('rect')
      .data(binsData)
      .enter()
      .append('rect')
      .attr('x', d => x(d.x0))
      .attr('y', d => y(d.length))
      .attr('width', d => Math.max(0, x(d.x1) - x(d.x0)))
      .attr('height', d => height - y(d.length))
      .style('fill', 'linear-gradient(to right, #B0BEC5, #CFD8DC)')
      .style('stroke', '#CFD8DC')
      .style('stroke-width', 0.5);
    container.append('g')
      .attr('transform', `translate(0,${height})`)
      .call(d3.axisBottom(x)
        .ticks(2)
        .tickValues([minV, maxV])
        .tickFormat(d3.format('.2f')));
    svg.selectAll('.axis text').style('fill', '#333').style('font-size', '10px');
    svg.selectAll('.axis path, .axis line').style('stroke', '#e0e0e0');
    var brush = d3.brushX()
      .extent([[0, 0], [width - margin.left - margin.right, height]])
      .on('brush end', function(event) {
        var selection = d3.brushSelection(this);
        if (selection) {
          var [x0, x1] = selection.map(x.invert);
          document.getElementById(`minThresh${rowId}`).value = x0.toFixed(2);
          document.getElementById(`maxThresh${rowId}`).value = x1.toFixed(2);
        }
      });
    var brushGroup = container.append('g')
      .attr('class', 'brush')
      .call(brush)
      .call(brush.move, [x(minV), x(minV + binSize)]);
    syncBrushWithInputs(rowId);
  }

  function syncBrushWithInputs(rowId) {
    var svg = d3.select(`#hist${rowId}`);
    var minInput = document.getElementById(`minThresh${rowId}`);
    var maxInput = document.getElementById(`maxThresh${rowId}`);
    var minV = parseFloat(minInput.value) || 0;
    var maxV = parseFloat(maxInput.value) || 0;
    var brush = d3.select(`#hist${rowId} .brush`).call(d3.brushX().move, null);
    var x = d3.scaleLinear()
      .domain([d3.min(varData[document.getElementById(`varSelect${rowId}`).value].map(v => parseFloat(v)).filter(v => !isNaN(v))),
               d3.max(varData[document.getElementById(`varSelect${rowId}`).value].map(v => parseFloat(v)).filter(v => !isNaN(v)))])
      .range([0, 115]);
    brush.call(d3.brushX().move, [x(minV), x(maxV)]);
  }

  function applyFilters() {
    console.log('\n🔍 ==================== APPLY FILTERS ====================');
    console.log('🔍 Checking', filterRows.length, 'filter rows');
    
    var filters = [];
    var filtersAdded = 0;
    
    // Iterate through all active filter rows
    for (var i = 0; i < filterRows.length; i++) {
      var rowId = filterRows[i];
      
      var selectEl = document.getElementById(`varSelect${rowId}`);
      var minEl = document.getElementById(`minThresh${rowId}`);
      var maxEl = document.getElementById(`maxThresh${rowId}`);
      
      if (!selectEl || !minEl || !maxEl) {
        console.warn(`⚠️ Row ${rowId}: DOM elements missing`);
        continue;
      }
      
      var col = selectEl.value;
      var minThresh = parseFloat(minEl.value);
      var maxThresh = parseFloat(maxEl.value);
      
      console.log(`📋 Row ${rowId}: col="${col}", min=${minThresh}, max=${maxThresh}`);
      
      // Only add if all parameters are valid
      if (col && col !== "" && !isNaN(minThresh) && !isNaN(maxThresh) && minThresh <= maxThresh) {
        filters.push({col: col, minThresh: minThresh, maxThresh: maxThresh});
        filtersAdded++;
        console.log(`✅ Filter ${filtersAdded} added: ${col} ∈ [${minThresh}, ${maxThresh}]`);
      } else {
        console.warn(`⚠️ Row ${rowId} skipped - incomplete data`);
      }
    }
    
    console.log(`\n📊 Summary: ${filtersAdded} valid filters`);
    filters.forEach((f, idx) => {
      console.log(`   ${idx + 1}. ${f.col} in [${f.minThresh}, ${f.maxThresh}]`);
    });
    
    currentFilters = filters;
    
    console.log('\n📤 Sending to Shiny as JSON...');
    console.log('JSON:', JSON.stringify(filters));
    
    // ⭐ Send as JSON string
    Shiny.setInputValue('applyFilters', JSON.stringify(filters), {priority: 'event'});
    console.log('🔍 =====================================================\n');
  }

  window.clearFilterUI = function() {
    console.log('🗑️ Clearing filter UI');
    var container = document.getElementById('filtersContainer');
    if (container) {
      container.innerHTML = '';
    }
    filterRows = [];
    nextRowId = 1;
    varData = {};
    currentFilters = [];
    addFilterRow();
    console.log('✅ Filter UI cleared and reset');
  };

  window.filterDebug = function() {
    console.log('\n🔍 ===== FILTER DEBUG =====');
    console.log('Active filter rows:', filterRows);
    console.log('Current filters:', currentFilters);
    console.log('Available variables:', Object.keys(varData));
    filterRows.forEach(rowId => {
      const selectEl = document.getElementById(`varSelect${rowId}`);
      const minEl = document.getElementById(`minThresh${rowId}`);
      const maxEl = document.getElementById(`maxThresh${rowId}`);
      console.log(`  Row ${rowId}:`, {
        exists: !!(selectEl && minEl && maxEl),
        variable: selectEl?.value || "NONE",
        min: minEl?.value || "EMPTY",
        max: maxEl?.value || "EMPTY"
      });
    });
    console.log('==========================\n');
  };

  // ⭐ Don't export removeRow - it's internal
  window.applyNumericFilters = function(filters) {
    if (!filters || filters.length === 0) {
      console.log('📤 Resetting filters');
      Shiny.setInputValue('applyFilters', JSON.stringify([]), {priority: 'event'});
      return;
    }
    console.log('📤 Applying external filters:', filters);
    Shiny.setInputValue('applyFilters', JSON.stringify(filters), {priority: 'event'});
  };
})();