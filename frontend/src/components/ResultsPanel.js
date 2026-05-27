// Results Panel for displaying analysis results
export default class ResultsPanel {
    constructor() {
        this.nodeActivitiesContainer = document.getElementById('node-activities');
        this.influenceScoresContainer = document.getElementById('influence-scores');
        this.pathwayStatsContainer = document.getElementById('pathway-stats');
        
        this.init();
    }
    
    init() {
        console.log('ResultsPanel initialized');
    }
    
    displayResults(results) {
        if (!results) {
            console.warn('No results to display');
            return;
        }
        
        // Display node activities
        this.displayNodeActivities(results.node_activities || {});
        
        // Display influence scores
        this.displayInfluenceScores(results.influence_scores || {});
        
        // Display pathway statistics
        this.displayPathwayStats(results);
        
        console.log('Results displayed successfully');
    }
    
    displayNodeActivities(nodeActivities) {
        this.nodeActivitiesContainer.innerHTML = '';
        
        if (Object.keys(nodeActivities).length === 0) {
            this.nodeActivitiesContainer.innerHTML = '<p class="no-data">No activity data available</p>';
            return;
        }
        
        // Create table
        const table = document.createElement('table');
        table.className = 'results-table';
        
        // Table header
        const thead = document.createElement('thead');
        thead.innerHTML = `
            <tr>
                <th title="Unique identifier for each network node">Node ID</th>
                <th title="Activity level from 0-100. Higher = more active, lower = less active">
                    Activity Level
                    <span class="info-icon" onclick="showInfoPopup('nodeActivity')" style="margin-left: 5px; cursor: pointer; color: #3b82f6; font-weight: bold;">ℹ️</span>
                </th>
                <th title="Classification based on activity level: High (>2), Medium (0.5-2), Low (<0.5)">Status</th>
            </tr>
        `;
        table.appendChild(thead);

        // Table body
        const tbody = document.createElement('tbody');
        
        // Filter to show only nodes with meaningful activity changes (not baseline ~1.0 activity)
        // In 0-1 scale: baseline is ~0.01 (1%), significant changes are >1.5% or <0.5%
        const filteredActivities = Object.entries(nodeActivities)
            .filter(([nodeId, activity]) => {
                const activityLevel = activity * 100;
                return activityLevel < 0.5 || activityLevel > 1.5; // Show downregulated (<0.5%) or upregulated (>1.5%)
            });
        
        // Sort by activity level (highest first)
        const sortedActivities = filteredActivities.sort(([,a], [,b]) => b - a);
        
        sortedActivities.forEach(([nodeId, activity]) => {
            const row = document.createElement('tr');
            
            // Convert from 0-1 scale to 0-100 activity level
            const activityLevel = activity * 100;
            const displayValue = activityLevel.toFixed(1);
            const status = this.getActivityStatus(activity);
            
            // Bar width should be clamped to 100% max for display
            const barWidth = Math.min(activityLevel, 100);
            
            row.innerHTML = `
                <td>
                    <span class="node-id" title="${nodeId}">
                        ${this.truncateNodeId(nodeId)}
                    </span>
                </td>
                <td>
                    <div class="activity-bar">
                        <div class="activity-fill ${status}" style="width: ${barWidth}%"></div>
                        <span class="activity-text">${displayValue}</span>
                    </div>
                </td>
                <td><span class="status-badge ${status}">${status}</span></td>
            `;
            
            tbody.appendChild(row);
        });
        
        table.appendChild(tbody);
        
        // Create search box
        const searchContainer = document.createElement('div');
        searchContainer.className = 'table-search-container';
        searchContainer.innerHTML = `
            <input type="text" 
                   class="table-search-input" 
                   id="node-activities-search"
                   placeholder="Search by node ID, activity level, or status..."
                   style="width: 100%; padding: 8px 12px; margin-bottom: 10px; border: 1px solid var(--border-color); border-radius: 4px; font-size: 0.875rem;">
        `;
        
        // Create scrollable wrapper for the table
        const tableWrapper = document.createElement('div');
        tableWrapper.className = 'table-wrapper';
        table.id = 'node-activities-table';
        
        tableWrapper.appendChild(table);
        this.nodeActivitiesContainer.appendChild(searchContainer);
        this.nodeActivitiesContainer.appendChild(tableWrapper);
        
        // Add search functionality
        const searchInput = document.getElementById('node-activities-search');
        searchInput.addEventListener('input', (e) => {
            this.filterTable('node-activities-table', e.target.value);
        });
        
        // Add styles if not already present
        this.addTableStyles();
    }
    
    displayInfluenceScores(influenceScores) {
        this.influenceScoresContainer.innerHTML = '';
        
        if (Object.keys(influenceScores).length === 0) {
            this.influenceScoresContainer.innerHTML = '<p class="no-data">No influence data available</p>';
            return;
        }
        
        // Create table
        const table = document.createElement('table');
        table.className = 'results-table';
        
        // Table header
        const thead = document.createElement('thead');
        thead.innerHTML = `
            <tr>
                <th title="Unique identifier for each network node">Node ID</th>
                <th title="Mathematical measure of how much this node influences the entire network">
                    Influence Score
                    <span class="info-icon" onclick="showInfoPopup('influenceScore')" style="margin-left: 5px; cursor: pointer; color: #3b82f6; font-weight: bold;">ℹ️</span>
                </th>
                <th title="Ranking by influence (1st = most influential)">Rank</th>
            </tr>
        `;
        table.appendChild(thead);
        
        // Table body
        const tbody = document.createElement('tbody');
        
        // Sort by influence score (highest first)
        const sortedInfluence = Object.entries(influenceScores)
            .sort(([,a], [,b]) => b - a);
        
        sortedInfluence.forEach(([nodeId, score], index) => {
            const row = document.createElement('tr');
            
            row.innerHTML = `
                <td>
                    <span class="node-id" title="${nodeId}">
                        ${this.truncateNodeId(nodeId)}
                    </span>
                </td>
                <td>
                    <span class="influence-score">${score.toFixed(3)}</span>
                </td>
                <td>
                    <span class="rank-badge rank-${this.getRankClass(index)}">
                        ${index + 1}
                    </span>
                </td>
            `;
            
            tbody.appendChild(row);
        });
        
        table.appendChild(tbody);
        
        // Create search box
        const searchContainer = document.createElement('div');
        searchContainer.className = 'table-search-container';
        searchContainer.innerHTML = `
            <input type="text" 
                   class="table-search-input" 
                   id="influence-scores-search"
                   placeholder="Search by node ID or influence score..."
                   style="width: 100%; padding: 8px 12px; margin-bottom: 10px; border: 1px solid var(--border-color); border-radius: 4px; font-size: 0.875rem;">
        `;
        
        // Create scrollable wrapper for the table
        const tableWrapper = document.createElement('div');
        tableWrapper.className = 'table-wrapper';
        table.id = 'influence-scores-table';
        
        tableWrapper.appendChild(table);
        this.influenceScoresContainer.appendChild(searchContainer);
        this.influenceScoresContainer.appendChild(tableWrapper);
        
        // Add search functionality
        const searchInput = document.getElementById('influence-scores-search');
        searchInput.addEventListener('input', (e) => {
            this.filterTable('influence-scores-table', e.target.value);
        });
    }
    
    displayPathwayStats(results) {
        this.pathwayStatsContainer.innerHTML = '';
        
        const stats = this.calculateStats(results);
        
        stats.forEach(stat => {
            const statElement = document.createElement('div');
            statElement.className = 'stat-item';
            statElement.innerHTML = `
                <span class="stat-label">${stat.label}</span>
                <span class="stat-value">${stat.value}</span>
            `;
            this.pathwayStatsContainer.appendChild(statElement);
        });
    }
    
    calculateStats(results) {
        const stats = [];
        
        // Node count
        const nodeCount = results.node_activities ? Object.keys(results.node_activities).length : 0;
        stats.push({ label: 'Total Nodes', value: nodeCount });
        
        // Active nodes (>1.5x baseline)
        if (results.node_activities) {
            const activeNodes = Object.values(results.node_activities).filter(activity => activity > 0.015).length;
            stats.push({ label: 'Active Nodes', value: activeNodes });
        }
        
        // Observed nodes
        if (results.observations) {
            const observedCount = Object.keys(results.observations).length;
            stats.push({ label: 'Observed', value: observedCount });
        }
        
        // Solution quality (if available)
        if (results.solution_quality) {
            const quality = (results.solution_quality * 100).toFixed(1);
            stats.push({ label: 'Solution Quality', value: `${quality}%` });
        }
        
        // Convergence status
        if (results.converged !== undefined) {
            stats.push({ label: 'Converged', value: results.converged ? 'Yes' : 'No' });
        }
        
        // Analysis time
        if (results.analysis_time) {
            const time = results.analysis_time < 1 ? 
                `${Math.round(results.analysis_time * 1000)}ms` : 
                `${results.analysis_time.toFixed(2)}s`;
            stats.push({ label: 'Analysis Time', value: time });
        }
        
        return stats;
    }
    
    getActivityStatus(activity) {
        // Based on 0-1 scale where 0.01 = baseline (1%)
        const activityPercent = activity * 100;
        if (activityPercent >= 2.0) return 'high';     // >= 2% (upregulated)
        if (activityPercent >= 0.5) return 'medium';   // 0.5-2% (near baseline)
        return 'low';                                   // < 0.5% (downregulated)
    }
    
    getRankClass(index) {
        if (index === 0) return 'first';
        if (index === 1) return 'second';
        if (index === 2) return 'third';
        return 'other';
    }
    
    truncateNodeId(nodeId) {
        if (nodeId && nodeId.length > 12) {
            return nodeId.substring(0, 12) + '...';
        }
        return nodeId || '';
    }
    
    addTableStyles() {
        // Add styles for results tables if not already present
        if (document.querySelector('#results-table-styles')) {
            return;
        }
        
        const style = document.createElement('style');
        style.id = 'results-table-styles';
        style.textContent = `
            .table-wrapper {
                max-height: 600px;
                overflow-y: auto;
                border: 1px solid var(--border-color);
                border-radius: 6px;
                background: white;
            }
            
            .table-wrapper::-webkit-scrollbar {
                width: 8px;
            }
            
            .table-wrapper::-webkit-scrollbar-track {
                background: #f1f1f1;
                border-radius: 4px;
            }
            
            .table-wrapper::-webkit-scrollbar-thumb {
                background: #c1c1c1;
                border-radius: 4px;
            }
            
            .table-wrapper::-webkit-scrollbar-thumb:hover {
                background: #a8a8a8;
            }
            
            .results-table {
                width: 100%;
                border-collapse: collapse;
                font-size: 0.875rem;
                background: white;
                margin: 0;
            }
            
            .results-table th {
                background: var(--bg-tertiary);
                padding: 0.75rem 0.5rem;
                text-align: left;
                font-weight: 600;
                border-bottom: 2px solid var(--border-color);
                color: var(--text-primary);
                position: sticky;
                top: 0;
                z-index: 10;
            }
            
            .results-table td {
                padding: 0.625rem 0.5rem;
                border-bottom: 1px solid var(--border-color);
                vertical-align: middle;
            }
            
            .results-table tr:nth-child(even) {
                background: #f8fafc;
            }
            
            .results-table tr:hover {
                background: #e8f4fd;
            }
            
            .node-id {
                font-family: 'Consolas', monospace;
                font-size: 0.8rem;
                color: var(--text-secondary);
                cursor: help;
            }
            
            .activity-bar {
                position: relative;
                width: 80px;
                height: 20px;
                background: var(--bg-tertiary);
                border-radius: 4px;
                overflow: hidden;
            }
            
            .activity-fill {
                height: 100%;
                border-radius: 4px;
                transition: width 0.3s ease;
            }
            
            .activity-fill.high {
                background: linear-gradient(90deg, #10b981, #059669);
            }
            
            .activity-fill.medium {
                background: linear-gradient(90deg, #fbbf24, #f59e0b);
            }
            
            .activity-fill.low {
                background: linear-gradient(90deg, #ef4444, #dc2626);
            }
            
            .activity-text {
                position: absolute;
                top: 50%;
                left: 50%;
                transform: translate(-50%, -50%);
                font-size: 0.75rem;
                font-weight: 600;
                color: var(--text-primary);
                text-shadow: 0 0 2px rgba(255,255,255,0.8);
            }
            
            .status-badge {
                padding: 0.25rem 0.5rem;
                border-radius: 4px;
                font-size: 0.75rem;
                font-weight: 600;
                text-transform: uppercase;
                letter-spacing: 0.025em;
            }
            
            .status-badge.high {
                background: #dcfce7;
                color: #166534;
            }
            
            .status-badge.medium {
                background: #fef3c7;
                color: #92400e;
            }
            
            .status-badge.low {
                background: #fee2e2;
                color: #991b1b;
            }
            
            .influence-score {
                font-family: 'Consolas', monospace;
                font-weight: 600;
                color: var(--primary-color);
            }
            
            .rank-badge {
                display: inline-flex;
                align-items: center;
                justify-content: center;
                width: 24px;
                height: 24px;
                border-radius: 50%;
                font-size: 0.75rem;
                font-weight: bold;
                color: white;
            }
            
            .rank-badge.rank-first {
                background: linear-gradient(135deg, #ffd700, #ffed4e);
                color: #92400e;
            }
            
            .rank-badge.rank-second {
                background: linear-gradient(135deg, #c0c0c0, #e5e5e5);
                color: #374151;
            }
            
            .rank-badge.rank-third {
                background: linear-gradient(135deg, #cd7f32, #d97706);
                color: white;
            }
            
            .rank-badge.rank-other {
                background: var(--secondary-color);
            }
            
            .no-data {
                text-align: center;
                color: var(--text-muted);
                font-style: italic;
                padding: 2rem;
            }
        `;
        
        document.head.appendChild(style);
    }
    
    // Method to export results as CSV
    exportToCSV(results) {
        const csvData = [];
        
        // Add header
        csvData.push(['Node ID', 'Fold Change', 'Influence Score', 'Status']);
        
        // Add data rows
        const nodeActivities = results.node_activities || {};
        const influenceScores = results.influence_scores || {};
        
        Object.keys(nodeActivities).forEach(nodeId => {
            const activity = nodeActivities[nodeId] * 100; // Convert to fold change
            const influence = influenceScores[nodeId] || 0;
            const status = this.getActivityStatus(nodeActivities[nodeId]);
            
            csvData.push([nodeId, activity.toFixed(2), influence.toFixed(6), status]);
        });
        
        // Convert to CSV string
        const csvContent = csvData.map(row => row.join(',')).join('\n');
        
        // Download
        const blob = new Blob([csvContent], { type: 'text/csv' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `deltasignal-results-${Date.now()}.csv`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }
    
    // Filter table rows based on search query
    filterTable(tableId, query) {
        const table = document.getElementById(tableId);
        if (!table) return;
        
        const rows = table.querySelectorAll('tbody tr');
        const searchQuery = query.toLowerCase().trim();
        
        rows.forEach(row => {
            if (searchQuery === '') {
                row.style.display = '';
                return;
            }
            
            // Get all cell text content
            const cells = Array.from(row.querySelectorAll('td'));
            const rowText = cells.map(cell => {
                // For cells with complex HTML, get just the text content
                return cell.textContent || cell.innerText || '';
            }).join(' ').toLowerCase();
            
            // Show row if search query matches any cell content
            if (rowText.includes(searchQuery)) {
                row.style.display = '';
            } else {
                row.style.display = 'none';
            }
        });
        
        // Update row count or show "no results" message
        const visibleRows = Array.from(rows).filter(row => row.style.display !== 'none');
        
        // Add/update search results indicator
        let indicator = table.parentElement.querySelector('.search-results-indicator');
        if (!indicator) {
            indicator = document.createElement('div');
            indicator.className = 'search-results-indicator';
            indicator.style.cssText = `
                font-size: 0.75rem;
                color: var(--text-secondary);
                margin-top: 8px;
                text-align: right;
            `;
            table.parentElement.appendChild(indicator);
        }
        
        if (searchQuery === '') {
            indicator.textContent = '';
        } else {
            indicator.textContent = `${visibleRows.length} of ${rows.length} rows shown`;
        }
    }
    
    // Clear all displayed results
    clear() {
        this.nodeActivitiesContainer.innerHTML = '';
        this.influenceScoresContainer.innerHTML = '';
        this.pathwayStatsContainer.innerHTML = '';
    }
    
    // Show info popup for table explanations
    showInfoPopup(type) {
        const popups = {
            nodeActivity: {
                title: 'Node Activity Level',
                content: `
                    <p><strong>What it shows:</strong> How active each node is compared to baseline activity.</p>
                    <p><strong>Scale:</strong> 0-100 where:</p>
                    <ul>
                        <li><strong>1</strong> = Baseline activity (normal)</li>
                        <li><strong>&gt;2</strong> = Upregulated (more active)</li>
                        <li><strong>&lt;0.5</strong> = Downregulated (less active)</li>
                    </ul>
                    <p><strong>Color coding:</strong></p>
                    <ul>
                        <li><span style="color: #10b981; font-weight: bold;">Green</span> = High activity (&gt;2)</li>
                        <li><span style="color: #f59e0b; font-weight: bold;">Yellow</span> = Medium activity (0.5-2)</li>
                        <li><span style="color: #ef4444; font-weight: bold;">Red</span> = Low activity (&lt;0.5)</li>
                    </ul>
                `
            },
            influenceScore: {
                title: 'Influence Score',
                content: `
                    <p><strong>What it measures:</strong> How much each node influences the overall network behavior when perturbed.</p>
                    <p><strong>How it's calculated:</strong> Uses mathematical sensitivity analysis to measure how changes in one node ripple through the entire network.</p>
                    <p><strong>Interpretation:</strong></p>
                    <ul>
                        <li><strong>High scores</strong> = Key control points - small changes cause big network effects</li>
                        <li><strong>Low scores</strong> = Less influential - changes have minimal network impact</li>
                    </ul>
                    <p><strong>Use case:</strong> Target high-influence nodes for maximum therapeutic effect.</p>
                `
            }
        };
        
        const popup = popups[type];
        if (!popup) return;
        
        // Create modal backdrop
        const backdrop = document.createElement('div');
        backdrop.className = 'info-popup-backdrop';
        backdrop.style.cssText = `
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(0, 0, 0, 0.5);
            z-index: 1000;
            display: flex;
            align-items: center;
            justify-content: center;
        `;
        
        // Create modal
        const modal = document.createElement('div');
        modal.className = 'info-popup-modal';
        modal.style.cssText = `
            background: white;
            border-radius: 8px;
            padding: 1.5rem;
            max-width: 500px;
            width: 90%;
            box-shadow: 0 20px 25px -5px rgba(0, 0, 0, 0.1);
            position: relative;
        `;
        
        modal.innerHTML = `
            <button class="info-popup-close" style="
                position: absolute;
                top: 1rem;
                right: 1rem;
                background: none;
                border: none;
                font-size: 1.5rem;
                cursor: pointer;
                color: #6b7280;
                line-height: 1;
            ">&times;</button>
            <h3 style="margin: 0 0 1rem 0; color: #1f2937; font-size: 1.25rem;">${popup.title}</h3>
            <div style="color: #4b5563; line-height: 1.6;">${popup.content}</div>
        `;
        
        backdrop.appendChild(modal);
        document.body.appendChild(backdrop);
        
        // Close handlers
        const closePopup = () => {
            document.body.removeChild(backdrop);
        };
        
        backdrop.addEventListener('click', (e) => {
            if (e.target === backdrop) closePopup();
        });
        
        modal.querySelector('.info-popup-close').addEventListener('click', closePopup);
        
        // Close on escape key
        const escapeHandler = (e) => {
            if (e.key === 'Escape') {
                closePopup();
                document.removeEventListener('keydown', escapeHandler);
            }
        };
        document.addEventListener('keydown', escapeHandler);
    }
}

// Make showInfoPopup global so it can be called from onclick
window.showInfoPopup = function(type) {
    // Find the ResultsPanel instance and call its method
    // This is a bit of a hack, but necessary for onclick handlers in innerHTML
    const container = document.querySelector('#node-activities, #influence-scores');
    if (container && window.resultsPanel) {
        window.resultsPanel.showInfoPopup(type);
    }
};