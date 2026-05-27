// Cytoscape.js Network Visualizer with Interactive Nodes
export default class NetworkVisualizer {
    constructor(containerId) {
        this.containerId = containerId;
        this.cy = null;
        this.currentLayout = 'dagre';
        this.showLabels = true;
        this.onNodeClick = null;
        this.networkData = null;
        this.perturbations = {}; 
        this.baselineResults = null; 
        
        console.log('Cytoscape NetworkVisualizer created');
    }
    
    async init() {
        console.log('NetworkVisualizer init - creating Cytoscape container');
        const container = document.getElementById(this.containerId);
        if (container) {
            container.innerHTML = `
                <div style="display: flex; align-items: center; justify-content: center; height: 100%; color: #666;">
                    <div style="text-align: center;">
                        <div style="font-size: 3rem; margin-bottom: 1rem;">=</div>
                        <p>Network visualization will appear here</p>
                        <p style="font-size: 0.8rem; margin-top: 0.5rem;">Click nodes to edit activity levels</p>
                    </div>
                </div>
            `;
        }
        console.log('NetworkVisualizer initialized');
    }
    
    async loadNetwork(networkData) {
        console.log('Loading network data into Cytoscape:', networkData);
        this.networkData = networkData;
        this.currentView = 'network';
        
        if (!networkData || !networkData.nodes) {
            console.error('No network data provided');
            return;
        }
        
        // Show network view by default
        this.showNetworkView();
        return Promise.resolve();
    }
    
    showNetworkView() {
        const container = document.getElementById(this.containerId);
        if (!container) {
            console.error('Container not found');
            return;
        }
        
        // Clear container and add structure
        container.innerHTML = `
            <div class="network-view-wrapper" style="position: relative; width: 100%; height: 100%;">
                <div id="cy-container" style="width: 100%; height: 100%;"></div>
                <div class="activity-legend" style="position: absolute; bottom: 20px; right: 20px; background: white; padding: 12px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); font-size: 12px; max-width: 250px;">
                    <div style="display: flex; align-items: center; gap: 8px; margin-bottom: 8px;">
                        <div style="font-weight: bold;">Network Legend</div>
                        <button id="legend-info-btn" style="background: #3b82f6; color: white; border: none; width: 18px; height: 18px; border-radius: 50%; font-size: 11px; cursor: pointer; display: flex; align-items: center; justify-content: center;" title="Click for detailed explanation">?</button>
                    </div>
                    
                    <div style="font-weight: 600; margin-bottom: 4px; color: #374151;">Node Colors (Activity):</div>
                    <div style="display: flex; flex-direction: column; gap: 3px; margin-bottom: 8px;">
                        <div style="display: flex; align-items: center; gap: 6px;">
                            <div style="width: 16px; height: 16px; background: #10b981; border-radius: 3px;"></div>
                            <span style="font-size: 11px;">Upregulated (>2x)</span>
                        </div>
                        <div style="display: flex; align-items: center; gap: 6px;">
                            <div style="width: 16px; height: 16px; background: #86efac; border-radius: 3px;"></div>
                            <span style="font-size: 11px;">Slightly up (1.5-2x)</span>
                        </div>
                        <div style="display: flex; align-items: center; gap: 6px;">
                            <div style="width: 16px; height: 16px; background: #6b7280; border-radius: 3px;"></div>
                            <span style="font-size: 11px;">Normal (0.75-1.5x)</span>
                        </div>
                        <div style="display: flex; align-items: center; gap: 6px;">
                            <div style="width: 16px; height: 16px; background: #fca5a5; border-radius: 3px;"></div>
                            <span style="font-size: 11px;">Slightly down (0.5-0.75x)</span>
                        </div>
                        <div style="display: flex; align-items: center; gap: 6px;">
                            <div style="width: 16px; height: 16px; background: #ef4444; border-radius: 3px;"></div>
                            <span style="font-size: 11px;">Downregulated (<0.5x)</span>
                        </div>
                    </div>
                    
                    <div style="font-weight: 600; margin-bottom: 4px; color: #374151;">Node Borders:</div>
                    <div style="display: flex; flex-direction: column; gap: 3px;">
                        <div style="display: flex; align-items: center; gap: 6px;">
                            <div style="width: 16px; height: 16px; background: #f9fafb; border: 2px dashed #f59e0b; border-radius: 3px;"></div>
                            <span style="font-size: 11px;">Perturbed (your input)</span>
                        </div>
                        <div style="display: flex; align-items: center; gap: 6px;">
                            <div style="width: 16px; height: 16px; background: #f9fafb; border: 2px solid #9ca3af; border-radius: 3px;"></div>
                            <span style="font-size: 11px;">Calculated (system response)</span>
                        </div>
                    </div>
                </div>
            </div>
        `;
        
        const cyContainer = container.querySelector('#cy-container');
        
        // Set up legend info button
        const legendBtn = container.querySelector('#legend-info-btn');
        if (legendBtn) {
            legendBtn.addEventListener('click', () => this.showDetailedLegend());
        }
        
        // Create Cytoscape elements
        const elements = this.createCytoscapeElements(this.networkData);
        
        // Initialize Cytoscape with fallback layout
        const layoutConfig = this.getLayoutConfig();
        
        this.cy = cytoscape({
            container: cyContainer,
            elements: elements,
            style: this.getCytoscapeStyle(),
            layout: layoutConfig,
            userZoomingEnabled: true,
            userPanningEnabled: true,
            boxSelectionEnabled: false,
            selectionType: 'single',
            touchTapThreshold: 8,
            desktopTapThreshold: 4
        });
        
        // Add click handler for nodes
        this.cy.on('tap', 'node', (evt) => {
            const node = evt.target;
            const nodeData = node.data();
            console.log('Node clicked:', nodeData);
            
            if (this.onNodeClick) {
                const fullNodeData = this.networkData.nodes.find(n => n.uuid === nodeData.id);
                if (fullNodeData) {
                    this.onNodeClick(fullNodeData);
                }
            }
        });
        
        // Add hover effects
        this.cy.on('mouseover', 'node', (evt) => {
            const node = evt.target;
            node.addClass('highlighted');
        });
        
        this.cy.on('mouseout', 'node', (evt) => {
            const node = evt.target;
            node.removeClass('highlighted');
        });
        
        console.log('Network view loaded with', elements.length, 'elements');
    }
    
    showNodesListView() {
        const container = document.getElementById(this.containerId);
        if (!container || !this.networkData?.nodes) return;
        
        // Clear cytoscape instance
        if (this.cy) {
            this.cy.destroy();
            this.cy = null;
        }
        
        // Create nodes list HTML
        const nodesHtml = this.networkData.nodes.map(node => {
            const activity = this.perturbations[node.uuid] || 
                            (this.baselineResults?.node_activities?.[node.uuid] * 100) || 1; // Default to 1x (neutral)
            const isPerturbed = node.uuid in this.perturbations;
            
            return `
                <div class="node-list-item ${isPerturbed ? 'perturbed' : ''}" data-node-id="${node.uuid}">
                    <div class="node-info">
                        <div class="node-name">${node.name || node.uuid}</div>
                        <div class="node-type">${node.entity_type || 'unknown'}</div>
                        <div class="node-id">${node.uuid}</div>
                    </div>
                    <div class="node-activity">
                        <div class="activity-value">${Math.round(activity)}</div>
                        <div class="activity-bar">
                            <div class="activity-fill" style="width: ${activity}%; background: ${this.getActivityColor(activity/100)};"></div>
                        </div>
                        ${isPerturbed ? '<div class="perturbation-marker">●</div>' : ''}
                    </div>
                </div>
            `;
        }).join('');
        
        container.innerHTML = `
            <div class="nodes-list-view">
                <div class="list-header">
                    <h3>Network Nodes (${this.networkData.nodes.length})</h3>
                    <div class="list-stats">
                        <span>Perturbed: ${Object.keys(this.perturbations).length}</span>
                    </div>
                </div>
                <div class="nodes-list">
                    ${nodesHtml}
                </div>
            </div>
        `;
        
        // Add click listeners for nodes
        container.querySelectorAll('.node-list-item').forEach(item => {
            item.addEventListener('click', (e) => {
                const nodeId = item.dataset.nodeId;
                const node = this.networkData.nodes.find(n => n.uuid === nodeId);
                if (node && this.onNodeClick) {
                    this.onNodeClick(node);
                }
            });
        });
        
        console.log('Nodes list view loaded');
    }
    
    showEdgesListView() {
        const container = document.getElementById(this.containerId);
        if (!container || !this.networkData?.edges) return;
        
        // Clear cytoscape instance
        if (this.cy) {
            this.cy.destroy();
            this.cy = null;
        }
        
        // Create edges list HTML
        const edgesHtml = this.networkData.edges.map((edge, index) => {
            const parentNode = this.networkData.nodes.find(n => n.uuid === edge.parent_uuid);
            const childNode = this.networkData.nodes.find(n => n.uuid === edge.child_uuid);
            const parentName = parentNode?.name || edge.parent_uuid;
            const childName = childNode?.name || edge.child_uuid;
            
            return `
                <div class="edge-list-item ${edge.is_positive ? 'activation' : 'inhibition'}">
                    <div class="edge-info">
                        <div class="edge-connection">
                            <span class="source-node">${this.truncateText(parentName, 15)}</span>
                            <span class="edge-arrow ${edge.is_positive ? 'positive' : 'negative'}">
                                ${edge.is_positive ? '→' : '⊣'}
                            </span>
                            <span class="target-node">${this.truncateText(childName, 15)}</span>
                        </div>
                        <div class="edge-type">${edge.is_positive ? 'Activation' : 'Inhibition'}</div>
                        <div class="edge-details">
                            <span>AND: ${edge.is_and ? 'Yes' : 'No'}</span>
                            <span>Weight: ${edge.stoichiometry || 1}</span>
                        </div>
                    </div>
                </div>
            `;
        }).join('');
        
        const activationCount = this.networkData.edges.filter(e => e.is_positive).length;
        const inhibitionCount = this.networkData.edges.filter(e => !e.is_positive).length;
        
        container.innerHTML = `
            <div class="edges-list-view">
                <div class="list-header">
                    <h3>Network Edges (${this.networkData.edges.length})</h3>
                    <div class="list-stats">
                        <span class="activation">Activation: ${activationCount}</span>
                        <span class="inhibition">Inhibition: ${inhibitionCount}</span>
                    </div>
                </div>
                <div class="edges-list">
                    ${edgesHtml}
                </div>
            </div>
        `;
        
        console.log('Edges list view loaded');
    }
    
    switchView(viewType) {
        this.currentView = viewType;
        
        switch(viewType) {
            case 'network':
                this.showNetworkView();
                break;
            case 'nodes':
                this.showNodesListView();
                break;
            case 'edges':
                this.showEdgesListView();
                break;
            case 'pathway':
                this.showPathwayView();
                break;
            default:
                this.showNetworkView();
        }
    }
    
    showPathwayView() {
        const container = document.getElementById(this.containerId);
        if (!container) return;
        
        // Placeholder for pathway view - could integrate with PathwayBrowser
        container.innerHTML = `
            <div class="pathway-view">
                <div class="placeholder">
                    <div class="placeholder-icon">🧬</div>
                    <p>Pathway visualization coming soon</p>
                </div>
            </div>
        `;
    }
    
    getLayoutConfig() {
        // Try dagre first, fallback to breadthfirst if dagre is not available
        try {
            if (typeof dagre !== 'undefined' && this.cy && this.cy.layout) {
                return {
                    name: 'dagre',
                    directed: true,
                    padding: 20,
                    spacingFactor: 1.2,
                    nodeDimensionsIncludeLabels: true,
                    rankDir: 'TB'
                };
            }
        } catch (e) {
            console.warn('Dagre layout not available, using breadthfirst');
        }
        
        // Fallback to breadthfirst layout
        return {
            name: 'breadthfirst',
            directed: true,
            padding: 20,
            spacingFactor: 1.5,
            avoidOverlap: true,
            nodeDimensionsIncludeLabels: true
        };
    }
    
    createCytoscapeElements(networkData) {
        const elements = [];
        
        // Add nodes
        if (networkData.nodes) {
            networkData.nodes.forEach(node => {
                const activity = this.perturbations[node.uuid] || 0.01; // Default to 1x baseline (neutral)
                
                elements.push({
                    data: {
                        id: node.uuid,
                        label: node.name || node.uuid.substring(0, 8),
                        name: node.name || node.uuid,
                        entity_type: node.entity_type,
                        activity: activity,
                        isPerturbed: false  // Explicitly initialize as false
                    }
                });
            });
        }
        
        // Add edges
        if (networkData.edges) {
            networkData.edges.forEach((edge, index) => {
                elements.push({
                    data: {
                        id: `edge-${index}`,
                        source: edge.parent_uuid,
                        target: edge.child_uuid,
                        is_positive: edge.is_positive.toString(),
                        interaction_type: edge.is_positive ? 'activation' : 'inhibition'
                    }
                });
            });
        }
        
        return elements;
    }
    
    getCytoscapeStyle() {
        return [
            // Node styles
            {
                selector: 'node',
                style: {
                    'width': 60,
                    'height': 60,
                    'background-color': (ele) => this.getActivityColor(ele.data('activity')),
                    'border-width': 2,
                    'border-color': '#9ca3af',
                    'border-style': 'solid',
                    'label': 'data(label)',
                    'text-valign': 'center',
                    'text-halign': 'center',
                    'font-size': 10,
                    'font-weight': 'bold',
                    'color': '#1e293b',
                    'text-outline-width': 2,
                    'text-outline-color': '#ffffff',
                    'shape': 'ellipse',
                    'cursor': 'pointer'
                }
            },
            // Perturbed node style
            {
                selector: 'node[isPerturbed]',
                style: {
                    'border-width': 4,
                    'border-color': '#f59e0b',
                    'border-style': 'dashed'
                }
            },
            // Highlighted node style
            {
                selector: 'node.highlighted',
                style: {
                    'border-width': 4,
                    'border-color': '#3b82f6',
                    'transform': 'scale(1.1)'
                }
            },
            // Edge styles
            {
                selector: 'edge',
                style: {
                    'width': 3,
                    'line-color': '#6b7280',
                    'target-arrow-color': '#6b7280',
                    'target-arrow-shape': 'triangle',
                    'curve-style': 'bezier',
                    'arrow-scale': 1.2
                }
            },
            // Activation edges
            {
                selector: 'edge[is_positive = "true"]',
                style: {
                    'line-color': '#6b7280',
                    'target-arrow-color': '#6b7280',
                    'line-style': 'solid',
                    'target-arrow-shape': 'triangle'
                }
            },
            // Inhibition edges  
            {
                selector: 'edge[is_positive = "false"]',
                style: {
                    'line-color': '#6b7280',
                    'target-arrow-color': '#6b7280',
                    'line-style': 'solid',
                    'target-arrow-shape': 'tee'
                }
            }
        ];
    }
    
    async updateWithResults(results) {
        console.log('Updating Cytoscape with results:', results);
        this.baselineResults = results;
        
        if (!this.cy || !results.node_activities) return;
        
        // Update node colors and labels based on activity
        Object.entries(results.node_activities).forEach(([nodeId, activity]) => {
            const node = this.cy.getElementById(nodeId);
            if (node.length > 0) {
                node.data('activity', activity);
                
                // Always update background color to reflect calculated activity
                node.style('background-color', this.getActivityColor(activity));
                
                // Only update border style for non-perturbed nodes
                // Perturbed nodes should keep their dashed yellow borders
                const isPerturbed = node.data('isPerturbed');
                if (!isPerturbed) {
                    node.style('border-style', 'solid');
                    node.style('border-color', '#9ca3af');
                    node.style('border-width', 2);
                }
                
                // Update label - just show the node name for clarity
                const originalLabel = node.data('name') || node.data('label').split('\n')[0];
                node.data('label', originalLabel);
            }
        });
        
        // Add results indicator
        const container = document.getElementById(this.containerId);
        if (container) {
            const existingIndicator = container.querySelector('.results-indicator');
            if (existingIndicator) {
                existingIndicator.remove();
            }
            
            const indicator = document.createElement('div');
            indicator.className = 'results-indicator';
            indicator.style.cssText = 'position: absolute; top: 10px; right: 10px; background: #10b981; color: white; padding: 0.5rem 0.75rem; border-radius: 20px; font-size: 0.75rem; font-weight: 600; z-index: 1000;';
            indicator.textContent = ` Analysis Complete (${Object.keys(results.node_activities || {}).length} nodes)`;
            
            container.style.position = 'relative';
            container.appendChild(indicator);
        }
    }
    
    updateWithPerturbationResults(results, perturbations) {
        console.log('Updating Cytoscape with perturbation results:', results);
        
        if (!this.cy || !results.node_activities) return;
        
        // Update node activities and highlight changes
        if (results.node_activities) {
            Object.entries(results.node_activities).forEach(([nodeId, newActivity]) => {
                const node = this.cy.getElementById(nodeId);
                if (node.length > 0) {
                    // Use baseline if available, otherwise assume 0.5 (normal) as baseline  
                    const baselineActivity = this.baselineResults?.node_activities?.[nodeId] ?? 0.5;
                    
                    node.data('activity', newActivity);
                    node.style('background-color', this.getActivityColor(newActivity));
                    
                    // Update label - keep it simple, just the node name
                    const originalLabel = node.data('name') || node.data('label').split('\n')[0];
                    node.data('label', originalLabel);
                    
                    // Only update border style for non-perturbed nodes
                    // Perturbed nodes should keep their dashed borders
                    const isPerturbed = node.data('isPerturbed');
                    if (!isPerturbed) {
                        node.style('border-style', 'solid');
                        node.style('border-color', '#9ca3af');
                        node.style('border-width', 2);
                    }
                }
            });
        }
        
        this.showPerturbationSummary(perturbations, results);
    }
    
    showPerturbationSummary(perturbations, results) {
        const container = document.getElementById(this.containerId);
        if (!container) return;
        
        const existingSummary = container.querySelector('.perturbation-summary');
        if (existingSummary) {
            existingSummary.remove();
        }
        
        const summary = document.createElement('div');
        summary.className = 'perturbation-summary';
        summary.style.cssText = `
            position: absolute; 
            top: 10px; 
            left: 10px; 
            background: #fff; 
            border: 1px solid #e5e7eb; 
            border-radius: 8px; 
            padding: 1rem; 
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
            max-width: 300px;
            font-size: 0.8rem;
            z-index: 1000;
        `;
        
        const perturbationCount = Object.keys(perturbations).length;
        const totalNodes = results.node_activities ? Object.keys(results.node_activities).length : 0;
        
        summary.innerHTML = `
            <h4 style="margin: 0 0 0.5rem 0; color: #374151; font-size: 0.9rem;"><h4 style="margin: 0 0 0.5rem 0; color: #374151; font-size: 0.9rem;">⚡ Perturbation Applied</h4>
            <p style="margin: 0 0 0.5rem 0; color: #6b7280;">${perturbationCount} nodes modified</p>
            <p style="margin: 0; color: #6b7280;">Impact on ${totalNodes} total nodes</p>
        `;
        
        container.appendChild(summary);
    }
    
    setPerturbation(nodeId, activity) {
        this.perturbations[nodeId] = activity;
        
        if (this.cy) {
            const node = this.cy.getElementById(nodeId);
            if (node.length > 0) {
                node.data('activity', activity);
                node.data('isPerturbed', true);
                node.style('background-color', this.getActivityColor(activity));
                node.style('border-color', '#f59e0b');
                node.style('border-width', 4);
                node.style('border-style', 'dashed');
                
                // Update label - just the node name, color indicates perturbation
                const originalLabel = node.data('name') || node.data('label').split('\n')[0];
                node.data('label', originalLabel);
            }
        }
        
        console.log('Set perturbation:', { nodeId, activity });
    }
    
    clearPerturbations() {
        this.perturbations = {};
        
        if (this.cy) {
            this.cy.nodes().forEach(node => {
                node.data('isPerturbed', false);
                node.style('border-color', '#9ca3af');
                node.style('border-width', 2);
                node.style('border-style', 'solid');
                
                // Restore baseline activity if available
                const nodeId = node.data('id');
                if (this.baselineResults?.node_activities && this.baselineResults.node_activities[nodeId] !== undefined) {
                    const baselineActivity = this.baselineResults.node_activities[nodeId];
                    node.data('activity', baselineActivity);
                    node.style('background-color', this.getActivityColor(baselineActivity));
                    
                    const originalLabel = node.data('name') || node.data('label').split('\n')[0];
                    node.data('label', originalLabel);
                } else {
                    node.data('activity', 0.01); // Default to 1x baseline (neutral)
                    node.style('background-color', this.getActivityColor(0.01));
                    const originalLabel = node.data('name') || node.data('label').split('\n')[0];
                    node.data('label', originalLabel);
                }
            });
        }
        
        console.log('Cleared all perturbations');
    }
    
    getPerturbations() {
        return { ...this.perturbations };
    }
    
    getActivityColor(activity) {
        // Activity is 0-1 where 0.01 = 1x (normal fold change)
        // Convert to fold change for better thresholds
        const foldChange = activity * 100;
        
        if (foldChange > 2) return '#10b981';      // >2x - green (upregulated)
        if (foldChange < 0.5) return '#ef4444';    // <0.5x - red (downregulated)
        if (foldChange > 1.5) return '#86efac';    // 1.5-2x - light green
        if (foldChange < 0.75) return '#fca5a5';   // 0.5-0.75x - light red
        return '#6b7280';                          // 0.75-1.5x - gray (near normal)
    }
    
    changeLayout(layoutName) {
        this.currentLayout = layoutName;
        if (this.cy) {
            let layoutOptions;
            
            if (layoutName === 'dagre' && typeof dagre !== 'undefined') {
                layoutOptions = {
                    name: 'dagre',
                    directed: true,
                    padding: 20,
                    spacingFactor: 1.2,
                    nodeDimensionsIncludeLabels: true,
                    rankDir: 'TB'
                };
            } else {
                // Fallback layout
                layoutOptions = {
                    name: 'breadthfirst',
                    directed: true,
                    padding: 20,
                    spacingFactor: 1.5,
                    avoidOverlap: true,
                    nodeDimensionsIncludeLabels: true
                };
            }
            
            this.cy.layout(layoutOptions).run();
        }
        console.log('Layout changed to:', layoutName);
    }
    
    fit() {
        if (this.cy) {
            this.cy.fit();
        }
        console.log('Fit called');
    }
    
    toggleLabels(show) {
        this.showLabels = show;
        if (this.cy) {
            if (show) {
                this.cy.style().selector('node').style('label', 'data(label)').update();
            } else {
                this.cy.style().selector('node').style('label', '').update();
            }
        }
        console.log('Labels toggled:', show);
    }
    
    filterByPathway(pathway) {
        console.log('Filter by pathway:', pathway);
    }
    
    highlightNodes(nodeIds) {
        if (this.cy) {
            // Reset all nodes
            this.cy.nodes().removeClass('highlighted');
            
            // Highlight specified nodes
            nodeIds.forEach(nodeId => {
                const node = this.cy.getElementById(nodeId);
                if (node.length > 0) {
                    node.addClass('highlighted');
                }
            });
        }
        console.log('Highlight nodes:', nodeIds);
    }
    
    getNetworkStats() {
        const stats = {
            nodeCount: this.networkData?.nodes?.length || 0,
            edgeCount: this.networkData?.edges?.length || 0,
            perturbationCount: Object.keys(this.perturbations).length,
            activationEdges: this.networkData?.edges?.filter(e => e.is_positive).length || 0,
            inhibitionEdges: this.networkData?.edges?.filter(e => !e.is_positive).length || 0
        };
        
        if (this.cy) {
            stats.visibleNodes = this.cy.nodes(':visible').length;
            stats.visibleEdges = this.cy.edges(':visible').length;
        }
        
        return stats;
    }
    
    showDetailedLegend() {
        const modal = document.createElement('div');
        modal.style.cssText = `
            position: fixed; top: 0; left: 0; width: 100%; height: 100%;
            background: rgba(0,0,0,0.5); z-index: 10000; display: flex; 
            align-items: center; justify-content: center;
        `;
        
        modal.innerHTML = `
            <div style="background: white; padding: 2rem; border-radius: 12px; max-width: 500px; max-height: 80vh; overflow-y: auto;">
                <div style="display: flex; justify-content: between; align-items: center; margin-bottom: 1.5rem;">
                    <h3 style="margin: 0; color: #1f2937;">Network Visualization Guide</h3>
                    <button id="close-legend" style="background: none; border: none; font-size: 1.5rem; cursor: pointer; color: #6b7280;">×</button>
                </div>
                
                <div style="margin-bottom: 1.5rem;">
                    <h4 style="color: #374151; margin-bottom: 0.5rem;">How to Use the Network:</h4>
                    <ol style="padding-left: 1.5rem; line-height: 1.6;">
                        <li><strong>Click on any node</strong> to set its activity level (perturbation)</li>
                        <li><strong>Click "Solve Network"</strong> to calculate downstream effects</li>
                        <li><strong>Use Node Editor</strong> to manage multiple perturbations</li>
                        <li><strong>View results</strong> in the Analysis Results panel below</li>
                    </ol>
                </div>
                
                <div style="margin-bottom: 1.5rem;">
                    <h4 style="color: #374151; margin-bottom: 0.5rem;">Node Colors (Activity Levels):</h4>
                    <div style="display: grid; grid-template-columns: auto 1fr; gap: 0.5rem 1rem; align-items: center;">
                        <div style="width: 20px; height: 20px; background: #ef4444; border-radius: 4px;"></div>
                        <span><strong>Red:</strong> Downregulated (activity reduced below normal)</span>
                        <div style="width: 20px; height: 20px; background: #fca5a5; border-radius: 4px;"></div>
                        <span><strong>Light Red:</strong> Slightly downregulated</span>
                        <div style="width: 20px; height: 20px; background: #6b7280; border-radius: 4px;"></div>
                        <span><strong>Gray:</strong> Normal activity (baseline, ~1x)</span>
                        <div style="width: 20px; height: 20px; background: #86efac; border-radius: 4px;"></div>
                        <span><strong>Light Green:</strong> Slightly upregulated</span>
                        <div style="width: 20px; height: 20px; background: #10b981; border-radius: 4px;"></div>
                        <span><strong>Green:</strong> Upregulated (activity increased above normal)</span>
                    </div>
                </div>
                
                <div style="margin-bottom: 1.5rem;">
                    <h4 style="color: #374151; margin-bottom: 0.5rem;">Node Borders:</h4>
                    <div style="display: grid; grid-template-columns: auto 1fr; gap: 0.5rem 1rem; align-items: center;">
                        <div style="width: 20px; height: 20px; background: #f9fafb; border: 2px dashed #f59e0b; border-radius: 4px;"></div>
                        <span><strong>Dashed Yellow:</strong> Perturbed nodes (your experimental input)</span>
                        <div style="width: 20px; height: 20px; background: #f9fafb; border: 2px solid #9ca3af; border-radius: 4px;"></div>
                        <span><strong>Solid Gray:</strong> Calculated nodes (system's computed response)</span>
                    </div>
                </div>
                
                <div>
                    <h4 style="color: #374151; margin-bottom: 0.5rem;">Activity Scale:</h4>
                    <p style="line-height: 1.6; margin: 0;">
                        Activity levels represent fold-change from normal expression:
                        <br>• <strong>0x</strong> = completely inactive
                        <br>• <strong>1x</strong> = normal/baseline activity  
                        <br>• <strong>2x</strong> = double normal activity
                        <br>• <strong>10x</strong> = ten times normal activity
                    </p>
                </div>
            </div>
        `;
        
        // Close modal when clicking outside or on close button
        modal.addEventListener('click', (e) => {
            if (e.target === modal || e.target.id === 'close-legend') {
                document.body.removeChild(modal);
            }
        });
        
        document.body.appendChild(modal);
    }
}