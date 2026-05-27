// Interactive Network Visualizer with Node Editing
export default class NetworkVisualizer {
    constructor(containerId) {
        this.containerId = containerId;
        this.cy = null;
        this.currentLayout = 'dagre';
        this.showLabels = true;
        this.onNodeClick = null;
        this.networkData = null;
        this.perturbations = {}; // Store node perturbations
        this.baselineResults = null; // Store baseline results for comparison
        
        console.log('Interactive NetworkVisualizer created');
    }
    
    async init() {
        console.log('NetworkVisualizer init - creating placeholder');
        const container = document.getElementById(this.containerId);
        if (container) {
            container.innerHTML = `
                <div style="display: flex; align-items: center; justify-content: center; height: 100%; color: #666;">
                    <div style="text-align: center;">
                        <div style="font-size: 3rem; margin-bottom: 1rem;">🔗</div>
                        <p>Network visualization will appear here</p>
                        <p style="font-size: 0.8rem; margin-top: 0.5rem;">Click nodes to edit activity levels</p>
                    </div>
                </div>
            `;
        }
        console.log('Interactive NetworkVisualizer initialized');
    }
    
    async loadNetwork(networkData) {
        console.log('Loading network data:', networkData);
        this.networkData = networkData;
        
        const container = document.getElementById(this.containerId);
        if (container && networkData) {
            // Create interactive network display
            container.innerHTML = `
                <div class="network-display">
                    <div class="network-nodes">
                        <h4 style="margin: 0 0 1rem 0; color: #374151; font-size: 1rem;">Network Nodes</h4>
                        <div id="node-grid" class="node-grid">
                            ${this.createNodeGrid(networkData.nodes)}
                        </div>
                    </div>
                    <div class="network-edges" style="margin-top: 1.5rem;">
                        <h4 style="margin: 0 0 1rem 0; color: #374151; font-size: 1rem;">Connections</h4>
                        <div id="edge-summary" class="edge-summary">
                            ${this.createEdgeSummary(networkData.edges)}
                        </div>
                    </div>
                    <div class="network-stats" style="margin-top: 1.5rem; padding: 1rem; background: #f8fafc; border-radius: 8px;">
                        <h4 style="margin: 0 0 0.5rem 0; color: #374151; font-size: 0.9rem;">Network Statistics</h4>
                        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; font-size: 0.8rem;">
                            <div>Nodes: <strong>${networkData.nodes?.length || 0}</strong></div>
                            <div>Edges: <strong>${networkData.edges?.length || 0}</strong></div>
                        </div>
                    </div>
                </div>
            `;
            
            // Add click listeners to nodes
            this.setupNodeClickListeners();
        }
        
        return Promise.resolve();
    }
    
    createNodeGrid(nodes) {
        if (!nodes || nodes.length === 0) {
            return '<p style="color: #94a3b8; font-style: italic;">No nodes available</p>';
        }
        
        return nodes.map(node => {
            const activity = this.perturbations[node.uuid] || 0;
            const isPerturbed = node.uuid in this.perturbations;
            
            return `
                <div class="network-node ${isPerturbed ? 'perturbed' : ''}" 
                     data-node-id="${node.uuid}" 
                     title="${node.name || node.uuid}"
                     style="cursor: pointer;">
                    <div class="node-header">
                        <span class="node-name">${this.truncateText(node.name || node.uuid, 12)}</span>
                        ${isPerturbed ? '<span class="perturbation-indicator">●</span>' : ''}
                    </div>
                    <div class="node-activity">
                        <div class="activity-bar">
                            <div class="activity-fill" style="width: ${activity * 100}%; background: ${this.getActivityColor(activity)};"></div>
                        </div>
                        <span class="activity-text">${(activity * 100).toFixed(0)}</span>
                    </div>
                    <div class="node-type">${node.entity_type || 'unknown'}</div>
                </div>
            `;
        }).join('');
    }
    
    createEdgeSummary(edges) {
        if (!edges || edges.length === 0) {
            return '<p style="color: #94a3b8; font-style: italic;">No connections available</p>';
        }
        
        const positiveCount = edges.filter(e => e.is_positive).length;
        const negativeCount = edges.filter(e => !e.is_positive).length;
        
        // Create edge list
        const edgeList = edges.slice(0, 10).map(edge => {
            const parentName = this.getNodeName(edge.parent_uuid);
            const childName = this.getNodeName(edge.child_uuid);
            const arrow = edge.is_positive ? '→' : '⊣';
            const color = edge.is_positive ? '#10b981' : '#ef4444';
            
            return `
                <div style="display: flex; align-items: center; gap: 0.5rem; padding: 0.25rem 0; border-bottom: 1px solid #f1f5f9;">
                    <span style="font-size: 0.75rem; color: #6b7280;">${this.truncateText(parentName, 8)}</span>
                    <span style="color: ${color}; font-weight: bold;">${arrow}</span>
                    <span style="font-size: 0.75rem; color: #6b7280;">${this.truncateText(childName, 8)}</span>
                </div>
            `;
        }).join('');
        
        return `
            <div style="margin-bottom: 1rem;">
                <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; font-size: 0.8rem; margin-bottom: 1rem;">
                    <div style="display: flex; align-items: center; gap: 0.5rem;">
                        <div style="width: 12px; height: 3px; background: #10b981; border-radius: 2px;"></div>
                        <span>Activation: <strong>${positiveCount}</strong></span>
                    </div>
                    <div style="display: flex; align-items: center; gap: 0.5rem;">
                        <div style="width: 12px; height: 3px; background: #ef4444; border-radius: 2px; border: 1px dashed #ef4444;"></div>
                        <span>Inhibition: <strong>${negativeCount}</strong></span>
                    </div>
                </div>
            </div>
            <div style="max-height: 200px; overflow-y: auto; border: 1px solid #e5e7eb; border-radius: 6px; padding: 0.5rem;">
                <h5 style="margin: 0 0 0.5rem 0; font-size: 0.75rem; color: #374151; font-weight: 600;">Connections:</h5>
                ${edgeList}
                ${edges.length > 10 ? `<div style="text-align: center; padding: 0.5rem; color: #6b7280; font-size: 0.75rem;">...and ${edges.length - 10} more</div>` : ''}
            </div>
        `;
    }
    
    setupNodeClickListeners() {
        const nodeElements = document.querySelectorAll('.network-node');
        nodeElements.forEach(nodeElement => {
            nodeElement.addEventListener('click', (e) => {
                const nodeId = nodeElement.dataset.nodeId;
                const node = this.networkData.nodes.find(n => n.uuid === nodeId);
                
                if (node && this.onNodeClick) {
                    this.onNodeClick(node);
                }
                
                // Visual feedback
                nodeElements.forEach(el => el.classList.remove('selected'));
                nodeElement.classList.add('selected');
            });
        });
    }
    
    async updateWithResults(results) {
        console.log('Updating with results:', results);
        this.baselineResults = results;
        
        // Update node display with activity levels
        if (results.node_activities) {
            Object.entries(results.node_activities).forEach(([nodeId, activity]) => {
                this.updateNodeDisplay(nodeId, activity, false);
            });
        }
        
        const container = document.getElementById(this.containerId);
        if (container) {
            // Add results indicator
            const existingIndicator = container.querySelector('.results-indicator');
            if (existingIndicator) {
                existingIndicator.remove();
            }
            
            const indicator = document.createElement('div');
            indicator.className = 'results-indicator';
            indicator.style.cssText = 'position: absolute; top: 10px; right: 10px; background: #10b981; color: white; padding: 0.5rem 0.75rem; border-radius: 20px; font-size: 0.75rem; font-weight: 600;';
            indicator.textContent = `✓ Analysis Complete (${Object.keys(results.node_activities || {}).length} nodes)`;
            
            container.style.position = 'relative';
            container.appendChild(indicator);
        }
    }
    
    updateWithPerturbationResults(results, perturbations) {
        console.log('Updating with perturbation results:', results);
        
        // Update node display with new activity levels and highlight changes
        if (results.node_activities && this.baselineResults?.node_activities) {
            Object.entries(results.node_activities).forEach(([nodeId, newActivity]) => {
                const baselineActivity = this.baselineResults.node_activities[nodeId] || 0;
                const change = newActivity - baselineActivity;
                
                this.updateNodeDisplay(nodeId, newActivity, true, change);
            });
        }
        
        // Show perturbation summary
        this.showPerturbationSummary(perturbations, results);
    }
    
    updateNodeDisplay(nodeId, activity, showChange = false, change = 0) {
        const nodeElement = document.querySelector(`[data-node-id="${nodeId}"]`);
        if (nodeElement) {
            const activityBar = nodeElement.querySelector('.activity-fill');
            const activityText = nodeElement.querySelector('.activity-text');
            
            if (activityBar) {
                activityBar.style.width = `${activity * 100}%`;
                activityBar.style.background = this.getActivityColor(activity);
            }
            
            if (activityText) {
                let text = `${(activity * 100).toFixed(0)}`;
                if (showChange && Math.abs(change) > 0.01) {
                    const changeAmount = (change * 100).toFixed(0);
                    const changeSign = change > 0 ? '+' : '';
                    text += ` (${changeSign}${changeAmount})`;
                    
                    // Add visual change indicator
                    nodeElement.style.borderLeft = `4px solid ${change > 0 ? '#10b981' : '#ef4444'}`;
                }
                activityText.textContent = text;
            }
        }
    }
    
    showPerturbationSummary(perturbations, results) {
        const container = document.getElementById(this.containerId);
        if (!container) return;
        
        // Remove existing summary
        const existingSummary = container.querySelector('.perturbation-summary');
        if (existingSummary) {
            existingSummary.remove();
        }
        
        // Create perturbation summary
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
        `;
        
        const perturbationCount = Object.keys(perturbations).length;
        const totalNodes = results.node_activities ? Object.keys(results.node_activities).length : 0;
        
        summary.innerHTML = `
            <h4 style="margin: 0 0 0.5rem 0; color: #374151; font-size: 0.9rem;">🎛️ Perturbation Applied</h4>
            <p style="margin: 0 0 0.5rem 0; color: #6b7280;">${perturbationCount} nodes modified</p>
            <p style="margin: 0; color: #6b7280;">Impact on ${totalNodes} total nodes</p>
        `;
        
        container.appendChild(summary);
    }
    
    setPerturbation(nodeId, activity) {
        this.perturbations[nodeId] = activity;
        this.updateNodeDisplay(nodeId, activity);
        
        // Update visual state
        const nodeElement = document.querySelector(`[data-node-id="${nodeId}"]`);
        if (nodeElement) {
            nodeElement.classList.add('perturbed');
        }
        
        console.log('Set perturbation:', { nodeId, activity });
    }
    
    clearPerturbations() {
        this.perturbations = {};
        
        // Reset all node visuals
        const nodeElements = document.querySelectorAll('.network-node');
        nodeElements.forEach(nodeElement => {
            nodeElement.classList.remove('perturbed');
            nodeElement.style.borderLeft = '';
        });
        
        // If we have baseline results, restore them
        if (this.baselineResults?.node_activities) {
            Object.entries(this.baselineResults.node_activities).forEach(([nodeId, activity]) => {
                this.updateNodeDisplay(nodeId, activity);
            });
        }
        
        console.log('Cleared all perturbations');
    }
    
    getPerturbations() {
        return { ...this.perturbations };
    }
    
    getActivityColor(activity) {
        if (activity > 0.7) return '#10b981'; // High - green
        if (activity > 0.3) return '#f59e0b'; // Medium - yellow
        return '#ef4444'; // Low - red
    }
    
    getNodeName(nodeId) {
        if (!this.networkData?.nodes) return nodeId;
        const node = this.networkData.nodes.find(n => n.uuid === nodeId);
        return node ? (node.name || node.uuid) : nodeId;
    }
    
    truncateText(text, maxLength) {
        if (text && text.length > maxLength) {
            return text.substring(0, maxLength - 3) + '...';
        }
        return text || '';
    }
    
    changeLayout(layoutName) {
        this.currentLayout = layoutName;
        console.log('Layout changed to:', layoutName);
    }
    
    fit() {
        console.log('Fit called');
    }
    
    toggleLabels(show) {
        this.showLabels = show;
        console.log('Labels toggled:', show);
    }
    
    filterByPathway(pathway) {
        console.log('Filter by pathway:', pathway);
    }
    
    highlightNodes(nodeIds) {
        console.log('Highlight nodes:', nodeIds);
    }
    
    getNetworkStats() {
        return {
            nodeCount: this.networkData?.nodes?.length || 0,
            edgeCount: this.networkData?.edges?.length || 0,
            perturbationCount: Object.keys(this.perturbations).length,
            activationEdges: this.networkData?.edges?.filter(e => e.is_positive).length || 0,
            inhibitionEdges: this.networkData?.edges?.filter(e => !e.is_positive).length || 0
        };
    }
}