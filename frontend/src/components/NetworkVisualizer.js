// Network Visualizer using Cytoscape.js
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
        
        // Initialize immediately
        this.init();
    }
    
    async init() {
        console.log('NetworkVisualizer init - creating Cytoscape container');
        
        // Import Cytoscape and extensions
        const cytoscape = (await import('cytoscape')).default;
        const dagre = (await import('cytoscape-dagre')).default;
        const cose = (await import('cytoscape-cose-bilkent')).default;
        
        // Register extensions
        cytoscape.use(dagre);
        cytoscape.use(cose);
        
        const container = document.getElementById(this.containerId);
        if (container) {
            container.innerHTML = `
                <div style="display: flex; align-items: center; justify-content: center; height: 100%; color: #666;">
                    <div style="text-align: center;">
                        <div style="font-size: 3rem; margin-bottom: 1rem;">🧬</div>
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

        // Detect network type
        const isReactionNetwork = networkData.network_type === 'reaction' ||
                                (networkData.entities && networkData.reactions);

        if (isReactionNetwork) {
            if (!networkData.entities && !networkData.reactions) {
                console.error('No reaction network data provided');
                return;
            }
        } else {
            if (!networkData.nodes) {
                console.error('No network data provided');
                return;
            }
        }

        // Store network type
        this.networkType = isReactionNetwork ? 'reaction' : 'regulatory';

        // Show network view by default
        await this.showNetworkView();
        return Promise.resolve();
    }
    
    async showNetworkView() {
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
        
        // Add legend info button listener
        const legendInfoBtn = document.getElementById('legend-info-btn');
        if (legendInfoBtn) {
            legendInfoBtn.addEventListener('click', () => {
                if (window.resultsPanel && window.resultsPanel.showNetworkLegendInfo) {
                    window.resultsPanel.showNetworkLegendInfo();
                }
            });
        }
        
        // Wait for DOM to be ready then create Cytoscape
        setTimeout(async () => {
            await this.createCytoscapeInstance();
        }, 100);
    }
    
    async createCytoscapeInstance() {
        // Import Cytoscape and extensions if not already imported
        const cytoscape = (await import('cytoscape')).default;
        const dagre = (await import('cytoscape-dagre')).default;
        const cose = (await import('cytoscape-cose-bilkent')).default;
        
        // Register extensions
        cytoscape.use(dagre);
        cytoscape.use(cose);
        
        // Initialize Cytoscape instance
        this.cy = cytoscape({
            container: document.getElementById('cy-container'),
            
            style: [
                {
                    selector: 'node',
                    style: {
                        'background-color': '#e2e8f0',
                        'border-color': '#64748b',
                        'border-width': 2,
                        'label': 'data(label)',
                        'text-valign': 'center',
                        'text-halign': 'center',
                        'font-size': '12px',
                        'font-family': 'Inter, sans-serif',
                        'font-weight': '500',
                        'color': '#1e293b',
                        'width': 'label',
                        'height': 'label',
                        'padding': '8px',
                        'shape': 'roundrectangle',
                        'text-wrap': 'wrap',
                        'text-max-width': '80px'
                    }
                },
                {
                    selector: 'node:selected',
                    style: {
                        'border-color': '#3b82f6',
                        'border-width': 3
                    }
                },
                {
                    selector: 'edge',
                    style: {
                        'width': 3,
                        'line-color': '#9ca3af',
                        'target-arrow-color': '#9ca3af',
                        'target-arrow-shape': 'triangle',
                        'curve-style': 'bezier'
                    }
                },
                {
                    selector: '.edge-positive',
                    style: {
                        'line-color': '#10b981',
                        'target-arrow-color': '#10b981'
                    }
                },
                {
                    selector: '.edge-negative',
                    style: {
                        'line-color': '#ef4444',
                        'target-arrow-color': '#ef4444',
                        'target-arrow-shape': 'tee'
                    }
                },
                {
                    selector: '.node-upregulated',
                    style: {
                        'background-color': '#10b981'
                    }
                },
                {
                    selector: '.node-downregulated',
                    style: {
                        'background-color': '#ef4444'
                    }
                },
                // Reaction network styles - Reaction nodes (diamonds)
                {
                    selector: '.reaction-node',
                    style: {
                        'shape': 'diamond',
                        'background-color': '#e2e8f0', // Default gray like other nodes
                        'border-color': '#64748b',
                        'border-width': 2,
                        'width': '40px',
                        'height': '40px',
                        'font-size': '10px',
                        'text-max-width': '60px',
                        'label': 'data(label)',
                        'text-valign': 'center',
                        'text-halign': 'center',
                        'font-family': 'Inter, sans-serif',
                        'font-weight': '500',
                        'color': '#1e293b',
                        'text-wrap': 'wrap'
                    }
                },
                // Entity nodes use same styling as regular nodes (roundrectangle)
                {
                    selector: '.entity-node',
                    style: {
                        'shape': 'roundrectangle', // Same as regulatory network
                        'background-color': '#e2e8f0',
                        'border-color': '#64748b',
                        'border-width': 2,
                        'label': 'data(label)',
                        'text-valign': 'center',
                        'text-halign': 'center',
                        'font-size': '12px',
                        'font-family': 'Inter, sans-serif',
                        'font-weight': '500',
                        'color': '#1e293b',
                        'width': 'label',
                        'height': 'label',
                        'padding': '8px',
                        'text-wrap': 'wrap',
                        'text-max-width': '80px'
                    }
                },
                // Activity-based colors (must come after entity/reaction styles to override them)
                {
                    selector: '.entity-node.node-upregulated, .reaction-node.node-upregulated, node.node-upregulated',
                    style: {
                        'background-color': '#10b981'
                    }
                },
                {
                    selector: '.entity-node.node-downregulated, .reaction-node.node-downregulated, node.node-downregulated',
                    style: {
                        'background-color': '#ef4444'
                    }
                },
                {
                    selector: '.edge-substrate',
                    style: {
                        'line-color': '#4ade80',
                        'target-arrow-color': '#4ade80',
                        'target-arrow-shape': 'triangle'
                    }
                },
                {
                    selector: '.edge-product',
                    style: {
                        'line-color': '#fb7185',
                        'target-arrow-color': '#fb7185',
                        'target-arrow-shape': 'triangle'
                    }
                },
                {
                    selector: '.edge-catalyst',
                    style: {
                        'line-color': '#a78bfa',
                        'target-arrow-color': '#a78bfa',
                        'target-arrow-shape': 'diamond',
                        'line-style': 'dashed'
                    }
                },
                {
                    selector: '.node-perturbed',
                    style: {
                        'border-style': 'dashed',
                        'border-color': '#f59e0b',
                        'border-width': 3
                    }
                }
            ],
            
            elements: this.networkType === 'reaction'
                ? this.convertReactionNetworkToCytoscape(this.networkData)
                : this.convertNetworkDataToCytoscape(this.networkData),
            
            layout: {
                name: 'dagre',
                rankDir: 'TB',
                nodeSep: 50,
                rankSep: 80
            },
            
            wheelSensitivity: 0.5,
            minZoom: 0.1,
            maxZoom: 3
        });
        
        // Add event listeners
        this.cy.on('tap', 'node', (evt) => {
            const node = evt.target;
            let nodeData = null;

            // Handle different network types
            if (this.networkType === 'reaction') {
                // Look in both entities and reactions
                if (this.networkData.entities) {
                    nodeData = this.networkData.entities.find(n => n.uuid === node.id());
                }
                if (!nodeData && this.networkData.reactions) {
                    nodeData = this.networkData.reactions.find(n => n.uuid === node.id());
                }
            } else {
                // Regular regulatory network
                if (this.networkData.nodes) {
                    nodeData = this.networkData.nodes.find(n => n.uuid === node.id());
                }
            }

            if (nodeData && this.onNodeClick) {
                this.onNodeClick(nodeData);
            }
        });
        
        // Sync zoom slider when user zooms with mouse wheel or other zoom events
        this.cy.on('zoom', () => {
            this.updateZoomSlider();
        });
        
        // Initialize zoom slider with current zoom level
        this.updateZoomSlider();

        // Initialize autocomplete functionality
        this.initializeAutocomplete();

        // Apply initial activity coloring for reaction networks
        if (this.networkType === 'reaction') {
            this.applyInitialActivityColoring();
        }

        console.log('Cytoscape instance created successfully');
    }
    
    convertNetworkDataToCytoscape(networkData) {
        if (!networkData) return { nodes: [], edges: [] };
        
        const elements = [];
        
        // Add nodes
        if (networkData.nodes) {
            networkData.nodes.forEach(node => {
                elements.push({
                    data: {
                        id: node.uuid,
                        label: node.name || node.uuid,
                        nodeData: node
                    }
                });
            });
        }
        
        // Add edges
        if (networkData.edges) {
            networkData.edges.forEach((edge, index) => {
                elements.push({
                    data: {
                        id: `edge_${index}`,
                        source: edge.parent_uuid,
                        target: edge.child_uuid,
                        edgeData: edge
                    },
                    classes: edge.is_positive ? 'edge-positive' : 'edge-negative'
                });
            });
        }
        
        return elements;
    }

    convertReactionNetworkToCytoscape(reactionNetwork) {
        if (!reactionNetwork) return { nodes: [], edges: [] };

        const elements = [];

        // Add entity nodes (circular)
        if (reactionNetwork.entities) {
            reactionNetwork.entities.forEach(entity => {
                const classes = ['entity-node', `entity-${entity.type}`];

                elements.push({
                    data: {
                        id: entity.uuid,
                        label: entity.name,
                        name: entity.name,
                        uuid: entity.uuid,
                        type: 'entity',
                        entityType: entity.type,
                        compartment: entity.compartment,
                        activity: entity.activity || 0,
                        nodeData: entity
                    },
                    classes: classes.join(' ')
                });
            });
        }

        // Add reaction nodes (diamond)
        if (reactionNetwork.reactions) {
            reactionNetwork.reactions.forEach(reaction => {
                const classes = ['reaction-node', `reaction-${reaction.type.replace('_', '')}`];

                elements.push({
                    data: {
                        id: reaction.uuid,
                        label: reaction.name,
                        name: reaction.name,
                        uuid: reaction.uuid,
                        type: 'reaction',
                        reactionType: reaction.type,
                        ecNumber: reaction.ec_number,
                        compartment: reaction.compartment,
                        activity: reaction.activity || 0,
                        nodeData: reaction
                    },
                    classes: classes.join(' ')
                });
            });
        }

        // Add connections
        if (reactionNetwork.connections) {
            reactionNetwork.connections.forEach((connection, index) => {
                const classes = [`edge-${connection.role}`];

                elements.push({
                    data: {
                        id: `edge_${index}`,
                        source: connection.source,
                        target: connection.target,
                        role: connection.role,
                        stoichiometry: connection.stoichiometry || 1,
                        edgeData: connection
                    },
                    classes: classes.join(' ')
                });
            });
        }

        return elements;
    }

    updateWithResults(results) {
        if (!this.cy || !results.node_activities) {
            console.error('updateWithResults: Missing cy or node_activities', {
                cy: !!this.cy,
                node_activities: !!results.node_activities
            });
            return;
        }

        console.log('=== updateWithResults DEBUG ===');
        console.log('Results structure:', Object.keys(results));
        console.log('Node activities keys:', Object.keys(results.node_activities));
        console.log('Sample activities:', Object.entries(results.node_activities).slice(0, 5));

        // Get all node IDs currently in the visualization
        const visualizationNodeIds = [];
        this.cy.nodes().forEach(node => {
            visualizationNodeIds.push(node.id());
        });
        console.log('Visualization node IDs:', visualizationNodeIds);

        // Check for matches between result keys and visualization node IDs
        const resultKeys = Object.keys(results.node_activities);
        const matches = resultKeys.filter(key => visualizationNodeIds.includes(key));
        const mismatches = resultKeys.filter(key => !visualizationNodeIds.includes(key));

        console.log('Matching node IDs:', matches);
        console.log('Non-matching result keys:', mismatches);


        // Update node activities and styles
        this.cy.nodes().forEach(node => {
            const nodeId = node.id();
            const activity = results.node_activities[nodeId];

            if (activity !== undefined) {
                // Convert from 0-1 scale to display scale
                const displayActivity = activity * 100;

                console.log(`Updating node ${nodeId}: activity ${activity} -> display ${displayActivity}`);

                // Apply color based on activity level
                let nodeClass = '';
                if (displayActivity > 150) {
                    nodeClass = 'node-upregulated';
                } else if (displayActivity < 75) {
                    nodeClass = 'node-downregulated';
                }

                // Remove old classes and add new ones
                node.removeClass('node-upregulated node-downregulated');
                if (nodeClass) {
                    node.addClass(nodeClass);
                }

                // Check if node is perturbed
                if (this.perturbations && this.perturbations[nodeId]) {
                    node.addClass('node-perturbed');
                }
            } else {
                console.log(`No activity found for node ${nodeId}`);
            }
        });

        console.log('=== updateWithResults COMPLETE ===');
    }

    simulateReactionNetworkResults() {
        console.log('Simulating reaction network results based on perturbations and initial activities');

        this.cy.nodes().forEach(node => {
            const nodeId = node.id();
            let activity = 50; // Default activity

            // Use perturbation activity if available
            if (this.perturbations && this.perturbations[nodeId]) {
                activity = this.perturbations[nodeId] * 100; // Convert to percentage
            } else {
                // Use initial activity from network data if available
                if (this.networkData.entities) {
                    const entity = this.networkData.entities.find(e => e.uuid === nodeId);
                    if (entity && entity.activity !== undefined) {
                        activity = entity.activity;
                    }
                }
                if (this.networkData.reactions) {
                    const reaction = this.networkData.reactions.find(r => r.uuid === nodeId);
                    if (reaction && reaction.activity !== undefined) {
                        activity = reaction.activity;
                    }
                }

                // Add some random variation for simulation
                activity = activity + (Math.random() - 0.5) * 20;
                activity = Math.max(0, Math.min(100, activity)); // Clamp 0-100
            }

            console.log(`Simulating node ${nodeId}: activity ${activity}`);

            // Apply color based on activity level
            let nodeClass = '';
            if (activity > 90) {
                nodeClass = 'node-upregulated';
            } else if (activity < 60) {
                nodeClass = 'node-downregulated';
            }

            // Remove old classes and add new ones
            node.removeClass('node-upregulated node-downregulated');
            if (nodeClass) {
                node.addClass(nodeClass);
                console.log(`Applied class ${nodeClass} to node ${nodeId}`);
            }

            // Mark perturbed nodes
            if (this.perturbations && this.perturbations[nodeId]) {
                node.addClass('node-perturbed');
            }
        });

        console.log('Reaction network simulation complete');
    }

    // Method to manually apply colors based on perturbations (for reaction networks)
    updateColorsFromPerturbations() {
        if (!this.cy || !this.perturbations) return;

        console.log('Updating colors from perturbations:', this.perturbations);

        this.cy.nodes().forEach(node => {
            const nodeId = node.id();
            const perturbationActivity = this.perturbations[nodeId];

            // Remove old classes first
            node.removeClass('node-upregulated node-downregulated node-perturbed');

            if (perturbationActivity !== undefined) {
                // Convert from 0-1 scale to display scale
                const displayActivity = perturbationActivity * 100;

                console.log(`Applying perturbation color to node ${nodeId}: activity ${perturbationActivity} -> display ${displayActivity}`);

                // Apply color based on activity level
                let nodeClass = '';
                if (displayActivity > 2) { // 2x upregulation
                    nodeClass = 'node-upregulated';
                } else if (displayActivity < 0.5) { // 0.5x downregulation
                    nodeClass = 'node-downregulated';
                }

                if (nodeClass) {
                    node.addClass(nodeClass);
                }

                // Mark as perturbed
                node.addClass('node-perturbed');
            }
        });

        console.log('Colors updated from perturbations');
    }

    applyInitialActivityColoring() {
        if (!this.cy || this.networkType !== 'reaction') return;

        // Apply activity coloring based on initial activity values
        this.cy.nodes().forEach(node => {
            const nodeId = node.id();
            let activity = null;

            // Find activity in entities or reactions
            if (this.networkData.entities) {
                const entity = this.networkData.entities.find(e => e.uuid === nodeId);
                if (entity && entity.activity !== undefined) {
                    activity = entity.activity;
                }
            }
            if (activity === null && this.networkData.reactions) {
                const reaction = this.networkData.reactions.find(r => r.uuid === nodeId);
                if (reaction && reaction.activity !== undefined) {
                    activity = reaction.activity;
                }
            }

            if (activity !== null) {
                // Apply color based on activity level (activities are 0-100 in sample data)
                let nodeClass = '';
                if (activity > 95) {
                    nodeClass = 'node-upregulated';
                } else if (activity < 80) {
                    nodeClass = 'node-downregulated';
                }

                // Remove old classes and add new ones
                node.removeClass('node-upregulated node-downregulated');
                if (nodeClass) {
                    node.addClass(nodeClass);
                }
            }
        });
    }
    
    async changeLayout(layoutName) {
        this.currentLayout = layoutName;
        await this.applyLayout(layoutName);
    }
    
    async applyLayout(layoutName) {
        if (!this.cy) return;
        
        let layoutOptions;
        
        switch (layoutName) {
            case 'dagre':
                layoutOptions = {
                    name: 'dagre',
                    rankDir: 'TB',
                    nodeSep: 50,
                    rankSep: 80,
                    animate: true,
                    animationDuration: 500
                };
                break;
                
            case 'cose':
                layoutOptions = {
                    name: 'cose-bilkent',
                    idealEdgeLength: 100,
                    nodeOverlap: 20,
                    refresh: 20,
                    fit: true,
                    padding: 30,
                    randomize: false,
                    componentSpacing: 100,
                    nodeRepulsion: 400000,
                    edgeElasticity: 100,
                    nestingFactor: 5,
                    gravity: 80,
                    numIter: 1000,
                    initialTemp: 200,
                    coolingFactor: 0.95,
                    minTemp: 1.0,
                    animate: true,
                    animationDuration: 500
                };
                break;
                
            case 'grid':
                layoutOptions = {
                    name: 'grid',
                    fit: true,
                    padding: 30,
                    avoidOverlap: true,
                    animate: true,
                    animationDuration: 500
                };
                break;
                
            default:
                layoutOptions = { name: 'random' };
        }
        
        const layout = this.cy.layout(layoutOptions);
        layout.run();
        
        return new Promise((resolve) => {
            layout.on('layoutstop', () => {
                resolve();
            });
        });
    }
    
    fit() {
        if (this.cy) {
            this.cy.fit(null, 50);
        }
    }
    
    toggleLabels(show) {
        this.showLabels = show;
        
        if (show) {
            this.cy.style().selector('node').style('label', 'data(label)').update();
        } else {
            this.cy.style().selector('node').style('label', '').update();
        }
    }
    
    searchNodes(searchTerm) {
        if (!this.cy || !searchTerm) {
            this.clearSearch();
            return;
        }
        
        const term = searchTerm.toLowerCase();
        const matchingNodes = [];
        
        // Search through all nodes
        this.cy.nodes().forEach(node => {
            const nodeData = node.data();
            const label = (nodeData.label || nodeData.id || '').toLowerCase();
            const id = (nodeData.id || '').toLowerCase();
            const name = (nodeData.name || '').toLowerCase();
            
            // Check if search term matches ID, label, or name
            const matches = label.includes(term) || id.includes(term) || name.includes(term);
            
            if (matches) {
                matchingNodes.push(node);
            }
        });
        
        if (matchingNodes.length > 0) {
            // Reset all nodes to normal style
            this.cy.nodes().removeClass('search-highlight search-match');
            
            // Highlight matching nodes
            matchingNodes.forEach(node => {
                node.addClass('search-match');
            });
            
            // Update styles for search results
            this.cy.style()
                .selector('.search-match')
                .style({
                    'border-width': '3px',
                    'border-color': '#ef4444',
                    'border-style': 'solid'
                })
                .update();
            
            // Fit to show the first match
            this.cy.fit(matchingNodes[0], 100);
            
            console.log(`Found ${matchingNodes.length} node(s) matching "${searchTerm}"`);
        } else {
            console.log(`No nodes found matching "${searchTerm}"`);
            this.clearSearch();
        }
    }
    
    clearSearch() {
        if (!this.cy) return;

        // Remove search highlighting
        this.cy.nodes().removeClass('search-highlight search-match');

        // Reset search-specific styles
        this.cy.style()
            .selector('.search-match')
            .style({
                'border-width': '2px',
                'border-color': 'data(borderColor)',
                'border-style': 'data(borderStyle)'
            })
            .update();

        // Hide autocomplete
        this.hideAutocomplete();
    }

    // Autocomplete functionality
    initializeAutocomplete() {
        this.autocompleteContainer = document.getElementById('search-autocomplete');
        this.searchInput = document.getElementById('node-search');
        this.selectedAutocompleteIndex = -1;
        this.autocompleteItems = [];

        if (!this.autocompleteContainer || !this.searchInput) return;

        // Add event listeners for autocomplete
        this.searchInput.addEventListener('input', (e) => {
            this.handleAutocompleteInput(e.target.value);
        });

        this.searchInput.addEventListener('keydown', (e) => {
            this.handleAutocompleteKeydown(e);
        });

        this.searchInput.addEventListener('blur', () => {
            // Delay hiding to allow for item clicks
            setTimeout(() => this.hideAutocomplete(), 150);
        });

        this.searchInput.addEventListener('focus', () => {
            if (this.searchInput.value.trim()) {
                this.handleAutocompleteInput(this.searchInput.value);
            }
        });

        // Handle window resize to reposition dropdown
        window.addEventListener('resize', () => {
            if (this.autocompleteContainer && this.autocompleteContainer.style.display === 'block') {
                this.positionAutocomplete();
            }
        });
    }

    getSearchableNodes() {
        if (!this.cy) return [];

        const searchableNodes = [];
        this.cy.nodes().forEach(node => {
            const nodeData = node.data();
            searchableNodes.push({
                id: nodeData.id,
                uuid: nodeData.uuid,
                name: nodeData.name || nodeData.label || nodeData.id,
                label: nodeData.label,
                node: node
            });
        });

        return searchableNodes;
    }

    handleAutocompleteInput(searchTerm) {
        if (!searchTerm.trim()) {
            this.hideAutocomplete();
            this.clearSearch();
            return;
        }

        const term = searchTerm.toLowerCase();
        const searchableNodes = this.getSearchableNodes();

        // Filter nodes that match the search term
        const matches = searchableNodes.filter(nodeData => {
            const name = (nodeData.name || '').toLowerCase();
            const id = (nodeData.id || '').toLowerCase();
            const uuid = (nodeData.uuid || '').toLowerCase();
            const label = (nodeData.label || '').toLowerCase();

            return name.includes(term) || id.includes(term) || uuid.includes(term) || label.includes(term);
        });

        // Limit to top 10 matches
        this.autocompleteItems = matches.slice(0, 10);
        this.selectedAutocompleteIndex = -1;

        this.showAutocomplete();
    }

    showAutocomplete() {
        if (!this.autocompleteContainer || this.autocompleteItems.length === 0) {
            this.hideAutocomplete();
            return;
        }

        let html = '';
        if (this.autocompleteItems.length === 0) {
            html = '<div class="search-autocomplete-empty">No nodes found</div>';
        } else {
            html = this.autocompleteItems.map((item, index) => {
                const isSelected = index === this.selectedAutocompleteIndex;
                return `
                    <div class="search-autocomplete-item ${isSelected ? 'selected' : ''}" data-index="${index}">
                        <div class="search-autocomplete-name">${this.escapeHtml(item.name)}</div>
                        <div class="search-autocomplete-id">${this.escapeHtml(item.uuid || item.id)}</div>
                    </div>
                `;
            }).join('');
        }

        this.autocompleteContainer.innerHTML = html;

        // Position the dropdown using fixed positioning
        this.positionAutocomplete();

        this.autocompleteContainer.style.display = 'block';

        // Add click event listeners to autocomplete items
        this.autocompleteContainer.querySelectorAll('.search-autocomplete-item').forEach(item => {
            item.addEventListener('mousedown', (e) => {
                e.preventDefault(); // Prevent blur event
                const index = parseInt(e.currentTarget.dataset.index);
                this.selectAutocompleteItem(index);
            });

            item.addEventListener('mouseenter', (e) => {
                // Update visual selection on hover
                this.updateAutocompleteSelection(parseInt(e.currentTarget.dataset.index));
            });
        });
    }

    positionAutocomplete() {
        if (!this.searchInput || !this.autocompleteContainer) return;

        const inputRect = this.searchInput.getBoundingClientRect();
        this.autocompleteContainer.style.top = (inputRect.bottom + window.scrollY) + 'px';
        this.autocompleteContainer.style.left = (inputRect.left + window.scrollX) + 'px';
        this.autocompleteContainer.style.width = inputRect.width + 'px';
    }

    hideAutocomplete() {
        if (this.autocompleteContainer) {
            this.autocompleteContainer.style.display = 'none';
            this.autocompleteContainer.innerHTML = '';
        }
        this.selectedAutocompleteIndex = -1;
        this.autocompleteItems = [];
    }

    handleAutocompleteKeydown(e) {
        if (!this.autocompleteContainer || this.autocompleteContainer.style.display === 'none') {
            // Handle escape to clear search even when autocomplete is hidden
            if (e.key === 'Escape') {
                e.preventDefault();
                this.searchInput.value = '';
                this.clearSearch();
                this.hideAutocomplete();
            }
            return;
        }

        switch (e.key) {
            case 'ArrowDown':
                e.preventDefault();
                this.selectedAutocompleteIndex = Math.min(
                    this.selectedAutocompleteIndex + 1,
                    this.autocompleteItems.length - 1
                );
                this.updateAutocompleteSelection(this.selectedAutocompleteIndex);
                break;

            case 'ArrowUp':
                e.preventDefault();
                this.selectedAutocompleteIndex = Math.max(this.selectedAutocompleteIndex - 1, -1);
                this.updateAutocompleteSelection(this.selectedAutocompleteIndex);
                break;

            case 'Enter':
                e.preventDefault();
                if (this.selectedAutocompleteIndex >= 0) {
                    this.selectAutocompleteItem(this.selectedAutocompleteIndex);
                } else if (this.autocompleteItems.length > 0) {
                    // If no item is selected, select the first one
                    this.selectAutocompleteItem(0);
                }
                break;

            case 'Escape':
                e.preventDefault();
                this.searchInput.value = '';
                this.hideAutocomplete();
                this.clearSearch();
                break;
        }
    }

    updateAutocompleteSelection(index) {
        this.selectedAutocompleteIndex = index;

        // Update visual selection
        this.autocompleteContainer.querySelectorAll('.search-autocomplete-item').forEach((item, i) => {
            if (i === index) {
                item.classList.add('selected');
            } else {
                item.classList.remove('selected');
            }
        });
    }

    selectAutocompleteItem(index) {
        if (index < 0 || index >= this.autocompleteItems.length) return;

        const selectedItem = this.autocompleteItems[index];
        this.searchInput.value = selectedItem.name;
        this.hideAutocomplete();

        // Perform search and focus on the selected node
        this.searchNodes(selectedItem.name);

        // Center the view on the selected node
        if (selectedItem.node) {
            this.cy.animate({
                center: { eles: selectedItem.node },
                zoom: Math.max(this.cy.zoom(), 1.5)
            }, {
                duration: 500
            });
        }
    }

    escapeHtml(unsafe) {
        if (!unsafe) return '';
        return unsafe
            .replace(/&/g, "&amp;")
            .replace(/</g, "&lt;")
            .replace(/>/g, "&gt;")
            .replace(/"/g, "&quot;")
            .replace(/'/g, "&#039;");
    }
    
    setZoom(zoomLevel) {
        if (!this.cy) return;
        
        // Set zoom level while keeping current pan position
        this.cy.zoom({
            level: zoomLevel,
            renderedPosition: {
                x: this.cy.width() / 2,
                y: this.cy.height() / 2
            }
        });
    }
    
    updateZoomSlider() {
        if (!this.cy) return;
        
        const currentZoom = this.cy.zoom();
        const slider = document.getElementById('zoom-slider');
        const label = document.getElementById('zoom-level');
        
        if (slider && label) {
            slider.value = currentZoom;
            label.textContent = Math.round(currentZoom * 100) + '%';
        }
    }
    
    filterByPathway(pathway) {
        if (!this.cy || !pathway) {
            return;
        }
        
        // Show all elements first
        this.cy.elements().show();
        
        if (pathway.nodeIds && pathway.nodeIds.length > 0) {
            // Hide elements not in pathway
            const pathwayNodeIds = new Set(pathway.nodeIds);
            
            this.cy.nodes().forEach(node => {
                if (!pathwayNodeIds.has(node.id())) {
                    node.hide();
                }
            });
            
            // Hide edges that don't connect visible nodes
            this.cy.edges().forEach(edge => {
                const source = edge.source();
                const target = edge.target();
                
                if (source.hidden() || target.hidden()) {
                    edge.hide();
                }
            });
            
            // Fit to visible elements
            this.cy.fit(this.cy.elements(':visible'), 50);
        }
    }
    
    highlightNodes(nodeIds) {
        if (!this.cy) return;
        
        // Clear existing highlights
        this.cy.elements().removeClass('highlighted');
        
        // Highlight specified nodes
        nodeIds.forEach(nodeId => {
            const node = this.cy.$(`#${nodeId}`);
            if (node.length > 0) {
                node.addClass('highlighted');
            }
        });
    }
    
    getNetworkStats() {
        if (!this.cy) return null;
        
        const nodes = this.cy.nodes();
        const edges = this.cy.edges();
        
        return {
            nodeCount: nodes.length,
            edgeCount: edges.length,
            activationEdges: edges.filter('.edge-positive').length,
            inhibitionEdges: edges.filter('.edge-negative').length,
            connectedComponents: this.cy.elements().components().length
        };
    }
    
    // Perturbation management methods
    getPerturbations() {
        return this.perturbations;
    }
    
    setPerturbation(nodeId, activity) {
        this.perturbations[nodeId] = activity;

        // Update visual styling if node exists
        if (this.cy) {
            const node = this.cy.$(`#${nodeId}`);
            if (node.length > 0) {
                // Remove old classes
                node.removeClass('node-upregulated node-downregulated node-perturbed');

                // Apply activity-based color
                const displayActivity = activity * 100;
                let nodeClass = '';
                if (displayActivity > 2) { // 2x upregulation
                    nodeClass = 'node-upregulated';
                } else if (displayActivity < 0.5) { // 0.5x downregulation
                    nodeClass = 'node-downregulated';
                }

                if (nodeClass) {
                    node.addClass(nodeClass);
                }

                // Mark as perturbed
                node.addClass('node-perturbed');
            }
        }

        console.log('Perturbation set:', nodeId, activity);
    }
    
    removePerturbation(nodeId) {
        delete this.perturbations[nodeId];

        // Remove visual styling if node exists
        if (this.cy) {
            const node = this.cy.$(`#${nodeId}`);
            if (node.length > 0) {
                node.removeClass('node-perturbed node-upregulated node-downregulated');

                // For reaction networks, reapply initial activity coloring
                if (this.networkType === 'reaction') {
                    let activity = null;

                    // Find activity in entities or reactions
                    if (this.networkData.entities) {
                        const entity = this.networkData.entities.find(e => e.uuid === nodeId);
                        if (entity && entity.activity !== undefined) {
                            activity = entity.activity;
                        }
                    }
                    if (activity === null && this.networkData.reactions) {
                        const reaction = this.networkData.reactions.find(r => r.uuid === nodeId);
                        if (reaction && reaction.activity !== undefined) {
                            activity = reaction.activity;
                        }
                    }

                    if (activity !== null) {
                        // Apply initial color based on activity level
                        let nodeClass = '';
                        if (activity > 95) {
                            nodeClass = 'node-upregulated';
                        } else if (activity < 80) {
                            nodeClass = 'node-downregulated';
                        }

                        if (nodeClass) {
                            node.addClass(nodeClass);
                        }
                    }
                }
            }
        }

        console.log('Perturbation removed:', nodeId);
    }
    
    clearPerturbations() {
        this.perturbations = {};
        
        // Remove all perturbation styling
        if (this.cy) {
            this.cy.nodes().removeClass('node-perturbed');
        }
        
        console.log('All perturbations cleared');
    }
    
    updateWithPerturbationResults(results, perturbations) {
        console.log('=== updateWithPerturbationResults DEBUG ===');
        console.log('Results:', results);
        console.log('Perturbations:', perturbations);

        if (!this.cy || !results.node_activities) {
            console.error('updateWithPerturbationResults: Missing cy or node_activities', {
                cy: !!this.cy,
                node_activities: !!results.node_activities
            });
            return;
        }
        
        // Debug: Log first few activity values to check the scale
        const entries = Object.entries(results.node_activities);
        console.log('Sample activity values:', entries.slice(0, 5).map(([id, val]) => `${id}: ${val}`));
        
        // Update node activities and highlight changes
        if (results.node_activities) {
            Object.entries(results.node_activities).forEach(([nodeId, newActivity]) => {
                const node = this.cy.$(`#${nodeId}`);
                if (node.length > 0) {
                    // Convert from 0-1 scale to display scale
                    const displayActivity = newActivity * 100;

                    console.log(`Updating perturbation result for node ${nodeId}: activity ${newActivity} -> display ${displayActivity}`);

                    // Apply color based on activity level (same as regular updateWithResults)
                    let nodeClass = '';
                    if (displayActivity > 150) {
                        nodeClass = 'node-upregulated';
                    } else if (displayActivity < 75) {
                        nodeClass = 'node-downregulated';
                    }

                    // Remove old activity classes and add new ones
                    node.removeClass('node-upregulated node-downregulated');
                    if (nodeClass) {
                        node.addClass(nodeClass);
                    }

                    // Keep perturbation styling for perturbed nodes
                    const isPerturbed = perturbations && perturbations[nodeId];
                    if (isPerturbed) {
                        node.addClass('node-perturbed');
                    }

                    // Update node data
                    node.data('activity', newActivity);

                    // Update label - keep it simple, just the node name
                    const originalLabel = node.data('label');
                    const nodeName = originalLabel ? originalLabel.split('\n')[0] : nodeId;
                    node.data('label', nodeName);
                }
            });
        }

        console.log('=== updateWithPerturbationResults COMPLETE ===');
    }
    
    getActivityColor(activity) {
        // Activity comes from backend as 0-1 scale
        // Convert to percentage to match table logic
        const activityPercent = activity * 100;
        
        // Use same thresholds as ResultsPanel for consistency
        if (activityPercent >= 2.0) return '#10b981';      // >= 2% - green (high/upregulated)
        if (activityPercent >= 0.5) return '#6b7280';      // 0.5-2% - gray (medium/near baseline)
        return '#ef4444';                                   // < 0.5% - red (low/downregulated)
    }
}