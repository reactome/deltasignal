// DeltaSignal Frontend Application
import NetworkVisualizer from './components/NetworkVisualizer.js';
import PathwayBrowser from './components/PathwayBrowser-simple.js';
import FileUploadHandler from './components/FileUploadHandler.js';
import ResultsPanel from './components/ResultsPanel.js';
import APIClient from './utils/APIClient-simple.js';

class DeltaSignalApp {
    constructor() {
        this.networkVisualizer = null;
        this.pathwayBrowser = null;
        this.fileUploadHandler = null;
        this.resultsPanel = null;
        this.apiClient = new APIClient();
        
        this.currentNetwork = null;
        this.currentResults = null;
        this.selectedNode = null;
        this.nodeEditor = null;
        this.savedScenarios = {};
        
        this.init();
    }
    
    init() {
        this.setupComponents();
        this.setupEventListeners();
        this.setupUI();
        
        console.log('🧬 DeltaSignal frontend initialized');
    }
    
    setupComponents() {
        // Initialize network visualizer (Cytoscape)
        this.networkVisualizer = new NetworkVisualizer('cytoscape-container');
        
        // Initialize pathway browser
        this.pathwayBrowser = new PathwayBrowser('pathway-browser');
        
        // Initialize file upload handler
        this.fileUploadHandler = new FileUploadHandler();
        
        // Initialize results panel
        this.resultsPanel = new ResultsPanel();
        
        // Make results panel globally available for info popups
        window.resultsPanel = this.resultsPanel;
        
        // Setup component event listeners
        this.fileUploadHandler.onFilesProcessed = (files) => {
            this.processNetworkFiles(files);
        };
        
        this.networkVisualizer.onNodeClick = (nodeData) => {
            this.handleNodeClick(nodeData);
        };
        
        this.pathwayBrowser.onPathwaySelect = (pathway) => {
            this.handlePathwaySelect(pathway);
        };
    }
    
    setupEventListeners() {
        // Upload button
        document.getElementById('upload-btn').addEventListener('click', () => {
            console.log('Upload button clicked!');
            console.log('FileUploadHandler:', this.fileUploadHandler);
            this.fileUploadHandler.showModal();
        });
        
        // Solve button
        document.getElementById('solve-btn').addEventListener('click', () => {
            this.solveSteadyState();
        });
        
        // Perturbation/Node Editor button
        document.getElementById('perturbation-btn').addEventListener('click', () => {
            this.showNodeEditor();
        });
        
        // Export button
        document.getElementById('export-btn').addEventListener('click', () => {
            this.exportResults();
        });
        
        // View tabs switching
        const viewTabs = document.querySelectorAll('.view-tab');
        const networkControls = document.getElementById('network-controls');
        const pathwayControls = document.getElementById('pathway-controls');
        const cytoscapeContainer = document.getElementById('cytoscape-container');
        const pathwayBrowser = document.getElementById('pathway-browser');
        
        viewTabs.forEach(tab => {
            tab.addEventListener('click', (e) => {
                const view = e.target.dataset.view;
                
                // Update active tab
                viewTabs.forEach(t => t.classList.remove('active'));
                e.target.classList.add('active');
                
                // Switch views
                if (view === 'network') {
                    networkControls.style.display = 'flex';
                    pathwayControls.style.display = 'none';
                    cytoscapeContainer.style.display = 'block';
                    pathwayBrowser.style.display = 'none';
                } else if (view === 'pathway') {
                    networkControls.style.display = 'none';
                    pathwayControls.style.display = 'flex';
                    cytoscapeContainer.style.display = 'none';
                    pathwayBrowser.style.display = 'block';
                }
            });
        });
        
        // View selector (within Network View)
        document.getElementById('view-selector').addEventListener('change', (e) => {
            this.networkVisualizer.switchView(e.target.value);
        });
        
        // Layout selector
        document.getElementById('layout-selector').addEventListener('change', (e) => {
            this.networkVisualizer.changeLayout(e.target.value);
        });
        
        // Fit buttons
        document.getElementById('network-fit').addEventListener('click', () => {
            this.networkVisualizer.fit();
        });
        
        document.getElementById('pathway-fit').addEventListener('click', () => {
            this.pathwayBrowser.fit();
        });
        
        // Label toggle
        document.getElementById('show-labels').addEventListener('click', (e) => {
            const showLabels = e.target.classList.toggle('active');
            this.networkVisualizer.toggleLabels(showLabels);
        });
        
        // Note: Node search with autocomplete is now handled directly in NetworkVisualizer
        // The search input events are managed by the NetworkVisualizer.initializeAutocomplete() method

        // Zoom slider
        document.getElementById('zoom-slider').addEventListener('input', (e) => {
            const zoomLevel = parseFloat(e.target.value);
            this.networkVisualizer.setZoom(zoomLevel);
            document.getElementById('zoom-level').textContent = Math.round(zoomLevel * 100) + '%';
        });
        
        // Results panel toggle
        document.getElementById('results-toggle').addEventListener('click', () => {
            const resultsPanel = document.querySelector('.results-panel');
            const resizeHandle = document.getElementById('resize-handle');
            const app = document.getElementById('app');
            resultsPanel.classList.toggle('collapsed');
            
            // Show/hide resize handle and adjust grid
            if (!resultsPanel.classList.contains('collapsed')) {
                resizeHandle.classList.add('active');
                const currentHeight = resultsPanel.offsetHeight || 300;
                app.style.gridTemplateRows = `auto 1fr auto ${currentHeight}px`;
            } else {
                resizeHandle.classList.remove('active');
                app.style.gridTemplateRows = `auto 1fr auto auto`;
            }
            
            // Toggle body class for browsers without :has() support
            document.body.classList.toggle('results-expanded', !resultsPanel.classList.contains('collapsed'));
        });
        
        // Add resize functionality
        const resizeHandle = document.getElementById('resize-handle');
        const resultsPanel = document.getElementById('results-panel');
        const mainContent = document.querySelector('.main-content');
        
        let isResizing = false;
        let startY = 0;
        let startHeight = 0;
        
        resizeHandle.addEventListener('mousedown', (e) => {
            isResizing = true;
            startY = e.clientY;
            startHeight = resultsPanel.offsetHeight;
            document.body.style.cursor = 'ns-resize';
            e.preventDefault();
        });
        
        document.addEventListener('mousemove', (e) => {
            if (!isResizing) return;
            
            const deltaY = startY - e.clientY;
            const newHeight = Math.min(Math.max(150, startHeight + deltaY), window.innerHeight - 200);
            
            // Set explicit height and ensure it stays at bottom
            resultsPanel.style.height = newHeight + 'px';
            resultsPanel.style.minHeight = newHeight + 'px';
            resultsPanel.style.maxHeight = newHeight + 'px';
            
            // Update grid to accommodate new size
            const app = document.getElementById('app');
            app.style.gridTemplateRows = `auto 1fr auto ${newHeight}px`;
        });
        
        document.addEventListener('mouseup', () => {
            if (isResizing) {
                isResizing = false;
                document.body.style.cursor = '';
            }
        });
        
        // Pathway selector
        document.getElementById('pathway-selector').addEventListener('change', (e) => {
            const pathwayId = e.target.value;
            if (pathwayId) {
                this.pathwayBrowser.loadPathway(pathwayId);
            }
        });
        
        // Node editor modal controls
        document.getElementById('node-editor-close').addEventListener('click', () => {
            this.hideNodeEditor();
        });
        
        document.getElementById('clear-perturbations-btn').addEventListener('click', () => {
            this.clearPerturbations();
        });
        
        // Close button (was apply-perturbations-btn)
        document.getElementById('close-editor-btn').addEventListener('click', () => {
            this.hideNodeEditor();
        });
        
        // Panel solve button
        document.getElementById('panel-solve-btn').addEventListener('click', () => {
            this.solveSteadyState();
        });
        
        document.getElementById('save-scenario-btn').addEventListener('click', () => {
            this.saveScenario();
        });
    }
    
    setupUI() {
        // Initially disable solve and export buttons
        document.getElementById('solve-btn').disabled = true;
        document.getElementById('export-btn').disabled = true;
        document.getElementById('pathway-selector').disabled = true;
    }
    
    async processNetworkFiles(files) {
        this.showLoading('Processing network files...');

        try {
            let networkData;

            // Check if this is a reaction network (loaded directly)
            if (files.reaction_network) {
                networkData = files.reaction_network;
                console.log('Loading reaction network directly:', networkData);
            } else {
                // Send files to backend for processing (regulatory network)
                networkData = await this.apiClient.parseNetwork(files);
            }

            this.currentNetwork = networkData;

            // Update visualizations
            await this.loadNetworkVisualization(networkData);

            // Enable controls
            document.getElementById('solve-btn').disabled = false;
            document.getElementById('panel-solve-btn').disabled = false;
            document.getElementById('perturbation-btn').disabled = false;
            document.getElementById('pathway-selector').disabled = false;

            // Populate pathway selector if pathways are available
            this.populatePathwaySelector(networkData.pathways || []);

            const networkType = networkData.network_type === 'reaction' ? 'reaction' : 'regulatory';
            this.showSuccess(`${networkType.charAt(0).toUpperCase() + networkType.slice(1)} network loaded successfully!`);

        } catch (error) {
            console.error('Error processing network files:', error);
            this.showError('Failed to process network files: ' + error.message);
        } finally {
            this.hideLoading();
        }
    }

    convertReactionToRegulatoryNetwork(reactionNetwork) {
        console.log('Converting reaction network to regulatory format...');

        const nodes = [];
        const edges = [];

        // Convert entities and reactions to nodes
        if (reactionNetwork.entities) {
            reactionNetwork.entities.forEach(entity => {
                nodes.push({
                    uuid: entity.uuid,
                    name: entity.name,
                    reactome_id: null,
                    entity_type: entity.type,
                    set_id: null
                });
            });
        }

        if (reactionNetwork.reactions) {
            reactionNetwork.reactions.forEach(reaction => {
                nodes.push({
                    uuid: reaction.uuid,
                    name: reaction.name,
                    reactome_id: null,
                    entity_type: 'reaction',
                    set_id: null
                });
            });
        }

        // Convert connections to edges
        if (reactionNetwork.connections) {
            reactionNetwork.connections.forEach(connection => {
                // Convert reaction network connections to regulatory edges
                let isPositive = true; // Most biochemical connections are positive

                // Determine edge type based on role
                if (connection.role === 'substrate' || connection.role === 'product') {
                    isPositive = true;
                } else if (connection.role === 'catalyst') {
                    isPositive = true;
                } else if (connection.role === 'modifier') {
                    isPositive = connection.type !== 'inhibition';
                }

                edges.push({
                    parent_uuid: connection.source,
                    child_uuid: connection.target,
                    is_and: false, // Most connections are OR logic
                    is_positive: isPositive,
                    stoichiometry: connection.stoichiometry || 1
                });
            });
        }

        const convertedNetwork = {
            nodes: nodes,
            edges: edges,
            pathways: reactionNetwork.pathways || [],
            metadata: {
                ...reactionNetwork.metadata,
                converted_from: 'reaction_network',
                original_type: reactionNetwork.network_type
            },
            stats: {
                nodeCount: nodes.length,
                edgeCount: edges.length,
                pathwayCount: (reactionNetwork.pathways || []).length
            }
        };

        console.log(`Converted ${nodes.length} nodes and ${edges.length} edges`);
        return convertedNetwork;
    }

    async loadNetworkVisualization(networkData) {
        // Load network topology in Cytoscape
        await this.networkVisualizer.loadNetwork(networkData);
        
        // If pathway mappings exist, show pathway browser
        if (networkData.pathways && networkData.pathways.length > 0) {
            await this.pathwayBrowser.initialize(networkData.pathways);
        }
    }
    
    async solveSteadyState() {
        if (!this.currentNetwork) {
            this.showError('No network loaded');
            return;
        }

        this.showLoading('Solving steady-state...');

        try {
            console.log('=== ANALYSIS DEBUG START ===');
            console.log('Current network type:', this.currentNetwork.network_type);
            console.log('Current network structure:', Object.keys(this.currentNetwork));

            // Get current perturbations
            const perturbations = this.networkVisualizer.getPerturbations();
            console.log('Perturbations to send:', perturbations);
            let results;

            if (Object.keys(perturbations).length > 0) {
                // Convert perturbations to observations format
                const observations = {};
                Object.entries(perturbations).forEach(([nodeId, activity]) => {
                    observations[nodeId] = {
                        activity: activity,
                        confidence: 1.0
                    };
                });

                console.log('Sending to API - Network:', {
                    type: this.currentNetwork.network_type,
                    nodeCount: this.currentNetwork.nodes?.length || 'N/A',
                    entityCount: this.currentNetwork.entities?.length || 'N/A',
                    reactionCount: this.currentNetwork.reactions?.length || 'N/A'
                });
                console.log('Sending to API - Observations:', observations);

                // Convert reaction network to regulatory format if needed
                let networkToSend = this.currentNetwork;
                if (this.currentNetwork.network_type === 'reaction') {
                    console.log('Converting reaction network to regulatory format for backend');
                    networkToSend = this.convertReactionToRegulatoryNetwork(this.currentNetwork);
                    console.log('Converted network:', networkToSend);
                }

                // Solve with perturbations
                results = await this.apiClient.solveSteadyState(networkToSend, {
                    method: 'perturbation',
                    observations: observations
                });

                console.log('API Response (with perturbations):', results);

                // Update with perturbation results
                console.log('About to call updateWithPerturbationResults');
                console.log('NetworkVisualizer exists:', !!this.networkVisualizer);
                try {
                    this.networkVisualizer.updateWithPerturbationResults(results, perturbations);
                    console.log('updateWithPerturbationResults completed successfully');
                } catch (error) {
                    console.error('updateWithPerturbationResults failed:', error);
                }
            } else {
                console.log('Sending to API - Network (no perturbations):', {
                    type: this.currentNetwork.network_type,
                    nodeCount: this.currentNetwork.nodes?.length || 'N/A',
                    entityCount: this.currentNetwork.entities?.length || 'N/A',
                    reactionCount: this.currentNetwork.reactions?.length || 'N/A'
                });

                // Convert reaction network to regulatory format if needed
                let networkToSend = this.currentNetwork;
                if (this.currentNetwork.network_type === 'reaction') {
                    console.log('Converting reaction network to regulatory format for backend');
                    networkToSend = this.convertReactionToRegulatoryNetwork(this.currentNetwork);
                    console.log('Converted network:', networkToSend);
                }

                // Solve without perturbations
                results = await this.apiClient.solveSteadyState(networkToSend);

                console.log('API Response (no perturbations):', results);

                // Update with normal results
                console.log('About to call updateWithResults');
                console.log('NetworkVisualizer exists:', !!this.networkVisualizer);
                try {
                    this.networkVisualizer.updateWithResults(results);
                    console.log('updateWithResults completed successfully');
                } catch (error) {
                    console.error('updateWithResults failed:', error);
                }
            }

            console.log('=== ANALYSIS DEBUG END ===');
            
            this.currentResults = results;
            
            // Update other visualizations
            this.pathwayBrowser.updateWithResults(results);

            // Show results panel
            console.log('About to call displayResults with:', results);
            console.log('Results panel exists:', !!this.resultsPanel);
            try {
                this.resultsPanel.displayResults(results);
                console.log('displayResults completed successfully');
            } catch (error) {
                console.error('displayResults failed:', error);
            }
            const resultsPanel = document.querySelector('.results-panel');
            const resizeHandle = document.getElementById('resize-handle');
            resultsPanel.classList.remove('collapsed');
            resizeHandle.classList.add('active');
            document.body.classList.add('results-expanded');
            
            // Set initial grid layout with results panel
            const app = document.getElementById('app');
            const initialHeight = 300; // Default height for results panel
            resultsPanel.style.height = initialHeight + 'px';
            app.style.gridTemplateRows = `auto 1fr auto ${initialHeight}px`;
            
            // Enable export
            document.getElementById('export-btn').disabled = false;
            
            this.showSuccess('Analysis completed!');
            
        } catch (error) {
            console.error('Error solving steady-state:', error);
            this.showError('Analysis failed: ' + error.message);
        } finally {
            this.hideLoading();
        }
    }
    
    exportResults() {
        if (!this.currentResults) {
            this.showError('No results to export');
            return;
        }
        
        // Create export data
        const exportData = {
            network: this.currentNetwork,
            results: this.currentResults,
            timestamp: new Date().toISOString(),
            version: '0.1.0'
        };
        
        // Download as JSON
        const blob = new Blob([JSON.stringify(exportData, null, 2)], {
            type: 'application/json'
        });
        
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `deltasignal-results-${Date.now()}.json`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
        
        this.showSuccess('Results exported successfully!');
    }
    
    populatePathwaySelector(pathways) {
        const selector = document.getElementById('pathway-selector');
        
        // Clear existing options except first
        while (selector.children.length > 1) {
            selector.removeChild(selector.lastChild);
        }
        
        // Add pathway options
        pathways.forEach(pathway => {
            const option = document.createElement('option');
            option.value = pathway.id;
            option.textContent = pathway.name || pathway.id;
            selector.appendChild(option);
        });
    }
    
    handleNodeClick(nodeData) {
        console.log('Node clicked:', nodeData);
        
        // Set this node as selected for perturbation
        this.selectedNode = nodeData;
        
        // Check if node already has a perturbation, otherwise use default
        const existingPerturbations = this.networkVisualizer.getPerturbations();
        if (!existingPerturbations[nodeData.uuid]) {
            // Add node to perturbation list with default activity (1x normal)
            const defaultActivity = 0.01; // Default activity level (1x fold change = normal)
            this.networkVisualizer.setPerturbation(nodeData.uuid, defaultActivity);
        }
        // If already has a perturbation, keep the existing value
        
        // Show node editor modal
        this.showNodeEditor();
        
        // Highlight corresponding elements in pathway browser
        if (this.pathwayBrowser && nodeData.reactomeId) {
            this.pathwayBrowser.highlightElement(nodeData.reactomeId);
        }
        
        this.showToast(`Node ${nodeData.name || nodeData.uuid} selected for editing`, 'info');
    }
    
    handlePathwaySelect(pathway) {
        console.log('Pathway selected:', pathway);
        
        // Filter network view to show only nodes in this pathway
        this.networkVisualizer.filterByPathway(pathway);
    }
    
    showNodeDetails(nodeData) {
        // This could be expanded to show a detailed modal or sidebar
        console.log('Node details:', nodeData);
    }
    
    // UI state management
    showLoading(message = 'Loading...') {
        const overlay = document.getElementById('loading-overlay');
        const text = document.getElementById('loading-text');
        text.textContent = message;
        overlay.classList.add('show');
    }
    
    hideLoading() {
        document.getElementById('loading-overlay').classList.remove('show');
    }
    
    showSuccess(message) {
        this.showToast(message, 'success');
    }
    
    showError(message) {
        this.showToast(message, 'error');
    }
    
    showToast(message, type = 'info') {
        // Simple toast notification
        const toast = document.createElement('div');
        toast.className = `toast toast-${type}`;
        toast.textContent = message;
        
        // Add toast styles if not already added
        if (!document.querySelector('#toast-styles')) {
            const style = document.createElement('style');
            style.id = 'toast-styles';
            style.textContent = `
                .toast {
                    position: fixed;
                    top: 20px;
                    right: 20px;
                    padding: 12px 20px;
                    border-radius: 6px;
                    color: white;
                    font-weight: 500;
                    z-index: 3000;
                    animation: toast-slide-in 0.3s ease;
                }
                .toast-success { background: var(--success-color); }
                .toast-error { background: var(--error-color); }
                .toast-info { background: var(--primary-color); }
                @keyframes toast-slide-in {
                    from { transform: translateX(100%); opacity: 0; }
                    to { transform: translateX(0); opacity: 1; }
                }
            `;
            document.head.appendChild(style);
        }
        
        document.body.appendChild(toast);
        
        // Auto remove after 3 seconds
        setTimeout(() => {
            if (toast.parentNode) {
                toast.parentNode.removeChild(toast);
            }
        }, 3000);
    }
    
    // Node Editor Methods
    showNodeEditor() {
        const modal = document.getElementById('node-editor-modal');
        modal.classList.add('show');
        this.updateNodeEditorList();
    }
    
    hideNodeEditor() {
        const modal = document.getElementById('node-editor-modal');
        modal.classList.remove('show');
    }
    
    updateNodeEditorList() {
        const container = document.getElementById('node-editor-list');
        const perturbations = this.networkVisualizer.getPerturbations();
        
        if (Object.keys(perturbations).length === 0) {
            container.innerHTML = '<p style="text-align: center; color: #94a3b8; font-style: italic;">No nodes selected. Click on a node in the network view.</p>';
            return;
        }
        
        container.innerHTML = Object.entries(perturbations)
            .map(([nodeId, activity]) => {
                let node = null;

                // Handle different network types
                if (this.currentNetwork.network_type === 'reaction') {
                    // Look in entities and reactions
                    if (this.currentNetwork.entities) {
                        node = this.currentNetwork.entities.find(n => n.uuid === nodeId);
                    }
                    if (!node && this.currentNetwork.reactions) {
                        node = this.currentNetwork.reactions.find(n => n.uuid === nodeId);
                    }
                } else {
                    // Regular regulatory network
                    if (this.currentNetwork.nodes) {
                        node = this.currentNetwork.nodes.find(n => n.uuid === nodeId);
                    }
                }

                const nodeName = node ? (node.name || nodeId) : nodeId;
                // Convert internal activity (0-1) to fold change (0-100x)
                const foldChange = activity * 100;
                
                return `
                    <div class="node-editor-item" style="display: flex; justify-content: space-between; align-items: center; padding: 0.75rem; border: 1px solid #e5e7eb; border-radius: 6px; margin-bottom: 0.5rem;">
                        <div>
                            <div style="font-weight: 600; color: #374151; font-size: 0.9rem;">${this.truncateText(nodeName, 20)}</div>
                            <div style="font-size: 0.75rem; color: #6b7280; font-family: monospace;">${nodeId}</div>
                        </div>
                        <div style="display: flex; align-items: center; gap: 0.5rem;">
                            <input 
                                type="range" 
                                min="0" 
                                max="100" 
                                step="0.5" 
                                value="${foldChange}"
                                data-node-id="${nodeId}"
                                class="node-activity-slider"
                                style="width: 100px;"
                            >
                            <input
                                type="number"
                                min="0"
                                max="100"
                                step="0.1"
                                value="${foldChange.toFixed(1)}"
                                data-node-id="${nodeId}"
                                class="node-activity-input"
                                style="width: 55px; padding: 2px 4px; font-size: 0.8rem; border: 1px solid #d1d5db; border-radius: 4px; text-align: right;"
                            >
                            <span style="font-weight: 600; color: ${foldChange < 0.5 ? '#ef4444' : foldChange > 2 ? '#10b981' : '#3b82f6'}; font-size: 0.8rem;">x</span>
                            <button class="remove-perturbation-btn btn-icon" data-node-id="${nodeId}" style="width: 24px; height: 24px; font-size: 0.7rem;">✕</button>
                        </div>
                    </div>
                `;
            }).join('');
        
        // Add event listeners to sliders
        container.querySelectorAll('.node-activity-slider').forEach(slider => {
            slider.addEventListener('input', (e) => {
                const nodeId = e.target.dataset.nodeId;
                const foldChange = parseFloat(e.target.value);
                // Convert fold change (0-100x) to internal 0-1 scale
                const activity = foldChange / 100;
                
                this.networkVisualizer.setPerturbation(nodeId, activity);
                
                // Update the paired input field
                const parentItem = e.target.closest('.node-editor-item');
                const input = parentItem.querySelector('.node-activity-input');
                input.value = foldChange.toFixed(1);
                
                // Update color based on value
                const xLabel = parentItem.querySelector('span');
                xLabel.style.color = foldChange < 0.5 ? '#ef4444' : foldChange > 2 ? '#10b981' : '#3b82f6';
            });
        });
        
        // Add event listeners to input fields
        container.querySelectorAll('.node-activity-input').forEach(input => {
            input.addEventListener('change', (e) => {
                const nodeId = e.target.dataset.nodeId;
                let foldChange = parseFloat(e.target.value);
                
                // Clamp to valid range
                foldChange = Math.max(0, Math.min(100, foldChange));
                e.target.value = foldChange.toFixed(1);
                
                // Convert fold change (0-100x) to internal 0-1 scale
                const activity = foldChange / 100;
                
                this.networkVisualizer.setPerturbation(nodeId, activity);
                
                // Update the paired slider
                const parentItem = e.target.closest('.node-editor-item');
                const slider = parentItem.querySelector('.node-activity-slider');
                slider.value = foldChange;
                
                // Update color based on value
                const xLabel = parentItem.querySelector('span');
                xLabel.style.color = foldChange < 0.5 ? '#ef4444' : foldChange > 2 ? '#10b981' : '#3b82f6';
            });
        });
        
        container.querySelectorAll('.remove-perturbation-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                e.preventDefault();
                e.stopPropagation();
                // Get nodeId from the button itself, not e.target (which might be the text)
                const nodeId = btn.dataset.nodeId;
                this.removePerturbation(nodeId);
                this.updateNodeEditorList();
            });
        });
    }
    
    clearPerturbations() {
        this.networkVisualizer.clearPerturbations();
        this.updateNodeEditorList();
        this.showToast('Perturbations cleared', 'info');
    }
    
    removePerturbation(nodeId) {
        const perturbations = this.networkVisualizer.getPerturbations();
        delete perturbations[nodeId];
        
        // Update the visualizer
        this.networkVisualizer.perturbations = perturbations;
        this.networkVisualizer.updateNodeDisplay(nodeId, 0);
        
        // Remove visual indicators
        const nodeElement = document.querySelector(`[data-node-id="${nodeId}"]`);
        if (nodeElement) {
            nodeElement.classList.remove('perturbed');
            nodeElement.style.borderLeft = '';
        }
    }
    
    async applyPerturbations() {
        const perturbations = this.networkVisualizer.getPerturbations();
        
        if (Object.keys(perturbations).length === 0) {
            this.showToast('No perturbations to apply', 'warning');
            return;
        }
        
        this.showLoading('Analyzing perturbation effects...');
        this.hideNodeEditor();
        
        try {
            // Create modified network with perturbations as observations
            const modifiedNetwork = { ...this.currentNetwork };
            const observations = {};
            
            // Convert perturbations to observations format
            Object.entries(perturbations).forEach(([nodeId, activity]) => {
                observations[nodeId] = {
                    activity: activity * 100, // Convert to percentage
                    confidence: 1.0 // High confidence for user-set values
                };
            });
            
            console.log('Sending perturbations to API:', observations);
            
            // Solve with perturbations
            const results = await this.apiClient.solveSteadyState(modifiedNetwork, {
                method: 'perturbation',
                observations: observations
            });
            
            // Update visualizations with perturbation results
            this.networkVisualizer.updateWithPerturbationResults(results, perturbations);
            this.pathwayBrowser.updateWithResults(results);
            this.resultsPanel.displayResults(results);
            
            // Show results panel
            const resultsPanel = document.querySelector('.results-panel');
            resultsPanel.classList.remove('collapsed');
            document.body.classList.add('results-expanded');
            
            this.showSuccess(`Perturbation analysis complete! ${Object.keys(perturbations).length} nodes modified.`);
            
        } catch (error) {
            console.error('Error applying perturbations:', error);
            this.showError('Failed to analyze perturbations: ' + error.message);
        } finally {
            this.hideLoading();
        }
    }
    
    saveScenario() {
        const scenarioName = document.getElementById('scenario-name').value.trim();
        if (!scenarioName) {
            this.showToast('Please enter a scenario name', 'warning');
            return;
        }
        
        const perturbations = this.networkVisualizer.getPerturbations();
        if (Object.keys(perturbations).length === 0) {
            this.showToast('No perturbations to save', 'warning');
            return;
        }
        
        this.savedScenarios[scenarioName] = {
            perturbations: { ...perturbations },
            timestamp: new Date().toISOString(),
            nodeCount: Object.keys(perturbations).length
        };
        
        this.updateSavedScenarios();
        document.getElementById('scenario-name').value = '';
        this.showToast(`Scenario "${scenarioName}" saved`, 'success');
    }
    
    updateSavedScenarios() {
        const container = document.getElementById('saved-scenarios');
        
        if (Object.keys(this.savedScenarios).length === 0) {
            container.innerHTML = '';
            return;
        }
        
        container.innerHTML = Object.entries(this.savedScenarios)
            .map(([name, scenario]) => `
                <div style="display: flex; justify-content: space-between; align-items: center; padding: 0.5rem; background: white; border: 1px solid #e5e7eb; border-radius: 4px; margin-bottom: 0.25rem;">
                    <div>
                        <span style="font-weight: 500; font-size: 0.8rem;">${name}</span>
                        <span style="color: #6b7280; font-size: 0.7rem; margin-left: 0.5rem;">(${scenario.nodeCount} nodes)</span>
                    </div>
                    <button class="load-scenario-btn btn-secondary" data-scenario="${name}" style="padding: 0.25rem 0.5rem; font-size: 0.7rem;">Load</button>
                </div>
            `).join('');
        
        // Add event listeners to load buttons
        container.querySelectorAll('.load-scenario-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const scenarioName = e.target.dataset.scenario;
                this.loadScenario(scenarioName);
            });
        });
    }
    
    loadScenario(scenarioName) {
        const scenario = this.savedScenarios[scenarioName];
        if (!scenario) return;
        
        // Clear current perturbations
        this.networkVisualizer.clearPerturbations();
        
        // Apply saved perturbations
        Object.entries(scenario.perturbations).forEach(([nodeId, activity]) => {
            this.networkVisualizer.setPerturbation(nodeId, activity);
        });
        
        this.updateNodeEditorList();
        this.showToast(`Scenario "${scenarioName}" loaded`, 'success');
    }
    
    truncateText(text, maxLength) {
        if (text && text.length > maxLength) {
            return text.substring(0, maxLength - 3) + '...';
        }
        return text || '';
    }
}

// Initialize app when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
    window.deltaSignalApp = new DeltaSignalApp();
});