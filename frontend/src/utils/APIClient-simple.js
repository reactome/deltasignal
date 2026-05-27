// Simple API Client - Working Version
export default class APIClient {
    constructor(baseURL = '') {
        this.baseURL = baseURL;
    }
    
    async parseNetwork(files) {
        // Convert files to FormData
        const formData = new FormData();
        
        // Add files to FormData
        if (files['logic-file']) {
            formData.append('logic_network', files['logic-file']);
        }
        
        if (files['uuid-file']) {
            formData.append('uuid_mapping', files['uuid-file']);
        }
        
        if (files['set-file']) {
            formData.append('set_mappings', files['set-file']);
        }
        
        if (files['observations-file']) {
            formData.append('observations', files['observations-file']);
        }
        
        try {
            console.log('Sending files to API...');
            const response = await fetch('/api/parse', {
                method: 'POST',
                body: formData
            });
            
            if (!response.ok) {
                throw new Error(`API returned ${response.status}: ${response.statusText}`);
            }
            
            const result = await response.json();
            console.log('Parse API response:', result);
            
            return {
                nodes: result.nodes || [],
                edges: result.edges || [],
                pathways: result.pathways || [],
                metadata: result.metadata || {},
                stats: {
                    nodeCount: result.nodes?.length || 0,
                    edgeCount: result.edges?.length || 0,
                    pathwayCount: result.pathways?.length || 0
                }
            };
            
        } catch (error) {
            console.error('Parse API error:', error);
            throw new Error(`Failed to parse network: ${error.message}`);
        }
    }
    
    async solveSteadyState(networkData, options = {}) {
        const payload = {
            network: networkData,
            options: options
        };
        
        try {
            console.log('Sending solve request to API...');
            const response = await fetch('/api/solve', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify(payload)
            });
            
            if (!response.ok) {
                throw new Error(`API returned ${response.status}: ${response.statusText}`);
            }
            
            const result = await response.json();
            console.log('Solve API response:', result);
            
            return {
                node_activities: result.node_activities || {},
                influence_scores: result.influence_scores || {},
                observations: result.observations || {},
                converged: result.converged || false,
                iterations: result.iterations || 0,
                solution_quality: result.solution_quality || 0,
                analysis_time: result.analysis_time || 0,
                metadata: result.metadata || {}
            };
            
        } catch (error) {
            console.error('Solve API error:', error);
            throw new Error(`Failed to solve steady-state: ${error.message}`);
        }
    }
    
    async checkHealth() {
        try {
            const response = await fetch('/api/health');
            return response.ok;
        } catch (error) {
            console.error('Health check failed:', error);
            return false;
        }
    }
}