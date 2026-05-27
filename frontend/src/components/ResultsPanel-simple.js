// Simple Results Panel - Working Version
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
        
        console.log('Displaying results:', results);
        
        // Display node activities
        this.displayNodeActivities(results.node_activities || {});
        
        // Display influence scores
        this.displayInfluenceScores(results.influence_scores || {});
        
        // Display pathway statistics
        this.displayPathwayStats(results);
        
        console.log('Results displayed successfully');
    }
    
    displayNodeActivities(nodeActivities) {
        if (!this.nodeActivitiesContainer) return;
        
        this.nodeActivitiesContainer.innerHTML = '';
        
        if (Object.keys(nodeActivities).length === 0) {
            this.nodeActivitiesContainer.innerHTML = '<p style="text-align: center; color: #666; font-style: italic; padding: 2rem;">No activity data available</p>';
            return;
        }
        
        // Create simple table
        const table = document.createElement('table');
        table.style.width = '100%';
        table.style.borderCollapse = 'collapse';
        table.style.fontSize = '0.875rem';
        
        // Table header
        const thead = document.createElement('thead');
        thead.innerHTML = `
            <tr>
                <th style="background: #f1f5f9; padding: 0.75rem 0.5rem; text-align: left; font-weight: 600; border-bottom: 2px solid #e2e8f0;">Node ID</th>
                <th style="background: #f1f5f9; padding: 0.75rem 0.5rem; text-align: left; font-weight: 600; border-bottom: 2px solid #e2e8f0;">Activity (%)</th>
                <th style="background: #f1f5f9; padding: 0.75rem 0.5rem; text-align: left; font-weight: 600; border-bottom: 2px solid #e2e8f0;">Status</th>
            </tr>
        `;
        table.appendChild(thead);
        
        // Table body
        const tbody = document.createElement('tbody');
        
        // Sort by activity level (highest first)
        const sortedActivities = Object.entries(nodeActivities)
            .sort(([,a], [,b]) => b - a)
            .slice(0, 10); // Show only top 10
        
        sortedActivities.forEach(([nodeId, activity], index) => {
            const row = document.createElement('tr');
            if (index % 2 === 0) {
                row.style.background = '#f8fafc';
            }
            
            const activityPercent = (activity * 100).toFixed(1);
            const status = this.getActivityStatus(activity);
            const statusColor = this.getStatusColor(status);
            
            row.innerHTML = `
                <td style="padding: 0.625rem 0.5rem; border-bottom: 1px solid #e2e8f0;">
                    <span style="font-family: monospace; font-size: 0.8rem; color: #64748b;" title="${nodeId}">
                        ${this.truncateNodeId(nodeId)}
                    </span>
                </td>
                <td style="padding: 0.625rem 0.5rem; border-bottom: 1px solid #e2e8f0;">
                    <div style="display: flex; align-items: center; gap: 0.5rem;">
                        <div style="width: 60px; height: 16px; background: #e2e8f0; border-radius: 8px; overflow: hidden;">
                            <div style="height: 100%; background: ${statusColor}; width: ${activityPercent}%; transition: width 0.3s ease;"></div>
                        </div>
                        <span style="font-weight: 600; color: #333;">${activityPercent}%</span>
                    </div>
                </td>
                <td style="padding: 0.625rem 0.5rem; border-bottom: 1px solid #e2e8f0;">
                    <span style="padding: 0.25rem 0.5rem; background: ${statusColor}20; color: ${statusColor}; border-radius: 4px; font-size: 0.75rem; font-weight: 600; text-transform: uppercase;">
                        ${status}
                    </span>
                </td>
            `;
            
            tbody.appendChild(row);
        });
        
        table.appendChild(tbody);
        this.nodeActivitiesContainer.appendChild(table);
    }
    
    displayInfluenceScores(influenceScores) {
        if (!this.influenceScoresContainer) return;
        
        this.influenceScoresContainer.innerHTML = '';
        
        if (Object.keys(influenceScores).length === 0) {
            this.influenceScoresContainer.innerHTML = '<p style="text-align: center; color: #666; font-style: italic; padding: 2rem;">No influence data available</p>';
            return;
        }
        
        // Sort by influence score (highest first)
        const sortedInfluence = Object.entries(influenceScores)
            .sort(([,a], [,b]) => b - a)
            .slice(0, 10); // Show only top 10
        
        const list = document.createElement('div');
        list.style.display = 'flex';
        list.style.flexDirection = 'column';
        list.style.gap = '0.5rem';
        
        sortedInfluence.forEach(([nodeId, score], index) => {
            const item = document.createElement('div');
            item.style.display = 'flex';
            item.style.justifyContent = 'space-between';
            item.style.alignItems = 'center';
            item.style.padding = '0.5rem 0.75rem';
            item.style.background = index % 2 === 0 ? '#f8fafc' : 'white';
            item.style.borderRadius = '6px';
            item.style.border = '1px solid #e2e8f0';
            
            const rankBadge = index < 3 ? this.getRankEmoji(index) : `${index + 1}.`;
            
            item.innerHTML = `
                <div style="display: flex; align-items: center; gap: 0.5rem;">
                    <span style="font-weight: bold; color: #3b82f6;">${rankBadge}</span>
                    <span style="font-family: monospace; font-size: 0.8rem; color: #64748b;" title="${nodeId}">
                        ${this.truncateNodeId(nodeId)}
                    </span>
                </div>
                <span style="font-family: monospace; font-weight: 600; color: #3b82f6;">
                    ${score.toFixed(3)}
                </span>
            `;
            
            list.appendChild(item);
        });
        
        this.influenceScoresContainer.appendChild(list);
    }
    
    displayPathwayStats(results) {
        if (!this.pathwayStatsContainer) return;
        
        this.pathwayStatsContainer.innerHTML = '';
        
        const stats = this.calculateStats(results);
        
        stats.forEach(stat => {
            const statElement = document.createElement('div');
            statElement.style.display = 'flex';
            statElement.style.justifyContent = 'space-between';
            statElement.style.padding = '0.5rem 0.75rem';
            statElement.style.background = '#f1f5f9';
            statElement.style.borderRadius = '6px';
            statElement.style.marginBottom = '0.5rem';
            
            statElement.innerHTML = `
                <span style="color: #64748b; font-size: 0.875rem;">${stat.label}</span>
                <span style="font-weight: 600; color: #3b82f6;">${stat.value}</span>
            `;
            
            this.pathwayStatsContainer.appendChild(statElement);
        });
    }
    
    calculateStats(results) {
        const stats = [];
        
        // Node count
        const nodeCount = results.node_activities ? Object.keys(results.node_activities).length : 0;
        stats.push({ label: 'Total Nodes', value: nodeCount });
        
        // Active nodes (>50% activity)
        if (results.node_activities) {
            const activeNodes = Object.values(results.node_activities).filter(activity => activity > 0.5).length;
            stats.push({ label: 'Active Nodes', value: activeNodes });
        }
        
        // Observed nodes
        if (results.observations) {
            const observedCount = Object.keys(results.observations).length;
            stats.push({ label: 'Observed', value: observedCount });
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
        if (activity > 0.7) return 'high';
        if (activity > 0.3) return 'medium';
        return 'low';
    }
    
    getStatusColor(status) {
        switch(status) {
            case 'high': return '#10b981';
            case 'medium': return '#f59e0b';
            case 'low': return '#ef4444';
            default: return '#64748b';
        }
    }
    
    getRankEmoji(index) {
        switch(index) {
            case 0: return '🥇';
            case 1: return '🥈';
            case 2: return '🥉';
            default: return `${index + 1}.`;
        }
    }
    
    truncateNodeId(nodeId) {
        if (nodeId && nodeId.length > 12) {
            return nodeId.substring(0, 12) + '...';
        }
        return nodeId || '';
    }
    
    // Clear all displayed results
    clear() {
        if (this.nodeActivitiesContainer) this.nodeActivitiesContainer.innerHTML = '';
        if (this.influenceScoresContainer) this.influenceScoresContainer.innerHTML = '';
        if (this.pathwayStatsContainer) this.pathwayStatsContainer.innerHTML = '';
    }
}