// Simple Express server for development API endpoints
const express = require('express');
const cors = require('cors');
const multer = require('multer');
const path = require('path');
const swaggerJsdoc = require('swagger-jsdoc');
const swaggerUi = require('swagger-ui-express');

const app = express();
const PORT = 3001;

// Swagger definition
const swaggerOptions = {
  definition: {
    openapi: '3.0.0',
    info: {
      title: 'DeltaSignal API',
      version: '0.1.0',
      description: 'REST API for DeltaSignal pathway perturbation analysis',
      contact: {
        name: 'DeltaSignal Team'
      }
    },
    servers: [
      {
        url: `http://localhost:${PORT}`,
        description: 'Development server'
      }
    ]
  },
  apis: ['./server.js'] // Files containing annotations
};

const swaggerSpec = swaggerJsdoc(swaggerOptions);

// Simulate network propagation for perturbation analysis
function simulateNetworkPropagation(network, nodeActivities, observations) {
  // Create adjacency list for network with proper edge grouping
  const adjacency = {};
  network.nodes.forEach(node => {
    adjacency[node.uuid] = { 
      activators: [],  // Positive inputs
      inhibitors: []   // Negative inputs
    };
  });
  
  network.edges.forEach(edge => {
    // Ensure both nodes exist in adjacency list
    if (!adjacency[edge.parent_uuid]) {
      adjacency[edge.parent_uuid] = { activators: [], inhibitors: [] };
    }
    if (!adjacency[edge.child_uuid]) {
      adjacency[edge.child_uuid] = { activators: [], inhibitors: [] };
    }
    
    // Group edges by type for the target node
    if (edge.is_positive) {
      adjacency[edge.child_uuid].activators.push({
        source: edge.parent_uuid,
        weight: edge.stoichiometry || 1,
        is_and: edge.is_and || false
      });
    } else {
      adjacency[edge.child_uuid].inhibitors.push({
        source: edge.parent_uuid,
        weight: edge.stoichiometry || 1
      });
    }
  });
  
  // Parameters for Hill functions
  const hillCoeff = 2.0;   // Steepness of response  
  const halfMax = 0.1;     // K_m value (10% threshold) - matches Julia implementation
  const epsilon = 0.0001;  // Small value to avoid log(0)
  
  // Iteratively propagate activity changes
  const maxIterations = 10;
  const dampingFactor = 0.4; // Balanced damping for stable linear propagation
  
  for (let iter = 0; iter < maxIterations; iter++) {
    const updates = {};
    let maxChange = 0;
    
    // Calculate new activity for each non-observed node
    network.nodes.forEach(node => {
      if (observations[node.uuid]) {
        // Don't update observed/perturbed nodes
        return;
      }
      
      const activators = adjacency[node.uuid].activators;
      const inhibitors = adjacency[node.uuid].inhibitors;
      
      // Skip if no incoming edges (keep baseline)
      if (activators.length === 0 && inhibitors.length === 0) {
        return;
      }
      
      // Check if any parent nodes have changed from baseline (0.01)
      let hasChangedParents = false;
      const baseline = 0.01;
      const threshold = 0.001; // Small threshold for detecting changes
      
      activators.forEach(edge => {
        if (Math.abs(nodeActivities[edge.source] - baseline) > threshold) {
          hasChangedParents = true;
        }
      });
      
      inhibitors.forEach(edge => {
        if (Math.abs(nodeActivities[edge.source] - baseline) > threshold) {
          hasChangedParents = true;
        }
      });
      
      // Only update if at least one parent has changed from baseline
      if (!hasChangedParents) {
        return;
      }
      
      let targetActivity = 0.01; // Default baseline
      
      // Calculate activator contribution (geometric mean for AND, max for OR)
      let A = 1.0;
      if (activators.length > 0) {
        // Check if all activators use AND logic
        const useAND = activators.every(edge => edge.is_and);
        
        if (useAND && activators.length > 1) {
          // Geometric mean for AND logic (all inputs needed)
          let logSum = 0;
          let weightSum = 0;
          activators.forEach(edge => {
            const sourceActivity = nodeActivities[edge.source];
            // Transform input with Hill function
            const transformed = Math.pow(sourceActivity, hillCoeff) / 
                              (Math.pow(sourceActivity, hillCoeff) + Math.pow(halfMax, hillCoeff));
            logSum += edge.weight * Math.log(transformed + epsilon);
            weightSum += edge.weight;
          });
          A = Math.exp(logSum / weightSum);
        } else {
          // Maximum for OR logic (any input sufficient)
          let maxActivation = 0;
          activators.forEach(edge => {
            const sourceActivity = nodeActivities[edge.source];
            // Transform input with Hill function
            const transformed = Math.pow(sourceActivity, hillCoeff) / 
                              (Math.pow(sourceActivity, hillCoeff) + Math.pow(halfMax, hillCoeff));
            maxActivation = Math.max(maxActivation, transformed * edge.weight);
          });
          A = maxActivation;
        }
      }
      
      // Calculate inhibitor contribution (multiplicative suppression)
      let H = 1.0;
      if (inhibitors.length > 0) {
        inhibitors.forEach(edge => {
          const sourceActivity = nodeActivities[edge.source];
          // Hill repression: H = 1 / (1 + (source/K)^n)
          const ratio = sourceActivity / halfMax;
          const hillTerm = 1 / (1 + Math.pow(ratio, hillCoeff));
          H *= hillTerm;
        });
      }
      
      // Combined effect: activation modulated by inhibition
      const preActivation = A * H;
      
      // Final output with Hill transformation
      // For single activator with no inhibitors, OR when only one input is significantly perturbed
      const baselineLevel = 0.01;
      const perturbationThreshold = 0.001; // 0.1% change threshold
      
      // Count how many inputs are significantly different from baseline
      let perturbedInputCount = 0;
      activators.forEach(edge => {
        if (Math.abs(nodeActivities[edge.source] - baselineLevel) > perturbationThreshold) {
          perturbedInputCount++;
        }
      });
      inhibitors.forEach(edge => {
        if (Math.abs(nodeActivities[edge.source] - baselineLevel) > perturbationThreshold) {
          perturbedInputCount++;
        }
      });
      
      if ((activators.length === 1 && inhibitors.length === 0) || perturbedInputCount <= 1) {
        // True linear propagation: output = input (with decay factor for biological realism)
        targetActivity = preActivation * 0.85; // 85% efficiency - 15% signal decay per step
      } else {
        // Use hill function only when multiple inputs are actually perturbed
        targetActivity = Math.pow(preActivation, hillCoeff) / 
                        (Math.pow(preActivation, hillCoeff) + Math.pow(halfMax, hillCoeff));
      }
      
      // Apply damping for stability
      const currentActivity = nodeActivities[node.uuid];
      const newActivity = currentActivity * (1 - dampingFactor) + targetActivity * dampingFactor;
      const clampedActivity = Math.min(1.0, Math.max(0.0001, newActivity));
      
      // Track maximum change for convergence check
      const change = Math.abs(clampedActivity - currentActivity);
      maxChange = Math.max(maxChange, change);
      
      // Only update if there's a meaningful change
      if (change > 0.0001) {
        const oldFold = currentActivity * 100;
        const newFold = clampedActivity * 100;
        console.log(`  ${node.uuid}: ${oldFold.toFixed(1)}x -> ${newFold.toFixed(1)}x`);
        updates[node.uuid] = clampedActivity;
      }
    });
    
    // Apply updates
    Object.entries(updates).forEach(([nodeId, activity]) => {
      nodeActivities[nodeId] = activity;
    });
    
    console.log(`Iteration ${iter + 1}: Updated ${Object.keys(updates).length} nodes, max change: ${maxChange.toFixed(4)}`);
    
    // Check for convergence
    if (maxChange < 0.001) {
      console.log('Converged!');
      break;
    }
  }
}

// Middleware
app.use(cors());
app.use(express.json({ limit: '50mb' }));
app.use(express.static('.'));
app.use('/assets', express.static(path.join(__dirname, '../assets')));
app.use('/examples', express.static(path.join(__dirname, '../examples')));

// Swagger UI
app.use('/api-docs', swaggerUi.serve, swaggerUi.setup(swaggerSpec));

// Configure multer for file uploads
const storage = multer.memoryStorage();
const upload = multer({ storage: storage, limits: { fileSize: 10 * 1024 * 1024 } }); // 10MB limit

/**
 * @swagger
 * /api/health:
 *   get:
 *     summary: Health check endpoint
 *     description: Returns the current status of the API server
 *     tags:
 *       - System
 *     responses:
 *       200:
 *         description: Server is healthy
 *         content:
 *           application/json:
 *             schema:
 *               type: object
 *               properties:
 *                 status:
 *                   type: string
 *                   example: ok
 *                 timestamp:
 *                   type: string
 *                   format: date-time
 *                 version:
 *                   type: string
 *                   example: 0.1.0
 *                 service:
 *                   type: string
 *                   example: deltasignal-frontend-dev
 */
app.get('/api/health', (req, res) => {
  res.json({
    status: 'ok',
    timestamp: new Date().toISOString(),
    version: '0.1.0',
    service: 'deltasignal-frontend-dev'
  });
});

/**
 * @swagger
 * /api/parse:
 *   post:
 *     summary: Parse network files
 *     description: Parse TSV network files and return structured network data
 *     tags:
 *       - Network
 *     requestBody:
 *       content:
 *         multipart/form-data:
 *           schema:
 *             type: object
 *             required:
 *               - logic_network
 *               - uuid_mapping
 *             properties:
 *               logic_network:
 *                 type: string
 *                 format: binary
 *                 description: TSV file containing network logic (parent, child, is_and, is_positive, stoichiometry)
 *               uuid_mapping:
 *                 type: string
 *                 format: binary
 *                 description: TSV file mapping UUIDs to node names
 *               set_mappings:
 *                 type: string
 *                 format: binary
 *                 description: TSV file with set/group mappings (optional)
 *               observations:
 *                 type: string
 *                 format: binary
 *                 description: CSV file with observation data (optional)
 *     responses:
 *       200:
 *         description: Successfully parsed network files
 *         content:
 *           application/json:
 *             schema:
 *               type: object
 *               properties:
 *                 nodes:
 *                   type: array
 *                   items:
 *                     type: object
 *                     properties:
 *                       uuid:
 *                         type: string
 *                       name:
 *                         type: string
 *                       reactome_id:
 *                         type: string
 *                       entity_type:
 *                         type: string
 *                       set_id:
 *                         type: string
 *                         nullable: true
 *                 edges:
 *                   type: array
 *                   items:
 *                     type: object
 *                     properties:
 *                       parent_uuid:
 *                         type: string
 *                       child_uuid:
 *                         type: string
 *                       is_and:
 *                         type: boolean
 *                       is_positive:
 *                         type: boolean
 *                       stoichiometry:
 *                         type: integer
 *                 metadata:
 *                   type: object
 *       400:
 *         description: Bad request - missing required files
 *       500:
 *         description: Server error during parsing
 */
app.post('/api/parse', upload.fields([
  { name: 'logic_network', maxCount: 1 },
  { name: 'uuid_mapping', maxCount: 1 },
  { name: 'set_mappings', maxCount: 1 },
  { name: 'observations', maxCount: 1 }
]), async (req, res) => {
  try {
    console.log('Parsing network files...');
    console.log('Files received:', Object.keys(req.files || {}));
    
    // Simulate processing delay
    await new Promise(resolve => setTimeout(resolve, 1000));
    
    // Parse logic network file
    const logicFile = req.files['logic_network']?.[0];
    if (!logicFile) {
      return res.status(400).json({ error: 'Logic network file is required' });
    }
    
    const logicContent = logicFile.buffer.toString('utf-8');
    const lines = logicContent.split('\n').filter(line => line.trim());
    
    const nodes = new Set();
    const edges = [];
    
    // Parse TSV format
    lines.forEach((line, index) => {
      if (index === 0 && line.toLowerCase().includes('parent')) return; // Skip header
      
      const parts = line.split('\t');
      if (parts.length >= 5) {
        const [parent, child, isAnd, isPositive, stoichiometry] = parts;
        
        nodes.add(parent.trim());
        nodes.add(child.trim());
        
        edges.push({
          parent_uuid: parent.trim(),
          child_uuid: child.trim(),
          is_and: parseInt(isAnd) === 1,
          is_positive: parseInt(isPositive) === 1,
          stoichiometry: parseInt(stoichiometry) || 1
        });
      }
    });
    
    // Create node objects with simulated metadata
    const nodeObjects = Array.from(nodes).map(nodeId => ({
      uuid: nodeId,
      name: nodeId,  // Use the actual node ID as the name (e.g., "APC", "CTNNB1")
      reactome_id: `REACT:R-HSA-${Math.floor(Math.random() * 999999)}`,
      entity_type: ['protein', 'small_molecule', 'complex'][Math.floor(Math.random() * 3)],
      set_id: null
    }));
    
    const result = {
      nodes: nodeObjects,
      edges: edges,
      metadata: {
        parsed_at: new Date().toISOString(),
        parser_version: '0.1.0',
        total_files: Object.keys(req.files).length
      }
    };
    
    console.log(`Parsed network: ${nodeObjects.length} nodes, ${edges.length} edges`);
    res.json(result);
    
  } catch (error) {
    console.error('Error parsing network:', error);
    res.status(500).json({ error: error.message });
  }
});

/**
 * @swagger
 * /api/solve:
 *   post:
 *     summary: Solve steady-state network
 *     description: Calculate steady-state activities for network nodes with optional perturbations
 *     tags:
 *       - Analysis
 *     requestBody:
 *       required: true
 *       content:
 *         application/json:
 *           schema:
 *             type: object
 *             required:
 *               - network
 *             properties:
 *               network:
 *                 type: object
 *                 properties:
 *                   nodes:
 *                     type: array
 *                     items:
 *                       type: object
 *                       properties:
 *                         uuid:
 *                           type: string
 *                         name:
 *                           type: string
 *                   edges:
 *                     type: array
 *                     items:
 *                       type: object
 *                       properties:
 *                         parent_uuid:
 *                           type: string
 *                         child_uuid:
 *                           type: string
 *                         is_positive:
 *                           type: boolean
 *                         is_and:
 *                           type: boolean
 *                         stoichiometry:
 *                           type: number
 *               options:
 *                 type: object
 *                 properties:
 *                   method:
 *                     type: string
 *                     enum: [enhanced, perturbation]
 *                     description: Solving method to use
 *                   observations:
 *                     type: object
 *                     description: Node perturbations (uuid -> {activity, is_observed})
 *                     additionalProperties:
 *                       type: object
 *                       properties:
 *                         activity:
 *                           type: number
 *                           minimum: 0
 *                           maximum: 1
 *                           description: Activity level (0-1 scale, where 0.01 = 1x fold change)
 *                         is_observed:
 *                           type: boolean
 *     responses:
 *       200:
 *         description: Successfully solved network
 *         content:
 *           application/json:
 *             schema:
 *               type: object
 *               properties:
 *                 node_activities:
 *                   type: object
 *                   description: Final activity levels for each node (0-1 scale)
 *                   additionalProperties:
 *                     type: number
 *                 influence_scores:
 *                   type: object
 *                   description: Influence scores for each node
 *                   additionalProperties:
 *                     type: number
 *                 converged:
 *                   type: boolean
 *                   description: Whether the solution converged
 *                 iterations:
 *                   type: integer
 *                   description: Number of iterations performed
 *                 solution_quality:
 *                   type: number
 *                   description: Quality metric of the solution (0-1)
 *                 analysis_time:
 *                   type: number
 *                   description: Time taken for analysis in seconds
 *                 metadata:
 *                   type: object
 *                   properties:
 *                     solved_at:
 *                       type: string
 *                       format: date-time
 *                     method:
 *                       type: string
 *                     solver_version:
 *                       type: string
 *       400:
 *         description: Bad request - invalid network data
 *       500:
 *         description: Server error during solving
 */
app.post('/api/solve', async (req, res) => {
  try {
    console.log('Solving steady-state...');
    const { network, options = {} } = req.body;
    
    if (!network || !network.nodes) {
      return res.status(400).json({ error: 'Network data is required' });
    }
    
    // Simulate solving delay
    const solvingTime = 0.5 + Math.random() * 2;
    await new Promise(resolve => setTimeout(resolve, solvingTime * 1000));
    
    const nodeActivities = {};
    const influenceScores = {};
    
    // Check if we have observations/perturbations
    const observations = options.observations || {};
    const isPerturbationAnalysis = options.method === 'perturbation';
    
    console.log('Observations:', observations);
    console.log('Is perturbation analysis:', isPerturbationAnalysis);
    
    // Generate activities considering perturbations
    network.nodes.forEach(node => {
      let activity;
      
      // If this node has an observation/perturbation, use that value
      if (observations[node.uuid]) {
        // Activity is already in 0-1 scale from frontend (0 = 0x, 0.01 = 1x, 1 = 100x)
        activity = observations[node.uuid].activity;
        const foldChange = activity * 100;
        console.log(`Using observed activity for ${node.uuid}: ${activity} (${foldChange.toFixed(1)}x)`);
      } else {
        // Set baseline activity to 1x (0.01 in our 0-1 scale)
        activity = 0.01; // 1x baseline for all non-perturbed nodes
      }
      
      nodeActivities[node.uuid] = Math.min(1.0, Math.max(0.0, activity));
      
      // Influence scores based on connectivity (simulated)
      const baseInfluence = Math.random() * 50;
      const connectivityBonus = Math.random() * 100;
      influenceScores[node.uuid] = baseInfluence + connectivityBonus;
    });
    
    // If this is a perturbation analysis, simulate network propagation
    if (isPerturbationAnalysis && Object.keys(observations).length > 0) {
      console.log('Simulating network propagation...');
      simulateNetworkPropagation(network, nodeActivities, observations);
    }
    
    const result = {
      node_activities: nodeActivities,
      influence_scores: influenceScores,
      converged: true,
      iterations: 100 + Math.floor(Math.random() * 300),
      solution_quality: 0.8 + Math.random() * 0.15,
      analysis_time: solvingTime,
      metadata: {
        solved_at: new Date().toISOString(),
        method: options.method || 'enhanced',
        solver_version: '0.1.0-dev'
      }
    };
    
    console.log(`Solved network: ${Object.keys(nodeActivities).length} nodes`);
    res.json(result);
    
  } catch (error) {
    console.error('Error solving steady-state:', error);
    res.status(500).json({ error: error.message });
  }
});

// Serve main page
app.get('/', (req, res) => {
  res.sendFile(path.join(__dirname, 'index.html'));
});

// Start server
app.listen(PORT, () => {
  console.log(`🧬 DeltaSignal frontend server running on http://localhost:${PORT}`);
  console.log('API endpoints available:');
  console.log('  GET  /api/health');
  console.log('  POST /api/parse');  
  console.log('  POST /api/solve');
  console.log('');
  console.log(`📚 Swagger API Documentation: http://localhost:${PORT}/api-docs`);
});