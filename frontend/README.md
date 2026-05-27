# 🧬 DeltaSignal Frontend

Web-based visualization interface for DeltaSignal pathway analysis.

## Features

- **Dual Visualization**: PathwayBrowser-style biological view + Cytoscape.js network topology
- **Interactive Analysis**: Upload networks, solve steady-state, visualize results
- **Real-time Updates**: Activity levels and influence scores displayed on both views
- **Export Capabilities**: Download results in JSON/CSV formats
- **Responsive Design**: Works on desktop and tablet devices

## Technology Stack

- **Frontend**: Vanilla JavaScript (ES6+), HTML5, CSS3
- **Visualization**: 
  - Cytoscape.js for network topology
  - D3.js for pathway browser
- **Build System**: Vite
- **Backend API**: Express.js (development server)

## Quick Start

### Prerequisites

- Node.js 16+ and npm
- Modern web browser with ES6 support

### Installation

```bash
# Install dependencies
cd frontend
npm install

# Start development server
npm run dev

# Or start with backend API simulation
npm run start
```

The application will be available at http://localhost:3000

### Development vs Production

- **Development mode** (`npm run dev`): Uses Vite dev server with hot reload
- **API simulation mode** (`npm start`): Includes Express backend for testing without Julia

## Usage

1. **Upload Network Files**
   - Click "Upload Network" button
   - Provide required files:
     - Logic Network (TSV): Network topology
     - UUID Mapping (TSV): Maps UUIDs to Reactome IDs
   - Optional files:
     - Set Mappings (TSV): Complex/set expansions
     - Observations (CSV): Experimental data

2. **Analyze Network**
   - Click "Solve" to run steady-state analysis
   - View results in both pathway and network views
   - Activities shown as color coding and percentages
   - Observed nodes highlighted with special borders

3. **Export Results**
   - Click "Export" to download analysis results
   - Available formats: JSON, CSV

## File Formats

### Logic Network (TSV)
```
parent-001	child-001	1	1	1
parent-002	child-001	1	1	1
parent-003	child-002	0	-1	2
```
Columns: Parent UUID | Child UUID | AND/OR (1/0) | Pos/Neg (1/-1) | Stoichiometry

### UUID Mapping (TSV)
```
parent-001	REACT:R-HSA-123456	protein	set-001
parent-002	REACT:R-HSA-123457	protein	set-001
parent-003	REACT:R-HSA-123458	small_molecule	
```
Columns: Network UUID | Reactome DB ID | Entity Type | Set ID (optional)

### Observations (CSV)
```
node_uuid,activity,confidence
parent-001,75.0,0.9
parent-003,50.0,0.8
```

## API Endpoints

When running with backend simulation (`npm start`):

- `GET /api/health` - Health check
- `POST /api/parse` - Parse network files
- `POST /api/solve` - Solve steady-state
- `POST /api/export` - Export results
- `GET /api/pathways` - Get available pathways

## Project Structure

```
frontend/
├── src/
│   ├── components/           # UI components
│   │   ├── NetworkVisualizer.js    # Cytoscape.js network view
│   │   ├── PathwayBrowser.js       # D3.js pathway view
│   │   ├── FileUploadHandler.js    # File upload modal
│   │   └── ResultsPanel.js         # Analysis results display
│   ├── styles/
│   │   └── main.css         # Main stylesheet
│   ├── utils/
│   │   └── APIClient.js     # Backend communication
│   └── app.js               # Main application
├── index.html               # Entry point
├── vite.config.js          # Build configuration
├── server.js               # Development API server
└── package.json            # Dependencies
```

## Customization

### Adding New Layouts

To add a new network layout:

1. Install Cytoscape.js extension:
   ```bash
   npm install cytoscape-your-layout
   ```

2. Import and register in `NetworkVisualizer.js`:
   ```javascript
   const yourLayout = (await import('cytoscape-your-layout')).default;
   cytoscape.use(yourLayout);
   ```

3. Add layout option in `applyLayout()` method

### Styling Nodes and Edges

Modify the Cytoscape styles in `NetworkVisualizer.js`:

```javascript
{
    selector: '.your-class',
    style: {
        'background-color': '#custom-color',
        'border-width': 3
    }
}
```

### Custom Pathway Rendering

Extend `PathwayBrowser.js` to support additional entity types or rendering styles.

## Build and Deploy

### Development Build
```bash
npm run build
```

### Production Deployment

1. Build the frontend:
   ```bash
   npm run build
   ```

2. Deploy `dist/` folder to your web server

3. Configure backend API endpoints in `APIClient.js`

## Integration with DeltaSignal CLI

To connect with the Julia backend:

1. Start the DeltaSignal API server:
   ```bash
   julia src/api/server.jl --port 8080
   ```

2. Update `APIClient.js` baseURL to point to Julia server

3. Remove API simulation and use real endpoints

## Browser Support

- Chrome 80+
- Firefox 75+
- Safari 13+
- Edge 80+

ES6+ features used: modules, async/await, destructuring, template literals.

## Contributing

1. Follow existing code style and structure
2. Add JSDoc comments for new functions
3. Test with different network sizes and formats
4. Ensure responsive design on various screen sizes

## License

Same as parent DeltaSignal project.

---

For backend integration and Julia CLI usage, see the main [DeltaSignal README](../README.md).