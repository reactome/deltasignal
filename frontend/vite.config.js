import { defineConfig } from 'vite';
import legacy from '@vitejs/plugin-legacy';

export default defineConfig({
  plugins: [
    legacy({
      targets: ['defaults', 'not IE 11']
    })
  ],
  publicDir: 'public',
  server: {
    port: 3000,
    host: '0.0.0.0',
    open: false,  // Don't auto-open browser in Docker
    proxy: {
      '/api': {
        target: 'http://julia-api:8080',  // Use Docker service name
        changeOrigin: true,
        secure: false
      }
    }
  },
  build: {
    outDir: 'dist',
    assetsDir: 'assets',
    sourcemap: true,
    rollupOptions: {
      output: {
        manualChunks: {
          cytoscape: ['cytoscape', 'cytoscape-cose-bilkent', 'cytoscape-dagre'],
          d3: ['d3']
        }
      }
    }
  },
  optimizeDeps: {
    include: ['cytoscape', 'cytoscape-cose-bilkent', 'cytoscape-dagre', 'd3']
  }
});