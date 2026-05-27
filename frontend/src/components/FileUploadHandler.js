// File Upload Handler for network files
export default class FileUploadHandler {
    constructor() {
        this.modal = document.getElementById('upload-modal');
        this.files = {};
        this.onFilesProcessed = null;
        
        this.init();
    }
    
    init() {
        this.setupEventListeners();
        console.log('FileUploadHandler initialized');
    }
    
    setupEventListeners() {
        // Modal controls
        document.getElementById('modal-close').addEventListener('click', () => {
            this.hideModal();
        });
        
        document.getElementById('upload-cancel').addEventListener('click', () => {
            this.hideModal();
        });
        
        document.getElementById('upload-confirm').addEventListener('click', () => {
            this.processFiles();
        });
        
        // Example files buttons
        document.getElementById('load-examples-btn').addEventListener('click', () => {
            this.loadExampleFiles();
        });

        document.getElementById('load-reaction-example-btn').addEventListener('click', () => {
            this.loadReactionExample();
        });
        
        // File input listeners
        const fileInputs = ['logic-file', 'uuid-file', 'set-file', 'observations-file'];
        fileInputs.forEach(inputId => {
            document.getElementById(inputId).addEventListener('change', (e) => {
                this.handleFileSelect(inputId, e.target.files[0]);
            });
        });
        
        // Close modal on backdrop click
        this.modal.addEventListener('click', (e) => {
            if (e.target === this.modal) {
                this.hideModal();
            }
        });
    }
    
    showModal() {
        console.log('showModal called!');
        console.log('Modal element:', this.modal);
        this.files = {};
        this.updateConfirmButton();
        this.clearFileInputs();
        this.modal.classList.add('show');
        console.log('Modal classes after show:', this.modal.classList.toString());
    }
    
    hideModal() {
        this.modal.classList.remove('show');
    }
    
    handleFileSelect(inputType, file) {
        if (file) {
            this.files[inputType] = file;
            console.log(`File selected for ${inputType}:`, file.name);
        } else {
            delete this.files[inputType];
        }
        
        this.updateConfirmButton();
    }
    
    updateConfirmButton() {
        const confirmBtn = document.getElementById('upload-confirm');
        
        // Required files: logic-file and uuid-file
        const hasRequiredFiles = this.files['logic-file'] && this.files['uuid-file'];
        
        confirmBtn.disabled = !hasRequiredFiles;
    }
    
    clearFileInputs() {
        const fileInputs = ['logic-file', 'uuid-file', 'set-file', 'observations-file'];
        fileInputs.forEach(inputId => {
            document.getElementById(inputId).value = '';
        });
    }
    
    async processFiles() {
        if (!this.files['logic-file'] || !this.files['uuid-file']) {
            alert('Logic network and UUID mapping files are required');
            return;
        }
        
        try {
            this.hideModal();
            
            // Call callback with raw files (not content)
            if (this.onFilesProcessed) {
                this.onFilesProcessed(this.files);
            }
            
        } catch (error) {
            console.error('Error processing files:', error);
            alert('Error processing files: ' + error.message);
        }
    }
    
    readFileContent(file) {
        return new Promise((resolve, reject) => {
            const reader = new FileReader();
            
            reader.onload = (e) => {
                resolve({
                    name: file.name,
                    content: e.target.result,
                    type: file.type,
                    size: file.size
                });
            };
            
            reader.onerror = () => {
                reject(new Error(`Failed to read file: ${file.name}`));
            };
            
            reader.readAsText(file);
        });
    }
    
    // Utility method to validate file formats
    validateFiles(files) {
        const errors = [];
        
        // Check logic network file
        if (files['logic-file']) {
            const logicFile = files['logic-file'];
            if (!logicFile.name.match(/\.(tsv|txt)$/i)) {
                errors.push('Logic network file must be TSV format');
            }
        }
        
        // Check UUID mapping file
        if (files['uuid-file']) {
            const uuidFile = files['uuid-file'];
            if (!uuidFile.name.match(/\.(tsv|txt)$/i)) {
                errors.push('UUID mapping file must be TSV format');
            }
        }
        
        // Check observations file
        if (files['observations-file']) {
            const obsFile = files['observations-file'];
            if (!obsFile.name.match(/\.csv$/i)) {
                errors.push('Observations file must be CSV format');
            }
        }
        
        return errors;
    }
    
    async loadExampleFiles() {
        console.log('Loading example files...');
        
        try {
            // Fetch example files from server
            const exampleFiles = {
                'logic-file': { url: '/examples/sample_logic_network.tsv', name: 'sample_logic_network.tsv' },
                'uuid-file': { url: '/examples/sample_uuid_mapping.tsv', name: 'sample_uuid_mapping.tsv' },
                'set-file': { url: '/examples/sample_set_mappings.tsv', name: 'sample_set_mappings.tsv' },
                'observations-file': { url: '/examples/sample_observations.csv', name: 'sample_observations.csv' }
            };
            
            for (const [inputType, fileInfo] of Object.entries(exampleFiles)) {
                try {
                    const response = await fetch(fileInfo.url);
                    if (response.ok) {
                        const content = await response.text();
                        
                        // Create a File object from the content
                        const blob = new Blob([content], { 
                            type: inputType === 'observations-file' ? 'text/csv' : 'text/tab-separated-values' 
                        });
                        const file = new File([blob], fileInfo.name, { 
                            type: blob.type,
                            lastModified: Date.now()
                        });
                        
                        this.files[inputType] = file;
                        this.updateFileStatus(inputType, `✅ ${fileInfo.name} loaded`);
                        console.log(`Loaded example file: ${fileInfo.name}`);
                    } else {
                        this.updateFileStatus(inputType, '⚠️ Example not available');
                    }
                } catch (error) {
                    console.warn(`Could not load example file ${fileInfo.name}:`, error);
                    this.updateFileStatus(inputType, '⚠️ Example not available');
                }
            }
            
            this.updateConfirmButton();
            
        } catch (error) {
            console.error('Error loading example files:', error);
            alert('Could not load example files. Please upload your own files.');
        }
    }
    
    async loadReactionExample() {
        console.log('Loading reaction network example...');

        try {
            // Fetch the reaction network example
            const response = await fetch('/examples/glycolysis_reaction_network.json');
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }

            const reactionNetwork = await response.json();

            // Close the modal
            this.hideModal();

            // Process the reaction network directly
            if (this.onFilesProcessed) {
                // Pass the reaction network as a special file type
                this.onFilesProcessed({
                    reaction_network: reactionNetwork
                });
            }

        } catch (error) {
            console.error('Error loading reaction network example:', error);
            alert('Could not load reaction network example. Please check if the file exists.');
        }
    }

    updateFileStatus(inputType, message) {
        const statusElement = document.getElementById(`${inputType}-status`);
        if (statusElement) {
            statusElement.textContent = message;
            statusElement.style.fontSize = '0.75rem';
            statusElement.style.color = message.includes('✅') ? '#059669' : '#d97706';
            statusElement.style.marginTop = '0.25rem';
        }
    }

    // Get processed files in format expected by API
    getProcessedFiles() {
        return this.files;
    }
}