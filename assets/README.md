# DeltaSignal Brand Assets

## Logo Files

- **`logo.svg`** - Main logo with full branding (light theme)
- **`logo-dark.svg`** - Dark theme version of the main logo
- **`logo-text.svg`** - Text-focused logo for headers and banners
- **`favicon.svg`** - Square favicon/icon version (32x32)

## Logo Design Elements

### Color Palette
- **Primary Blue**: `#3b82f6` → `#1d4ed8` (gradient)
- **Delta Green**: `#10b981` → `#059669` (gradient)  
- **Activation Green**: `#10b981`
- **Inhibition Red**: `#ef4444`
- **Network Gray**: `#64748b`
- **Warning Yellow**: `#fbbf24`

### Symbolism
- **Delta (Δ)**: Represents change/perturbation analysis
- **Network Nodes**: Biological pathway components
- **Arrows**: Activation/positive regulation
- **T-bars**: Inhibition/negative regulation
- **Connections**: Pathway relationships

### Typography
- **Font Family**: Inter (primary), system fallbacks
- **Main Title**: 24px, weight 700
- **Tagline**: 11px, weight 500
- **Version Badge**: 8px, weight 600

## Usage Guidelines

### ✅ Correct Usage
- Use original SVG files when possible for scalability
- Maintain aspect ratios when resizing
- Use on backgrounds that provide sufficient contrast
- Use favicon.svg for browser icons and small applications

### ❌ Avoid
- Don't stretch or distort the logo
- Don't change colors unless adapting for dark theme
- Don't separate the delta symbol from the network elements
- Don't use on backgrounds that reduce readability

### File Sizes
- **Main Logo**: ~3KB SVG
- **Favicon**: ~1KB SVG
- **Text Logo**: ~2KB SVG

## Integration Examples

### HTML
```html
<!-- Main logo in header -->
<img src="assets/logo.svg" alt="DeltaSignal" height="50">

<!-- Favicon -->
<link rel="icon" type="image/svg+xml" href="assets/favicon.svg">
```

### Markdown
```markdown
<!-- README header -->
<img src="assets/logo.svg" alt="DeltaSignal Logo" width="400">
```

### CSS
```css
.logo-image {
  height: 50px;
  width: auto;
  max-width: 200px;
}
```

## Dark Theme Support

The logo includes a dark theme variant (`logo-dark.svg`) with:
- Lighter text colors (`#e2e8f0`)
- Enhanced node contrast
- Adjusted opacity for better visibility on dark backgrounds

Use the appropriate version based on your theme context.