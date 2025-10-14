# MultiQC Sketchy Template

A child theme of the default MultiQC template using the [Bootswatch Sketchy](https://bootswatch.com/sketchy/) theme.

## Theme Files

The Sketchy theme SCSS files (`_variables.scss` and `_bootswatch.scss`) are sourced from:

- **Source**: https://bootswatch.com/5/sketchy/
- **Repository**: https://github.com/thomaspark/bootswatch
- **License**: MIT

## Usage

```bash
multiqc . --template sketchy
```

Or in your MultiQC configuration file:

```yaml
template: sketchy
```

## Building

```bash
npm install
npm run build
```
