# MultiQC Brite Template

A child theme of the default MultiQC template using the [Bootswatch Brite](https://bootswatch.com/brite/) theme.

## Theme Files

The Brite theme SCSS files (`_variables.scss` and `_bootswatch.scss`) are sourced from:

- **Source**: https://bootswatch.com/5/brite/
- **Repository**: https://github.com/thomaspark/bootswatch
- **License**: MIT

## Usage

```bash
multiqc . --template brite
```

Or in your MultiQC configuration file:

```yaml
template: brite
```

## Building

```bash
npm install
npm run build
```
