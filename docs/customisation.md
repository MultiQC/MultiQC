# Customising Reports
MultiQC offers a few ways to customise reports to easily add your own
branding and some additional report-level information. These features
are primarily designed for core genomics facilities.

## Report Logo
To add your own custom logo to reports, you can add the following
three lines to your MultiQC configuration file:

```yaml
custom_logo: '/abs/path/to/logo.png'
custom_logo_url: 'https://www.example.com'
custom_logo_title: 'Our Institute Name'
```
Only `custom_logo` is needed. The URL will make the logo open up
a new web browser tab with your address and the title sets the mouse
hover title text.