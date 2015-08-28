from setuptools import setup, find_packages

version = '0.1'

setup(
    name = 'multiqc',
    version = version,
    author = 'Phil Ewels',
    author_email = 'phil.ewels@scilifelab.se',
    description = "Create aggregate analyses reports for many samples",
    long_description = 'A modular tool to aggregate results from bioinformatics analyses across many samples into a single report',
    keywords = 'bioinformatics',
    url = 'https://github.com/ewels/MultiQC',
    license = 'MIT',
    packages = find_packages(exclude=['ez_setup', 'examples', 'tests']),
    package_data = {'': ['*.*']},
    include_package_data = True,
    zip_safe = False,
    scripts = ['scripts/multiqc'],
    install_requires = [
        'jinja2',
        'simplejson',
        'pyyaml'
    ]
)
