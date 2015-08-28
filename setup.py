from setuptools import setup, find_packages
import os

version = '0.1'

template_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk('templates') for f in filenames]

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
    # data_files = template_files,
    included_files="git",
    zip_safe = False,
    scripts = ['scripts/multiqc'],
    install_requires = [
        'jinja2',
        'simplejson',
        'pyyaml'
    ]
)