# Installing MultiQC

## Installing Python
To see if you have python installed, run `python --version` on the command line.
If you see version 2.7+, 3.4+ or 3.5+ then you can skip this step.

We recommend using virtual environments to manage your Python installation.
Our favourite is Anaconda, a cross-platform tool to manage Python environments.
You can installation instructions for Anaconda
[here](http://conda.pydata.org/docs/install/quick.html).

Once Anaconda is installed, you can create an environment with the
following commands:

```
conda create --name py2.7 python=2.7
source activate py2.7
# Windows: activate py2.7
```

You'll want to add the `source activate py2.7` line to your `.bashrc` file so
that the environment is loaded every time you load the terminal.

## Installation with `pip`
This is the easiest way to install MultiQC. `pip` is the package manager for
the Python Package Manager. It comes bundled with recent versions of Python,
otherwise you can find installation instructions [here](http://pip.readthedocs.org/en/stable/installing/).

You can now install MultiQC from
[PyPI](https://pypi.python.org/pypi/multiqc) as follows:
```
pip install multiqc
```

If you would like the development version, the command is:
```
pip install git+https://github.com/ewels/MultiQC.git
```

## Installing with `conda`
If you're using `conda` as described above, you can install MultiQC from the `bioconda`
channel as follows:
```
conda install -c bioconda multiqc
```

## Manual installation
If you'd rather not use either of these tools, you can clone the code and install the code yourself:
```
git clone https://github.com/ewels/MultiQC.git
cd MultiQC
python setup.py install
```

`git` not installed? No problem - just download the flat files:
```
curl -LOk https://github.com/ewels/MultiQC/archive/master.zip
unzip master.zip
cd MultiQC-master
python setup.py install
```

## Updating MultiQC
You can update MultiQC from [PyPI](https://pypi.python.org/pypi/multiqc)
at any time by running the following command:
```
pip --update multiqc
```

To update the development version, use:
```
pip install --force git+https://github.com/ewels/MultiQC.git
```

If you cloned the `git` repo, just pull the latest changes and install:
```
cd MultiQC
git pull
python setup.py install
```

If you downloaded the flat files, just repeat the installation procedure.

