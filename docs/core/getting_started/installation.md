---
title: Installation
description: How to install MultiQC on your system
---

# Installing MultiQC

MultiQC is written in Python and can be installed in a number of ways.
Which method you should use depends on how you're using MultiQC and how familiar you are with the Python ecosystem.

If you're new to software packaging, this page can be a little overwhelming.
If in doubt, a general rule is:

- _Running MultiQC in a pipeline?_ &nbsp; Use [Docker](#docker) or [Singularity](#singularity).
- _Running MultiQC locally?_ &nbsp; Use [Pip](#pip) or [Conda](#conda).

:::tip{title="Installation cheat sheet"}

Know what you're doing with this kind of thing? Here's a quick reference:

<table><thead>
<tr><th>Method</th><th>Command</th></tr>
</thead><tbody><tr>
<td>Pip</td>
<td>

```bash
pip install multiqc
```

</td></tr><tr>
<td>Pip (dev version)</td>
<td>

```bash
pip install --upgrade --force-reinstall git+https://github.com/MultiQC/MultiQC.git
```

</td></tr><tr>
<td>Conda</td>
<td>

```bash
conda install multiqc
```

</td></tr><tr>
<td>Docker</td>
<td>

```bash
docker run -t -v `pwd`:`pwd` -w `pwd` multiqc/multiqc multiqc .
```

</td></tr></tbody></table>

:::

## Installing Python

MultiQC is written in Python and needs a Python installation to run.

To run MultiQC manually install, you'll typically install it into a local Python environment.
MultiQC requires Python version 3.8 or above.

### System Python

Python comes installed on most operating systems. You can install MultiQC directly here, but it is _not_ recommended. This often causes problems and it's a little risky to mess with it.

:::danger
If you find yourself prepending `sudo` to any MultiQC commands, take a step back and think about Python virtual environments / conda instead.
:::

### Python with Conda

To see if you have python installed, run `python --version` on the command line.
MultiQC needs Python version 3.8+.

We recommend using virtual environments to manage your Python installation.
Our favourite is _conda_, a cross-platform tool to manage Python environments.
You can installation instructions for Miniconda
[here](https://docs.conda.io/en/latest/miniconda.html).

Once conda is installed, you can create a Python environment with the following commands:

```bash
conda create --name py3.11 python=3.11
conda activate py3.11
```

You'll want to add the `conda activate py3.11` line to your `.bashrc` file,
so that the environment is loaded every time you load the terminal.

### Using a specific python interpreter

If you prefer, you can also run MultiQC with a specific python interpreter.
The command line usage and flags are then exactly the same as if you ran just `multiqc`.

For example:

```bash
python -m multiqc .
python3 -m multiqc .
~/my_env/bin/python -m multiqc .
```

## Installing MultiQC locally

There are a few different ways to install MultiQC into your local Python environment:

### Conda

MultiQC is available on [BioConda](https://bioconda.github.io/).

```bash
conda install multiqc
```

:::info

The order of conda channels is important!
Please make sure that you have [configured your conda channels](https://bioconda.github.io/#usage) prior to installing anything with BioConda:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

:::

:::warning

In the past we used `-c bioconda` in the installation command, but this is no longer the correct usage.
Doing so will likely cause weird stuff to happen (such as only being able to install very old versions).

:::

### Pip / PyPI

This is the easiest way to install MultiQC. `pip` is the package manager for the Python Package Manager. It comes bundled with recent versions of Python,
otherwise you can find installation instructions [here](https://pip.pypa.io/en/stable/installation/).

You can install MultiQC from [PyPI](https://pypi.python.org/pypi/multiqc) as follows:

```bash
pip install multiqc
```

Use the `--upgrade` flag to update to the latest version.

If you would like the development version, the command is:

```bash
pip install git+https://github.com/MultiQC/MultiQC.git
```

To update the dev version between releases, use `--upgrade --force-reinstall`. This is needed as the version number isn't changing.

If you have problems with read-only directories, you can install to
your home directory with the `--user` parameter:

```bash
pip install --user multiqc
```

### Spack

MultiQC [is available on spack](https://packages.spack.io/package.html?name=py-multiqc) as `py-multiqc`:

```bash
spack install py-multiqc
```

### FreeBSD

If you're using the [FreeBSD](https://www.freebsd.org/) operating system, you can install MultiQC via [FreeBSD ports](https://www.freebsd.org/ports/):

```bash
pkg install py39-multiqc
```

This will install a prebuilt binary using only highly-portable
optimizations.

FreeBSD ports can also be built and installed from source:

```bash
cd /usr/ports/biology/py-multiqc
make install
```

To report issues with a FreeBSD port, please submit a PR on the
[FreeBSD bug reports page](https://www.freebsd.org/support/bugreports.html).

### Cloning the repository

If you'd rather not use either of these tools, you can clone the code and install the code yourself:

```bash
git clone https://github.com/MultiQC/MultiQC.git
cd MultiQC
pip install .
```

This will fetch the latest development code. To update to the latest changes, use `git pull`.

`git` not installed? No problem - just download the flat files:

```bash
curl -LOk https://github.com/MultiQC/MultiQC/archive/main.zip
unzip main.zip
cd MultiQC-main
pip install .
```

### Nix

If you're using the [nix package manager](https://nixos.org/download.html#download-nixm) with [flakes](https://nixos.wiki/wiki/Flakes) enabled, you can
run `nix develop`in the cloned MultiQC repository to enter a shell
with required dependencies. To build MultiQC, run `nix build`.

## MultiQC container images

### Docker

A Docker container is provided on Docker Hub called [`multiqc/multiqc`](https://hub.docker.com/r/ewels/multiqc/).
It's based on an `python-slim` base image to give the smallest image size possible.

To use, call the `docker run` with your current working directory mounted as a volume and working directory. Then just specify the MultiQC command at the end as usual:

```bash
docker run -t -v `pwd`:`pwd` -w `pwd` multiqc/multiqc multiqc .
```

- `-t`: Runs docker with a pseudo-tty, for nice terminal colours
- `-v`: Mounts the current working directory into the container
- `-w`: Sets the working directory in the container as your local working directory

You can specify additional MultiQC parameters as normal at the end of the command:

```bash
docker run -t -v `pwd`:`pwd` -w `pwd` multiqc/multiqc multiqc . --title "My amazing report" -b "This was made with docker"
```

By default, docker will use the `:latest` tag. For MultiQC, this is set to be the most recent release.
To use the most recent development code, use `multiqc/multiqc:dev`.
You can also specify specific versions, eg: `multiqc/multiqc:v1.20`.

Note that all files on the command line (eg. config files) must also be mounted in the docker container to be accessible.
For more help, look into [the Docker documentation](https://docs.docker.com/engine/reference/commandline/run/).

:::warning{title="Docker image name change"}
The docker image used to be called `ewels/multiqc`.
All releases prior to MultiQC v1.19 can be found at [ewels/multiqc](https://hub.docker.com/r/ewels/multiqc/)
and everything from v1.20 onwards can be found at [multiqc/multiqc](https://hub.docker.com/r/multiqc/multiqc/).
:::

:::tip{title="Tip: Docker bash alias"}

The docker command above is a little verbose, so if you are using this a lot it may be worth adding the following bash alias to your `~/.bashrc` file:

```bash
alias multiqc="docker run -tv `pwd`:`pwd` -w `pwd` multiqc/multiqc multiqc"
```

Once applied (first log out and in again) you can then just use the `multiqc` command instead:

```bash
multiqc .
```

:::

:::note{title="Compute architectures"}

These docker images are [multi-platform images](https://docs.docker.com/build/building/multi-platform/) – each build contains two digests, one for `linux/amd64` and one for `linux/arm64`.

Generally, the Docker client should be clever enough to pull the digest appropriate for your local compute architecture.
However, if you wish you can force it with the `--platform` flag.

```bash
docker pull --platform linux/arm64 multiqc/multiqc:latest
```

:::

### GitHub Packages

If you prefer, the Docker image above is also available from [GitHub packages](https://github.com/MultiQC/MultiQC/pkgs/container/multiqc).
Usage is identical, the only difference is that the URI has a `ghcr.io/` prefix:

```bash
docker pull ghcr.io/multiqc/multiqc
```

This image was also renamed, versions up to v1.19 can be found at [`ghcr.io/ewels/multiqc`](https://github.com/users/ewels/packages/container/package/multiqc).

### Singularity

To build a singularity container image from the docker image, use the following command: _(change `1.20` to the current MultiQC version)_

```bash
singularity build multiqc-1.20.sif docker://multiqc/multiqc:v1.20
```

Then, use `singularity run` to run the image with the normal MultiQC arguments:

```bash
singularity run multiqc-1.20.sif my_results/ --title "Report made using Singularity"
```

:::info{title="Import errors with Singularity"}

Sometimes, Singularity can be over-ambitious with sharing file paths which can result in the Python environment in your local system interacting with Python inside the image. This can give rise to `ImportError` errors for `numpy` and other packages.

The giveaway for when this is the problem is that traceback will list python package paths which are on your system and look different that of MultiQC inside the container (eg. `/usr/lib/python3.8/site-packages/multiqc/`).

To fix this, run the command `export PYTHONNOUSERSITE=1` before running MultiQC. This variable [tells Python](https://docs.python.org/3/using/cmdline.html#envvar-PYTHONNOUSERSITE) not to add site-packages to the system path when loading, which should avoid the conflicts.
:::

:::tip

If you prefer, you can download a pre-built Singularity image from BioContainers, see below.

:::

### BioContainers

[BioContainers](https://biocontainers.pro/) is a project that automatically builds Docker and Singularity container images from [BioConda](https://bioconda.github.io/). The images are less fine-tuned for MultiQC so tend to have a larger filesize, but they should work well and are convenient.

To see available images, visit the BioContainers [registry page for MultiQC](https://biocontainers.pro/tools/multiqc).

## Using MultiQC in a Python script

You can import and run MultiQC from within a Python script, using
the `multiqc.run()` function as follows:

```python
import multiqc
multiqc.run("/path/to/dir")
```

More development of interactive usage is planned for the future.
Currently you can't do a lot more than just running MultiQC.

## Galaxy

### On the main Galaxy instance

The easiest and fast manner to use MultiQC is to use the [usegalaxy.org](https://usegalaxy.org/) main Galaxy instance where you will find [MultiQC Galaxy tool](https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fengineson%2Fmultiqc%2Fmultiqc%2F1.0.0.0&version=1.0.0.0&__identifer=2sjdq8d9r3l) under the _NGS: QC and manipualtion_ tool panel section.

### On your instance

You can install MultiQC on your own Galaxy instance through your Galaxy admin space, searching on the [main Toolshed](https://toolshed.g2.bx.psu.edu/) for the [MultiQC repository](https://toolshed.g2.bx.psu.edu/view/iuc/multiqc/3bad335ccea9) available under the _visualization_, _statistics_ and _Fastq Manipulation_ sections.

## Environment modules

Many people using MultiQC will be working on a HPC environment.
Every server / cluster is different, and you're probably best off asking
your friendly sysadmin to install MultiQC for you. However, with that
in mind, here are a few general tips for installing MultiQC into an
environment module system:

MultiQC comes in two parts - the `multiqc` python package and the
`multiqc` executable script. The former must be available in `$PYTHONPATH`
and the script must be available on the `$PATH`.

A typical installation procedure with an environment module Python install
might look like this: _(Note that `$PYTHONPATH` must be defined before `pip` installation.)_

```bash
VERSION=0.7
INST=/path/to/software/multiqc/$VERSION
module load python/3.11
mkdir $INST
export PYTHONPATH=$INST/lib/python2.7/site-packages
pip install --install-option="--prefix=$INST" multiqc
```

Once installed, you'll need to create an environment module file.
Again, these vary between systems a lot, but here's an example:

```bash
#%Module1.0#####################################################################
##
## MultiQC
##

set components [ file split [ module-info name ] ]
set version [ lindex $components 1 ]
set modroot /path/to/software/multiqc/$version

proc ModulesHelp { } {
    global version modroot
    puts stderr "\tMultiQC - use MultiQC $version"
    puts stderr "\n\tVersion $version\n"
}
module-whatis   "Loads MultiQC environment."

# load required modules
module load python/3.11

# only one version at a time
conflict multiqc

# Make the directories available
prepend-path    PATH        $modroot/bin
prepend-path	PYTHONPATH	$modroot/lib/python3.11/site-packages
```
