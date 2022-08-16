# ASQSQIS: Simulating open quantum systems with QuTiP

Presenters:
- **Alex Pitchford**
- **Simon Cross**

Institutions:
- **[RIKEN](https://www.riken.jp/en/)**,
  **[Theoretical Quantum Physics Laboratory](https://dml.riken.jp/)**
- **[Aberystwyth University](https://www.aber.ac.uk/en/about-us/faculties/business-physical-sciences/)**
- **[QuTiP](https://qutip.org)**

## Hello!

If you're attending the African School on Quantum Simulation and Quantum
Information Science (ASQSQIS) in Kigali in September, and are looking for
instructions for setting up QuTiP, you have found them.

To get set up you'll need to do the following:

1. Obtain a copy of the contents of this GitHub repository.

2. Install the required Python packages by following the instruction in one of
   the options below.

3. Check that everything is working by running the `smoke-test.ipynb`
   Jupyter notebook.

There are instructions for what to do if you get stuck further down.

Don't panic, have fun and see you in Kigali!


## Download the repository

If you're familiar with GitHub, just clone this repository.

Otherwise, you can download the latest contents of the repository
[here](https://github.com/hodgestar/qutip-asqsqis-2022/archive/refs/heads/main.zip).


## Installing locally with pip

1. If you don't have Python 3.9 or later installed, install it however
   you usually do. If you are completely new to Python, follow the
   [Getting Started](https://www.python.org/about/gettingstarted/) guide.

2. If you don't have `pip` and `virtualenv` installed, install them by
   following the [instructions](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/). Once you've made it passed these
   first two steps, the hardest part is done.

3. Create a virtual environment by running
   `python3 -m venv -m qutip-summer-school-env`.

4. Activate the environment you have just created by running
   `source qutip-summer-school-env/bin/activate` (on Linux or Mac OS) or
   `.\qutip-summer-school-env\Scripts\activate` (on Windows).

5. Install QuTiP and the other required libraries by running
   `pip install -r requirements.txt`.

6. Open the Jupyter notebook by running `jupyter notebook`.


## Installing locally with conda

1. Install miniconda by following the [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). This is the hardest step.

2. Create a conda environment named `qutip-summer-school-env` by running
   `conda env create --file environment.yml`.

3. Activate the environment you have just created by running
   `conda activate qutip-summer-school-env`.

4. Open the Jupyter notebook by running `jupyter notebook`.


## Running remotely with Binder

I'd highly recommend installing and running Python locally for the summer
school, but if that really isn't possible, you are welcome to try running
Jupyter and the necessary dependencies online using either Binder, by clicking
the button below:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hodgestar/qutip-asqsqis-2022/main)

Or by using [Google Colab](https://colab.research.google.com/github/hodgestar/qutip-asqsqis-2022/blob/main/).

If you are using Google Colab, you'll need to first install QuTiP by running
the command `!pip install 'qutip[full]'` in a cell by itself.


## Testing your installation

Once you have the dependencies installed and Jupyter notebook open, you should
open the notebook named `smoke-test.ipynb` and run it.

The notebook tests the basic plotting and QuTiP features we will need for the
summer school and contains instructions for how to tell if it ran correctly.

If the test notebook ran correctly, you're good to go!


## What to do if you're stuck

If you're having trouble setting up, you're welcome to create an issue in this
repository [here](https://github.com/hodgestar/qutip-asqsqis-2022/issues/new)
and we'll do our best to help.

If you're looking for help with QuTiP, you can message the
[QuTiP Google Group](https://groups.google.com/g/qutip), and we or someone
else will reply there.

You're also welcome to contact us via the summer school organizers.
