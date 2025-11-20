Getting started with Astroclimes
=====

.. _installation:
Installation
-----
To install Astroclimes, simply clone the `GitHub <https://github.com/MarceloAron/Astroclimes>`_ repository using:

>>> git clone https://github.com/MarceloAron/Astroclimes.git

The GitHub installation comes with all of the necessary Python scripts and some exemplary auxiliary files, listed below: 

1. One model emission line spectrum from `SKYCALC <https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC>`_.
2. A few PHOENIX model stellar spectra from the `Göttingen Spectral Library <https://phoenix.astro.physik.uni-goettingen.de/>`_. 
3. Some GGG2020 atmospheric profiles computed with `ginput <https://ginput.readthedocs.io/en/latest/>`_. For details on how to generate these atmospheric profiles yourself, check :doc:`generateGGG2020profs`. 
4. Publicly available `CARMENES  <https://carmenes.caha.es/>`_ observations of Tau Boo on the night of March 26th 2018, for which the atmospheric profiles were computed. 
5. The files with the coefficients necessary to calculate the collision-induced absorption (CIA) for :math:`O_2` collisions with other atmospheric particles. If needed, other coefficients are available on the `HITRAN database <https://www.hitran.org/cia/>`_. 
6. Some model planetary emission spectra computed with GENESIS (`Gandhi & Madhusudhan 2017 <https://ui.adsabs.harvard.edu/abs/2017MNRAS.472.2334G/abstract>`_) with differing water abundances.

Due to their large file size, the molecular cross_sections for line-by-line absorption are not included in the GitHub repository, but can be downloaded from `Google Drive <https://drive.google.com/drive/folders/1Awtc4UlRM1GR_sOg3eOVMGlF5WyD3Wue?usp=sharing>`_. Cross-sections for other molecules not included may be provided upon request. 

.. _prerequisites:
Pre-requisites
-----
Before running Astroclimes, make sure you have all of the necessary Python packages installed, which are listed below. The minimum Python version recommended is 3.9. This can be done with ``conda`` environments. If you have Anaconda installed in your machine, you may do so by running:

>>> conda create --name <my_env> python=3.9

To install each package, activate the environment running ``conda activate <my_env>`` and then:

>>> conda install <package>

or

>>> conda install conda-forge::<package>

Python packages pre-requisites:

  - astropy>=5.3.4
  - corner>=2.2.1
  - emcee>=3.1.6
  - h5py>=3.12.1
  - matplotlib>=3.9.2
  - numpy>=1.26.4
  - pandas>=2.2.2
  - scipy>=1.13.1
  - tqdm>=4.67.1

Once you have installed all of the necessary scripts and auxiliary files, you may check out :doc:`setup`.