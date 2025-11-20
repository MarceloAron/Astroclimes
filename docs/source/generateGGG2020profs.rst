.. _genGGG2020profs:
Generating GGG2020 atmospheric profiles
=====
The GGG pipeline is used by the Total Carbon Column Observing Network (TCCON) to process their spectra and compute the total column abundances of the relevant gases. The current version of the pipeline is GGG2020 and one of its applications is an algorithm that allows the user to compute a priori profiles for several gases. This algorithm is called `ginput` and is publicly available on `GitHub <https://github.com/TCCON/py-ginput>`_.

.. _installginput:

Installing ginput
------------
To install and learn how to use ginput, see their `documentation <https://ginput.readthedocs.io/en/latest/>`_. 
Astroclimes comes with its own installation of ginput found in the directory ``py-ginput/``, which should work by simply copying it to your local repository.

.. _donwloadGEOSFPfiles:

Downloading the necessary auxiliary files
------------

ginput relies on meteorological data to compute atmospheric profiles. These come from the Goddard Earth Observing System Forward Product (GEOS FP), 
provided by the Global Modelling and Assimilation Office (GMAO) at NASA Goddard Space Flight Center through the `online data portal <https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/>`_ in the NASA Center for Climate Simulation (NCCS). The auxiliary files required are the 3D meteorological (``inst3_3d_asm_Nv``), 3D chemical (``inst3_3d_chm_Nv``) and 2D surface (``inst3_2d_asm_Nx``) information files. You can either download them manually from the online data portal, or download them automatically with the help of the ``get_GEOS_aux_files.py`` script which can be found in the ``py-ginput/`` directory. To do so, simply change the path assigned to the ``home_directory`` variable and the path to the data on the ``spec_file_list`` variable and run:

>>> python get_GEOS_aux_files.py

This script will create, if not already existent, the directory ``GEOSFP/`` inside ``Astroclimes/data/`` and the subdirectories ``GEOSFP/Nv/`` and ``GEOSFP/Nx/``, where the auxiliary files will be stored.
At this stage, the spectroscopic data is only used to get the observations dates. You may wish to run this script before obtaining the observational data, as long as you specify the dates in a different way.

.. _runginput:

Running ginput
----------------
Once you have the GEOS-FP auxiliary files, the GGG2020 atmospheric profiles can be generated with the help of the ``generate_GGG20202_atm_profs.py`` script, also in the ``py-ginput/`` directory.
Again, simply change the path assigned to the ``home_directory`` variable and run:

>>> python generate_GGG20202_atm_profs.py

.. note::
    This script will generate the ``.mod``, ``.vmr`` and ``map`` files (the latter in ``.txt`` format) and store them in their appropriate subdirectories inside ``Astroclimes/atmosphere_profiles/GGG2020/fp/al/``.
    For some reason, the subdirectory for the ``map`` files is not automatically created, so you must make sure the directory ``Astroclimes/atmosphere_profiles/GGG2020/fp/al/maps-vertical/`` already exists.

The ``generate_GGG20202_atm_profs.py`` script is by default tailored to Calar Alto. If you wish to calculate atmospheric profiles for a different site, you must change the ``site`` variable, provided the site is among the TCCON site list (if it is not, you can simply add a custom site to the ``tccon_sites.py`` file found in ``Astroclimes/py-ginput/mod_maker/``), or alternatively
change the options on the lines calling the ``run_ginput.py`` script from ``--site`` to ``--lat`` and ``--lon`` and specify the coordinates.
