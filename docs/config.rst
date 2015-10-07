Running Condor simulations
==========================

There are generally two ways of configuring and running simulations with Condor:

**A) Simulation with a configuration file** (`more info <#a-simulation-with-a-configuration-file>`_)

  Write a configuration file based on the examples below, name the file 'condor.conf' and run the Condor executable in the same folder::

     $ condor -n [number_of_patterns]
     
  Condor will create an HDF5 file named 'condor.cxi', which contains the simulated pattern(s) and additional output. The file follows the standards of CXI-DB files (`www.cxidb.org <http://www.cxidb.org>`_) and can be opened with the file viewer *Owl* (`github.com/FilipeMaia/owl/wiki/Owl-wiki <http://github.com/FilipeMaia/owl/wiki/Owl-wiki>`_).

  For help and more options run::

     $ condor -h
   
**B) Simulation in Python** (`more info <#b-simulation-directly-in-python>`_)

  Create a :class:`condor.experiment.Experiment` instance and call :meth:`condor.experiment.Experiment.propagate` to obtaining the simulation result::

    E = condor.Experiment(source, sample, detector)
    res = E.propagate()
    intensity_pattern = res["entry_1"]["data_1"]["data"]
    fourier_pattern = res["entry_1"]["data_1"]["data_fourier"]

  The simulation result is a dictionary that contains the simulated pattern and additional output.


A) Simulation with a configuration files
----------------------------------------

A Condor configuration file is composed of at least three sections:

`1) Source`_ ``[source]``
     
`2) Particle`_ - at least one particle section:

   `a) Uniform sphere`_ ``[particle_sphere]``

   `b) Uniform spheroid`_ ``[partice_spheroid]``

   `c) Refractive index map`_ ``[particle_map]``

   `d) Atom positions`_ ``[particle_molecule]``

`3) Detector`_ ``[detector]``

.. note:: All section titles have to be unique in a configuration file. If you want to specify more than one particle sections of the same particle model make the section title unique by appending an underscore and a number to the standard title (e.g. ``[particle_sphere_2]``).

1) Source
^^^^^^^^^

This section configures the :class:`condor.source.Source` class instance.

**Example:**

.. literalinclude:: ../examples/configfile/source.conf

2) Particle
^^^^^^^^^^^
		    
a) Uniform sphere
"""""""""""""""""

This section configures a :class:`condor.particle.particle_sphere.ParticleSphere` class instance.

**Example:**

.. literalinclude:: ../examples/configfile/particle_sphere.conf

b) Uniform spheroid
"""""""""""""""""""

This section configures a :class:`condor.particle.particle_spheroid.ParticleSpheroid` class instance.

**Example:**

.. literalinclude:: ../examples/configfile/particle_spheroid.conf
			  
c) Refractive index map
"""""""""""""""""""""""

This section configures a :class:`condor.particle.particle_map.ParticleMap` class instance.

**Example:**

.. literalinclude:: ../examples/configfile/particle_map.conf

d) Atom positions
"""""""""""""""""

This section configures a :class:`condor.particle.particle_molecule.ParticleMolecule` class instance.

**Example:**

.. literalinclude:: ../examples/configfile/particle_molecule.conf
		 

3) Detector
^^^^^^^^^^^

This section configures a :class:`condor.detector.Detector` class instance.

**Example:**

.. literalinclude:: ../examples/configfile/detector.conf


B) Simulation directly in Python
--------------------------------

Simulations are carried out with an instance of the :class:`condor.experiment.Experiment` class. Its constructor requires three arguments

  1) A Source instance :class:`condor.source.Source`

  2) A dictionary with at least one Particle instance:

     - Uniform sphere - :class:`condor.particle.particle_sphere.ParticleSphere` (the key has to start with ``'particle_sphere'``)
       
     - Uniform spheroid - :class:`condor.particle.particle_spheroid.ParticleSpheroid` (the key has to start with ``'particle_spheroid'``)

     - Refractive index map - :class:`condor.particle.particle_map.ParticleMap` (the key has to start with ``'particle_map'``)

     - Atom positions - :class:`condor.particle.particle_molecule.ParticleMolecule` (the key has to start with ``'particle_molecule'``)
     
  3) A Detector instance - :class:`condor.detector.Detector`

Calling the method :meth:`condor.experiment.propagate` starts the simulation of a single diffraction pattern. The method returns a dictionary with various outputs.

Example
^^^^^^^

.. literalinclude:: ../examples/scripts/simple/example.py
   :language: python

