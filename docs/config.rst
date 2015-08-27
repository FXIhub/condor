Running condor simulations
==========================

There are generally two ways of configuring and running simulations with Condor:

**A) Configuraion file** (`more info <#a-configuration-file>`_)

  Write a configuration file based on the examples below, name the file 'condor.conf' and run the Condor executable in the same folder::

     $ condor -n [number_of_patterns]
     
  Condor will create an HDF5 file named 'condor.cxi', which contains the simulated pattern(s) and additional output. The file follows the standards of CXI-DB files (`www.cxidb.org <http://www.cxidb.org>`_) and can be opened with the file viewer *Owl* (`github.com/FilipeMaia/owl/wiki/Owl-wiki <http://github.com/FilipeMaia/owl/wiki/Owl-wiki>`_).

  For help and more options run::

     $ condor -h
   
**B) Scripting in Python** (`more info <#b-scripting-in-python>`_)

  Create a :class:`condor.experiment.Experiment` instance and call :meth:`condor.experiment.Experiment.propagate` to obtaining the simulation result::

    E = condor.Experiment(source, sample, detector)
    res = E.propagate()
    intensity_pattern = res["entry_1"]["data_1"]["data"]
    fourier_pattern = res["entry_1"]["data_1"]["data_fourier"]

  The simulation result is a dictionary that contains the simulated pattern and additional output.


A) Configuration file
---------------------

A Condor configuration file is composed of at least four sections:

  `1) Source`_ ``[source]``
     
  `2) Sample`_ ``[sample]``

  `3) Particle`_
  
  At least one of the following particle sections:

     `a) Uniform sphere`_ ``[particle_sphere]``

     `b) Uniform spheroid`_ ``[partice_spheroid]``

     `c) Refractive index map`_ ``[particle_map]``

     `d) Atom positions`_ ``[particle_molecule]``

  `4) Detector`_ ``[detector]``

.. note:: All section titles have to be unique in one configuration file. If you want to specify many particle sections of the same particle model type just append a unique identifier to the standard header (e.g. ``[particle_sphere_2]``).

1) Source
^^^^^^^^^

This section configures the :class:`condor.source.Source` class.

**Example:**

.. literalinclude:: ../examples/configfile/source.conf

2) Sample
^^^^^^^^^

This section configures the :class:`condor.sample.Sample` class instance.

**Example:**

.. literalinclude:: ../examples/configfile/sample.conf

3) Particle
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
		 

4) Detector
^^^^^^^^^^^

This section configures a :class:`condor.detector.Detector` class instance.

**Example:**

.. literalinclude:: ../examples/configfile/detector.conf


B) Scripting in Python
----------------------

The class *Experiment* holds all information of a simulation. Its constructor requires three arguments

  1) Source instance - :class:`condor.source.Source`

  2) Sample instance - :class:`condor.sample.Sample`,

     which contains at least one Particle instance out of

     - Uniform sphere - :class:`condor.particle.particle_sphere.ParticleSphere`
       
     - Uniform spheroid - :class:`condor.particle.particle_spheroid.ParticleSpheroid`

     - Refractive index map - :class:`condor.particle.particle_map.ParticleMap`

     - Atom positions - :class:`condor.particle.particle_molecule.ParticleMolecule`
     
  3) Detector instance - :class:`condor.detector.Detector`

Example
^^^^^^^

.. literalinclude:: ../examples/simple/example.py
   :language: python

