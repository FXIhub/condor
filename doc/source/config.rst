.. _configuration-file:

The configuration file
======================

.. note:: Default configuration values are defined in the file *default.conf* located in the root directory of Condor. To write your own configuration file use *default.conf* as a starting point.

The configuration file is subdivided into three sections: *source*, *sample* and *detector*.

1. Source
---------

   The source is defined by the photon *wavelength*, the focal spot size (*focus_diameter*) and the *pulse_energy*:

   .. code-block:: none

      [source]
      # Wavelength [m]
      wavelength = 1.13E-09
      # Focal spot size (flat-top beam profile) [m]
      focus_diameter = 1.E-06
      # Pulse energy [J]
      pulse_energy = 1E-03

2. Sample 
----------

   The sample section looks different depending on the specified *sample_type*. For a description of the calculation of the diffraction patterns with the individual sample types please got to the section :doc:`sample_types`.

   For uniform bodies a material has to be specified. This can be done by passing a material_type from the following list:

   * *water*
   * *cell*
   * *virus*
   * *protein*
   * *latexball*

   But it can also be customized with a particular atomic composition and an estimate for the mass density. Here is a simple example for a 100-nm-sized water ball:

   .. code-block:: none     
	
      [sample]
      sample_type = uniform_sphere
      diameter = 400E-09
      material_type = custom
      cH = 2
      cO = 1
      massdensity = 1000

   The refractive index is calculated on the basis of scattering factors taken from the `Henke tables <http://henke.lbl.gov/optical_constants/>`_.

   Here is a list of all sample types with their parameters:

2.1 Uniform sphere
^^^^^^^^^^^^^^^^^^

   .. code-block:: none     

      sample_type = uniform_sphere
      # Sample diameter [m]
      diameter = 100E-09

2.2 Uniform spheroid
^^^^^^^^^^^^^^^^^^^^

   .. code-block:: none     

      sample_type = uniform_spheroid
      # Sample diameter [m] perpendicular to axis of rotation symmetry
      diameter_a = 100E-09
      # Sample diameter [m] along axis of rotation symmetry
      diameter_c = 200E-09

2.3 3-dimensional map of the refractive index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     **a) Generate a uniform body of a particular geometry:**

	You can choose from a set of different options for the particle *geometry*:
	
	* *sphere*

	  .. code-block:: none     

	     sample_type = map3d
	     geometry = sphere
	     # Sample diameter [m]
	     diameter = 450E-09

	* *spheroid*

	  .. code-block:: none     

	     sample_type = map3d
	     geometry = spheroid
	     # Sample diameter [m] perpendicular to axis of rotation symmetry
	     diameter_a = 100E-09
	     # Sample diameter [m] along axis of rotation symmetry
	     diameter_c = 200E-09

	* *cube*

	  .. code-block:: none     

	     sample_type = map3d
	     geometry = cube
	     # Edge length [m]
	     edge_length = 100E-09

	* *icosahedron*

	  .. code-block:: none     

	     sample_type = map3d
	     geometry = cube
	     # Sample diameter [m] (the volume of the icosahedron equals the volume of a sphere with the given diameter)
	     diameter = 100E-09

     **b) Read your own map from an HDF5 file**

      In order to read the map from a file specify its location and sampling step size (edge length of one voxel). The map has to have equal dimensions and is expected to be saved in the dataset named */data*. 

      .. code-block:: none     
	 
	 sample_type = map3d
	 # Location of the HDF5 file. The refractive index map has to be saved in the data set "/data". 
	 map3d_fine = /path/to/your/file.h5
	 # Sample diameter [m] (the volume of the icosahedron equals the volume of a sphere with the given diameter)
	 dx_fine = 100E-09

2.5 Atom positions
^^^^^^^^^^^^^^^^^^


3. Detector
-----------

The detector configuration is defining the pixel locations. With the parameters *x_gap_size_in_pixel*, *y_gap_size_in_pixel* and *hole_diameter_in_pixel* the pixel mask can be designed.

.. code-block:: none     

   # sample-detector distance [m]
   distance = 0.74
   # pixel width and height [m]
   pixel_size = 75E-06
   # pixels binned by 'binning' x 'binning'
   binning = 4
   # absolute number of pixels in x/y direction (unbinned)
   Nx = 1024
   Ny = 1024
   # Center position in pixel (pixel (0,0) has its center at x=0.0 y=0.0)
   # Make sure that border pixel is not existing twice! Center should be lying on a pixel
   cx = middle
   cy = middle
   # Central gap between detector halves in pixel
   x_gap_size_in_pixel = 23
   y_gap_size_in_pixel = 0
   # Central hole in detector
   hole_diameter_in_pixel = 70

