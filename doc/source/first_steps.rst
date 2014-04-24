Getting started with Penguin
============================

.. note:: For getting started it is recommended to take a look at the file *example.py* located in the root directory of Penguin.

The usage of Penguin is based on two objects:

* Input: :py:class:`penguin.Input`
* Output: :py:class:`penguin.Output`

1. Configure the simulation
---------------------------

For initialization of an instance of the input object a configuration file is passed::

  import penguin
  I = penguin.Input("location/of/configuration/file.conf")

Missing values are set to the default values as specified in *default.conf*.

To avoid unexpected behavior do NOT manipulate any variable within the input object. Only use it for retrieving information as shown in the following example::

  import penguin
  from matplotlib import pyplot
  I = penguin.Input()
  map3d = I.sample.map3d_fine
  map3d_proj = map3d.sum(0)
  pyplot.imshow(abs(map3d_proj)); pyplot.show()

The configuration after initialization is saved in the dictionary *configuration.confDict*. A dictionary of that form can be used instead of a configruation file to instantiate a new input object::

  import penguin
  I = penguin.Input()
  C = I.configuration.confDict
  C["sample"]["diameter"] = 300E-09
  I_new = penguin.Input(C)

2. Generate diffraction data
----------------------------

Instantiation of the output object starts the calculation of diffraction data on the basis of the parameters defined in the input object that is passed::

  import penguin
  I = penguin.Input()
  O = penguin.Output(I)
  intensities = O.get_intensity_pattern()

The complex-valued diffraction values are saved in the attribute *amplitudes* of the output object. A diffraction pattern can be retrieved by calling the method *get_intensity_pattern()*.

  

 
