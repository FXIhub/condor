Installation
============

1. Make sure that you have the followin python packages installed:

   * `numpy <www.numpy.org>`_
   * `scipy <www.scipy.org>`_
   * `h5py <www.h5py.org>`_
 
  If any of these is missing install it with *easy_install* or *pip*::

    pip install numpy scipy h5py

2. Clone the `python_tools <https://bitbucket.org/maxhantke/python_tools>`_ repository and install it::

     git clone https://maxhantke@bitbucket.org/maxhantke/python_tools.git
     cd python_tools
     python setup.py install
     cd ..

3. Clone the `Penguin <https://github.com/mhantke/penguin>`_ repository and install it::

     git clone https://github.com/mhantke/penguin.git
     cd penguin
     python setup.py install
     cd ..

4. Run the example script::

     python example.py

   A set of simulated images should appear in a folder called *example_out*.
