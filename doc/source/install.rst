Installation
============

1. Make sure that you have the followin python packages installed:

   * `numpy <www.numpy.org>`_
   * `scipy <www.scipy.org>`_
   * `h5py <www.h5py.org>`_
 
  If any of these is missing install it with *easy_install* or *pip*::

    pip install numpy scipy h5py

2. Install the `NFFT library <https://www-user.tu-chemnitz.de/~potts/nfft/>`_::

     wget https://www-user.tu-chemnitz.de/~potts/nfft/download/nfft-3.2.3.tar.gz
     cd nfft-3.2.3
     ./configure
     make
     sudo make install
     cd ..

3. Clone the `Condor <https://github.com/mhantke/condor>`_ repository recursively and install it::

     git clone --recursive https://github.com/mhantke/condor.git
     cd condor
     python setup.py install
     cd ..

4. Run the example script::

     cd tests/
     python example.py

   A set of simulated images should appear in a folder called *example_out*.
