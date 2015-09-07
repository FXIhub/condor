Installation
============

1) Dependencies
---------------

a) Python packages
^^^^^^^^^^^^^^^^^^

Make sure that you have the followin python packages installed:

  - *numpy*
  - *scipy*
  - *h5py*

If any of these packages is missing simply install it with *pip* (or *easy_install*):

.. code::
   
   $ pip install numpy scipy h5py

You might need to prepend ``sudo`` to obtain root privileges.
   
b) NFFT
^^^^^^^

Install the `NFFT library <https://www-user.tu-chemnitz.de/~potts/nfft/>`_:

.. code::
   
   wget https://www-user.tu-chemnitz.de/~potts/nfft/download/nfft-3.2.3.tar.gz
   tar xvzf nfft-3.2.3.tar.gz
   cd nfft-3.2.3
   ./configure
   make
   sudo make install
   cd ..

c) libspimage and spsim (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install the libraries `libspimage <https://github.com/FilipeMaia/libspimage>`_ and the `spsim <https://github.com/FilipeMaia/libspimage>`_:

i) libspimage
,,,,,,,,,,,,,

.. code:: bash

   git clone https://github.com/FilipeMaia/libspimage.git
   cd libspimage
   mkdir build
   cd build
   ccmake ..

Hit the **C** key for automatic configuration.
   
Make sure everything is set up correctly, then press the **G** key for generation of the Makefile.

.. code:: bash

   make
   sudo make install
   cd ../..

ii) spsim
,,,,,,,,,

.. code:: bash

   git clone https://github.com/FilipeMaia/spsim.git
   cd spsim
   mkdir build
   cd build
   ccmake ..

Hit the **C** key for automatic configuration.
   
Make sure everything is set up correctly, then press the **G** key for generation of the Makefile.

.. code::

   make
   sudo make install
   cd ../..   

2) Install the Condor package
-----------------------------

.. code:: bash

   git clone --recursive https://github.com/mhantke/condor.git
   cd condor
   python setup.py install

3) Run the tests (optional)
---------------------------

.. code:: bash

   python tests.py

If all tests were successful you will read:

  ``=> SUCCESS: All tests passed successfully``
   
