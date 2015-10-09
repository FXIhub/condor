Reading CXI files
=================

CXI files are really HDF5 files with a defined structure and naming convention. The files follow the standards of CXI-DB files (`www.cxidb.org <http://www.cxidb.org>`_).

The CXI files that are produced by Condor have the following structure:

+-------------------------------+--------------------------------------------+---------------------------------------+
| Dataset path                  | Array dimensions                           | Content                               |
                                | (from fastest to slowest changing dim.)    |                                       |
+===============================+============================================+=======================================+
| /entry_1/data_1/data          | (# patterns, # pixels in Y, # pixels in X) | Intensity patterns                    |
+-------------------------------+--------------------------------------------+---------------------------------------+
| /entry_1/data_1/data_fourier  | (# patterns, # pixels in Y, # pixels in X) | Complex valued diffraction amplitudes |
+-------------------------------+--------------------------------------------+---------------------------------------+


CXI files can be opened with any HDF5 file reader. Here three examples shall demonstrate the procedure.

A) Python
---------

Check first whether you have installed the python package *h5py*. In case of a working *h5py* installation you should be able to execute the following command without any error.

.. code::

   python -c "import h5py"

If *h5py* is not installed just do:

.. code::

   pip install h5py

The following Python script demonstrates how to read the simulated intensity patterns.
 

   
B) Matlab
---------

C) Owl
------

 and can be opened with the file viewer *Owl* (`github.com/FilipeMaia/owl/wiki/Owl-wiki <http://github.com/FilipeMaia/owl/wiki/Owl-wiki>`_).
