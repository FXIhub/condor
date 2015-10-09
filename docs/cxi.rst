Reading CXI files
=================

1) File structure
-----------------

CXI files are really HDF5 files with a defined structure and naming convention. The files follow the standards of CXI-DB files (`www.cxidb.org <http://www.cxidb.org>`_).

The CXI files that are produced by Condor have the following structure:


+---------------------------------------------+----------------------------------------------------------+---------------------------------+
| Path                                        | Content                                                  | Shape                           |
+=============================================+==========================================================+=================================+
| **/entry_1/data_1/**                        | **Image data group**                                     |                                 |
+---------------------------------------------+----------------------------------------------------------+---------------------------------+
| /entry_1/data_1/data                        | Intensity patterns                                       | (N,X,Y)                         |
+---------------------------------------------+----------------------------------------------------------+---------------------------------+
| /entry_1/data_1/data_fourier                | Diffraction amplitudes (complex valued)                  | (N,X,Y)                         |
+---------------------------------------------+----------------------------------------------------------+---------------------------------+
| /entry_1/data_1/mask                        | Pixel masks (see                                         | (N,X,Y)                         |
|                                             | :class:`condor.utils.pixelmask.PixelMask <Pixel masks>`) |                                 | 
+---------------------------------------------+----------------------------------------------------------+---------------------------------+
| **/source/**                                | **Source data group**                                    |                                 |
+---------------------------------------------+----------------------------------------------------------+---------------------------------+
| /source/pulse_energy                        | pulse_energy                                             | N                               |
+---------------------------------------------+----------------------------------------------------------+---------------------------------+
| /source/. . .                               | . . .                                                    | . . .                           |
+---------------------------------------------+----------------------------------------------------------+---------------------------------+
| **/particles/particle_0i**                  | **Particle #i data group**                               |                                 |
+---------------------------------------------+----------------------------------------------------------+---------------------------------+
| /particles/particle_00/extrinsic_quaternion | Quaternions defining particle orientation                | (N,4)                           |
+---------------------------------------------+----------------------------------------------------------+---------------------------------+
| /particles/particle_00/. . .                | . . .                                                    | . . .                           |
+---------------------------------------------+----------------------------------------------------------+---------------------------------+
| **/detector**                               | **Detector data group**                                  |                                 |
+---------------------------------------------+----------------------------------------------------------+---------------------------------+
| /detector/cx                                | X-coord. for center of the diffraction patterns          | N                               |
+---------------------------------------------+----------------------------------------------------------+---------------------------------+
| /detector/. . .                             | . . .                                                    | . . .                           |
+---------------------------------------------+----------------------------------------------------------+---------------------------------+

CXI files can be opened with any HDF5 file reader. For browsing and inspecting diffraction patterns produced with *Condor* the viewer *Owl* can come handy (`github.com/FilipeMaia/owl/wiki/Owl-wiki <http://github.com/FilipeMaia/owl/wiki/Owl-wiki>`_).



2) Examples
-----------

Here are examples for Python and Matlab demonstrating how to access the data stored in a CXI file.

a) Python
^^^^^^^^^
.. code::

   import h5py, numpy

   with h5py.File('condor.cxi', 'r') as f:
     intensity_pattern0 = numpy.asarray(f["entry_1/data_1/data"])

   print "Maximum intensity value in first pattern: %f photons" % intensity_pattern[0].max()
   
b) Matlab
^^^^^^^^^
.. code::
   
   fid = H5F.open('condor.cxi');
   dset_id = H5D.open(fid,'/entry_1/data_1/data');

   intensity_patterns = H5D.read(dset_id);

   H5D.close(dset_id);
   H5F.close(fid); 

