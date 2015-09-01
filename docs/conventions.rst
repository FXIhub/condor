Conventions
===========

Coodrinate system
-----------------

Positions in three dimensional space are described in a right-handed cartesian coordinate system with :math:`x`, :math:`y` and :math:`z` axis. In the lab reference frame the :math:`x` and :math:`z` corresponds typically to *horizontal* coordinate and :math:`y` to the *vertical* coordinate.  

Rotations
---------

Rotations follow the right hand rule.


Representations of non-scalar variables
---------------------------------------

Vectors
^^^^^^^

Vectors :math:`\vec{v}=(x, y, z)` are repesented by numpy arrays :math:`[x, y, z]`.

Matrices
^^^^^^^^

Two-dimensional matrices :math:`M_{i, j}` (row index :math:`i`; column index :math:`j`) are represented by two-dimensional arrays :math:`M[i, j]`.

Grid data
^^^^^^^^^

Two-dimensional (e.g. image) data :math:`D(x,y)` are represented by two-dimensional numpy arrays :math:`D[y, x]` with :math:`x` being the fastest changing dimension.

Three-dimensional volume data :math:`D(x,y,z)` are represented by numpy arrays :math:`D[z, y, x]` with :math:`x` being the fastest and :math:`z` being the slowest changing dimension.

Quaternions
^^^^^^^^^^^

Quaternions :math:`q = w + ix + jy + kz` are represented by numpy arrays of length 4 :math:`[w, x, y, z]` and can represent rotations around the unit vector :math:`\vec{u}` by the angle :math:`\theta` if :math:`q = \cos(\theta/2) + i u_x \sin(\theta/2) + j u_y \sin(\theta/2) + k u_z \sin(\theta/2)`

Physical units
--------------

All physical variables are in (derived) SI units if not stated otherwise in the variable name.
