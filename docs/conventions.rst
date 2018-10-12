Conventions
===========

Coordinate system
-----------------

Positions in three dimensional space are described in a right-handed cartesian coordinate system with :math:`x`, :math:`y` and :math:`z` axis. In the lab reference frame :math:`x` and :math:`y` are often referred to by *horizontal* and *vertical* coordinate, respectively. The :math:`z` axis is codirectional with the primary beam propagation.

Rotations
---------

Rotations follow the right hand rule and can be expressed by quaternions, Euler angles or rotation matrices. For more details about the implementation of rotations see :mod:`condor.utils.rotation`.


Representations of non-scalar variables
---------------------------------------

Vectors
^^^^^^^

Vectors :math:`\vec{v}=(x, y, z)` are repesented by numpy arrays :math:`[x, y, z]`.

Matrices
^^^^^^^^

Two-dimensional matrices :math:`M_{i, j}` (row index :math:`i`; column index :math:`j`) are represented by two-dimensional numpy arrays :math:`M[i, j]` (:math:`j` being the fastest changing dimension).

Grid data
^^^^^^^^^

Two-dimensional (e.g. image) data :math:`D(x,y)` are represented by two-dimensional numpy arrays :math:`D[y, x]` (:math:`x` being the fastest changing dimension).

Three-dimensional volume data :math:`D(x,y,z)` are represented by numpy arrays :math:`D[z, y, x]` (:math:`x` being the fastest and :math:`z` being the slowest changing dimension).

Quaternions
^^^^^^^^^^^

Quaternions :math:`q = w + ix + jy + kz` are represented by numpy arrays :math:`[w, x, y, z]` and represent rotations around the unit vector :math:`\vec{u} = (u_x,u_y,u_z)` by the angle :math:`\theta`:

.. math::

   q = \cos(\theta/2) + i \, u_x \sin(\theta/2) + j \, u_y \sin(\theta/2) + k \, u_z \sin(\theta/2)

Physical units
--------------

All physical variables are in (derived) SI units if not stated otherwise in the variable name.
