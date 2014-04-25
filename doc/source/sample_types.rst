Sample types
============

The sample type defines how the sample is being modeled. The corresponding parameter is *sample_type* and it can be set to:

1. *uniform_sphere*
2. *uniform_spheroid*
3. *map3d*
4. *atom_positions*

For details about the configuration see :doc:`here <config>`.

The implementation of each sample type is shortly described below. 

1. Uniform sphere
-----------------

     For a *uniform_sphere* the refractive index map is not calculated explicitly. Instead an analytical formula is used to calculate the diffraction pattern directly [#f1]_ :

     .. math:: F = \sqrt{I_0} \, \rho_e \, r_0 \, \frac{p}{D} \, \frac{4}{3} \pi R^3 \, \frac{ 3 \left[\sin(qR) - qR \cos(qR) \right] } { (qR)^3 }

     * :math:`I_0`: Intensity at the interaction point
     * :math:`\rho_e`: Electron density of the sample (uniform)
     * :math:`r_0`: Classical electron radius
     * :math:`p`: Pixel edge length
     * :math:`D`: Detector distance
     * :math:`R`: Sample radius
     * :math:`q`: Absolute value of the scattering vector

2. Uniform spheroid (ellipsoid of revolution)
---------------------------------------------

     Also for the *uniform_spheroid* an analytical formula can be used to calculate the diffraction pattern directly without generating a refractive index map [#f2]_ :

     .. math:: F = \sqrt{I_0} \, \rho_e \, r_0 \, \frac{p}{D} \, \frac{4}{3} \pi \, a^2 c \frac{ 3 \left[ \sin(qH) - qH \cos(qH) \right] }{(qH)^3}

     .. math:: H = \sqrt{a^2 \sin^2(g)+c^2 \cos^2(g)}

     .. math:: g = \arccos \left( \frac{ -q_x \cos(\theta) \sin(\phi) + q_y \cos(\theta) \cos(\phi) }{ \sqrt{q_x^2+q_y^2} } \right)

     * :math:`I_0`: Intensity at the interaction point
     * :math:`\rho_e`: Electron density of the sample (uniform)
     * :math:`r_0`: Classical electron radius
     * :math:`p`: Pixel edge length
     * :math:`D`: Detector distance
     * :math:`a`: Half-diameter measured perpendicular to the axis of rotational symmetry
     * :math:`c`: Half-diameter measured along the axis of rotational symmetry
     * :math:`q`: Absolute value of the scattering vector
     * :math:`q_x`: x-component of the scattering vector (perpendicular to beam axis)
     * :math:`q_y`: y-component of the scattering vector (perpendicular to beam axis)
     * :math:`\phi`: Angle defining the orientation of the sample (in-plane rotation)
     * :math:`\theta`: Angle defining the orientation of the sample (out-of-plane rotation)       

3. 3-dimensional map of the refractive index
--------------------------------------------

     For *map3d* a 3-dimensional map of the refractive index which is sampled on a regular cartesian grid is used to calculate the diffraction data. A discrete Fourier transform (DFT) is calculated using the `NFFT package <https://www-user.tu-chemnitz.de/~potts/nfft/>`_. The sampled points in diffraction space lie on the Ewald sphere which means curvature is taken into account. Three Euler angles define the orientation of the refractive index map in space.

     .. math:: F =  \sqrt{I_0} \, \rho_e \, r_0 \, \frac{p}{D} \, \mbox{DFT}\left(n-1\right) \Delta x^3 

     * :math:`I_0`: Intensity at the interaction point
     * :math:`\rho_e`: Electron density of the sample (uniform)
     * :math:`r_0`: Classical electron radius
     * :math:`p`: Pixel edge length
     * :math:`D`: Detector distance
     * :math:`n`: 3-dimensional refractive index map (complex-valued)
     * :math:`\Delta x`: Sampling step size of the refractive index map


4. Discrete positions of all atoms
----------------------------------


.. rubric:: References

.. [#f1] Feigin, L. A. & Svergun, D. I. Structure Analysis by Small-Angle X-Ray and Neutron Scattering, Plenum Press / Springer (1987).
.. [#f2] Hamzeh F. M. & Bragg, R. H. Small angle scattering of x-rays from groups of nonrandomly oriented ellipsoids of revolution of low concentration. J. Appl. Phys. 45, 3189-3195 (1974).

