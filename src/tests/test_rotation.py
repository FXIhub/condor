import unittest

import numpy
import condor.utils.rotation as rotation

class TestCaseRotation(unittest.TestCase):
    def test_R_x(self):
        # v =  [x, y, z]
        v_a = numpy.array([0, 1, 0])
        a = numpy.pi/2.
        v_b = rotation.R_x(a).dot(v_a)
        # Rotation matrix around x-axis - observing the right hand rule
        v_b_expected = numpy.array([0, 0, 1])
        for v_b_i, v_b_expected_i in zip(v_b, v_b_expected):
            self.assertAlmostEqual(v_b_i, v_b_expected_i)

    def test_R_y(self):
        # v =  [x, y, z]
        v_a = numpy.array([1, 0, 0])
        a = -numpy.pi/2.
        v_b = rotation.R_y(a).dot(v_a)
        # Rotation matrix around y-axis - observing the right hand rule 
        v_b_expected = numpy.array([0, 0, 1])
        for v_b_i, v_b_expected_i in zip(v_b, v_b_expected):
            self.assertAlmostEqual(v_b_i, v_b_expected_i)

    def test_R_z(self):
        # v =  [x, y, z]
        v_a = numpy.array([1, 0, 0])
        a = numpy.pi/2.
        v_b = rotation.R_z(a).dot(v_a)
        # Rotation matrix around z-axis - observing the right hand rule 
        v_b_expected = numpy.array([0, 1, 0])
        for v_b_i, v_b_expected_i in zip(v_b, v_b_expected):
            self.assertAlmostEqual(v_b_i, v_b_expected_i)

    def test_quat(self):
        angle = numpy.pi/3.
        u = numpy.array([1., 2., 3.])
        u /= numpy.sqrt((u**2).sum())
        q = rotation.quat(angle, u[0], u[1], u[2])
        q_expected = numpy.array([numpy.cos(angle/2.),
                                  numpy.sin(angle/2.)*u[0],
                                  numpy.sin(angle/2.)*u[1],
                                  numpy.sin(angle/2.)*u[2]])
        for q_i, q_expected_i in zip(q, q_expected):
            self.assertAlmostEqual(q_i, q_expected_i)
            
    def test_rotate_quat(self):
        angle = numpy.pi/4.
        for i_rot in range(3):
            v0 = numpy.array([0.,0.,0.])
            v0[(i_rot+1)%3] = 1.
            rotation_axis = numpy.array([0.,0.,0.])
            rotation_axis[i_rot] = 1.
            # q = [w, x, y, z] = [cos(a/2), ux sin(a/2), uy sin(a/2), uz sin(a/2)]
            q = numpy.array([numpy.cos(angle/2.),
                             numpy.sin(angle/2.)*rotation_axis[0],
                             numpy.sin(angle/2.)*rotation_axis[1],
                             numpy.sin(angle/2.)*rotation_axis[2]])
            v1 = rotation.rotate_quat(v0, q)
            v1_expected = numpy.array([0.,0.,0.])
            v1_expected[(i_rot+1)%3] = numpy.sqrt(2.)/2.
            v1_expected[(i_rot+2)%3] = numpy.sqrt(2.)/2.
            for v1_i, v1_expected_i in zip(v1, v1_expected):
                self.assertAlmostEqual(v1_i, v1_expected_i)

    def test_quat_from_rotmx(self):
        a = 1.2345
        R = rotation.R_x(a)
        q = rotation.quat_from_rotmx(R)
        # q = [w, x, y, z] = [cos(a/2), ux sin(a/2), uy sin(a/2), uz sin(a/2)]
        q_expected = numpy.array([numpy.cos(a/2.), numpy.sin(a/2.), 0., 0.])
        for q_i, q_expected_i in zip(q, q_expected):
            self.assertAlmostEqual(q_i, q_expected_i)

    def test_unique_representation_quat(self):
        a = 1.2345
        # These quaternions represent equivalent roatations
        l = 1/numpy.sqrt(3)
        q_same = [rotation.unique_representation_quat(rotation.quat(a, l, l, l)),
                  rotation.unique_representation_quat(rotation.quat(-a, -l, -l, -l)),
                  rotation.unique_representation_quat(rotation.quat(a + 2*numpy.pi, l, l, l)),
                  rotation.unique_representation_quat(rotation.quat(a - 2*numpy.pi, l, l, l))]
        q_expected = q_same[0]
        for q in q_same[1:]:
            for q_i, q_expected_i in zip(q, q_expected):
                self.assertAlmostEqual(q_i, q_expected_i)

    def test_rotmx_from_quat(self, N=1000):
        err = numpy.ones(N)
        for i in range(N):
            theta = numpy.random.rand()*2*numpy.pi
            for R,quat in [(rotation.R_x, rotation.quat_x), (rotation.R_y, rotation.quat_y), (rotation.R_z, rotation.quat_z)]:
                R_x0 = R(theta)
                R_x1 = rotation.rotmx_from_quat(quat(theta))
                for i in range(3):
                    for j in range(3):
                        self.assertAlmostEqual(R_x0[i][j], R_x1[i][j])
    
    def test_rotation_conversions_betw_formalisms(self, N=1000):
        """
        Create random euler angles with different axes formalisms and test conversion in a circular fashion.
        Euler angles => rotation matrix => quaternion => unique representation of quaternion (q_a) => Euler angles => rotation materix => quaterinion => unique representation of quaternion (q_b).
        Finally compare the two quaternions q_a and q_b.
        Repeat N times.
        """
        rotation_matrices = {
            'xyz': (rotation.R_x, rotation.R_y, rotation.R_z),
            'xzy': (rotation.R_x, rotation.R_z, rotation.R_y),
            'yxz': (rotation.R_y, rotation.R_x, rotation.R_z),
            'yzx': (rotation.R_y, rotation.R_z, rotation.R_x),
            'zxy': (rotation.R_z, rotation.R_x, rotation.R_y),
            'zyx': (rotation.R_z, rotation.R_y, rotation.R_x),
            'yxy': (rotation.R_y, rotation.R_x, rotation.R_y),
            'zxz': (rotation.R_z, rotation.R_x, rotation.R_z),
            'xyx': (rotation.R_x, rotation.R_y, rotation.R_x),
            'zyz': (rotation.R_z, rotation.R_y, rotation.R_z),
            'xzx': (rotation.R_x, rotation.R_z, rotation.R_x),
            'yzy': (rotation.R_y, rotation.R_z, rotation.R_y),
        }

        succ = 0
        NN = len(rotation_matrices)*N
        err = numpy.ones(NN)
        k = 0
        for ax,(R1,R2,R3) in rotation_matrices.items():
            for i in range(N):
                euler_a = numpy.random.random(3) * numpy.pi * 2
                euler_a[1] /= 2.
                # Euler angles => rotaion matrix
                R_a = R1(euler_a[0]).dot(R2(euler_a[1]).dot(R3(euler_a[2])))
                # Rotation matrix => unique quaternion
                q_a = rotation.quat_from_rotmx(R_a)
                q_a = rotation.unique_representation_quat(q_a)
                # Unique quaternion => Euler angles
                euler_b = rotation.euler_from_quat(q_a, ax)
                # Euler angles => quaternion
                R_b = R1(euler_b[0]).dot(R2(euler_b[1]).dot(R3(euler_b[2])))
                # Rotation matrix => unique quaternion
                q_b = rotation.quat_from_rotmx(R_b)
                q_b = rotation.unique_representation_quat(q_b)
                # Calculate discrepancy
                for q_a_i, q_b_i in zip(q_a, q_b):
                    self.assertAlmostEqual(q_a_i, q_b_i)

    def test_Rotation_compare_rotation_formalisms(self):
        # Test vector (before and after rotation)
        v0          = numpy.array([0., 1., 0.])
        v1_expected = numpy.array([0., 0., 1.])
        # Positive 90 degree rotation around x-axis
        a = numpy.pi/2.
        
        kwargs_list = [
            {"formalism" : "quaternion",       "values" : rotation.quat(a, 1., 0., 0.)},
            {"formalism" : "rotation_matrix",  "values" : rotation.R_x(a)},
            {"formalism" : "euler_angles_xyz", "values" : numpy.array([a, 0., 0.])},
            {"formalism" : "euler_angles_zxz", "values" : numpy.array([0., a, 0.])},
            {"formalism" : "euler_angles_zyx", "values" : numpy.array([0., 0., a])},
        ]
        for kwargs in kwargs_list:
            R = rotation.Rotation(**kwargs)
            v1 = R.rotate_vector(v0)
            for v1_i, v1_expected_i in zip(v1, v1_expected):
                self.assertAlmostEqual(v1_i, v1_expected_i)

    def test_Rotation_euler_angles_axes(self):
        angle = 0.25*numpy.pi
        for i_rot,rotation_formalism in zip(range(3),["euler_angles_xzx","euler_angles_yxy","euler_angles_zxz"]):
            rotation_values = numpy.array([angle,0.,0.])
            rotation_mode = "extrinsic"
            R = rotation.Rotation(values=rotation_values, formalism=rotation_formalism)
            v0 = numpy.array([0.,0.,0.])
            v0[(i_rot+1)%3] = 1.
            v1 = R.rotate_vector(v0)
            v1_expected = numpy.array([0.,0.,0.])
            v1_expected[(i_rot+1)%3] = numpy.sqrt(2.)/2.
            v1_expected[(i_rot+2)%3] = numpy.sqrt(2.)/2.
            for v1_i, v1_expected_i in zip(v1, v1_expected):
                self.assertAlmostEqual(v1_i, v1_expected_i)
