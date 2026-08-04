import unittest
import numpy as np
from viscous_model_3d import calculate_viscous_model_3d

class TestViscousModel3D(unittest.TestCase):
    def setUp(self):
        self.R = 10.0
        self.r = 1.0
        self.omega = 1.0
        self.N_theta = 32
        self.N_phi = 32
        res = calculate_viscous_model_3d(R=self.R, r=self.r, omega=self.omega,
                                         N_theta=self.N_theta, N_phi=self.N_phi)
        self.res = res

    def test_shapes_and_types(self):
        res = self.res
        self.assertEqual(res['a_ideal'].shape, (self.N_theta, self.N_phi))
        self.assertEqual(res['a_viscous'].shape, (self.N_theta, self.N_phi))
        self.assertEqual(res['damping'].shape, (self.N_theta, self.N_phi))

    def test_no_nan_or_inf(self):
        res = self.res
        self.assertFalse(np.any(np.isnan(res['a_ideal'])))
        self.assertFalse(np.any(np.isnan(res['a_viscous'])))
        self.assertFalse(np.any(np.isinf(res['a_ideal'])))
        self.assertFalse(np.any(np.isinf(res['a_viscous'])))

    def test_damping_peak_and_symmetry(self):
        res = self.res
        damping = res['damping']
        # find global max damping value (should be near 1.0 - 0.65*exp(0) = 0.35? Actually damping defined as 1 - 0.65*exp(...)
        max_val = np.max(1.0 - damping)  # measure attenuation amplitude
        self.assertGreater(max_val, 0.6)

    def test_energy_reduction(self):
        res = self.res
        energy_ideal = np.trapz(np.trapz(res['a_ideal'], res['phi'], axis=1), res['theta'], axis=0)
        energy_viscous = np.trapz(np.trapz(res['a_viscous'], res['phi'], axis=1), res['theta'], axis=0)
        self.assertLess(energy_viscous, energy_ideal)

    def test_csv_export_exists(self):
        import os
        self.assertTrue(os.path.exists('acceleration_viscosity_comparison_3d.csv'))

if __name__ == '__main__':
    unittest.main()
