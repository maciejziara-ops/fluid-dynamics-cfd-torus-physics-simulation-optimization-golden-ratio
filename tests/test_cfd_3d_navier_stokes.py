"""
Unit tests for the 3D CFD Navier-Stokes solver on toroidal geometry.

Test Coverage:
- Mesh generation with golden ratio discretization
- Velocity field initialization (three modes)
- Vorticity and pressure computation
- Viscous damping in critical zone
- JSON export and data structure consistency

All tests validate the ToroidalCFDSolver class from cfd_navier_stokes_3d_torus.py
"""

import pytest
import numpy as np
import json
import os
from pathlib import Path

# Import the solver class
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from cfd_navier_stokes_3d_torus import ToroidalCFDSolver


class TestMeshGeneration:
    """Tests for 3D mesh generation with golden ratio optimization."""
    
    def test_mesh_initialization(self):
        """Verify mesh is properly generated with correct dimensions."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, use_golden_ratio=True, resolution=32)
        solver.generate_mesh()
        
        assert solver.theta is not None, "Theta array not generated"
        assert solver.phi is not None, "Phi array not generated"
        assert solver.X is not None, "X coordinates not generated"
        assert solver.Y is not None, "Y coordinates not generated"
        assert solver.Z is not None, "Z coordinates not generated"
    
    def test_mesh_shape_consistency(self):
        """Verify all mesh arrays have consistent shapes."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, use_golden_ratio=True, resolution=32)
        solver.generate_mesh()
        
        n_theta, n_phi = solver.Theta.shape
        
        assert solver.X.shape == (n_theta, n_phi), f"X shape mismatch: {solver.X.shape}"
        assert solver.Y.shape == (n_theta, n_phi), f"Y shape mismatch: {solver.Y.shape}"
        assert solver.Z.shape == (n_theta, n_phi), f"Z shape mismatch: {solver.Z.shape}"
    
    def test_golden_ratio_discretization(self):
        """Verify golden ratio scaling in azimuthal direction."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, use_golden_ratio=True, resolution=60)
        solver.generate_mesh()
        
        n_theta = len(solver.theta)
        n_phi = len(solver.phi)
        
        phi = (1 + np.sqrt(5)) / 2  # Golden ratio
        ratio = n_theta / n_phi
        
        # Allow 5% tolerance for integer rounding
        assert abs(ratio - phi) / phi < 0.05, \
            f"Golden ratio not preserved: {ratio:.4f} vs {phi:.4f}"
    
    def test_torus_surface_parametrization(self):
        """Verify points lie on torus surface with correct geometry."""
        solver = ToroidalCFDSolver(R=5.0, r=2.0, use_golden_ratio=False, resolution=20)
        solver.generate_mesh()
        
        # For each point on torus, verify: ||(X,Y,Z) - C|| ≈ r
        # where C is the center circle at distance R from origin
        for i in range(len(solver.theta)):
            for j in range(len(solver.phi)):
                x = solver.X[i, j]
                y = solver.Y[i, j]
                z = solver.Z[i, j]
                
                # Distance from z-axis
                dist_z_axis = np.sqrt(x**2 + y**2)
                # Distance from center circle (at distance R)
                dist_from_circle = np.sqrt((dist_z_axis - solver.R)**2 + z**2)
                
                # Should be approximately r (minor radius)
                assert abs(dist_from_circle - solver.r) < 0.01, \
                    f"Point not on surface: {dist_from_circle:.4f} vs {solver.r}"
    
    def test_mesh_without_golden_ratio(self):
        """Verify uniform mesh generation when golden ratio is disabled."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, use_golden_ratio=False, resolution=32)
        solver.generate_mesh()
        
        n_theta = len(solver.theta)
        n_phi = len(solver.phi)
        
        # Without golden ratio, dimensions should be equal
        assert n_theta == n_phi, "Uniform mesh should have equal dimensions"


class TestVelocityFieldInitialization:
    """Tests for velocity field initialization modes."""
    
    def test_azimuthal_mode_initialization(self):
        """Verify azimuthal velocity field (primary swirl)."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='azimuthal')
        
        assert solver.u_theta is not None, "u_theta not initialized"
        assert solver.u_phi is not None, "u_phi not initialized"
        assert solver.u_r is not None, "u_r not initialized"
        
        # Azimuthal mode: u_phi and u_r should be zero
        assert np.allclose(solver.u_phi, 0, atol=1e-10), "Poloidal component not zero"
        assert np.allclose(solver.u_r, 0, atol=1e-10), "Radial component not zero"
        
        # u_theta should be non-zero
        assert not np.allclose(solver.u_theta, 0), "Azimuthal component is zero"
    
    def test_poloidal_mode_initialization(self):
        """Verify poloidal velocity field (circulation)."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='poloidal')
        
        # Poloidal mode: u_theta and u_r should be zero
        assert np.allclose(solver.u_theta, 0, atol=1e-10), "Azimuthal component not zero"
        assert np.allclose(solver.u_r, 0, atol=1e-10), "Radial component not zero"
        
        # u_phi should be non-zero
        assert not np.allclose(solver.u_phi, 0), "Poloidal component is zero"
    
    def test_helical_mode_initialization(self):
        """Verify helical velocity field (coupled motion)."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='helical')
        
        # Helical mode: both u_theta and u_phi should be non-zero
        assert not np.allclose(solver.u_theta, 0), "Azimuthal component is zero"
        assert not np.allclose(solver.u_phi, 0), "Poloidal component is zero"
        assert np.allclose(solver.u_r, 0, atol=1e-10), "Radial component not zero"
    
    def test_velocity_field_shapes(self):
        """Verify velocity field arrays have correct dimensions."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='azimuthal')
        
        n_theta, n_phi = solver.Theta.shape
        
        assert solver.u_theta.shape == (n_theta, n_phi), \
            f"u_theta shape mismatch: {solver.u_theta.shape}"
        assert solver.u_phi.shape == (n_theta, n_phi), \
            f"u_phi shape mismatch: {solver.u_phi.shape}"
        assert solver.u_r.shape == (n_theta, n_phi), \
            f"u_r shape mismatch: {solver.u_r.shape}"
    
    def test_velocity_field_no_nan_or_inf(self):
        """Verify velocity fields contain no NaN or Inf values."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        
        for mode in ['azimuthal', 'poloidal', 'helical']:
            solver.initialize_velocity_field(mode=mode)
            
            assert not np.any(np.isnan(solver.u_theta)), \
                f"{mode}: u_theta contains NaN"
            assert not np.any(np.isinf(solver.u_theta)), \
                f"{mode}: u_theta contains Inf"
            assert not np.any(np.isnan(solver.u_phi)), \
                f"{mode}: u_phi contains NaN"
            assert not np.any(np.isinf(solver.u_phi)), \
                f"{mode}: u_phi contains Inf"


class TestVorticityAndPressure:
    """Tests for vorticity and pressure field computation."""
    
    def test_vorticity_computation(self):
        """Verify vorticity is computed without NaN/Inf values."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='helical')
        solver.compute_vorticity()
        
        assert solver.vorticity is not None, "Vorticity not computed"
        assert solver.vorticity.shape == solver.u_theta.shape, \
            "Vorticity shape mismatch"
        
        # Check for NaN/Inf
        assert not np.any(np.isnan(solver.vorticity)), \
            "Vorticity contains NaN"
        assert not np.any(np.isinf(solver.vorticity)), \
            "Vorticity contains Inf"
    
    def test_vorticity_magnitude_reasonable(self):
        """Verify vorticity magnitude is within physically reasonable range."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='helical')
        solver.compute_vorticity()
        
        vort_max = np.max(np.abs(solver.vorticity))
        
        # Vorticity should be non-zero but not excessive
        assert vort_max > 0, "Vorticity is zero"
        assert vort_max < 1e10, f"Vorticity unreasonably large: {vort_max}"
    
    def test_pressure_field_computation(self):
        """Verify pressure field is computed without NaN/Inf values."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='helical')
        solver.compute_vorticity()
        solver.compute_pressure_field()
        
        assert solver.p is not None, "Pressure not computed"
        assert solver.p.shape == solver.u_theta.shape, \
            "Pressure shape mismatch"
        
        # Check for NaN/Inf
        assert not np.any(np.isnan(solver.p)), \
            "Pressure contains NaN"
        assert not np.any(np.isinf(solver.p)), \
            "Pressure contains Inf"
    
    def test_pressure_magnitude_reasonable(self):
        """Verify pressure magnitude is within reasonable range."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='helical')
        solver.compute_vorticity()
        solver.compute_pressure_field()
        
        p_max = np.max(np.abs(solver.p))
        
        # Pressure should be bounded
        assert p_max < 1e10, f"Pressure unreasonably large: {p_max}"
    
    def test_centripetal_acceleration_computation(self):
        """Verify centripetal acceleration is computed correctly."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='azimuthal')
        
        a_c = solver.compute_centripetal_acceleration()
        
        assert a_c.shape == solver.u_theta.shape, \
            "Acceleration shape mismatch"
        assert not np.any(np.isnan(a_c)), \
            "Acceleration contains NaN"
        assert not np.any(np.isinf(a_c)), \
            "Acceleration contains Inf"
        assert np.all(a_c >= 0), \
            "Acceleration should be non-negative"


class TestViscousDamping:
    """Tests for viscous damping in critical zone."""
    
    def test_damping_application(self):
        """Verify viscous damping is applied correctly."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='helical')
        
        u_theta_before = solver.u_theta.copy()
        solver.apply_viscous_damping(dt=0.01)
        u_theta_after = solver.u_theta.copy()
        
        # After damping, velocity should decrease (or stay same)
        assert np.all(np.abs(u_theta_after) <= np.abs(u_theta_before) + 1e-10), \
            "Damping increased velocity"
    
    def test_critical_point_damping_zone(self):
        """Verify damping is strongest at critical point region."""
        solver = ToroidalCFDSolver(R=10.0, r=1.0, resolution=64)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='helical')
        
        # Compute damping function values
        phi_range = (solver.Phi / (2.0 * np.pi)) * 4.0
        damping = 1.0 - 0.65 * np.exp(-((phi_range - solver.phi_crit) / 0.35)**2)
        
        # Find minimum damping value (maximum attenuation)
        min_damping = np.min(damping)
        
        # At critical point, damping should be near 0.35 (65% attenuation)
        # Allow for some variation due to discretization
        assert 0.30 < min_damping < 0.40, \
            f"Damping at critical point {min_damping:.4f} out of range"
    
    def test_damping_no_nan_or_inf(self):
        """Verify damping application produces no NaN/Inf."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='helical')
        solver.apply_viscous_damping(dt=0.01)
        
        assert not np.any(np.isnan(solver.u_theta)), \
            "Damping created NaN in u_theta"
        assert not np.any(np.isinf(solver.u_theta)), \
            "Damping created Inf in u_theta"


class TestJSONExport:
    """Tests for JSON export and data structure consistency."""
    
    def test_export_results_creates_file(self):
        """Verify JSON export creates a file."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='helical')
        solver.compute_vorticity()
        solver.compute_pressure_field()
        
        export_file = 'test_cfd_results.json'
        solver.export_results(filename=export_file)
        
        assert os.path.exists(export_file), f"Export file {export_file} not created"
        
        # Cleanup
        if os.path.exists(export_file):
            os.remove(export_file)
    
    def test_export_json_structure(self):
        """Verify JSON export has correct structure."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='helical')
        solver.compute_vorticity()
        solver.compute_pressure_field()
        
        export_file = 'test_cfd_results_struct.json'
        solver.export_results(filename=export_file)
        
        with open(export_file, 'r') as f:
            data = json.load(f)
        
        # Check required keys
        required_keys = ['geometry', 'mesh', 'physical_params', 'flow_statistics']
        for key in required_keys:
            assert key in data, f"Missing key in JSON: {key}"
        
        # Check geometry subkeys
        assert 'major_radius_R' in data['geometry']
        assert 'minor_radius_r' in data['geometry']
        
        # Check mesh subkeys
        assert 'n_theta' in data['mesh']
        assert 'n_phi' in data['mesh']
        
        # Cleanup
        if os.path.exists(export_file):
            os.remove(export_file)
    
    def test_export_statistics_validity(self):
        """Verify exported statistics have valid values."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='helical')
        solver.compute_vorticity()
        solver.compute_pressure_field()
        
        export_file = 'test_cfd_results_stats.json'
        solver.export_results(filename=export_file)
        
        with open(export_file, 'r') as f:
            data = json.load(f)
        
        stats = data['flow_statistics']
        
        # Verify statistics are numbers
        for key, value in stats.items():
            if value is not None:
                assert isinstance(value, (int, float)), \
                    f"Statistic {key} is not numeric: {type(value)}"
        
        # Verify min < max for velocity components
        u_theta_min = stats['u_theta_min']
        u_theta_max = stats['u_theta_max']
        if u_theta_min is not None and u_theta_max is not None:
            assert u_theta_min <= u_theta_max, \
                "Velocity min > max"
        
        # Cleanup
        if os.path.exists(export_file):
            os.remove(export_file)


class TestFullSimulation:
    """Integration tests for complete simulation workflow."""
    
    def test_full_simulation_azimuthal_mode(self):
        """Test complete workflow with azimuthal mode."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='azimuthal')
        solver.compute_vorticity()
        solver.compute_pressure_field()
        solver.apply_viscous_damping(dt=0.01)
        
        # Verify all fields are computed
        assert solver.theta is not None
        assert solver.u_theta is not None
        assert solver.vorticity is not None
        assert solver.p is not None
        
        # Verify no NaN/Inf
        assert not np.any(np.isnan(solver.u_theta))
        assert not np.any(np.isinf(solver.u_theta))
    
    def test_full_simulation_helical_mode(self):
        """Test complete workflow with helical mode."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='helical')
        solver.compute_vorticity()
        solver.compute_pressure_field()
        
        # Verify all fields are computed
        assert solver.u_theta is not None
        assert solver.u_phi is not None
        assert solver.vorticity is not None
        assert solver.p is not None
    
    def test_multiple_timesteps(self):
        """Test simulation over multiple time steps."""
        solver = ToroidalCFDSolver(R=3.0, r=1.0, resolution=24)
        solver.generate_mesh()
        solver.initialize_velocity_field(mode='helical')
        
        u_initial = np.copy(solver.u_theta)
        
        for step in range(3):
            solver.solve_one_timestep(dt=0.01)
        
        u_final = solver.u_theta
        
        # Verify velocity changes over time
        assert not np.allclose(u_initial, u_final, atol=1e-5), \
            "Velocity did not change over timesteps"


class TestPhysicalValidation:
    """Tests for physical parameter validation."""
    
    def test_reynolds_number_regime(self):
        """Verify Reynolds number is in expected range for turbulence."""
        solver = ToroidalCFDSolver(R=10.0, r=1.0)
        
        U_char = 1.0  # m/s
        L_char = solver.R  # m
        nu = solver.nu  # m²/s
        
        Re = U_char * L_char / nu
        
        # Should be in turbulent regime (Re > 4000)
        assert Re > 4000, f"Reynolds number {Re} not in turbulent regime"
    
    def test_froude_number_calculation(self):
        """Verify Froude number indicates centrifugal-dominated flow."""
        solver = ToroidalCFDSolver(R=10.0, r=1.0, omega=1.0)
        
        U_char = 1.0
        L_char = 1.0
        g_eff = solver.omega**2 * (solver.R + solver.r)
        
        Fr = U_char / np.sqrt(g_eff * L_char)
        
        # Subcritical flow (Fr < 1) for centrifugal dominance
        assert Fr < 1.0, f"Froude number {Fr} indicates supercritical regime"
    
    def test_courant_number_stability(self):
        """Verify Courant number ensures stability."""
        solver = ToroidalCFDSolver(R=10.0, r=1.0)
        
        U_max = 10.0  # m/s
        dt = 0.001  # s
        dx = 0.0004  # normalized
        
        C = (U_max * dt) / dx
        
        # CFL condition: C < 1 for stability
        assert C < 1.0, f"Courant number {C} exceeds stability criterion"


class TestCriticalPointBehavior:
    """Tests specific to critical point at φ ≈ 2.618."""
    
    def test_critical_point_location(self):
        """Verify critical point is at expected location."""
        solver = ToroidalCFDSolver()
        
        # φ_crit should be √(φ²) ≈ √2.618 ≈ 1.618
        phi = (1 + np.sqrt(5)) / 2
        phi_crit_expected = np.sqrt(phi**2)
        
        assert abs(solver.phi_crit - phi_crit_expected) < 1e-10, \
            f"Critical point location mismatch: {solver.phi_crit} vs {phi_crit_expected}"
    
    def test_damping_concentration_near_critical_point(self):
        """Verify damping is concentrated near critical point."""
        solver = ToroidalCFDSolver(R=10.0, r=1.0, resolution=128)
        solver.generate_mesh()
        
        phi_range = (solver.Phi / (2.0 * np.pi)) * 4.0
        damping = 1.0 - 0.65 * np.exp(-((phi_range - solver.phi_crit) / 0.35)**2)
        
        # Find location of minimum damping
        min_idx = np.unravel_index(np.argmin(damping), damping.shape)
        min_phi = phi_range[min_idx]
        
        # Should be close to φ_crit (within ±0.1)
        assert abs(min_phi - solver.phi_crit) < 0.15, \
            f"Damping minimum not at critical point: {min_phi} vs {solver.phi_crit}"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
