import numpy as np
import pandas as pd

def generate_torus_grid(R=10.0, r=1.0, N_theta=64, N_phi=64):
    """Generate a regular torus parameter grid.

    Parameters
    - R: major radius
    - r: minor radius
    - N_theta: number of samples around the major circle
    - N_phi: number of samples around the tube cross-section

    Returns
    - theta, phi: 1D arrays
    - Theta, Phi: 2D meshgrid arrays
    - coords: (N_theta, N_phi, 3) array of XYZ coordinates
    """
    theta = np.linspace(0.0, 2.0 * np.pi, N_theta, endpoint=False)
    phi = np.linspace(0.0, 2.0 * np.pi, N_phi, endpoint=False)
    Theta, Phi = np.meshgrid(theta, phi, indexing='ij')

    X = (R + r * np.cos(Phi)) * np.cos(Theta)
    Y = (R + r * np.cos(Phi)) * np.sin(Theta)
    Z = r * np.sin(Phi)

    coords = np.stack([X, Y, Z], axis=-1)
    return theta, phi, Theta, Phi, coords


def compute_acceleration_field_3d(R=10.0, r=1.0, omega=1.0, N_theta=64, N_phi=64):
    """Compute ideal and viscous-attenuated centripetal acceleration fields on a torus grid.

    The model is a separable extension of the 1D viscous_model into the toroidal surface
    by applying the same phi-dependent damping around a critical phi value.

    Returns:
    - Theta, Phi: meshgrid arrays
    - a_ideal: (N_theta, N_phi) array of ideal acceleration magnitudes
    - a_viscous: same shape, with damping applied
    - damping: damping field applied (same shape)
    """
    theta, phi, Theta, Phi, coords = generate_torus_grid(R=R, r=r, N_theta=N_theta, N_phi=N_phi)

    # Use phi coordinate for the cross-section dependency similar to 1D model
    # Map phi in [0, 2pi] to the earlier model domain [0,4] by a simple scaling
    # so that phi_range ~ [0,4] corresponds to [0, 2*pi]
    phi_range = (Phi / (2.0 * np.pi)) * 4.0

    phi_crit = np.sqrt(2.6180339887)

    # Ideal centripetal acceleration magnitude at each point (scalar)
    a_ideal = (omega**2) * np.abs(r * (phi_range**2) - R + r)

    # Damping field copied across theta (axisymmetric in theta)
    viscosity_damping = 1.0 - 0.65 * np.exp(-((phi_range - phi_crit) / 0.35)**2)

    a_viscous = a_ideal * viscosity_damping

    return Theta, Phi, a_ideal, a_viscous, viscosity_damping


def calculate_viscous_model_3d(R=10.0, r=1.0, omega=1.0, N_theta=64, N_phi=64):
    """High-level function to compute and optionally export 3D torus viscous model results.

    Returns:
    - dict with keys: theta, phi, Theta, Phi, a_ideal, a_viscous, damping, coords
    """
    theta, phi, Theta, Phi, coords = generate_torus_grid(R=R, r=r, N_theta=N_theta, N_phi=N_phi)
    Theta, Phi, a_ideal, a_viscous, damping = compute_acceleration_field_3d(R, r, omega, N_theta, N_phi)

    result = {
        'theta': theta,
        'phi': phi,
        'Theta': Theta,
        'Phi': Phi,
        'coords': coords,
        'a_ideal': a_ideal,
        'a_viscous': a_viscous,
        'damping': damping,
    }

    # Export a small CSV summary (coarse sampling) to keep file small
    df = pd.DataFrame({
        'theta_idx': Theta.reshape(-1),
        'phi_idx': Phi.reshape(-1),
        'a_ideal': a_ideal.reshape(-1),
        'a_viscous': a_viscous.reshape(-1),
        'damping': damping.reshape(-1),
    })
    df.to_csv('acceleration_viscosity_comparison_3d.csv', index=False)

    return result
