"""INTERSPACE Stage 1: dimensionless, planar classical particle mechanics.

Morse pairs are fixed at their initial nearest-neighbour connectivity. The
projectile uses a finite-range, conservative quadratic contact potential.
"""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Parameters:
    bond_energy: float = 1.0
    bond_width: float = 3.0
    equilibrium_distance: float = 1.0
    contact_stiffness: float = 120.0
    contact_radius: float = 0.65
    mass: float = 1.0
    projectile_mass: float = 4.0

@dataclass
class State:
    positions: np.ndarray
    velocities: np.ndarray
    masses: np.ndarray
    bonds: np.ndarray
    time: float = 0.0

    def copy(self):
        return State(self.positions.copy(), self.velocities.copy(),
                     self.masses.copy(), self.bonds.copy(), self.time)


def make_lattice(nx=8, ny=7, spacing=1.0, projectile_speed=1.0,
                 projectile_x=-1.5, projectile_y=None, params=Parameters()):
    """Create square lattice and one projectile (last particle).

    No supports/constraints: the complete system conserves momentum.
    """
    if nx < 2 or ny < 2 or spacing <= 0:
        raise ValueError('Require nx,ny >= 2 and positive spacing')
    xy = np.array([(i*spacing, j*spacing) for i in range(nx) for j in range(ny)], dtype=float)
    if projectile_y is None:
        projectile_y = (ny-1)*spacing/2
    xy = np.vstack([xy, [projectile_x, projectile_y]])
    v = np.zeros_like(xy)
    v[-1, 0] = projectile_speed
    masses = np.full(len(xy), params.mass)
    masses[-1] = params.projectile_mass
    bonds = []
    for i in range(nx):
        for j in range(ny):
            a = i*ny+j
            if i+1 < nx: bonds.append((a, (i+1)*ny+j))
            if j+1 < ny: bonds.append((a, i*ny+j+1))
    return State(xy, v, masses, np.asarray(bonds, dtype=int))


def forces_and_potential(state, params=Parameters()):
    p = state.positions
    f = np.zeros_like(p)
    potential = 0.0
    a, b = state.bonds[:, 0], state.bonds[:, 1]
    d = p[b] - p[a]
    r = np.linalg.norm(d, axis=1)
    if np.any(r < 1e-12): raise ValueError('Coincident bonded particles')
    exp_term = np.exp(np.clip(-params.bond_width*(r-params.equilibrium_distance), -60, 60))
    potential += float(np.sum(params.bond_energy*(1-exp_term)**2 - params.bond_energy))
    derivative = 2*params.bond_energy*params.bond_width*(1-exp_term)*exp_term
    pair_f = (derivative/r)[:, None]*d
    np.add.at(f, a, pair_f)
    np.add.at(f, b, -pair_f)

    # Last particle is projectile. Contact acts against all lattice particles.
    d = p[-1] - p[:-1]
    r = np.linalg.norm(d, axis=1)
    if np.any(r < 1e-12): raise ValueError('Coincident contact particles')
    overlap = np.maximum(params.contact_radius-r, 0.0)
    potential += float(0.5*params.contact_stiffness*np.dot(overlap, overlap))
    # F_projectile = k*overlap*(projectile - atom)/r
    contact_f = (params.contact_stiffness*overlap/r)[:, None]*d
    f[-1] += contact_f.sum(axis=0)
    f[:-1] -= contact_f
    return f, potential


def diagnostics(state, params=Parameters()):
    _, potential = forces_and_potential(state, params)
    kinetic = 0.5*np.sum(state.masses[:, None]*state.velocities**2)
    momentum = np.sum(state.masses[:, None]*state.velocities, axis=0)
    return {'time': float(state.time), 'kinetic': float(kinetic),
            'potential': float(potential), 'total_energy': float(kinetic+potential),
            'momentum_x': float(momentum[0]), 'momentum_y': float(momentum[1]),
            'projectile_x': float(state.positions[-1, 0]),
            'projectile_vx': float(state.velocities[-1, 0]),
            'projectile_kinetic': float(0.5*state.masses[-1]*np.dot(state.velocities[-1],state.velocities[-1]))}


def step(state, dt, params=Parameters()):
    """Velocity-Verlet for a conservative autonomous Hamiltonian."""
    if dt <= 0: raise ValueError('dt must be positive')
    f0, _ = forces_and_potential(state, params)
    accel0 = f0/state.masses[:, None]
    state.positions += state.velocities*dt + 0.5*accel0*dt**2
    f1, _ = forces_and_potential(state, params)
    state.velocities += 0.5*(accel0+f1/state.masses[:, None])*dt
    state.time += dt
    return state


def simulate(state, duration=3.0, dt=0.001, sample_every=10, params=Parameters()):
    """Return sampled diagnostics and final state; duration is rounded to steps."""
    if duration <= 0 or sample_every < 1: raise ValueError('Invalid duration/sample_every')
    nsteps = int(round(duration/dt))
    if nsteps < 1: raise ValueError('No integration steps')
    records = [diagnostics(state, params)]
    for k in range(1, nsteps+1):
        step(state, dt, params)
        if k % sample_every == 0 or k == nsteps:
            records.append(diagnostics(state, params))
    return records, state
