"""Stage 5: phenomenological interatomic viscous harvesting.

Pairwise dashpots remove relative radial kinetic energy; a fraction is stored,
remainder becomes heat. All quantities dimensionless. No microscopic actuator
or actual electrical conversion mechanism is implied.
"""
from dataclasses import dataclass
import numpy as np
from .core import Parameters, make_lattice, forces_and_potential

@dataclass(frozen=True)
class RecoveryConfig:
    damping: float = 1.0
    efficiency: float = 0.7
    capacity: float = 10.0
    dt: float = 0.0005
    duration: float = 4.0
    sample_every: int = 20


def dashpot_forces(state, damping):
    """Equal/opposite central damping forces for bonded lattice atoms."""
    f = np.zeros_like(state.positions)
    if damping == 0: return f, 0.0
    a, b = state.bonds.T
    delta = state.positions[b] - state.positions[a]
    r = np.linalg.norm(delta, axis=1)
    if np.any(r < 1e-12): raise ValueError('Coincident bonded particles')
    unit = delta / r[:, None]
    relative = np.sum((state.velocities[b]-state.velocities[a])*unit, axis=1)
    force = damping*relative[:, None]*unit
    np.add.at(f, a, force)
    np.add.at(f, b, -force)
    return f, float(damping*np.dot(relative,relative))


def run_recovery(config=RecoveryConfig(), params=Parameters(), speed=1.0):
    if config.dt <= 0 or config.duration <= 0 or config.sample_every < 1: raise ValueError('Invalid integration settings')
    if config.damping < 0 or not 0 <= config.efficiency <= 1 or config.capacity < 0: raise ValueError('Invalid recovery settings')
    state = make_lattice(projectile_speed=speed, params=params)
    f, pot = forces_and_potential(state, params)
    k0 = 0.5*np.sum(state.masses[:,None]*state.velocities**2)
    e0 = k0 + pot
    p0 = np.sum(state.masses[:,None]*state.velocities, axis=0)
    storage, heat = 0.0, 0.0
    records=[]
    nsteps = int(round(config.duration/config.dt))
    def record():
        _, u = forces_and_potential(state, params)
        k = 0.5*np.sum(state.masses[:,None]*state.velocities**2)
        proj = 0.5*state.masses[-1]*np.dot(state.velocities[-1],state.velocities[-1])
        momentum = np.sum(state.masses[:,None]*state.velocities, axis=0)
        records.append(dict(time=float(state.time), projectile_x=float(state.positions[-1,0]), projectile_vx=float(state.velocities[-1,0]), projectile_ke=float(proj), lattice_ke=float(k-proj), mechanical=float(k+u), stored=float(storage), heat=float(heat), energy_residual=float(k+u+storage+heat-e0), momentum_residual=float(np.linalg.norm(momentum-p0))))
    record()
    for i in range(nsteps):
        # Conservative velocity Verlet; explicit dissipative midpoint update.
        f0,_ = forces_and_potential(state, params)
        a0 = f0/state.masses[:,None]
        vhalf = state.velocities + 0.5*config.dt*a0
        state.positions += config.dt*vhalf
        f1,_ = forces_and_potential(state, params)
        vconservative = vhalf + 0.5*config.dt*f1/state.masses[:,None]
        state.velocities = vconservative
        # Pairwise symmetric damping kick, energy removed measured exactly
        # by the mechanical kinetic-energy difference (including cross terms).
        fd,_ = dashpot_forces(state,config.damping)
        vnew = state.velocities + config.dt*fd/state.masses[:,None]
        before = 0.5*np.sum(state.masses[:,None]*state.velocities**2)
        after = 0.5*np.sum(state.masses[:,None]*vnew**2)
        extracted = before-after
        if extracted < -1e-11: raise RuntimeError('Damping unstable: reduce dt')
        extracted = max(0.0,float(extracted))
        accepted = min(config.efficiency*extracted,config.capacity-storage)
        storage += accepted
        heat += extracted-accepted
        state.velocities = vnew
        state.time += config.dt
        if (i+1)%config.sample_every==0 or i==nsteps-1: record()
    return records, state
