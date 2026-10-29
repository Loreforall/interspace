"""Stage 3: prescribed, phenomenological interatomic coupling.

Dimensionless; not a microscopic electromagnetic interaction law. A bond's
reference midpoint defines its fixed actuation weight. External work is
tracked explicitly; the field source is an idealized actuator, not storage.
"""
from dataclasses import dataclass
import numpy as np
from .core import Parameters, forces_and_potential
from .impact import ImpactConfig, build, energy_partition

@dataclass(frozen=True)
class CouplingConfig:
    stiffness: float = 8.0
    center_x: float = 0.0
    width: float = 2.0
    ramp_start: float = 0.0
    ramp_duration: float = 0.5
    active_until: float = 20.0


def envelope(t, cfg):
    """Smooth C1 on/off envelope and its analytic time derivative."""
    if cfg.ramp_duration <= 0 or cfg.active_until < cfg.ramp_start + cfg.ramp_duration:
        raise ValueError('Invalid ramp timing')
    def smooth(z):
        if z <= 0: return 0., 0.
        if z >= 1: return 1., 0.
        return z*z*(3-2*z), 6*z*(1-z)
    a, da = smooth((t-cfg.ramp_start)/cfg.ramp_duration)
    b, db = smooth((cfg.active_until-t)/cfg.ramp_duration)
    return a*b, (da*b-a*db)/cfg.ramp_duration


def reference_weights(reference, bonds, cfg):
    if cfg.width <= 0 or cfg.stiffness < 0: raise ValueError('Invalid coupling parameters')
    mids = (reference[bonds[:,0],0]+reference[bonds[:,1],0])/2
    return np.exp(-0.5*((mids-cfg.center_x)/cfg.width)**2)


def coupling_terms(state, reference, cfg):
    a,b=state.bonds.T
    d=state.positions[b]-state.positions[a]
    r=np.linalg.norm(d,axis=1)
    if np.any(r < 1e-12): raise ValueError('Coincident bonded particles')
    r0=np.linalg.norm(reference[b]-reference[a],axis=1)
    w=reference_weights(reference,state.bonds,cfg)
    amplitude, derivative=envelope(state.time,cfg)
    delta=r-r0
    k=cfg.stiffness*w*amplitude
    pair=(k*delta/r)[:,None]*d
    f=np.zeros_like(state.positions)
    np.add.at(f,a,pair)
    np.add.at(f,b,-pair)
    u=0.5*np.dot(k,delta**2)
    power=0.5*cfg.stiffness*derivative*np.dot(w,delta**2)
    return f,float(u),float(power)


def combined_forces(state, reference, cfg, params):
    f,u=forces_and_potential(state,params)
    fc,uc,power=coupling_terms(state,reference,cfg)
    return f+fc,u+uc,power,uc


def coupled_step(state, reference, dt, cfg, params):
    """Nonautonomous velocity-Verlet; trapezoidal quadrature of dU/dt."""
    if dt <= 0: raise ValueError('dt must be positive')
    f0,_,p0,_=combined_forces(state,reference,cfg,params)
    a0=f0/state.masses[:,None]
    state.positions += dt*state.velocities + 0.5*dt**2*a0
    state.time += dt
    f1,_,p1,_=combined_forces(state,reference,cfg,params)
    state.velocities += 0.5*dt*(a0+f1/state.masses[:,None])
    return 0.5*dt*(p0+p1)


def run_coupled(impact=ImpactConfig(), coupling=CouplingConfig(), params=Parameters()):
    if impact.duration<=0 or impact.dt<=0 or impact.sample_every<1: raise ValueError('Invalid integration settings')
    state=build(impact,params)
    ref=state.positions.copy()
    initial=energy_partition(state,params)
    e0=initial['total_energy']
    p0=np.array([initial['momentum_x'],initial['momentum_y']])
    work=0.
    records=[]
    def sample():
        _,u,power,uc=combined_forces(state,ref,coupling,params)
        k=0.5*np.sum(state.masses[:,None]*state.velocities**2)
        momentum=np.sum(state.masses[:,None]*state.velocities,axis=0)
        total=k+u
        return {'time':float(state.time),'projectile_x':float(state.positions[-1,0]),
                'projectile_vx':float(state.velocities[-1,0]),
                'projectile_kinetic':float(0.5*state.masses[-1]*np.dot(state.velocities[-1],state.velocities[-1])),
                'lattice_kinetic':float(k-0.5*state.masses[-1]*np.dot(state.velocities[-1],state.velocities[-1])),
                'coupling_potential':uc,'total_energy':float(total),
                'actuator_work':float(work),'actuator_power':power,
                'energy_residual':float(total-e0-work),
                'momentum_x':float(momentum[0]),'momentum_y':float(momentum[1])}
    records.append(sample())
    n=int(round(impact.duration/impact.dt))
    for i in range(1,n+1):
        work+=coupled_step(state,ref,impact.dt,coupling,params)
        if i%impact.sample_every==0 or i==n: records.append(sample())
    residual=max(abs(r['energy_residual']) for r in records)
    momentum=max(np.linalg.norm(np.array([r['momentum_x'],r['momentum_y']])-p0) for r in records)
    summary={'initial_projectile_kinetic':initial['projectile_kinetic'],
             'final_projectile_kinetic':records[-1]['projectile_kinetic'],
             'net_actuator_work':float(work),'final_coupling_potential':records[-1]['coupling_potential'],
             'max_abs_energy_balance_residual':float(residual),
             'max_momentum_error':float(momentum),
             'final_projectile_x':float(state.positions[-1,0]),
             'final_projectile_vx':float(state.velocities[-1,0]),
             'duration':float(state.time)}
    return records,summary,state
