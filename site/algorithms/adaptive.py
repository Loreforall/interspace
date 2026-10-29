"""Stage 6: constrained feedback harvesting in a dimensionless lattice.

A mechanical dashpot is a phenomenological surrogate, not an EM field.
Storage saturation reduces electrical loading and therefore effective damping.
"""
from dataclasses import dataclass
import numpy as np
from .core import Parameters, make_lattice, forces_and_potential
from .recovery import dashpot_forces

@dataclass(frozen=True)
class AdaptiveConfig:
    mode: str = 'adaptive'  # passive, fixed, adaptive
    dt: float = 0.0005
    duration: float = 5.0
    sample_every: int = 20
    max_damping: float = 1.5
    fixed_damping: float = 0.75
    efficiency: float = 0.7
    capacity: float = 1.0
    load_floor: float = 0.0  # saturated storage -> zero harvesting load
    sensor_gain: float = 2.0
    response_time: float = 0.03
    max_slew: float = 10.0
    speed: float = 1.0


def run_adaptive(cfg=AdaptiveConfig(), params=Parameters()):
    if cfg.mode not in ('passive','fixed','adaptive'): raise ValueError('Invalid mode')
    if cfg.dt<=0 or cfg.duration<=0 or cfg.sample_every<1 or cfg.response_time<=0: raise ValueError('Invalid timing')
    if cfg.max_damping<0 or not 0<=cfg.fixed_damping<=cfg.max_damping or cfg.capacity<0: raise ValueError('Invalid damping/capacity')
    if not 0<=cfg.efficiency<=1 or not 0<=cfg.load_floor<=1 or cfg.sensor_gain<0 or cfg.max_slew<0: raise ValueError('Invalid controls')
    state=make_lattice(projectile_speed=cfg.speed,params=params)
    _,u0=forces_and_potential(state,params)
    initial=float(0.5*np.sum(state.masses[:,None]*state.velocities**2)+u0)
    p0=np.sum(state.masses[:,None]*state.velocities,axis=0)
    stored=heat=0.0
    command=0.0
    records=[]
    n=int(round(cfg.duration/cfg.dt))
    if n<1: raise ValueError('No integration steps')
    def record():
        _,u=forces_and_potential(state,params)
        ke=0.5*np.sum(state.masses[:,None]*state.velocities**2)
        pke=0.5*state.masses[-1]*np.dot(state.velocities[-1],state.velocities[-1])
        momentum=np.sum(state.masses[:,None]*state.velocities,axis=0)
        records.append(dict(time=float(state.time),projectile_x=float(state.positions[-1,0]),
          projectile_vx=float(state.velocities[-1,0]),projectile_ke=float(pke),lattice_ke=float(ke-pke),
          mechanical=float(ke+u),stored=float(stored),heat=float(heat),
          damping=float(command),energy_residual=float(ke+u+stored+heat-initial),
          momentum_residual=float(np.linalg.norm(momentum-p0))))
    record()
    for i in range(n):
        f0,_=forces_and_potential(state,params)
        vh=state.velocities+0.5*cfg.dt*f0/state.masses[:,None]
        state.positions+=cfg.dt*vh
        f1,_=forces_and_potential(state,params)
        state.velocities=vh+0.5*cfg.dt*f1/state.masses[:,None]
        # Feedback senses relative bond motion, rather than oracle knowledge of future impact.
        _,relative_power=dashpot_forces(state,1.0)
        if cfg.mode=='passive': target=0.0
        elif cfg.mode=='fixed': target=cfg.fixed_damping
        else: target=cfg.max_damping*(1-np.exp(-cfg.sensor_gain*np.sqrt(relative_power)))
        if cfg.mode=='adaptive':
            requested=command+cfg.dt*(target-command)/cfg.response_time
            command+=np.clip(requested-command,-cfg.max_slew*cfg.dt,cfg.max_slew*cfg.dt)
            command=float(np.clip(command,0,cfg.max_damping))
        else: command=target
        # Electrical load falls with storage fill. A saturated store cannot maintain
        # full harvesting damping without an external dump/load.
        fill=1.0 if cfg.capacity==0 else min(1.0,stored/cfg.capacity)
        load=cfg.load_floor+(1-cfg.load_floor)*(1-fill)
        effective=command*load
        if effective:
            fd,_=dashpot_forces(state,effective)
            before=0.5*np.sum(state.masses[:,None]*state.velocities**2)
            vnew=state.velocities+cfg.dt*fd/state.masses[:,None]
            after=0.5*np.sum(state.masses[:,None]*vnew**2)
            extracted=float(before-after)
            if extracted < -1e-11: raise RuntimeError('Unstable dissipative step; reduce dt')
            extracted=max(0.0,extracted)
            accepted=min(cfg.efficiency*extracted,max(0.0,cfg.capacity-stored))
            stored+=accepted
            heat+=extracted-accepted
            state.velocities=vnew
        state.time+=cfg.dt
        if (i+1)%cfg.sample_every==0 or i==n-1: record()
    return records,state
