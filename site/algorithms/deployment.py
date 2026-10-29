"""Stage 9: energy-accounted, prescribed-trajectory dome deployment reference.

A kinematic reference, not a validated deployable actuator. The shell's moving
nodes are commanded along a smooth path; ideal actuator forces are inferred
from Newton's law and the network force. Ring nodes remain stationary.
"""
from dataclasses import dataclass
import numpy as np
from .hemisphere import DomeConfig, mesh, internal

@dataclass(frozen=True)
class DeploymentConfig:
    dome: DomeConfig = DomeConfig()
    activation_time_s: float = 0.08
    duration_s: float = 0.10
    dt_s: float = 0.0001
    max_total_power_W: float = 250.0
    max_total_force_N: float = 1000.0
    stowed_height_fraction: float = 0.02

def quintic(t):
    """Position, first and second time derivatives on unit interval."""
    t=np.clip(t,0.,1.)
    return 10*t**3-15*t**4+6*t**5, 30*t**2-60*t**3+30*t**4, 60*t-180*t**2+120*t**3

def simulate(cfg=DeploymentConfig()):
    if cfg.activation_time_s<=0 or cfg.dt_s<=0 or cfg.duration_s<cfg.activation_time_s:
        raise ValueError('Invalid timing')
    if not (0<=cfg.stowed_height_fraction<1):
        raise ValueError('Invalid stowed height fraction')
    full,edges,_,fixed=mesh(cfg.dome)
    stowed=full.copy(); stowed[:,2]*=cfg.stowed_height_fraction
    # The stowed configuration defines the initial stress-free network.
    rest=np.linalg.norm(stowed[edges[:,1]]-stowed[edges[:,0]],axis=1)
    delta=full-stowed
    times=np.arange(0,cfg.duration_s+cfg.dt_s/2,cfg.dt_s)
    rows=[]; net_work=0.; supplied=0.; regenerated=0.; heat=0.
    last_t=0.;last_power=0.;max_power=0.;max_force=0.
    for t in times:
        q,v,a=quintic(t/cfg.activation_time_s)
        pos=stowed+q*delta
        vel=v*delta/cfg.activation_time_s
        acc=a*delta/cfg.activation_time_s**2
        vel[fixed]=0.;acc[fixed]=0.
        internal_force,elastic,heat_rate=internal(pos,vel,edges,rest,cfg.dome.stiffness,cfg.dome.damping)
        actuator=cfg.dome.node_mass*acc-internal_force
        actuator[fixed]=0.
        power=float(np.sum(actuator*vel))
        force_norm=float(np.linalg.norm(actuator,axis=1).sum())
        if t>0:
            dt=t-last_t
            net_work+=0.5*(last_power+power)*dt
            supplied+=0.5*(max(0,last_power)+max(0,power))*dt
            regenerated+=0.5*(max(0,-last_power)+max(0,-power))*dt
            heat+=heat_rate*dt # first-order heat quadrature; residual converges under refinement
        kinetic=0.5*cfg.dome.node_mass*np.sum(vel[~fixed]**2)
        residual=net_work-(kinetic+elastic+heat)
        rows.append((t,q,power,force_norm,net_work,supplied,regenerated,kinetic,elastic,heat,residual))
        max_power=max(max_power,max(0,power));max_force=max(max_force,force_norm)
        last_t=t;last_power=power
    return {'columns':['time_s','deployment_fraction','actuator_power_W','sum_actuator_force_N','net_work_J','positive_work_J','negative_work_J','kinetic_J','elastic_J','heat_J','energy_residual_J'],
            'history':np.array(rows), 'initial':stowed,'deployed':full,'edges':edges,
            'max_positive_power_W':max_power,'max_sum_force_N':max_force,
            'power_limit_pass':max_power<=cfg.max_total_power_W,
            'force_limit_pass':max_force<=cfg.max_total_force_N,
            'final_fraction':float(rows[-1][1]),
            'note':'Ideal inverse-dynamics actuator model. No motor geometry, material calibration, collision, or deployment stability proven.'}
