"""Stage 4: physically motivated *single-particle* induced-dipole field model.

SI units. Nonretarded, quasistatic, isotropic linear polarizability, prescribed
Gaussian electric-field amplitude E(x,t). This is NOT a Maxwell solver and
is not a field that can stop arbitrary macroscopic objects. The prescribed
field is an external source; its energy budget is not closed.
"""
from dataclasses import dataclass
import numpy as np

EPS0 = 8.8541878128e-12
KB = 1.380649e-23

@dataclass(frozen=True)
class FieldConfig:
    alpha_si: float = 1.648777274e-41  # SI polarizability of order one atomic volume, F m^2
    amplitude_v_m: float = 1e7
    width_m: float = 1e-9
    center_m: float = 0.0
    particle_mass_kg: float = 6.6335209e-26  # ~40 amu
    initial_speed_m_s: float = 100.0
    initial_x_m: float = -4e-9
    duration_s: float = 8e-11
    dt_s: float = 2e-15
    ramp_s: float = 0.0  # 0 = static; otherwise smooth turn-on


def validate(cfg):
    if cfg.alpha_si < 0 or cfg.amplitude_v_m < 0 or cfg.width_m <= 0 or cfg.particle_mass_kg <= 0 or cfg.dt_s <= 0 or cfg.duration_s <= 0 or cfg.ramp_s < 0:
        raise ValueError('Invalid field or integration parameters')


def ramp(t, cfg):
    if cfg.ramp_s == 0 or t >= cfg.ramp_s:
        return 1.0, 0.0
    if t <= 0:
        return 0.0, 0.0
    z = t/cfg.ramp_s
    return z*z*(3-2*z), 6*z*(1-z)/cfg.ramp_s


def interaction(x, t, cfg):
    """Return force [N], potential [J], explicit source power [W].

    U=-alpha*E^2/2, F=+alpha/2 grad(E^2).
    A static attractive dipole trap cannot form a repulsive barrier for alpha>0.
    """
    a, adot = ramp(t,cfg)
    z=(x-cfg.center_m)/cfg.width_m
    e2=cfg.amplitude_v_m**2 * np.exp(-z*z) * a*a
    u=-0.5*cfg.alpha_si*e2
    force=2*z/cfg.width_m*u  # F = -dU/dx = 2 z U / width
    power=-cfg.alpha_si*cfg.amplitude_v_m**2*np.exp(-z*z)*a*adot
    return float(force),float(u),float(power)


def simulate_field(cfg=FieldConfig(), sample_every=200):
    validate(cfg)
    if sample_every < 1: raise ValueError('sample_every must be >=1')
    n=int(round(cfg.duration_s/cfg.dt_s))
    if n < 1 or n > 2000000: raise ValueError('Invalid number of steps')
    x=cfg.initial_x_m; v=cfg.initial_speed_m_s; t=0.; work=0.
    f,u,p=interaction(x,t,cfg)
    e0=0.5*cfg.particle_mass_kg*v*v+u
    rows=[]
    def record():
        _,up,pp=interaction(x,t,cfg)
        energy=0.5*cfg.particle_mass_kg*v*v+up
        rows.append({'time_s':t,'x_m':x,'velocity_m_s':v,'kinetic_j':0.5*cfg.particle_mass_kg*v*v,
                     'potential_j':up,'source_work_j':work,'energy_residual_j':energy-e0-work,'source_power_w':pp})
    record()
    min_v=v; max_force=abs(f)
    for k in range(1,n+1):
        xnew=x+v*cfg.dt_s+0.5*f/cfg.particle_mass_kg*cfg.dt_s**2
        tnew=k*cfg.dt_s
        fnew,_,pnew=interaction(xnew,tnew,cfg)
        vnew=v+0.5*(f+fnew)/cfg.particle_mass_kg*cfg.dt_s
        work+=0.5*(p+pnew)*cfg.dt_s
        x,v,t,f,p=xnew,vnew,tnew,fnew,pnew
        min_v=min(min_v,v); max_force=max(max_force,abs(f))
        if k%sample_every==0 or k==n: record()
    kinetic_initial=0.5*cfg.particle_mass_kg*cfg.initial_speed_m_s**2
    well_depth=0.5*cfg.alpha_si*cfg.amplitude_v_m**2
    summary={'initial_kinetic_j':kinetic_initial,'potential_well_depth_j':well_depth,
             'well_depth_to_initial_kinetic':well_depth/kinetic_initial if kinetic_initial else None,
             'max_force_n':max_force,'min_velocity_m_s':min_v,'final_velocity_m_s':v,
             'net_source_work_j':work,'max_abs_energy_residual_j':max(abs(r['energy_residual_j']) for r in rows),
             'field_energy_density_peak_j_m3':0.5*EPS0*cfg.amplitude_v_m**2,
             'field_pressure_scale_pa':0.5*EPS0*cfg.amplitude_v_m**2,
             'interpretation':'Attractive induced-dipole interaction; no demonstrated barrier or energy harvesting'}
    return rows,summary


def required_field_for_energy(alpha_si, kinetic_j):
    if alpha_si <= 0 or kinetic_j < 0: raise ValueError('Invalid alpha or energy')
    return np.sqrt(2*kinetic_j/alpha_si)
