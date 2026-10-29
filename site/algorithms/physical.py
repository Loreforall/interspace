"""SI-unit electromechanical impact reference model.

Single effective compression mode with linear piezoelectric coupling.
This is a realizable class of *material-backed* harvesters, not a free-space shield.
State: compression x [m], velocity v [m/s], voltage V [V],
load energy [J], internal heat [J], released heat [J].
"""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class PhysicalParameters:
    mass_kg: float = 0.02
    stiffness_n_m: float = 8000.0
    damping_ns_m: float = 2.0
    coupling_n_v: float = 0.03
    capacitance_f: float = 1e-6
    resistance_ohm: float = 5000.0
    initial_velocity_m_s: float = 2.0
    heat_capacity_j_k: float = 20.0
    thermal_conductance_w_k: float = 0.5
    ambient_k: float = 293.15

    def __post_init__(self):
        for name in ('mass_kg','stiffness_n_m','capacitance_f','resistance_ohm','heat_capacity_j_k'):
            if getattr(self,name) <= 0: raise ValueError(f'{name} must be positive')
        for name in ('damping_ns_m','thermal_conductance_w_k'):
            if getattr(self,name) < 0: raise ValueError(f'{name} must be nonnegative')

def derivatives(y, p):
    x,v,V,Eload,Qinternal,Qreleased = y
    # Reciprocity: mechanical coupling force theta*V, electrical current theta*v.
    dx=v
    dv=(-p.stiffness_n_m*x-p.damping_ns_m*v-p.coupling_n_v*V)/p.mass_kg
    dV=(p.coupling_n_v*v-V/p.resistance_ohm)/p.capacitance_f
    dEload=V*V/p.resistance_ohm
    # thermal conduction to ambient with initially zero excess temperature
    dQreleased=p.thermal_conductance_w_k*(Qinternal/p.heat_capacity_j_k)
    dQinternal=p.damping_ns_m*v*v-dQreleased
    return np.array([dx,dv,dV,dEload,dQinternal,dQreleased],dtype=float)

def energy(y,p):
    x,v,V,Eload,Qinternal,Qreleased=np.asarray(y).T
    return (.5*p.mass_kg*v*v + .5*p.stiffness_n_m*x*x +
            .5*p.capacitance_f*V*V + Eload + Qinternal + Qreleased)

def run(p=None, duration_s=0.08, dt_s=2e-6):
    p=p or PhysicalParameters()
    if duration_s<=0 or dt_s<=0: raise ValueError('positive duration and timestep required')
    n=int(np.ceil(duration_s/dt_s)); h=duration_s/n
    t=np.linspace(0,duration_s,n+1)
    y=np.zeros((n+1,6)); y[0,1]=p.initial_velocity_m_s
    for i in range(n):
        a=derivatives(y[i],p); b=derivatives(y[i]+h*a/2,p)
        c=derivatives(y[i]+h*b/2,p); d=derivatives(y[i]+h*c,p)
        y[i+1]=y[i]+h*(a+2*b+2*c+d)/6
    e=energy(y,p); initial=e[0]
    return {'time_s':t,'states':y,'energy_j':e,'initial_energy_j':initial,
            'max_energy_residual_j':float(np.max(np.abs(e-initial))),
            'max_relative_energy_residual':float(np.max(np.abs(e-initial))/max(initial,1e-30)),
            'harvest_fraction':float(y[-1,3]/initial) if initial else 0.,
            'heat_fraction':float((y[-1,4]+y[-1,5])/initial) if initial else 0.}
