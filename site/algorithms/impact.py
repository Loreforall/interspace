"""Stage 2 passive impact mechanics; dimensionless and conservative.

Damage is *diagnostic only*: bond extension beyond a threshold is flagged,
not broken. No unsupported thermal or recoverable-energy claims.
"""
from dataclasses import dataclass
import numpy as np
from .core import Parameters, make_lattice, forces_and_potential, step

@dataclass(frozen=True)
class ImpactConfig:
    nx: int = 8
    ny: int = 7
    spacing: float = 1.0
    speed: float = 1.0
    projectile_x: float = -1.5
    offset_y: float = 0.0
    duration: float = 4.0
    dt: float = 0.001
    sample_every: int = 10
    strain_threshold: float = 0.15
    defect_fraction: float = 0.0
    seed: int = 0


def build(config=ImpactConfig(), params=Parameters()):
    if not (0 <= config.defect_fraction < 1):
        raise ValueError('defect_fraction must be in [0,1)')
    if config.strain_threshold <= 0:
        raise ValueError('strain_threshold must be positive')
    state = make_lattice(config.nx, config.ny, config.spacing, config.speed,
                         config.projectile_x,
                         (config.ny-1)*config.spacing/2+config.offset_y, params)
    if config.defect_fraction:
        rng = np.random.default_rng(config.seed)
        nremove = int(round(len(state.bonds)*config.defect_fraction))
        if nremove:
            retained = rng.choice(len(state.bonds), len(state.bonds)-nremove, replace=False)
            state.bonds = state.bonds[np.sort(retained)]
    return state


def bond_metrics(state, reference_positions, strain_threshold=0.15):
    a,b = state.bonds.T
    ref = np.linalg.norm(reference_positions[b]-reference_positions[a],axis=1)
    cur = np.linalg.norm(state.positions[b]-state.positions[a],axis=1)
    strain = (cur-ref)/ref
    return {'max_tensile_strain':float(np.max(strain)),
            'max_absolute_strain':float(np.max(np.abs(strain))),
            'bonds_over_threshold':int(np.count_nonzero(strain > strain_threshold)),
            'bond_count':len(strain)}


def energy_partition(state, params=Parameters()):
    _, potential = forces_and_potential(state, params)
    k = 0.5*state.masses*np.sum(state.velocities**2,axis=1)
    return {'projectile_kinetic':float(k[-1]),
            'lattice_kinetic':float(np.sum(k[:-1])),
            'potential':float(potential),
            'total_energy':float(np.sum(k)+potential),
            'momentum_x':float(np.dot(state.masses,state.velocities[:,0])),
            'momentum_y':float(np.dot(state.masses,state.velocities[:,1]))}


def run_impact(config=ImpactConfig(), params=Parameters()):
    if config.duration <= 0 or config.dt <= 0 or config.sample_every < 1:
        raise ValueError('Invalid integration settings')
    state=build(config,params)
    reference=state.positions.copy()
    initial=energy_partition(state,params)
    nsteps=int(round(config.duration/config.dt))
    if nsteps < 1: raise ValueError('No integration steps')
    records=[]
    peak_contact=0.0
    peak_displacement=0.0
    max_strain=0.0
    ever_over_threshold=np.zeros(len(state.bonds),dtype=bool)
    min_projectile_x=float(state.positions[-1,0])
    max_projectile_x=min_projectile_x
    def sample():
        en=energy_partition(state,params)
        bond=bond_metrics(state,reference,config.strain_threshold)
        displacement=np.linalg.norm(state.positions[:-1]-reference[:-1],axis=1)
        return {'time':float(state.time), 'projectile_x':float(state.positions[-1,0]),
                'projectile_vx':float(state.velocities[-1,0]),
                'max_lattice_displacement':float(displacement.max()),**en,**bond}
    records.append(sample())
    for i in range(1,nsteps+1):
        step(state,config.dt,params)
        sep=np.linalg.norm(state.positions[:-1]-state.positions[-1],axis=1)
        peak_contact=max(peak_contact,float(np.max(params.contact_stiffness*np.maximum(params.contact_radius-sep,0))))
        dis=np.linalg.norm(state.positions[:-1]-reference[:-1],axis=1)
        peak_displacement=max(peak_displacement,float(dis.max()))
        a,b=state.bonds.T
        r0=np.linalg.norm(reference[b]-reference[a],axis=1)
        strain=(np.linalg.norm(state.positions[b]-state.positions[a],axis=1)-r0)/r0
        max_strain=max(max_strain,float(np.max(strain)))
        ever_over_threshold |= strain>config.strain_threshold
        min_projectile_x=min(min_projectile_x,float(state.positions[-1,0]))
        max_projectile_x=max(max_projectile_x,float(state.positions[-1,0]))
        if i%config.sample_every==0 or i==nsteps: records.append(sample())
    final=energy_partition(state,params)
    e0=initial['total_energy']
    errors=np.array([r['total_energy']-e0 for r in records])
    mom=np.array([[r['momentum_x']-initial['momentum_x'],r['momentum_y']-initial['momentum_y']] for r in records])
    # Translational kinetic energy transfer is not dissipation or harvesting.
    summary={'initial_projectile_kinetic':initial['projectile_kinetic'],
             'final_projectile_kinetic':final['projectile_kinetic'],
             'final_lattice_kinetic':final['lattice_kinetic'],
             'projectile_kinetic_change':final['projectile_kinetic']-initial['projectile_kinetic'],
             'peak_pair_contact_force':peak_contact,
             'peak_lattice_displacement':peak_displacement,
             'peak_tensile_strain':max_strain,
             'unique_bonds_exceeding_strain_threshold':int(ever_over_threshold.sum()),
             'min_projectile_x':min_projectile_x,
             'max_projectile_x':max_projectile_x,
             'max_abs_energy_error':float(np.max(np.abs(errors))),
             'max_momentum_error':float(np.max(np.linalg.norm(mom,axis=1))),
             'total_bonds':len(state.bonds),
             'end_time':state.time}
    return records,summary,state
