"""Stage 8: reduced-order 3D hemispherical impact model (SI units).

This is a mesoscopic shell/network model, NOT an atomistic model or a
free-standing field. Energy moves as deformation waves through connected
nodes; viscous losses are heat, not automatically recoverable electricity.
"""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class DomeConfig:
    radius: float = 0.25
    rings: int = 9
    sectors: int = 24
    node_mass: float = 0.004
    stiffness: float = 1.5e4
    damping: float = 0.8
    contact_stiffness: float = 2e4
    contact_damping: float = 1.0
    projectile_mass: float = 0.05
    projectile_radius: float = 0.025
    impact_speed: float = 2.0
    dt: float = 1e-5
    duration: float = 0.012


def mesh(config=DomeConfig()):
    """Polar rings from pole to equator; fixed ring is last ring."""
    if config.rings < 3 or config.sectors < 6:
        raise ValueError('Need at least 3 rings and 6 sectors')
    pts = [[0., 0., config.radius]]
    levels = [[0]]
    for j in range(1, config.rings + 1):
        phi = (np.pi / 2) * j / config.rings
        ids = []
        for k in range(config.sectors):
            theta = 2 * np.pi * k / config.sectors
            ids.append(len(pts))
            pts.append([config.radius*np.sin(phi)*np.cos(theta),
                        config.radius*np.sin(phi)*np.sin(theta),
                        config.radius*np.cos(phi)])
        levels.append(ids)
    edges = set()
    for level in levels[1:]:
        for k in range(config.sectors):
            edges.add(tuple(sorted((level[k],level[(k+1)%config.sectors]))))
    for k in range(config.sectors):
        edges.add((0,levels[1][k]))
    for j in range(1,config.rings):
        for k in range(config.sectors):
            edges.add(tuple(sorted((levels[j][k],levels[j+1][k]))))
            edges.add(tuple(sorted((levels[j][k],levels[j+1][(k+1)%config.sectors]))))
    pos = np.array(pts, dtype=float)
    edge = np.array(sorted(edges), dtype=int)
    rest = np.linalg.norm(pos[edge[:,1]]-pos[edge[:,0]],axis=1)
    fixed = np.zeros(len(pos),bool); fixed[levels[-1]]=True
    return pos,edge,rest,fixed


def internal(pos,vel,edge,rest,k,c):
    a,b=edge.T
    delta=pos[b]-pos[a]
    length=np.linalg.norm(delta,axis=1)
    length=np.maximum(length,1e-12)
    direction=delta/length[:,None]
    relative=np.sum((vel[b]-vel[a])*direction,axis=1)
    tension=k*(length-rest)+c*relative
    f=tension[:,None]*direction
    force=np.zeros_like(pos)
    np.add.at(force,a,f); np.add.at(force,b,-f)
    elastic=0.5*k*np.sum((length-rest)**2)
    heat_rate=c*np.sum(relative**2)
    return force,elastic,heat_rate


def run(config=DomeConfig(), sample_every=10):
    """Integrate shell+projectile with semi-implicit Euler and explicit ledgers.

    Projectile is a sphere striking the pole from above. Contact acts only
    along z on the pole node; suitable for centered normal impacts only.
    """
    initial,edges,rest,fixed=mesh(config)
    pos=initial.copy(); vel=np.zeros_like(pos)
    proj_z=config.radius+config.projectile_radius+0.001
    proj_v=-config.impact_speed
    heat=0.; support_work=0.
    initial_energy=0.5*config.projectile_mass*proj_v**2
    history=[]
    steps=int(round(config.duration/config.dt))
    for n in range(steps+1):
        force,elastic,heat_rate=internal(pos,vel,edges,rest,config.stiffness,config.damping)
        overlap=max(0., pos[0,2]+config.projectile_radius-proj_z)
        approach=max(0.,vel[0,2]-proj_v)
        contact=(config.contact_stiffness*overlap+config.contact_damping*approach) if overlap>0 else 0.
        force[0,2]-=contact
        ke_nodes=0.5*config.node_mass*np.sum(vel[~fixed]**2)
        ke_proj=0.5*config.projectile_mass*proj_v**2
        contact_energy=0.5*config.contact_stiffness*overlap**2
        total=ke_nodes+ke_proj+elastic+contact_energy+heat
        if n%sample_every==0 or n==steps:
            history.append((n*config.dt,proj_z,proj_v,ke_proj,ke_nodes,elastic,contact_energy,heat,total-initial_energy,overlap))
        if n==steps: break
        # Anchor reactions have zero power because fixed-node velocity is zero.
        vel[~fixed]+=config.dt*force[~fixed]/config.node_mass
        pos[~fixed]+=config.dt*vel[~fixed]
        pos[fixed]=initial[fixed];vel[fixed]=0.
        proj_v+=config.dt*contact/config.projectile_mass
        proj_z+=config.dt*proj_v
        heat+=config.dt*(heat_rate+config.contact_damping*approach**2*(overlap>0))
    names=['time_s','projectile_z_m','projectile_v_m_s','projectile_ke_J','shell_ke_J','elastic_J','contact_J','heat_J','energy_residual_J','overlap_m']
    return {'history':np.array(history),'columns':names,'initial':initial,'positions':pos,
            'edges':edges,'fixed':fixed,'initial_energy_J':initial_energy,
            'max_abs_energy_residual_J':float(np.max(np.abs(np.array(history)[:,8]))),
            'note':'Idealized network; heat is not recovered energy; no deployment model.'}
