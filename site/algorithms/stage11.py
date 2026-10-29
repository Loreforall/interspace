"""Stage 11 reproducible reduced-order verification and bounded parameter study.
No claim of full material calibration, deployment feasibility, or field barrier.
"""
from dataclasses import replace
import numpy as np
from .transport3d import TransportConfig
from .stage10_extended import StorageConfig, run

def study(base=TransportConfig(),storage=StorageConfig(),angles=(0.0,0.38,0.7),
          couplings=(0.0,0.004,0.008),dt_factors=(1.0,0.5)):
    records=[]
    for phi in angles:
        for theta in couplings:
            for fac in dt_factors:
                cfg=replace(base,impact_phi=float(phi),coupling_N_per_V=float(theta),
                            dome=replace(base.dome,dt=base.dome.dt*fac))
                out=run(cfg,storage)
                records.append(dict(phi_rad=float(phi),coupling_N_per_V=float(theta),dt_s=cfg.dome.dt,
                    final_projectile_velocity_m_s=out['final_projectile_velocity_m_s'],
                    electrical_output_J=out['load_energy_J'],stored_energy_J=out['stored_energy_J'],
                    max_energy_residual_J=out['max_abs_energy_residual_J'],
                    storage_residual_J=out['storage_partition_residual_J']))
    return records

def convergence_pairs(records):
    result=[]
    groups={}
    for r in records: groups.setdefault((r['phi_rad'],r['coupling_N_per_V']),[]).append(r)
    for (phi,coupling),group in groups.items():
        group=sorted(group,key=lambda r:r['dt_s'],reverse=True)
        if len(group)>1:
            a,b=group[0],group[-1]
            result.append(dict(phi_rad=phi,coupling_N_per_V=coupling,
                velocity_difference_m_s=abs(a['final_projectile_velocity_m_s']-b['final_projectile_velocity_m_s']),
                electrical_energy_difference_J=abs(a['electrical_output_J']-b['electrical_output_J']),
                coarse_energy_residual_J=a['max_energy_residual_J'],fine_energy_residual_J=b['max_energy_residual_J']))
    return result
