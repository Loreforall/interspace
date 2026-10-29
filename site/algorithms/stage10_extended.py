"""Stage 10 completion candidate: storage ledger and radial wave arrival diagnostics.

Uses the Stage 10 reduced-order model. Electrical storage is an idealized
post-load converter; it does NOT feed back into transducer voltage. It cannot
be interpreted as a validated rectifier or battery model.
"""
from dataclasses import dataclass, replace
import numpy as np
from .transport3d import TransportConfig, run as run_transport
from .hemisphere import mesh

@dataclass(frozen=True)
class StorageConfig:
    efficiency: float = 0.75
    capacity_J: float = 1e-4
    initial_J: float = 0.0
    wave_threshold_fraction: float = 1e-5

def run(cfg=TransportConfig(), storage=StorageConfig(), sample_every=10):
    if not 0 <= storage.efficiency <= 1 or storage.capacity_J < 0 or not 0 <= storage.initial_J <= storage.capacity_J:
        raise ValueError('Invalid storage parameters')
    if storage.wave_threshold_fraction <= 0:
        raise ValueError('Wave threshold must be positive')
    base=run_transport(cfg,sample_every)
    h=base['history']; cols=base['columns']
    load=h[:,cols.index('resistor_energy_J')]
    generated=np.maximum(load-load[0],0)
    potential=storage.efficiency*generated
    captured=np.minimum(potential,storage.capacity_J-storage.initial_J)
    stored=storage.initial_J+captured
    converter_heat=(1-storage.efficiency)*generated
    overflow=np.maximum(potential-captured,0)
    # Existing resistor ledger reinterpreted as downstream electrical delivery:
    # storage + conversion loss + overflow = electrical output from the load.
    new=np.column_stack([h,stored,converter_heat,overflow])
    new_cols=cols+['stored_electrical_J','conversion_heat_J','storage_overflow_J']
    audit=generated-(captured+converter_heat+overflow)
    return {**base,'history':new,'columns':new_cols,'stored_energy_J':float(stored[-1]),
            'converter_heat_J':float(converter_heat[-1]),'overflow_J':float(overflow[-1]),
            'storage_partition_residual_J':float(np.max(np.abs(audit))),
            'storage_note':'Ideal downstream accounting only; no converter feedback, rectifier dynamics, or real battery model.'}
