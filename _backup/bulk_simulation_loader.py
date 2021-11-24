import pickle
import os
from pathlib import Path

def load_simulation(fn):
    f = Path(fn)
    with f.open('rb') as inf:
        return pickle.load( inf )

ms = dict()

for f in Path("simulation_results").glob("full_run_*.pickle"):
    if not os.path.getsize(f):
        continue
    ms_m = load_simulation(f)
    for k,v in ms_m.items():
        ms[k] = v
            
strats = [
    ('friend', 'nomination (random)'),
    ('friendClose', 'nomination (random, close)'),
    ('friendHighDeg', 'nomination (degree)'),
    ('friendHighDegClose', 'nomination (degree, close)'),
    ('nominate_local_betweenness', 'nomination (local betweenness)'),
    ('local_betweenness', 'local betweenness'),
    ('targeted', 'degree'),
    ('rand', 'random'),
    ('none', 'no intervention'),
]