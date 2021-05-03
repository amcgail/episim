from common_imports import *
from epi_model import *
from epi_model import *

__all__ = [
    'targeted','none'
]

def targeted(a, vaccinateN=100):
    deg = nx.degree(a.G)

    for pers in sorted(range(a.Nnodes), key=lambda x: -deg[x])[:vaccinateN]:
        a.vaccinate(pers)
    a.infect(choice(range(a.Nnodes)))
    a.run()
    return a.meas

def none(a, vaccinateN=None, infectN=1):
    from random import choice

    for n in sample(range(a.Nnodes), infectN):
        a.infect(n)
    a.run()
    return a.meas