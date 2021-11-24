# YAY
from frozendict import frozendict
from epi_model import *
import math

def iter_args_it(kvpairs):
    N = 1
    Ns = [N]
    
    psort = sorted(kvpairs)
    
    for param in psort:
        kvs = kvpairs[param]
        try:
            if type(kvs) in {dict, set, list}:
                N *= len(kvs)
        except TypeError:
            pass
            
        Ns.append(N)
        
        
    for i in range(N):
        res_v, res_n = {}, {}
        for pi, p in enumerate(psort):
            kvs = kvpairs[p]
            if type(kvs) == list:
                kvs = {z:z for z in kvs}
            elif type(kvs) != dict:
                kvs = {kvs:kvs}
                
            vali = math.floor(i / Ns[pi]) % len(kvs)
            val = kvs[ sorted(kvs)[ vali ] ]
                
            res_v[p] = val
            res_n[p] = sorted(kvs)[ vali ]
            
        yield ( frozendict(res_n), res_v )
        
def iter_args(kvpairs):
    return list(iter_args_it(kvpairs))


def friendHighDegRandTop3(sim, vaccinateN):
    return sampling.friendHighDegRandTopN(sim, vaccinateN, N=3)
def friendHighDegRandTop5(sim, vaccinateN):
    return sampling.friendHighDegRandTopN(sim, vaccinateN, N=5)

def friendHighDegNormalErr10(sim, vaccinateN):
    return sampling.friendHighDegNormalErr(sim, vaccinateN, sigma=10)
def friendHighDegNormalErr20(sim, vaccinateN):
    return sampling.friendHighDegNormalErr(sim, vaccinateN, sigma=20)

strats = [
    #sampling.friendClose,
    #sampling.friendHighDegClose,
    sampling.friendHighDeg,
    sampling.targeted,
    sampling.none,
    sampling.rand,
    #sampling.friendWeightedChain,
]

strats += 1*[
    sampling.friend,
    sampling.friendHighDegChain,
    friendHighDegRandTop3,
    friendHighDegRandTop5,
    friendHighDegNormalErr10,
    friendHighDegNormalErr20,
]

argset1 = iter_args({
    'strat': {x.__name__: x for x in strats},
    'R0_mean': 2.5,
    'VACCINATE_P': [0.05, 0.10, 0.20, 0.50],
    'INITIAL_INFECT': 20
})

argset2 = iter_args({
    'strat': {x.__name__: x for x in strats},
    'R0_mean': [1, 2.5, 4],
    'VACCINATE_P': 0.20,
    'INITIAL_INFECT': 20
})

argset3 = iter_args({
    'strat': {x.__name__: x for x in strats},
    'R0_mean': 2.5,
    'VACCINATE_P': 0.20,
    'INITIAL_INFECT': [5, 10, 20]
})



"""
old way to do it, just in case:

from time import time

for strat in strats:
    print('Starting', strat.__name__)
    st_full = time()
    first_it = True
    for VACCINATE_P in [0.05, 0.10, 0.20, 0.50]:
        for INITIAL_INFECT in [5, 10, 20]:
            
            sname = (strat.__name__, f"vacc={VACCINATE_P:0.2f}", f"inf={INITIAL_INFECT:d}")

            for i in range(N_PER - len(models[sname])):

                ns_to_remove = strat(sim, int(Nn*VACCINATE_P))
                ns_to_remove = [n2id[x] for x in ns_to_remove]

                ALPHA = [0.5 if i not in ns_to_remove else 0 for i in range(Nn)]

                args = dict(base_args)
                args['alpha'] = ALPHA
                #args['initI_asym'] = INITIAL_INFECT

                model = ExtSEIRSNetworkModel(**args)
                
                # reconstruct X myself
                ns_to_infect = sample( [i for i in range(Nn) if i not in ns_to_remove], INITIAL_INFECT )
                
                for i in range(Nn):
                    if i in ns_to_infect:
                        model.X[i] = model.I_sym
                
                # hard-code this to avoid ridiculousness
                model.A_deltabeta = model.A
                
                # sanity checks
                assert all( 0<=i<=Nn for i in ns_to_remove )
                assert all( 0<=i<=Nn for i in ns_to_infect )

                st = time()
                model.run(100,verbose=False)
                if first_it:
                    t_iter = time()-st
                    print(f'\t{t_iter:0.1f} seconds for first simulation. At that rate it\'ll take {N_iters*t_iter/3600:0.1f} hours.')
                    first_it = False
                    
                models[sname].append({
                    k: getattr(model, k)
                    for k in attrs_to_keep
                })
                
    print( f'{sname} finished after {time() - st_full:0.1f} seconds' )
"""