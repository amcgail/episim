from epi_model import *

tnets = [
    temporalNetwork.load('high school').range(1,2),
    temporalNetwork.load('high school2')
]

from threading import Thread

ms = defaultdict(list)

from time import time

def start_worker(S2E, tneti, 
                 VACC_Ps=[0.05, 0.10, 0.20],
                 strats=[sampling.friend],
                 N_SIMULATIONS=200
                ):
    
    tnet = tnets[tneti]
    
    print('loading network...')
    
    p = dict(params.covid_estimate)
    p['s2e'] = params.daily_to_momentary(tnet, S2E)
    sim = simulations.SEIR_daily(tnet, p)

    for VACC_P in VACC_Ps:
        N_T_VACC = int(tnet.Nnodes * VACC_P)

        for i, strat in enumerate(strats):
            NAME = strat.__name__#+("_%s" % i)
            NAME = (tneti, S2E, VACC_P, NAME)
            
            st = time()
            print("Starting on ", NAME, N_SIMULATIONS, 'simulations')

            for i in range( N_SIMULATIONS - len(ms[NAME]) ):
                if (i+1)%20 == 0:
                    print("simulation %s" % (i+1))

                sim.init_attributes()

                to_vacc = strat(sim, vaccinateN=N_T_VACC)

                for x in to_vacc:
                    xi = sim.tnet.nodes.index(x)
                    sim.state_change(xi, 'vacc')

                for who in sample([x for x in range(tnet.Nnodes) if not sim.vacc[x]], 20):
                    sim.state_change(who, 'inf')

                sim.run(100)
                ms[NAME].append( sim.meas )

            print(f'finished!', NAME, f'{time() - st:0.1f} seconds')
     

def _cart_prod(L1, L2):
    if L1 == []:
        return L2
    if L2 == []:
        return L1
    
    for x in L1:
        for y in L2:
            if type(x) != list:
                x = [x]
            if type(y) != list:
                y = [y]
            
            yield x+y
                
def cart_prod(L1, L2):
    return list(_cart_prod(L1,L2))

#def cp(L1, L2):
#    return cart_prod(L1, L2)

def cp(*args):
    res = args[0]
    for a in args[1:]:
        res = list(_cart_prod(res,a))
    return res


def make_file_name( name ):
    f = Path(__file__).parent.joinpath('simulation_results', f"{name}.pickle")
    if f.exists():
        raise Exception("given output file exists. delete it or choose another name to continue")
    return str(f)