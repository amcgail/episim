from networkx.algorithms.shortest_paths import weighted
from seirsplus.models import *
from seirsplus.networks import *
from seirsplus.sim_loops import *
from seirsplus.utilities import *
import networkx
import matplotlib.pyplot as pyplot

from generate_SEIRS import *
from simulation_helper import * # iter_args in particular

import json
from time import time
from frozendict import frozendict as fzd

__all__ = ['simulation_manager']

# top level parameters

T = 100 # time of the simulations
INIT_EXPOSED_PCT = 20/330 # initially exposed... this is modified later anyways now...

## Specifying contact networks
_tnets = {}
def get_tnet(name):
    if name not in _tnets:
        if name == 'high school 1':
            tnet = temporalNetwork.load('high school').range(1,2)
        elif name == 'high school 2':
            tnet = temporalNetwork.load('high school2')
        else:
            raise Exception(f"{name} is not a valid network...")

        _tnets[name] = tnet

    return _tnets[name]

def bootstrap(inf1, inf2, nboots=1000, CI=0.05):
    n1 = inf1['ninf']
    n2 = inf2['ninf']
    
    #assert( inf1['init_sus'] == inf2['init_sus'] ) # this is a simplifying assumption I use throughout... same sized simulations, so ratio of effectiveness == ratio of susceptible infected
    # lol literally can't make this assumption if I'm comparing any strategy with doing nothing!!
    
    samps1 = np.random.choice(n1, size=(nboots, n1.shape[0]), replace=True)
    samps2 = np.random.choice(n2, size=(nboots, n2.shape[0]), replace=True)
    
    means1 = samps1.mean(axis=1)
    means2 = samps2.mean(axis=1)
    
    rats = (means1/inf1['init_sus']) / (means2/inf2['init_sus'])
    meanrat = rats.mean()
    lowrat = np.quantile(rats, CI/2)
    highrat = np.quantile(rats, 1-CI/2)
    
    return lowrat, meanrat, highrat

class simulation_manager:

    def __init__(self, generate_network=False, N=None) -> None:
        self.models = defaultdict(list)
        self.most_recent_model = None
        
        self.base_arg_sets = {}
        self.base_args = None

        #self.default_args = {} # defaults for the defaults
        #self.default_args.update(default_args)

        self.generate_network = generate_network
        self.Nn = N

    def generate_params(self, verbose=False,
        R0_mean=2.5, R0_coeffvar = 0.2,
        presymptomaticPeriod_mean = 2.2, presymptomaticPeriod_coeffvar = 0.5,
        latentPeriod_mean = 5.5, latentPeriod_coeffvar = 0.6,
        symptomaticPeriod_mean = 4.0, symptomaticPeriod_coeffvar = 0.4,
        onsetToHospitalizationPeriod_mean = 1e6, onsetToHospitalizationPeriod_coeffvar = 0.45,
        hospitalizationToDischargePeriod_mean = 1e6, hospitalizationToDischargePeriod_coeffvar = 0.45,
        hospitalizationToDeathPeriod_mean = 1e6, hospitalizationToDeathPeriod_coeffvar = 0.45
    ):
        
        Nn = len(self.G.nodes)
        if verbose:
            print(Nn, "nodes")
        
        # Generate a distribution of expected latent periods (time in Exposed state) and presymptomatic periods (time in Pre-symptomatic infectious state). 
        # The `sigma` and `lamda` rates are calculated as the inverse of the expected exposed and pre-symptomatic periods assigned to each individual, respectively.

        SIGMA   = 1 / gamma_dist(latentPeriod_mean, latentPeriod_coeffvar, Nn)
        LAMDA   = 1 / gamma_dist(presymptomaticPeriod_mean, presymptomaticPeriod_coeffvar, Nn)

        # Generate a distribution of expected (a)symptomatic periods (time in symptomatic or asymptomatic state). The `gamma` rates are calculated as the inverse of the expected (a)symptomatic periods assigned to each individual. 
        # The expected total infectious period for each individual is the sum of their expected pre-symptomatic and (a)symptomatic periods.

        GAMMA   = 1 / gamma_dist(symptomaticPeriod_mean, symptomaticPeriod_coeffvar, Nn)

        infectiousPeriod = 1/LAMDA + 1/GAMMA

        # Generate a distribution of expected onset-to-hospitalization periods (time in symptomatic state before entering hospitalized state for those with severe cases) 
        #    and hospitalization-to-discharge periods (time in hospitalized state for those with non-fatal cases). 
        # The `eta` and `gamma_H` rates are calculated as the inverse of the expected onset-to-hospitalization periods and hospitalization-to-discharge periods 
        #    assigned to each individual, respectively.

        ETA     = 1 / gamma_dist(onsetToHospitalizationPeriod_mean, onsetToHospitalizationPeriod_coeffvar, Nn)
        GAMMA_H = 1 / gamma_dist(hospitalizationToDischargePeriod_mean, hospitalizationToDischargePeriod_coeffvar, Nn)

        # Generate a distribution of hospitalization-to-death periods (time in hospitalized state for those with fatal cases). The `mu_H` rates are calculated as the inverse of the expected hospitalization-to-death periods.
        MU_H    = 1 / gamma_dist(hospitalizationToDeathPeriod_mean, hospitalizationToDeathPeriod_coeffvar, Nn)

        # Specify the percentage of cases that are asymptomatic. 

        PCT_ASYMPTOMATIC = 0 # 0.8
        PCT_HOSPITALIZED = 0 # 0.0004

        # Here we specify fatality rates for hospitalized cases, 
        #    again using rates taken from [Verity et al. (2020)](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext).

        PCT_FATALITY = 0 # 0.3627

        # The transmissibility parameter *β* can be related to the basic reproduction number *R<sub>0</sub>* 
        #    (i.e., the expected number of new infections generated by a single infectious individual in a completely susceptible population) 
        #    by the standard formula: *β = R<sub>0</sub>𝛾*. 
        # *R<sub>0</sub>* is a more interpretable parameter, so we specify transmissibility in terms of *R<sub>0</sub>* and then calculate the corresponding *β* values.

        # First, we generate a distribution of individual *R<sub>0</sub>* values 
        #    (i.e., the expected number of new infections generated by a single *particular* infectious individual in a completely susceptible population).

        R0 = gamma_dist(R0_mean, R0_coeffvar, Nn)

        if verbose:
            dist_info(R0, "Individual R0", bin_size=0.1, plot=True, colors='crimson')
            plt.show()

            # information about the distributions just generated
            dist_info([1/LAMDA, 1/SIGMA, 1/LAMDA+1/SIGMA], ["latent period", "pre-symptomatic period", "total incubation period"], plot=True, colors=['gold', 'darkorange', 'black'], reverse_plot=True)
            plt.show()
            dist_info([1/LAMDA, 1/GAMMA, 1/LAMDA+1/GAMMA], ["pre-symptomatic period", "(a)symptomatic period", "total infectious period"], plot=True, colors=['darkorange', 'crimson', 'black'], reverse_plot=True)
            plt.show()
            #dist_info([1/ETA, 1/GAMMA_H, 1/ETA+1/GAMMA_H], ["onset-to-hospitalization period", "hospitalization-to-discharge period", "onset-to-discharge period"], plot=True, colors=['crimson', 'violet', 'black'], reverse_plot=True)
            #dist_info([1/ETA, 1/MU_H, 1/ETA+1/MU_H], ["onset-to-hospitalization period", "hospitalization-to-death period", "onset-to-death period"], plot=True, colors=['crimson', 'darkgray', 'black'], reverse_plot=True)

        # Individuals are ultimately assigned an [*Individual Transmissibility Value*](https://github.com/ryansmcgee/seirsplus/wiki/ExtSEIRSNetworkModel-Class#transmissibility-parameters) (*β<sup>i</sup>*), which are stored in the `beta` attribute of the model object. 
        # The means of the Individual Transmissibility Values for infectious subpopulations are used to calculate the 
        #    [global transmission](https://github.com/ryansmcgee/seirsplus/wiki/Extended-SEIRS-Model-Description#global-transmission) terms. 
        #    Individual Transmissibility Values may also be used to generate the Pairwise Transmissibility Values used for 
        #    [local transmission](https://github.com/ryansmcgee/seirsplus/wiki/Extended-SEIRS-Model-Description#local-transmission) terms, as we will specify in a few steps.

        BETA = 1/infectiousPeriod * R0

        # Individuals can be assigned different Individual Transmissibility Values for use when they are asymptomatic and quarantine states. 
        # Here we set the transmissibility of quarantined individuals such that the mean effective reproduction number for quarantined individuals is about 0.3. 
        #    This supposes that individuals have different behavior, hygiene, etc. when they are quarantining relative to baseline.

        BETA_Q = BETA * (0.3/R0_mean)

        # Here we set individual susceptibilities (default susceptibility is 1). 
        ALPHA = 1

        P_GLOBALINTXN = 0 # no casual contact outside of the network
        Q_GLOBALINTXN = 0 # same for "quarantine" (even though that's not happening)
        
        DELTA = np.zeros(Nn)+1

        base_args = {
            'G': self.G,
            'p': P_GLOBALINTXN,
            'beta': BETA,
            #'beta_local': wmat_scaled,
            'sigma': SIGMA,
            'lamda': LAMDA,
            'gamma': GAMMA,
            'gamma_asym': GAMMA,
            'eta':ETA,
            'gamma_H':GAMMA_H,
            'mu_H':MU_H,
            'a':PCT_ASYMPTOMATIC,
            'h':PCT_HOSPITALIZED,
            'f':PCT_FATALITY,
            'alpha':ALPHA,
            'beta_pairwise_mode':None, # doesn't matter. we set deltabeta later
            'delta_pairwise_mode':None, # doesn't matter. we set deltabeta later
            'delta': DELTA,
            'G_Q':self.G, # doesn't matter
            'q':Q_GLOBALINTXN,
            'beta_Q':BETA_Q,
            'isolation_time':14,
            #'initI_asym':INIT_EXPOSED # note that I'm putting them into the asymptomatic infectious group, to be similar to the other...
        }

        return base_args

    def generate_net_params(self, verbose=False, 
        R0_mean=2.5, R0_coeffvar = 0.2,
        presymptomaticPeriod_mean = 2.2, presymptomaticPeriod_coeffvar = 0.5,
        latentPeriod_mean = 5.5, latentPeriod_coeffvar = 0.6,
        symptomaticPeriod_mean = 4.0, symptomaticPeriod_coeffvar = 0.4,
        onsetToHospitalizationPeriod_mean = 11.0, onsetToHospitalizationPeriod_coeffvar = 0.45,
        hospitalizationToDischargePeriod_mean = 11.0, hospitalizationToDischargePeriod_coeffvar = 0.45,
        hospitalizationToDeathPeriod_mean = 7.0, hospitalizationToDeathPeriod_coeffvar = 0.45
    ):

        assert self.Nn is not None

        demographic_graphs, individual_ageGroups, households = generate_demographic_contact_network(
                                                                    N=self.Nn, demographic_data=household_country_data('US'), 
                                                                    distancing_scales=[0.7], isolation_groups=[])

        G_baseline   = demographic_graphs['baseline']
        G_quarantine = demographic_graphs['distancingScale0.7']

        households_indices = [household['indices'] for household in households]

        SIGMA   = 1 / gamma_dist(latentPeriod_mean, latentPeriod_coeffvar, self.Nn)
        LAMDA   = 1 / gamma_dist(presymptomaticPeriod_mean, presymptomaticPeriod_coeffvar, self.Nn)

        GAMMA   = 1 / gamma_dist(symptomaticPeriod_mean, symptomaticPeriod_coeffvar, self.Nn)

        infectiousPeriod = 1/LAMDA + 1/GAMMA

        ETA     = 1 / gamma_dist(onsetToHospitalizationPeriod_mean, onsetToHospitalizationPeriod_coeffvar, self.Nn)
        GAMMA_H = 1 / gamma_dist(hospitalizationToDischargePeriod_mean, hospitalizationToDischargePeriod_coeffvar, self.Nn)

        # Generate a distribution of hospitalization-to-death periods (time in hospitalized state for those with fatal cases). The `mu_H` rates are calculated as the inverse of the expected hospitalization-to-death periods.

        MU_H    = 1 / gamma_dist(hospitalizationToDeathPeriod_mean, hospitalizationToDeathPeriod_coeffvar, self.Nn)

        PCT_ASYMPTOMATIC = 0.25
        PCT_ASYMPTOMATIC = [0.8 if age in ['0-9', '10-19'] else PCT_ASYMPTOMATIC for age in individual_ageGroups]

        ageGroup_pctHospitalized = {'0-9':      0.0000,
                                    '10-19':    0.0004,
                                    '20-29':    0.0104,
                                    '30-39':    0.0343,
                                    '40-49':    0.0425,
                                    '50-59':    0.0816,
                                    '60-69':    0.118,
                                    '70-79':    0.166,
                                    '80+':      0.184 }
        PCT_HOSPITALIZED = [ageGroup_pctHospitalized[ageGroup] for ageGroup in individual_ageGroups]

        ageGroup_hospitalFatalityRate = {'0-9':     0.0000,
                                        '10-19':   0.3627,
                                        '20-29':   0.0577,
                                        '30-39':   0.0426,
                                        '40-49':   0.0694,
                                        '50-59':   0.1532,
                                        '60-69':   0.3381,
                                        '70-79':   0.5187,
                                        '80+':     0.7283 }

        PCT_FATALITY = [ageGroup_hospitalFatalityRate[ageGroup] for ageGroup in individual_ageGroups]

        print(R0_mean, R0_coeffvar, self.Nn)
        R0 = gamma_dist(R0_mean, R0_coeffvar, self.Nn)

        # Individuals are ultimately assigned an [*Individual Transmissibility Value*](https://github.com/ryansmcgee/seirsplus/wiki/ExtSEIRSNetworkModel-Class#transmissibility-parameters) (*β<sup>i</sup>*), which are stored in the `beta` attribute of the model object. 
        # The means of the Individual Transmissibility Values for infectious subpopulations are used to calculate the [global transmission](https://github.com/ryansmcgee/seirsplus/wiki/Extended-SEIRS-Model-Description#global-transmission) terms. Individual Transmissibility Values may also be used to generate the Pairwise Transmissibility Values used for [local transmission](https://github.com/ryansmcgee/seirsplus/wiki/Extended-SEIRS-Model-Description#local-transmission) terms, as we will specify in a few steps.

        BETA = 1/infectiousPeriod * R0

        # Individuals can be assigned different Individual Transmissibility Values for use when they are asymptomatic and quarantine states. 
        # Here we set the transmissibility of quarantined individuals such that the mean effective reproduction number for quarantined individuals is about 0.3. This supposes that individuals have different behavior, hygiene, etc. when they are quarantining relative to baseline.

        BETA_Q = BETA * (0.3/R0_mean)

        # Now we specify how we would like the [*Pairwise Transmissibility Values*](https://github.com/ryansmcgee/seirsplus/wiki/ExtSEIRSNetworkModel-Class#transmissibility-parameters) (*β<sup>i</sup>*), which define the local transmissibility for each pair of close contacts, to be calculated. There are [multiple ways to specify these pairwise transmissibilities](https://github.com/ryansmcgee/seirsplus/wiki/ExtSEIRSNetworkModel-Class#pairwise-transmissibility-values) (such as providing a matrix), but here we will set the `beta_pairwise_mode` argument of the `ExtSEIRSNetworkModel` constructor to `'infected'`, which will direct the model to automatically generate a matrix of Pairwise Transmissibility Values such that the transmissibility of each infectious-susceptible interaction is equal to the infected individual's transmissiblity.
        BETA_PAIRWISE_MODE  = 'infected'

        # Here we designate that we would like the model to automatically calculate [Connectivity Correction Factors](https://github.com/ryansmcgee/seirsplus/wiki/Extended-SEIRS-Model-Description#connectivity-correction-factor) for each pair of interacting close contacts. This pairwise factor is optional, but it can be used to weight the transmissibility of interactions according to the connectivity of the interacting individuals. Here we choose to weight interactions according to a ratio of the pair's mean degree to the population's mean degree (see [Specifying connectivity Correction Factors](https://github.com/ryansmcgee/seirsplus/wiki/ExtSEIRSNetworkModel-Class#connectivity-correction-factors) for more information). 
        # Using this definition, when two individuals whose average degree is an order of magnitude greater than the average degree of the population overall, then the propensity of exposure in their interaction is weighted to be twice that of two averagely-connected individuals. 
        DELTA_PAIRWISE_MODE = 'mean'

        # Here we set individual susceptibilities (default susceptibility is 1). 
        # In particular, we specify that children are half as susceptible as adults.
        ALPHA = [0.5 if age in ['0-9', '10-19'] else 1.0 for age in individual_ageGroups]

        # In the stochastic network model, an individual comes into contact with a random individual from the population at large (e.g., in a public space) with probability *p* or with an individual from their set of close contacts with probability *(1-p)*. Transmission that occurs between an individual and the population at large is referred to as [global transmission](https://github.com/ryansmcgee/seirsplus/wiki/Extended-SEIRS-Model-Description#global-transmission), and transmission between an individual and one of their close contacts (network neighbors) is referred to as [local transmission](https://github.com/ryansmcgee/seirsplus/wiki/Extended-SEIRS-Model-Description#local-transmission). The parameter *p* defines the locality of the network: for *p=0* an individual only interacts with their close contacts, while *p=1* represents a uniformly mixed population.
        # Here we set *p* to reflect 20% of interactions being with incidental or casual contacts outside their set of close contacts.
        P_GLOBALINTXN = 0.2 #AHHH DO I WANT THIS?!

        # The parameter *q* (down)weights the rate of interactions with the population at large while one is in quarantine relative to baseline.
        # Here we set *q* to 0.05, which supposes that global interactions are quite rare (but nonzero) for quarantined individuals.
        Q_GLOBALINTXN = 0.05

        if verbose:
            dist_info(R0, "Individual R0", bin_size=0.1, plot=True, colors='crimson')

            # information about the distributions just generated
            dist_info([1/LAMDA, 1/SIGMA, 1/LAMDA+1/SIGMA], ["latent period", "pre-symptomatic period", "total incubation period"], plot=True, colors=['gold', 'darkorange', 'black'], reverse_plot=True)
            dist_info([1/LAMDA, 1/GAMMA, 1/LAMDA+1/GAMMA], ["pre-symptomatic period", "(a)symptomatic period", "total infectious period"], plot=True, colors=['darkorange', 'crimson', 'black'], reverse_plot=True)
            #dist_info([1/ETA, 1/GAMMA_H, 1/ETA+1/GAMMA_H], ["onset-to-hospitalization period", "hospitalization-to-discharge period", "onset-to-discharge period"], plot=True, colors=['crimson', 'violet', 'black'], reverse_plot=True)
            #dist_info([1/ETA, 1/MU_H, 1/ETA+1/MU_H], ["onset-to-hospitalization period", "hospitalization-to-death period", "onset-to-death period"], plot=True, colors=['crimson', 'darkgray', 'black'], reverse_plot=True)

        base_args = {
            'G': G_baseline,
            'p': P_GLOBALINTXN,
            'beta': BETA,
            'sigma': SIGMA,
            'lamda': LAMDA,
            'gamma': GAMMA,
            'gamma_asym': GAMMA,
            'eta':ETA,
            'gamma_H':GAMMA_H,
            'mu_H':MU_H,
            'a':PCT_ASYMPTOMATIC,
            'h':PCT_HOSPITALIZED,
            'f':PCT_FATALITY,
            'alpha':ALPHA,
            'beta_pairwise_mode':BETA_PAIRWISE_MODE,
            'delta_pairwise_mode':DELTA_PAIRWISE_MODE,
            'G_Q':G_quarantine,
            'q':0,
            'beta_Q':BETA_Q,
            'isolation_time':14,
            'initI_asym':0, # note that I'm putting them into the asymptomatic infectious group, to be similar to the other...
        }

        return base_args

    def param_info(self):
        if self.base_args is None:
            self.base_args = list(self.base_arg_sets.values())[0]

        network_info(self.base_args['G'], "Baseline", plot=True)
        network_info(self.base_args['G_Q'], "Quarantine", plot=True)
        dist_info(
            [1/np.array( self.base_args['lamda'] ), 1/np.array( self.base_args['gamma'] ), 1/np.array( self.base_args['lamda'] )+1/np.array( self.base_args['gamma'] )], 
            ["pre-symptomatic period", "(a)symptomatic period", "total infectious period"], 
            plot=True,
            colors=['darkorange', 'crimson', 'black'], 
            reverse_plot=True
        )
        dist_info([1/np.array( self.base_args['eta'] ), 1/np.array( self.base_args['gamma_H'] ), 1/np.array( self.base_args['eta'] )+1/np.array( self.base_args['gamma_H'] )], 
            ["onset-to-hospitalization period", "hospitalization-to-discharge period", "onset-to-discharge period"], plot=True, 
            colors=['crimson', 'violet', 'black'], reverse_plot=True)

        R0 = ( 1/np.array( self.base_args['lamda'] ) + 1/np.array( self.base_args['gamma'] ) ) * np.array( self.base_args['beta'] )
        dist_info(np.array( self.base_args['R0'] ), "Individual R0", bin_size=0.1, plot=True, colors='crimson')
        dist_info([1/np.array( self.base_args['eta'] ), 1/np.array( self.base_args['mu_H'] ), 1/np.array( self.base_args['eta'] )+1/np.array( self.base_args['mu_H'] )],
         ["onset-to-hospitalization period", "hospitalization-to-death period", "onset-to-death period"], plot=True, 
         colors=['crimson', 'darkgray', 'black'], reverse_plot=True)



    def load_net(self, name):
        self.tnet = get_tnet(name)

        self.n2id = {n:i for i,n in enumerate(self.tnet.node_ids)}

        counts = Counter([(i,j) for _,i,j in self.tnet.edgelist])
        wts = np.array([c/(30*60/20) for (i,j),c in counts.items()]) # 30 minutes ~= 1X transmissibility
        self.G = nx.Graph()
        self.G.add_weighted_edges_from( [(self.n2id[i], self.n2id[j], wts[ei]) for ei,((i,j),c) in enumerate(counts.items())] )

        self.init_network()

        # these should be reset, for a new network...
        self.base_arg_sets = {}
        self.base_args = None

    def get_params(self, verbose=False, generate_new=False, **kwargs):
        # I used to just use the value of R0_mean, not a whole dictionary
        # R0_mean=None, 
        #if R0_mean is not None:
        #    kwargs = {'R0_mean': R0_mean}

        if not generate_new:
            # keep a network for each set of network generation parameters
            kwargs = fzd(kwargs)
            #print(kwargs)
            if kwargs not in self.base_arg_sets:
                if self.generate_network:
                    self.base_arg_sets[kwargs] = self.generate_net_params(verbose=verbose, **kwargs)
                else:
                    self.base_arg_sets[kwargs] = self.generate_params(verbose=verbose, **kwargs)

            return self.base_arg_sets[kwargs]
        else:
            if self.generate_network:
                return self.generate_net_params(verbose=False, **kwargs)
            else:
                return self.generate_params(verbose=False, **kwargs)

    def init_network(self, draw=False):
        # ================================================
        # ======   Reformulate temporal network    =======
        # ================================================

        # in case the graph was loaded from something other than epi_model
        # really should get rid of this...

        nodes = sorted(self.G.nodes)
        edges = list(self.G.edges(data='weight'))
        self.tnet = weightedNetwork(nodes, edges)

        self.Nn = len(self.G.nodes)

        if False:
            plt.hist( wts, bins=50, log=True );
            plt.vlines(1,0,50000);
            print( ", ".join([f"{x:0.2f}" for x in [np.quantile(wts,0.25), np.quantile(wts,0.50), np.quantile(wts,0.75), np.quantile(wts,0.95)]]) )
            plt.show()

            

            nx.draw(self.G, node_size=0, edge_color=(0.2,0.2,0.2,0.05))
            plt.show()
            network_info(self.G, "Baseline", plot=True)
            plt.show()

    def run_and_checkpoint(self, argsets, N_PER=20, verbose=True, interval=20, dump_name=None, generate_new=False):
        import math

        assert dump_name is not None

        STEP = interval
        MAX = N_PER
        START = max(STEP-1, min(len( self.models[x] ) for (x,_) in argsets) if len(self.models) else 0)
        START = math.ceil((START+1)/STEP)*STEP

        print(f"starting at {START}")

        for i in range(START,MAX+STEP,STEP):
            #print("round {i}")
            self.run(argsets, N_PER=i, verbose=verbose, generate_new=generate_new)
            self.dump_models(dump_name)
            print(f"Checkpointing after {i}/{MAX} simulations")

    def run(self, argsets, N_PER=20, generate_new=False, verbose=True):

        # this is just a dummy, to pass to the sampling algorithm
        # ideally this would all be migrated back, to just taking a Graph,
        # but I'm lazy right now

        p = dict(params.covid_estimate)

        # determine which attributes to keep.
        # this could (should) be hard-coded...
        dummy_args = self.get_params() # dummy attrs
        model = ExtSEIRSNetworkModel(**dummy_args) # dummy model
        attrs_to_keep = [k for k in model.__dict__ if 'num' == k[:3]]
        attrs_to_keep += ['tseries']

        N_iters = len(argsets) * N_PER

        for argseti, (sname, argset) in enumerate(argsets):

            if verbose:
                print(f'Doing {sname}. {argseti}/{len(argsets)}. {N_PER - len(self.models[sname])} to do.')

            #print(sname, argset)

            # parameters are exploded from within `argset`
            strat = argset['strat']
            VACCINATE_P = argset['VACCINATE_P']
            INITIAL_INFECT = argset['INITIAL_INFECT']

            param_generation_args = {}
            for K in ['R0_mean','R0_coeffvar']: # you can add others here. I just am lazy
                if K in argset:
                    param_generation_args[K] = argset[K]

            # update parameters, depending on R0_mean, R0_coeffvar, etc.
            #print(param_generation_args)
            self.base_args = self.get_params(generate_new=generate_new,**param_generation_args)

            # bunch of bullshit hacking bullshit
            nodes = list(self.base_args['G'].nodes)
            edges = list(self.base_args['G'].edges(data='weight'))
            self.tnet = weightedNetwork(nodes, edges)
            sim = simulations.SEIR_daily(self.tnet, p)
            n2id = {n:i for i,n in enumerate(self.tnet.nodes)}
            
            # figure out what the default alpha should be...
            vaccinate_alpha = 0
            if 'vaccinate_alpha' in argset:
                vaccinate_alpha = argset['vaccinate_alpha']
            
            #if (argseti+1)%int(len(argsets)/100) == 0:
            first_it = True

            for i in range(N_PER - len(self.models[sname])):

                ns_to_remove = strat(sim, int(self.Nn*VACCINATE_P))
                ns_to_remove = [n2id[x] for x in ns_to_remove]

                ALPHA = [0.5 if i not in ns_to_remove else vaccinate_alpha for i in range(self.Nn)]

                args = dict(self.base_args)
                args['alpha'] = ALPHA
                #args['initI_asym'] = INITIAL_INFECT

                model = ExtSEIRSNetworkModel(**args)
                self.most_recent_model = model

                # reconstruct X myself
                ns_to_infect = sample( [i for i in range(self.Nn) if i not in ns_to_remove], INITIAL_INFECT )

                for i in range(self.Nn):
                    if i in ns_to_infect:
                        model.X[i] = model.I_sym

                # hard-code this to avoid ridiculousness
                # no documentation on all those options sooo.
                # ignore the type
                model.A_deltabeta = type(model.A)( np.multiply( model.A.toarray(), model.beta ) )

                # sanity checks
                assert all( 0<=i<=self.Nn for i in ns_to_remove )
                assert all( 0<=i<=self.Nn for i in ns_to_infect )

                st = time()
                model.run(100,verbose=False)
                if first_it:
                    t_iter = time()-st
                    if verbose:
                        print(f'{t_iter:0.1f} seconds for first simulation. At that rate it\'ll take {N_iters*t_iter/3600:0.1f} hours.')
                    first_it = False

                self.models[sname].append({
                    k: getattr(model, k)
                    for k in attrs_to_keep
                })

    def final_vals(self, typ, var):
        return [m[var][-1] for m in self.models[typ]]
    
    def summarize_params(self):
        # parameters sanity check
        names = "A_deltabeta,A_delta_pairwise,A_beta_pairwise,beta_local,beta".split(",")
        sq = int(np.ceil(np.sqrt(len(names))))
        plt.figure(figsize=(10,5))

        model = self.most_recent_model

        for xi,x in enumerate(names):
            plt.subplot(sq,sq,xi+1)
            attr = getattr(model, x)
            if type(attr) == scipy.sparse.csr_matrix:
                attr = attr.toarray()
            plt.hist(attr.flatten(), bins=20, log=True);
            plt.title(x)
        plt.show()

    def summarize_models(self):
        plt.figure(figsize=(20,20))

        nsqrs = int(np.ceil( np.sqrt(len(self.models)) ))

        for ki,k in enumerate(self.models):
            plt.subplot(nsqrs, nsqrs, ki+1)
            for i in range( min(30,len(self.models[k])) ):
                m = self.models[k][i]

                plt.plot(m['tseries'], m['numNodes'] - m['numS'], label=i, linewidth=1, color='black', alpha=0.5)
            #plt.legend()
            plt.title(k, fontsize=8)
            plt.ylim(0,m['numNodes'])

        plt.show()

    def load_models(self, NAME):

        def decode_argset(X):
            args_decode = dict(X)
            for K in ['G', 'G_Q']:
                Gb = args_decode[K]
                G = nx.Graph()
                
                # add nodes. you have to do this separately, as there are sometimes disconnected nodes
                G.add_nodes_from(Gb['nodes'])

                # add edges // there's some bug where the weight is sometimes None?? only observed on synthetic network
                edj = Gb['edges']
                edj = [(a,b,c if c is not None else 1) for (a,b,c) in edj]
                G.add_weighted_edges_from(edj)

                # replace the value in args
                args_decode[K] = G

            return args_decode
            
        with open(f'simulation_results/{NAME}.pickle', 'rb') as outf:
            self.models = pickle.load(outf)

        with open(f'simulation_results/{NAME}.args.json', 'r', encoding='utf8') as inf:
            argsets = json.load(inf)

            # making accomodations for my old style
            if type(argsets) == dict:
                argsets_new = {}
                for k,v in argsets.items():
                    if type(k) in [str, float]:
                        k = fzd({"R0_mean": float(k)})
                    argsets_new[k] = decode_argset(v)
                argsets = argsets_new
            # this is the new style
            elif type(argsets) == list:
                argsets = {fzd(k):decode_argset(v) for k,v in argsets}

            self.base_arg_sets = argsets

        self.G = list(self.base_arg_sets.values())[0]['G']
        self.init_network()


            

    def dump_models(self, NAME):
        with open(f'simulation_results/{NAME}.pickle', 'wb') as outf:
            pickle.dump(self.models, outf)

        # bleh...
        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if hasattr(obj, 'tolist'):
                    return obj.tolist()
                if type(obj) == nx.classes.graph.Graph:
                    return 'graph'
                return json.JSONEncoder.default(self, obj)
            
        def encode_argset(X):
            args_encode = dict(X)
            for K in ['G','G_Q']:
                Gb = args_encode[K]
                args_encode[K] = {
                    'nodes':list(Gb.nodes),
                    'edges':list( Gb.edges(data="weight") )
                }
            return args_encode

        with open(f'simulation_results/{NAME}.args.json', 'w', encoding='utf8') as outf:
            arg_dump = [
                (dict(k), encode_argset(v)) # keys are dictionaries (frozendicts) of parameter values now... 
                for k,v in self.base_arg_sets.items()
            ]
            #print(arg_dump)
            json.dump(arg_dump, outf, cls=NumpyEncoder)

    def _getQuantileBounds(self, l, q, z=1.96):
        l = sorted(l)
        n = len(l)

        i = int( np.ceil( n*q - z*np.sqrt(n*q*(1-q)) ) )
        j = int( np.ceil( n*q + z*np.sqrt(n*q*(1-q)) ) )

        return (l[i], l[j])

    def info(self, k):

        if type(k) != fzd:
            if type(k) != dict:
                raise Exception('get info for dict or frozendict!')
            k = fzd(k)

        ks = iter_args(k)
        if len(ks) > 1:
            return [ self.info(k) for k,_ in ks ]
        else:
            k = ks[0][0]


        row = {}

        row['Nn'] = list(self.models[k])[0]['numNodes'] # assuming constant
        
        if k['strat'] == 'none':
            row['N_vacc'] = 0            
        else:
            row['N_vacc'] = int(row['Nn']*k['VACCINATE_P'])

        row['init_sus'] = row['Nn'] - k['INITIAL_INFECT'] - row['N_vacc']
        
        ninf = [ (row['Nn'] - k['INITIAL_INFECT'] - m['numS'])[-1] for m in self.models[k] ]
        row['ninf'] = np.array(ninf)
        #print(k)
        #print(row['init_sus'], [m['numS'][-1] for m in self.models[k]], [m['numS'][1] for m in self.models[k]])
        #print(ninf)
        #print("")
        
        for kk,vv in k.items():
            row[kk] = vv

        row['mean'] = np.mean( ninf )
        row['std'] = np.std( ninf )
        row['stderr'] = np.std( ninf ) / np.sqrt( len(ninf) )
        row['Nsims'] = len(ninf)

        for Q in [25, 50, 75]:
            row[f'{Q}p'] = np.quantile( ninf, Q/100 )
            row[f'{Q}pL95'], row[f'{Q}pH95'] = self._getQuantileBounds( ninf, Q/100 )

        row[f'min'] = np.min(ninf)
        row[f'max'] = np.max(ninf)
        
        row['P_sus_inf'] = row['mean'] / row['init_sus']
        
        none_k = dict(k)
        none_k['strat'] = 'none'
        none_k = fzd(none_k)

        rand_k = dict(k)
        rand_k['strat'] = 'rand'
        rand_k = fzd(rand_k)

        # compare to rand for all except none and rand strategies :P
        if k['strat'] not in ['none', 'rand']:
            
            if rand_k in self.models and len(self.models[rand_k]):
                rand_info = self.info(rand_k)                

                row['P_sus_inf_rel_rand'] = row['P_sus_inf'] / rand_info['P_sus_inf']
                row['P_sus_inf_rel_rand_err'] = (row['stderr'] / row['init_sus']) / rand_info['P_sus_inf'] # should I bootstrap?
                
                row['P_sus_inf_rel_rand_BS'] = bootstrap( row, rand_info )
                
        # compare to none for all except none strategy :P
        if k['strat'] not in ['none']:

            if none_k in self.models and len(self.models[none_k]):
                none_info = self.info(none_k)

                row['P_sus_inf_rel_none'] = row['P_sus_inf'] / none_info['P_sus_inf']
                row['P_sus_inf_rel_none_err'] = (row['stderr'] / row['init_sus']) / none_info['P_sus_inf'] # should I bootstrap? lol... months later, I say yes. 20m of work, discovering these error bars are totally wrong.
                
                row['P_sus_inf_rel_none_BS'] = bootstrap( row, none_info )

            
        #if np.sum( np.array(ninf) > row['init_sus'] ) != 0:
        #    print(np.array(ninf), row['init_sus'], row['Nn'], row['VACCINATE_P'])
        #    print(k)

        #assert( np.sum( np.array(ninf) > row['init_sus'] ) == 0 )

        #if np.sum( np.array(ninf) < 0 ) != 0:
        #    print(ninf)
        assert( np.sum( np.array(ninf) < 0 ) == 0 )

        return row


# ============================================
# example usage
# ============================================

if __name__ == '__main__':
    strats = [
        sampling.none,
        #sampling.rand
    ]

    argset = iter_args({
        'strat': {x.__name__: x for x in strats},
        'R0_mean': [0.1, 0.5, 1,2.5,3.5,5],
        'VACCINATE_P': 0.20,
        'INITIAL_INFECT': 20
    })

    manager = simulation_manager()
    manager.load_net('high school 1')

    if False:
        weights = [w for (u,v,w) in manager.G.edges.data('weight') ]
        plt.hist(weights, bins=30)
        plt.show()

    manager.run(argset, N_PER=50)

    manager.summarize_models()

    for k in manager.models:
        print(k, np.mean( [x['numS'][-1] for x in manager.models[k]] ))

    """




    def friendHighDegRandTop3(sim, vaccinateN):
        return sampling.friendHighDegRandTopN(sim, vaccinateN, N=3)
    def friendHighDegRandTop5(sim, vaccinateN):
        return sampling.friendHighDegRandTopN(sim, vaccinateN, N=5)

    strats = [
        sampling.friendClose,
        sampling.friendHighDegClose,
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
    ]

    # this is one way, but it's too many permutations.
    # no way I actually need this many experiments...

    argsets = iter_args({
        'strat': {x.__name__: x for x in strats},
        'R0_mean': [1, 2.5, 4],
        'VACCINATE_P': [0.05, 0.10, 0.20, 0.50],
        'INITIAL_INFECT': [5, 10, 20]
    })

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

    argsets = argset1 + argset2 + argset3
    len(argsets)
    """