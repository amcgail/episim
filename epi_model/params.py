from .common_imports import *

# for those not relying on contact, the situation is easier
N = 3600*24 / 20
B = 1 - np.power( 1 - 0.09, 1/N )
B

mean_rec = (3600*24/20) * 14
mean_e2i = (3600*24/20) * 10

covid_estimate = {
    's2e': 6e-04,
    'e2i': 1/mean_e2i,
    'i2r': 1/mean_rec, 
}



def daily_to_momentary( net, param ):

    contact_count = Counter([x[1] for x in net.edgelist if x[0] < net.day_breaks[1]])
    contact_amts = [x for x in contact_count.values() if x >= 15*60 / 20] # discard < 15 minutes

    def mean_inf_prob( beta ):
        probs = [1 - np.power(1-beta, c) for c in contact_amts]
        return np.mean(probs)

    from scipy.optimize import fsolve
    def to_solve( cutoff ):
        def woop(beta):
            return mean_inf_prob( beta ) - cutoff
        return woop

    guess = 5e-6
    sol = fsolve(to_solve(param), guess)
    assert(sol.shape[0]==1)

    return sol[0]

def e2i_i2r_correlated_a():
    # [here](https://www.medrxiv.org/content/10.1101/2020.09.04.20188516v1) for the relationship between symptoms arriving and infectiousness
        
    symp = lognormal(1.63, 0.5)
    
    # assume they're just constants. be nice
    inf = max(0,symp-2)
    inf_end = inf+10

def e2i_i2r_correlated_b():
    # [here](https://www.medrxiv.org/content/10.1101/2020.09.04.20188516v1) for the relationship between symptoms arriving and infectiousness
        
    sd = 2.8
    mu = -0.07
    df = 2*(sd**2)/(sd**2 - 1)
    TOT = tdist.rvs( df ) + mu
    
    symp = lognormal(1.63, 0.5)
    
    # assume they're super complicated. be mean
    inf = TOT+symp
    inf_end = inf+10

def e2i_lognormal():
    # From [this paper](https://bmjopen.bmj.com/content/10/8/e039652): 
    # Incubation period of COVID-19: a rapid systematic review and meta-analysis of observational research

    from numpy.random import lognormal
    return lognormal(1.63, 0.5)

def e2i_weibull():
    # this is only approximate, I can't figure out the real parameters.

    lamb_sig = 4.4
    k = 1.04
    return lamb_sig * weibull(k)


__all__ = [
    'covid_estimate',
    'daily_to_momentary',
]