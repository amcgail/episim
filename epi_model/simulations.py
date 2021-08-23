from .common_imports import *
from .networks import weightedNetwork, temporalNetwork

         


class SEIR_daily:

    states = [
        'inf', 'sus', 'rec', 'exp', 'vacc'
    ]

    def __init__(self, tnet, params):
        
        self.params = dict( params,
                           rescale_step_factor = 1 )
        
        self.tnet = tnet
        
        self.s2e = self.params['s2e']
        self.e2i = self.params['e2i']
        self.i2r = self.params['i2r']

        self.e2i_T = {}
        self.i2r_T = {}

        self.rescale_step_factor = self.params['rescale_step_factor']
        
        daily_counts = defaultdict(lambda:defaultdict(int))
        for time, f, t in self.tnet.edgelist:
            if f >= t:
                continue

            day = np.sum( time >= np.array(self.tnet.day_breaks) ) - 1
            f = int(f); t = int(t);
            daily_counts[ day ][ (f, t) ] += 1

        self.daily_counts = daily_counts
        
        self.init_attributes()
        
    def subc_init_attributes(self):
        pass
        
    def init_attributes(self):
        for s in self.states:
            if s == 'sus':
                val = np.ones(self.tnet.Nnodes)
            else:
                val = np.zeros(self.tnet.Nnodes)
            setattr(self, s, val)
        
        self.meas = defaultdict(list)
        self.T = 0

        self.e2i_T = {}
        self.i2r_T = {}
        
        self.subc_init_attributes()
        self.state_changes = []
        
    def state_change(self, ni, state):
        from numpy.random import weibull

        assert(state in self.states)
        for s in self.states:
            getattr(self, s)[ni] = int(s == state)
        
        if state == 'exp' and callable(self.e2i): # draw from distribution of e2i times to plan...
            self.e2i_T[ni] = self.e2i() + self.T
            
        if state == 'inf' and callable(self.i2r): # draw from distribution of i2r times to plan...
            self.i2r_T[ni] = self.i2r() + self.T

        self.state_changes.append( (self.T, ni, state) )
        

    def stateT(self, node, time):
        # for historical inspection of network evolution
        
        my_c = [ (t, s) for (t,n,s) in self.state_changes if n == node and t <= time ]
        if not len(my_c):
            return 'sus'

        return my_c[-1][1]
    
    def transform_prob(self, orig, scale):
        return (1-(1-orig)**scale)

    def relevant_pair(self, a, b):

        if self.rec[a] or self.rec[b]: # if either are recovered
            return False

        if self.vacc[a] or self.vacc[b]: # if either are vaccinated
            return False

        if not(self.inf[a] or self.inf[b]): # if neither are infected
            return False

        if self.inf[a] and self.inf[b]: # if both are infected
            return False

        if self.sus[a] + self.sus[b] == 0: # if neither are susceptible
            return False
        
        return True

    def run(self, duration=1): # only for weighted):
        assert(duration > 0)
        assert(type(duration) == int) # number of days
        
        for s in range(duration):
            self.step()
                    
            # once per time interval
            self.measure()
            self.T += self.rescale_step_factor

    def step(self):
        current_day_i = self.T % len(self.tnet.days)
        for (a, b), num_k in self.daily_counts[ current_day_i ].items():
            if a >= b: continue # only roll once for each edge
            if not self.relevant_pair(a,b): continue

            inf = a if self.inf[a] else b
            sus = a if self.sus[a] else b
            
            #print('here', 1 - np.exp( -self.alpha * num_k ))

            # now roll the dice
            #if random() < 1 - np.exp( -self.alpha * num_k ):
            if random() < 1 - np.power( 1-self.s2e, num_k ):
                # infect sus
                self.state_change(sus, 'exp')
        
        next_day = (current_day_i+1) % len(self.tnet.days)
        tstart_nextday = self.tnet.mindayT[next_day]
        tend_today = self.tnet.maxdayT[current_day_i]
        time_elapsing = (tstart_nextday - self.tnet.day_breaks[next_day]) + \
            self.tnet.day_breaks[next_day+1] - tend_today

        # 3600*24/20
        self.roll_for_recovery( time_elapsing ) # in 20s intervals
        self.roll_for_e2i( time_elapsing )
            
    def measure(self):
        for s in self.states:
            self.meas[ s ].append( getattr(self, s).sum() )

    def roll_for_e2i(self, time_elapsing):
        for pers in range(self.tnet.Nnodes):
            if self.exp[pers]:
                if callable(self.e2i):
                    if self.T >= self.e2i_T[ pers ]:
                        self.state_change(pers, 'inf');
                        del self.e2i_T[ pers ];
                else:
                    #if random() < 1 - np.exp( -self.e2i * time_elapsing ):
                    if random() < 1 - np.power(1-self.e2i, time_elapsing):
                        self.state_change(pers, 'inf')
                
                        
    def roll_for_recovery(self, time_elapsing):
        for pers in range(self.tnet.Nnodes):
            if self.inf[pers]:
                if callable(self.i2r):
                    if self.T >= self.i2r_T[ pers ]:
                        self.state_change(pers, 'rec');
                        del self.i2r_T[ pers ];
                else:
                    #if random() < 1 - np.exp( -self.beta * 3600*24/20 ):
                    if random() < 1 - np.power(1-self.i2r, time_elapsing):
                        # recover them!
                        self.state_change(pers, 'rec')
                    


