from .common_imports import *
from .networks import unweightedNetwork, weightedNetwork, temporalNetwork

#from seirsplus.models import SEIRSNetworkModel


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

        # construct daily counts...
        daily_counts = defaultdict(lambda:defaultdict(int))

        if isinstance(self.tnet, weightedNetwork):
            for f,net in self.tnet.ego_edges.items():
                for t,weight in net.items():
                    f = int(f); t = int(t);
                    daily_counts[0][ (f,t) ] = weight

        elif isinstance(self.tnet, temporalNetwork):
            # rescale shit
            
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
        
    def state_change(self, n, state):
        from numpy.random import weibull

        #ni = list(self.tnet.G.nodes).index(n)
        ni = n

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
        #a = self.tnet.nodes.index(a)
        #b = self.tnet.nodes.index(b)

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

            a = self.tnet.nodes.index(a)
            b = self.tnet.nodes.index(b)

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
        for pers, pname in enumerate(self.tnet.nodes):
            if self.exp[pers]:
                if callable(self.e2i):
                    if self.T >= self.e2i_T[ pers ]:
                        self.state_change(pers, 'inf');
                        del self.e2i_T[ pers ];
                else:
                    if random() < 1 - np.exp( -self.e2i * time_elapsing ):
                    #if random() < 1 - np.power(1-self.e2i, time_elapsing):
                        self.state_change(pers, 'inf')
                
                        
    def roll_for_recovery(self, time_elapsing):
        for pers, pname in enumerate(self.tnet.nodes):
            if self.inf[pers]:
                if callable(self.i2r):
                    if self.T >= self.i2r_T[ pers ]:
                        self.state_change(pers, 'rec');
                        del self.i2r_T[ pers ];
                else:
                    if random() < 1 - np.exp( -self.i2r * time_elapsing ):
                    #if random() < 1 - np.power(1-self.i2r, time_elapsing):
                        # recover them!
                        self.state_change(pers, 'rec')
                    


def SEIR_gillespie():

    def transform_files():
        first = None

        def transform(x):
            global first
            sp = x.split(" ")
            if not len(sp)>1:
                return ""
            
            if first is None:
                first = int(sp[0])
            
            sp[0] = str(int( sp[0] ) - first)
            
            return " ".join(sp)

        # zero_out_the_file...
        dta_start = Path("TemporalGillespieAlgorithm","High-School_data_2013.csv").open('r').read().split("\n")

        with Path("TemporalGillespieAlgorithm","High-School-transform_data_2013.csv").open('w') as outf:
            outf.write( "\n".join( map(transform, dta_start) ) )

    def get_results():
        I_p = list(Path("TemporalGillespieAlgorithm").glob("*I_t*High-school-transform*"))[0].open('r').read().split("\t")
        R_p = list(Path("TemporalGillespieAlgorithm").glob("*R_t*High-school-transform*"))[0].open('r').read().split("\t")

        I_p = map(float, I_p[:-1])
        R_p = map(float, R_p[:-1])

        I_p = np.array(list(I_p))
        R_p = np.array(list(R_p))

class SEIR:

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

        # construct daily counts...
        daily_counts = defaultdict(lambda:defaultdict(int))

        if isinstance(self.tnet, weightedNetwork):
            for f,net in self.tnet.ego_edges.items():
                for t,weight in net.items():
                    f = int(f); t = int(t);
                    daily_counts[0][ (f,t) ] = weight

        elif isinstance(self.tnet, temporalNetwork):
            # rescale shit
            
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
        
        dbrk = self.tnet.day_breaks
        
        # where are we now
        current_day_i = max([i for i in range(len(dbrk)) if dbrk[i] <= self.T]) + (self.T // self.tnet.Tmax) * 5
        new_day_i = current_day_i + duration
        
        new_steps = self.tnet.Tmax * (new_day_i // 5) + dbrk[new_day_i % 5]
        n_steps = new_steps - self.T

        
        if self.rescale_step_factor != 1:
            if self.rescale_step_factor < 0:
                raise Exception("You can only rescale by some integer multiple")
                
            self.rescale_step_factor = int(self.rescale_step_factor)
            n_steps = (n_steps // self.rescale_step_factor)
        
        for s in range(n_steps):
            self.step()

                    
            # once per time interval
            self.measure()
            self.T += self.rescale_step_factor

    def step(self):

        if isinstance(self.tnet, weightedNetwork):
            if self.T in self.tnet.times:
                for a, b, w in self.tnet.edges:
                    if a >= b: continue # only roll once for each edge

                    a = self.tnet.nodes.index(a)
                    b = self.tnet.nodes.index(b)

                    if not self.relevant_pair(a,b): continue

                    inf = a if self.inf[a] else b
                    sus = a if self.sus[a] else b

                    
                    if random() < 1 - np.exp( -self.s2e * w * self.rescale_step_factor ):
                    #if random() < 1 - np.power( 1-self.s2e, w * self.rescale_step_factor ):
                        # infect sus
                        self.state_change(sus, 'exp')


        elif isinstance(self.tnet, temporalNetwork):
            for a, b in self.tnet.t_edges[ self.T % self.tnet.Tmax ]:
                if a >= b: continue # only roll once for each edge

                a = self.tnet.nodes.index(a)
                b = self.tnet.nodes.index(b)

                if not self.relevant_pair(a,b): continue

                inf = a if self.inf[a] else b
                sus = a if self.sus[a] else b

                # now roll the dice
                if random() < self.s2e:
                    # infect sus
                    self.state_change(sus, 'exp')

                    
        self.roll_for_recovery()
        self.roll_for_e2i()
            
    def measure(self):
        for s in self.states:
            self.meas[ s ].append( getattr(self, s).sum() )

    def roll_for_e2i(self):
        for pers in range(self.tnet.Nnodes):
            if self.exp[pers]:
                if random() < self.transform_prob( self.e2i, self.rescale_step_factor ):
                    # infect them!
                    self.state_change(pers, 'inf')
                    
    def roll_for_recovery(self):
        for pers in range(self.tnet.Nnodes):
            if self.inf[pers]:
                #if random() < 1 - np.power(1-self.i2r, self.rescale_step_factor):
                if random() < 1 - np.power(1-self.i2r, self.rescale_step_factor):
                    # recover them!
                    self.state_change(pers, 'rec')
           