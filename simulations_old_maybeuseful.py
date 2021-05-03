class epiSim:
    
    edgelist_cache = {}
    
    def __init__(self, tnet, alpha, beta, rescale_step_factor=1):
        self.tnet = tnet
        
        self.alpha = alpha
        self.beta = beta
        self.rescale_step_factor = rescale_step_factor
        
        self.init_attributes()
        
        self.state_changes = []
        
    def subc_init_attributes(self):
        pass
        
    def init_attributes(self):
        self.susceptible = np.ones(self.tnet.Nnodes)  # susceptible = sus_over_time[ time ]
        self.infected = np.zeros(self.tnet.Nnodes) #infected = inf_over_time[ time ]
        self.recovered = np.zeros(self.tnet.Nnodes) # recovered = rec_over_time[ time ]
        self.vaccinated = np.zeros(self.tnet.Nnodes) # recovered = rec_over_time[ time ]
        
        self.meas = defaultdict(list)
        self.T = 0
        
        self.subc_init_attributes()
        
    def infect(self, ni):
        self.infected[ni] = 1
        self.susceptible[ni] = 0
        self.recovered[ni] = 0
        
        self.state_changes.append( (self.T, ni, 'inf') )
        
    def recover(self, ni):
        self.infected[ni] = 0
        self.susceptible[ni] = 0
        self.recovered[ni] = 1
        
        self.state_changes.append( (self.T, ni, 'rec') )
        
    def vaccinate(self, ni):
        self.infected[ni] = 0
        self.susceptible[ni] = 0
        self.recovered[ni] = 0
        self.vaccinated[ni] = 1
        
        self.state_changes.append( (self.T, ni, 'vacc') )
        
        
    def stateT(a, node, time):
        # for historical inspection of network evolution
        
        my_c = [ (t, s) for (t,n,s) in a.state_changes if n == node and t <= time ]
        if not len(my_c):
            return 'sus'

        return my_c[-1][1]
    
    def relevant_pair(self, a, b):

        if self.recovered[a] or self.recovered[b]: # if either are recovered
            return False

        if self.vaccinated[a] or self.vaccinated[b]: # if either are vaccinated
            return False

        if not(self.infected[a] or self.infected[b]): # if neither are infected
            return False

        if self.infected[a] and self.infected[b]: # if both are infected
            return False

        if self.susceptible[a] + self.susceptible[b] == 0: # if neither are susceptible
            return False
        
        return True
    
    def transform_prob(self, orig, scale):
        return (1-(1-orig)**scale)

    def run_gillespie(self, duration=1, NSIM=1e5):
        assert(duration > 0)
        assert(type(duration) == int) # number of days
        
        if not isinstance(self.tnet, temporalNetwork):
            raise Exception('Gillespie algo has only been considered for temporal networks...')

        dbrk = self.tnet.day_breaks
        current_day_i = max([i for i in range(len(dbrk)) if dbrk[i] <= self.T]) + (self.T // self.tnet.Tmax) * 5
        new_day_i = current_day_i + duration

        tmp_dir = Path('.','temp')
        tmp_dir.mkdir(exist_ok=True)
        tmp_edges = tmp_dir.joinpath('tmp.txt')
        
        new_steps = self.tnet.Tmax * (new_day_i // 5) + dbrk[new_day_i % 5]
        n_steps = new_steps - self.T
        
        with tmp_edges.open('w') as outf:

            for s in range(n_steps):

                for a, b in self.tnet.t_edges[ self.T % self.tnet.Tmax ]:

                    outf.write("%s\n"%( "\t".join([
                        self.T, a, b
                    ])))

                self.T += 1

        from subprocess import Popen, PIPE
        p = Popen("C:/cygwin64/bin/mintty.exe", stdin=PIPE, stdout=PIPE)

        cmds = [
            "export PATH=$PATH:/cygdrive/c/Users/amcga/win-builds/bin",
            'cd "G:/My Drive/2020 ORGANISATION/0. right now right now/disease contact spread/0 analysis/TemporalGillespieAlgorithm"',
            #'g++ ./SIR-Poisson-homogeneous.cpp -o "poo" -O2 -I"C:/Users/amcga/Downloads/boost_1_75_0/boost_1_75_0"',
            './poo2.exe "%s" 1 %0.10f %0.10f 20000 %d 15' % (tmp_edges, self.alpha, self.beta, NSIM)
        ]

        p.stdin.write()

                

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
                    if not self.relevant_pair(a,b): continue

                    inf = a if self.infected[a] else b
                    sus = a if self.susceptible[a] else b

                    # now roll the dice
                    if random() < self.transform_prob( self.alpha, w * self.rescale_step_factor ):
                        # infect sus
                        self.infect(sus)
                        self.infectT.append( (sus, self.T) )


        elif isinstance(self.tnet, temporalNetwork):
            for a, b in self.tnet.t_edges[ self.T % self.tnet.Tmax ]:
                if a >= b: continue # only roll once for each edge
                if not self.relevant_pair(a,b): continue

                inf = a if self.infected[a] else b
                sus = a if self.susceptible[a] else b

                # now roll the dice
                if random() < self.alpha:
                    # infect sus
                    self.infect(sus)
                    self.infectT.append( (sus, self.T) )

                    
        self.roll_for_recovery()
            
    def measure(self):
        self.meas['inf'].append( self.infected.sum() )
        self.meas['rec'].append( self.recovered.sum() )
        self.meas['sus'].append( self.susceptible.sum() )

                    
    def roll_for_recovery(self):
        for pers in range(self.tnet.Nnodes):
            if self.infected[pers]:
                if random() < self.transform_prob( self.beta, self.rescale_step_factor ):
                    # recover them!
                    self.recover(pers)
                    

class SEIR:

    states = [
        'inf', 'sus', 'rec', 'exp', 'vacc'
    ]

    def __init__(self, tnet, alpha, beta, e2i, rescale_step_factor=1):
        self.tnet = tnet
        
        self.alpha = alpha
        self.beta = beta
        self.e2i = e2i
        self.rescale_step_factor = rescale_step_factor
        
        self.init_attributes()
        
        self.state_changes = []
        
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
        
        self.subc_init_attributes()
        
    def state_change(self, ni, state):
        assert(state in self.states)
        for s in self.states:
            getattr(self, s)[ni] = int(s == state)

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
                    if not self.relevant_pair(a,b): continue

                    inf = a if self.inf[a] else b
                    sus = a if self.sus[a] else b

                    # now roll the dice
                    if random() < self.transform_prob( self.alpha, w * self.rescale_step_factor ):
                        # infect sus
                        self.state_change(sus, 'exp')


        elif isinstance(self.tnet, temporalNetwork):
            for a, b in self.tnet.t_edges[ self.T % self.tnet.Tmax ]:
                if a >= b: continue # only roll once for each edge
                if not self.relevant_pair(a,b): continue

                inf = a if self.inf[a] else b
                sus = a if self.sus[a] else b

                # now roll the dice
                if random() < self.alpha:
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
                    # recover them!
                    self.state_change(pers, 'inf')
                    
    def roll_for_recovery(self):
        for pers in range(self.tnet.Nnodes):
            if self.inf[pers]:
                if random() < self.transform_prob( self.beta, self.rescale_step_factor ):
                    # recover them!
                    self.state_change(pers, 'rec')
           