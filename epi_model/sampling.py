from .common_imports import *

__all__ = ['sampling_strats']

class early_infect_analysis:
    def __init__(self, sim):
        from epi_model import epiSim
        from random import choice
        
        # run 100 simulations...
        NN = 100

        for i in range(NN):
            sim.init_attributes()
            
            if (i+1)%(NN//10)==0:
                print(i,"/",NN)
            sim.state_change(choice(range(sim.tnet.Nnodes)), 'inf')
            sim.run(20)
            
        all_infect_t = inf_times = [
            ( sc[1], sc[0] )
            for sc in sim.state_changes 
            if sc[0] != 0 and\
            sc[2] == 'inf'
        ]
        
        self.means = {
            n: np.mean([t for nn,t in all_infect_t if nn==n])
            for n in range(sim.tnet.Nnodes)
        }
        self.counts = {
            n: np.sum([1 for nn,t in all_infect_t if nn==n])
            for n in range(sim.tnet.Nnodes)
        }

def early_infect(sim, vaccinateN=100):
    if not hasattr(sim, 'early_infect'):
        sim.early_infect = early_infect_analysis(sim)
    
    top100a = sorted(range(sim.tnet.Nnodes), key=lambda i:-sim.early_infect.counts[i])[:vaccinateN*2]
    top100b = sorted(top100a, key=lambda i:sim.early_infect.means[i])
    tovaccinate = top100b[:vaccinateN]

    return tovaccinate



def none(*args, **kwargs):
    return []

def rand(sim, vaccinateN=100):
    from random import choice, sample
    return sample(range(sim.tnet.Nnodes), vaccinateN)

def friendHighDegChain(sim, vaccinateN=100):
    #a = epiSim(edgelist, alpha=alpha, beta=beta)
    start = sample(range(sim.tnet.Nnodes), vaccinateN)
    deg = nx.degree(sim.tnet.G)
    
    chosen = []
    current = choice(range(sim.tnet.Nnodes))
    
    ii = 0
    while len(chosen) < vaccinateN:
        chosen.append( current )

        ii += 1
        if ii > 100000:
            raise Exception("Excessive recursion...")
    
        p1 = current
        
        p2s = list(nx.neighbors(sim.tnet.G,p1))
        p2s = [x for x in p2s if x not in chosen]
        
        if not len(p2s):
            #print('no choices, skipped')
            current = choice([i for i in range(sim.tnet.Nnodes) if i not in chosen])
        else:
            current = sorted(p2s, key=lambda x:deg[x])[-1]
     
    return chosen
    #sim.tnet.infect(choice([x for x in range(sim.tnet.Nnodes) if not sim.tnet.vaccinated[x]]))
    #sim.tnet.run()

def friendHighDegRandTopN(sim, vaccinateN=100, N=5):
    
    start = sample(range(sim.tnet.Nnodes), vaccinateN)
    deg = nx.degree(sim.tnet.G)
    
    chosen = []
    
    ii = 0
    while len(chosen) < vaccinateN:
        ii += 1
        if ii > 100000:
            raise Exception("Excessive recursion...")
    
        p1 = choice(range(sim.tnet.Nnodes))
        
        p2s = list(nx.neighbors(sim.tnet.G,p1))
        p2s = [x for x in p2s if x not in chosen]
        
        if not len(p2s):
            #print('no choices, skipped')
            continue
        
        p2 = sorted(p2s, key=lambda x:deg[x])[-1]
        chosen.append(p2)
     
    return chosen
    #sim.tnet.infect(choice([x for x in range(sim.tnet.Nnodes) if not sim.tnet.vaccinated[x]]))
    #sim.tnet.run()
    
def friendHighDeg(sim, vaccinateN=100):
    #a = epiSim(edgelist, alpha=alpha, beta=beta)
    start = sample(range(sim.tnet.Nnodes), vaccinateN)
    deg = nx.degree(sim.tnet.G)
    
    chosen = []
    
    ii = 0
    while len(chosen) < vaccinateN:
        ii += 1
        if ii > 100000:
            raise Exception("Excessive recursion...")
    
        p1 = choice(range(sim.tnet.Nnodes))
        
        p2s = list(nx.neighbors(sim.tnet.G,p1))
        p2s = [x for x in p2s if x not in chosen]
        
        if not len(p2s):
            #print('no choices, skipped')
            continue
        
        p2 = sorted(p2s, key=lambda x:deg[x])[-1]
        chosen.append(p2)
     
    return chosen
    #sim.tnet.infect(choice([x for x in range(sim.tnet.Nnodes) if not sim.tnet.vaccinated[x]]))
    #sim.tnet.run()

def friendHighDegCloseChain(sim, vaccinateN=100, CUTOFF=0.9):
    cutoff = np.quantile([x[2] for x in sim.tnet.weighted.edges], CUTOFF)
    deg = {
        n: len([_ for x,w in sim.tnet.weighted.ego_edges[n].items() if w >= CUTOFF])
        for n in range(sim.tnet.Nnodes)
    }
    
    chosen = []
    current = choice(range(sim.tnet.Nnodes))
    
    ii = 0
    while len(chosen) < vaccinateN:
        chosen.append(current)
        
        ii += 1
        if ii > 100000:
            raise Exception("Excessive recursion...")
            
        p2s = list(nx.neighbors(sim.tnet.G,current))
        p2s = [x for x in p2s if x not in chosen if sim.tnet.weighted.ego_edges[p1][x] >= cutoff]
        
        if not len(p2s):
            current = choice([i for i in range(sim.tnet.Nnodes) if i not in chosen])
        else:
            current = sorted(p2s, key=lambda x:deg[x])[-1]
     
    return chosen
    #sim.tnet.infect(choice([x for x in range(sim.tnet.Nnodes) if not sim.tnet.vaccinated[x]]))
    #sim.tnet.run()


def friendHighDegClose(sim, vaccinateN=100, CUTOFF=0.9):
    cutoff = np.quantile([x[2] for x in sim.tnet.weighted.edges], CUTOFF)
    deg = {
        n: len([_ for x,w in sim.tnet.weighted.ego_edges[n].items() if w >= CUTOFF])
        for n in range(sim.tnet.Nnodes)
    }
    
    chosen = []
    
    ii = 0
    while len(chosen) < vaccinateN:
        ii += 1
        if ii > 100000:
            raise Exception("Excessive recursion...")
    
        p1 = choice(range(sim.tnet.Nnodes))
        
        p2s = list(nx.neighbors(sim.tnet.G,p1))
        p2s = [x for x in p2s if x not in chosen if sim.tnet.weighted.ego_edges[p1][x] >= cutoff]
        
        if not len(p2s):
            #print('no choices, skipped')
            continue
        
        p2 = sorted(p2s, key=lambda x:deg[x])[-1]
        chosen.append(p2)
     
    return chosen
    #sim.tnet.infect(choice([x for x in range(sim.tnet.Nnodes) if not sim.tnet.vaccinated[x]]))
    #sim.tnet.run()
    
def targeted(sim, vaccinateN=100):
    deg = nx.degree(sim.tnet.G)
    return sorted(range(sim.tnet.Nnodes), key=lambda x: -deg[x])[:vaccinateN]

def friend(sim, vaccinateN=100):
    chosen = []
    ii = 0
    while len(chosen) < vaccinateN:
        ii += 1
        if ii > 100000:
            raise Exception("Excessive recursion...")
            
        start = choice(range(sim.tnet.Nnodes))
        
        p2s = list(nx.neighbors(sim.tnet.G,start))
        p2s = [x for x in p2s if x not in chosen]
        
        if not len(p2s):
            continue
        
        p2 = choice(p2s)
        chosen.append(p2)
        
    return chosen

def friendWeightedChain(sim, vaccinateN=100):
    chosen = []
    ii = 0

    current = choice(range(sim.tnet.Nnodes))

    while len(chosen) < vaccinateN:
        ii += 1
        if ii > 100000:
            raise Exception("Excessive recursion...")
            
        chosen.append(current)
        
        p2s = list(nx.neighbors(sim.tnet.G,current))
        p2s = [x for x in p2s if x not in chosen]
        
        p2w = np.array([ sim.tnet.weighted.ego_edges[ current ][ p ] for p in p2s ])
        
        if not len(p2s):
            choice([i for i in range(sim.tnet.Nnodes) if i not in chosen])
        else:
            p2 = np.choice(p2s, p=p2w)
        
    return chosen

def friendWeighted(sim, vaccinateN=100):
    chosen = []
    ii = 0
    while len(chosen) < vaccinateN:
        ii += 1
        if ii > 100000:
            raise Exception("Excessive recursion...")
            
        start = choice(range(sim.tnet.Nnodes))
        
        p2s = list(nx.neighbors(sim.tnet.G,start))
        p2s = [x for x in p2s if x not in chosen]
        
        p2w = np.array([ sim.tnet.weighted.ego_edges[ start ][ p ] for p in p2s ])
        
        if not len(p2s):
            continue
        
        p2 = np.choice(p2s, p=p2w)
        chosen.append(p2)
        
    return chosen


def friendChain(sim, vaccinateN=100):
    chosen = []
    ii = 0

    current = choice(range(sim.tnet.Nnodes))

    while len(chosen) < vaccinateN:
        chosen.append(current)

        ii += 1
        if ii > 100000:
            raise Exception("Excessive recursion...")
        
        p2s = list(nx.neighbors(sim.tnet.G,current))
        p2s = [x for x in p2s if x not in chosen]
        
        if not len(p2s):
            # choose a rando, if there are no options...
            p2 = choice([i for i in range(sim.tnet.Nnodes) if not i not in chosen])
        else:
            p2 = choice(p2s)

        current = p2
        
    return chosen

def friendCloseChain(sim, vaccinateN=100, CUTOFF=0.9):
    cutoff = np.quantile([x[2] for x in sim.tnet.weighted.edges], CUTOFF)
    
    chosen = []
    current = choice(range(sim.tnet.Nnodes))
    ii = 0
    while len(chosen) < vaccinateN:
        ii += 1
        if ii > 100000:
            raise Exception("Excessive recursion...")
            
        chosen.append(current)
                    
        p2s = list(nx.neighbors(sim.tnet.G,current))
        p2s = [x for x in p2s if x not in chosen and sim.tnet.weighted.ego_edges[current][x] >= cutoff]
        
        if not len(p2s):
            current = choice([i for i in range(sim.tnet.Nnodes) if i not in chosen])
        else:
            current = choice(p2s)
        
    return chosen

def friendClose(sim, vaccinateN=100, CUTOFF=0.9):
    cutoff = np.quantile([x[2] for x in sim.tnet.weighted.edges], CUTOFF)
    
    chosen = []
    ii = 0
    while len(chosen) < vaccinateN:
        ii += 1
        if ii > 100000:
            raise Exception("Excessive recursion...")
            
        start = choice(range(sim.tnet.Nnodes))
        
        p2s = list(nx.neighbors(sim.tnet.G,start))
        p2s = [x for x in p2s if x not in chosen and sim.tnet.weighted.ego_edges[start][x] >= cutoff]
        
        if not len(p2s):
            continue
        
        p2 = choice(p2s)
        chosen.append(p2)
        
    return chosen


def contact(sim, vaccinateN=100): # needs updated?
        
    tmin = sim.tnet.day_breaks[0]
    tmax = (sim.tnet.day_breaks[1] + sim.tnet.day_breaks[0])/2
    tspan = tmax-tmin
    
    dta = sim.tnet.eldf

    def take_step(samp1):
        from random import random
        samp2 = []
        for s in samp1:
            
            int_times = []
            while not len(int_times):
                t_look = random()*tspan + tmin
                int_times = set(dta[(dta[1] == s) & (dta[0] > t_look)][0])
                int_times = int_times.union(set(dta[(dta[2] == s) & (dta[0] > t_look)][0]))
            next_t = min(int_times)

            interaction = dta[(
                (dta[1] == s)|
                (dta[2] == s)
            ) & (dta[0] == next_t)]
            interaction = interaction.iloc[0]

            alter = (interaction[1]!=s)*interaction[1]
            alter += (interaction[2]!=s)*interaction[2]
            
            #print(next_t, s, alter)
            samp2.append(alter)
        return samp2

    from random import sample
    tovaccinate = take_step(sample(range(sim.tnet.Nnodes), vaccinateN))    

    return tovaccinate
    
        
        
def contact_survey(sim, vaccinateN=100, forward_N=50):
    
    tmin = first_day_S
    tmax = first_day_S + 0.5*(first_day_E-first_day_S) # first half of the first day
    tspan = tmax-tmin

    def take_step(samp1):
        from random import random
        samp2 = []
        for s in samp1:

            s = nodes[s]

            int_times = []

            loop_i = 0
            fail = False
            while 1:

                loop_i += 1
                if loop_i > 10:
                    #print('Failed to find a good time')
                    fail = True
                    break

                t_look = random()*tspan + tmin
                int_times = dta[((dta.a == s)|(dta.b == s)) & (dta.t > t_look)]
                if int_times.shape[0] < forward_N:
                    continue

                break

            if fail:
                continue

            int_times = int_times.head(forward_N)

            ppl = list(int_times.a) + list(int_times.b)
            ppl = Counter(ppl)

            # sort by amount of contact
            ppl = sorted(
                list(ppl.items()),
                key=lambda x: -x[1]
            )
            ppl = ppl[1:] # eliminate "me"
            ppl = [x[0] for x in ppl] # eliminate counts. who cares

            # first person that's not already gonna be vaccinated
            for p in ppl:
                if p in samp2:
                    continue
                samp2.append(p)
                break

        return samp2

    from random import sample
    tovaccinate = []
    while len(tovaccinate):
        tovaccinate = take_step(sample(range(sim.tnet.Nnodes), vaccinateN))

    return tovaccinate


"""
from random import choice, sample

def sim_targeted(vaccinateN=100):
    a = epiSim(edgelist, alpha=alpha, beta=beta)

    deg = nx.degree(a.G)

    for pers in sorted(range(a.Nnodes), key=lambda x: -deg[x])[:vaccinateN]:
        a.vaccinate(pers)
    a.infect(choice(range(a.Nnodes)))
    a.run()
    return a.meas

def sim_friendHighDeg(vaccinateN=100):
    a = epiSim(edgelist, alpha=alpha, beta=beta)
    start = sample(range(a.Nnodes), vaccinateN)
    deg = nx.degree(a.G)
    for p1 in start:
        p2s = list(nx.neighbors(a.G,p1))
        p2s = [x for x in p2s if not a.vaccinated[x]]
        
        if not len(p2s):
            print('no choices, skipped')
            continue
        
        p2 = sorted(p2s, key=lambda x:deg[x])[-1]
        a.vaccinate(p2)
        
    a.infect(choice([x for x in range(a.Nnodes) if not a.vaccinated[x]]))
    a.run()
    return a.meas

def sim_friend(vaccinateN=100):
    a = epiSim(edgelist, alpha=alpha, beta=beta)
    start = sample(range(a.Nnodes), vaccinateN)
    for p1 in start:
        p2s = list(nx.neighbors(a.G,p1))
        p2s = [x for x in p2s if not a.vaccinated[x]]
        
        if not len(p2s):
            #print('no choices, skipped')
            continue
        
        p2 = choice(p2s)
        a.vaccinate(p2)
        
    a.infect(choice([x for x in range(a.Nnodes) if not a.vaccinated[x]]))
    a.run()
    return a.meas
        

first_day_S = dta.t.min()
first_day_E = dta[dta.t < first_day_S+50000].t.max()

def sim_contact(vaccinateN=100):
    a = epiSim(edgelist, alpha=alpha, beta=beta)
    
    
    tmin = first_day_S
    tmax = first_day_S + 0.5*(first_day_E-first_day_S) # first half of the first day
    tspan = tmax-tmin

    def take_step(samp1):
        from random import random
        samp2 = []
        for s in samp1:
            
            s = nodes[s]
        
            int_times = []
            while not len(int_times):
                t_look = random()*tspan + tmin
                int_times = set(dta[(dta.a == s) & (dta.t > t_look)].t)
                int_times = int_times.union(set(dta[(dta.b == s) & (dta.t > t_look)].t))
            next_t = min(int_times)

            interaction = dta[(
                (dta.a == s)|
                (dta.b == s)
            ) & (dta.t == next_t)]
            interaction = interaction.iloc[0]

            alter = (interaction.a!=s)*interaction.a
            alter += (interaction.b!=s)*interaction.b
            samp2.append(alter)
        return samp2

    from random import sample
    tovaccinate = take_step(sample(range(a.Nnodes), 100))    

    for pers in tovaccinate:
        a.vaccinate(nodes.index(pers))
    a.infect(choice(range(a.Nnodes)))
    a.run()
    return a.meas
    

def sim_random(vaccinateN=100):
    from random import choice, sample

    a = epiSim(edgelist, alpha=alpha, beta=beta)

    for pers in sample(range(a.Nnodes), vaccinateN):
        a.vaccinate(pers)
    a.infect(choice(range(a.Nnodes)))
    a.run()
    return a.meas
    

def sim_none(vaccinateN=None):
    from random import choice

    a = epiSim(edgelist, alpha=alpha, beta=beta)
    a.infect(choice(range(a.Nnodes)))
    a.run()
    return a.meas

def sim_contact_survey(vaccinateN=100, forward_N=50):
    from collections import Counter
    a = epiSim(edgelist, alpha=alpha, beta=beta)
    
    
    tmin = first_day_S
    tmax = first_day_S + 0.5*(first_day_E-first_day_S) # first half of the first day
    tspan = tmax-tmin

    def take_step(samp1):
        from random import random
        samp2 = []
        for s in samp1:

            s = nodes[s]

            int_times = []

            loop_i = 0
            fail = False
            while 1:

                loop_i += 1
                if loop_i > 10:
                    #print('Failed to find a good time')
                    fail = True
                    break

                t_look = random()*tspan + tmin
                int_times = dta[((dta.a == s)|(dta.b == s)) & (dta.t > t_look)]
                if int_times.shape[0] < forward_N:
                    continue

                break

            if fail:
                continue

            int_times = int_times.head(forward_N)

            ppl = list(int_times.a) + list(int_times.b)
            ppl = Counter(ppl)

            # sort by amount of contact
            ppl = sorted(
                list(ppl.items()),
                key=lambda x: -x[1]
            )
            ppl = ppl[1:] # eliminate "me"
            ppl = [x[0] for x in ppl] # eliminate counts. who cares

            # first person that's not already gonna be vaccinated
            for p in ppl:
                if p in samp2:
                    continue
                samp2.append(p)
                break

        return samp2

    from random import sample
    tovaccinate = take_step(sample(range(a.Nnodes), vaccinateN))

    for pers in tovaccinate:
        a.vaccinate(nodes.index(pers))
    a.infect(choice(range(a.Nnodes)))
    a.run()
    return a.meas
"""



def target_between_classes(sim, vaccinateN=100):
    new_el = [x for x in sim.tnet.edgelist if sim.tnet.node_attr['class'][x[1]] != sim.tnet.node_attr['class'][x[2]]]
    edge_w = Counter( [tuple(x[1:]) for x in new_el if x[1] < x[2]] )

    CUTOFF = 15 * 60 / 20

    wts = {}
    for node in range(sim.tnet.Nnodes):
        wts[node] = len([ 1 for x, c in edge_w.items() if c>=CUTOFF and x[0]==node])

    to_vacc = sorted( wts, key=lambda x: -wts[x] )[:vaccinateN]
    return to_vacc

def local_betweenness(sim, vaccinateN=100):
    wts = {}
    for n in range(sim.tnet.Nnodes):
        subg = nx.subgraph( sim.tnet.G, nx.neighbors(sim.tnet.G, n) )
        nn = len(subg.nodes)
        pair_diff = nn*(nn-1) - len(subg.edges)
        
        wts[n] = pair_diff
    #print(wts)
        
    to_vacc = sorted( wts, key=lambda x: -wts[x] )[:vaccinateN]
    return to_vacc

def nominate_local_betweenness(sim, vaccinateN=100):
    wts = {}
    for n in range(sim.tnet.Nnodes):
        subg = nx.subgraph( sim.tnet.G, nx.neighbors(sim.tnet.G, n) )
        nn = len(subg.nodes)
        pair_diff = nn*(nn-1) - len(subg.edges)
        
        wts[n] = pair_diff
    
    
    chosen = []
    
    ii = 0
    while len(chosen) < vaccinateN:
        ii += 1
        if ii > 100000:
            raise Exception("Excessive recursion...")
    
        p1 = choice(range(sim.tnet.Nnodes))
        
        p2s = list(nx.neighbors(sim.tnet.G,p1))
        p2s = [x for x in p2s if x not in chosen]
        
        if not len(p2s):
            #print('no choices, skipped')
            continue
        
        p2 = sorted(p2s, key=lambda x:wts[x])[-1]
        chosen.append(p2)
     
    return chosen


def rotate_categories( fn, category ):
    def woot(sim, vaccinateN=100):
        order = fn(sim, vaccinateN=sim.tnet.Nnodes)
        
        to_vacc = []
        cv = sim.tnet.node_attr[category]
        ci = 0
        cs = sorted(set(cv.values()))
        csd = {
            c: [x for x,y in cv.items() if y==c]
            for c in cs
        }
        
        while len(to_vacc) <= min(vaccinateN, len(order)):
            myc = csd[ cs[ci] ]
            if len(myc):
                to_vacc.append( myc.pop() )
                
            ci += 1
            ci %= len(cs)
        
        return to_vacc
    
    return woot

def fill_categories( fn, category, CUTOFF=30 ):
    def woot(sim, vaccinateN=100):
        order = fn(sim, vaccinateN=sim.tnet.Nnodes)
        
        to_vacc = []
        cv = sim.tnet.node_attr[category]
        ci = 0
        cs = sorted(set(cv.values()))
        csd = {
            c: [x for x,y in cv.items() if y==c]
            for c in cs
        }
        
        cs_num = {c:0 for c in cs}
        
        while len(to_vacc) <= min(vaccinateN, len(order)):
            cs_sort = sorted( cs_num.items(), key=lambda x:(x[1]//CUTOFF, x[0]) )
            cs_top = cs_sort[0][0]

            myc = csd[ cs_top ]
            if len(myc):
                to_vacc.append( myc.pop() )
                cs_num[cs_top] += 1
            
            #print(cs_num)
        
        return to_vacc
    
    return woot



sampling_strats = [
    early_infect,
    rand,
    friend,
    friendHighDeg,
    targeted,
    contact,
    target_between_classes,
    local_betweenness,
    nominate_local_betweenness,
    none
]