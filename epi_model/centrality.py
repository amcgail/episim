from .common_imports import *



def joint_deg(whatG, ns):
    al = set()
    for n in ns:
        al.update(nx.neighbors(whatG, n))
    return len(al)

def cross_deg(whatG, x):
    nbs = list(nx.neighbors(whatG,x))
    subg = nx.subgraph(whatG, nbs)
    
    s = set()
    
    for i in nbs:
        for j in nbs:
            if not(i<j):continue
                
            if (i,j) not in subg.edges:
                s.add((i,j))
                
    return s

def joint_cross_deg(whatG, xs):
    if not hasattr(whatG, 'cd'): # these are slow, so we need to iteratively compute them and save for later
        whatG.cd = {
            x: cross_deg(whatG, x)
            for x in whatG.nodes
        }
        
    return set.union(*[
        whatG.cd[x] for x in xs
    ])

def cross_degwt(whatG, x):
    nbs = list(nx.neighbors(whatG,x))
    
    if not hasattr(whatG, 'mat2'):
        mat = nx.adj_matrix(whatG)
        whatG.mat2 = mat ** 2
        
    nds = list(whatG.nodes)
    
    s = defaultdict(int)
    
    for i in nbs:
        for j in nbs:
            if not(i<j):continue
                
            if (i,j) not in whatG.edges:
                s[(i,j)] = 1/whatG.mat2[nds.index(i), nds.index(j)]
                
    return s

def ddadd(x,y):
    new = defaultdict(int)
    keys = set(x.keys()).union(y.keys())
    for k in keys:
        new[k] = x[k] + y[k]
        
    return new

def joint_cross_degwt(whatG, xs):
    if not hasattr(whatG, 'cdwts'): # these are slow, so we need to iteratively compute them and save for later
        whatG.cdwts = {
            x:cross_degwt(whatG, x)
            for x in whatG.nodes
        }
    
    full = defaultdict(int)
    for x in xs:
        full = ddadd(full, whatG.cdwts[x])
    return full






def iterated_joint(whatG, whatCentr, max_len=20):
    # calculate the sorted list
    joint_group = []
    nodes_left = list(whatG.nodes)
    
    def get_centr( x ):
        res = whatCentr(whatG, joint_group+[x])
        if type(res) in {dict, defaultdict}:
            return sum(res.values())
        elif type(res) in {list, set}:
            return len(res)
        elif type(res) in {float, int}:
            return res
        else: raise Exception("What is a centrality of type %s" % (type(res)))

    for LP in range(max_len):
        next_n = max(nodes_left, key=get_centr)
        joint_group += [next_n]
        nodes_left.remove( next_n )
        
    return joint_group