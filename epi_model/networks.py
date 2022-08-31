from typing import Dict
from .common_imports import *

class weightedNetwork:
    
    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.edges = edges # (F, T, w)
        self.edgelist = edges # some compatibility with tnets
        
        self.ego_edges = defaultdict(dict)
        for f,t,w in self.edges:
            self.ego_edges[f][t] = w
            self.ego_edges[t][f] = w

        self.Nnodes = len(self.nodes)
        self.days = [dt.datetime.now().date()]

        self.mindayT = [0]
        self.maxdayT = [int( 0.25 * 3600*24/20 )]
        self.day_breaks = [0, 3600*24/20]

        G = nx.Graph()
        G.add_weighted_edges_from([
            [f,t,w]
            for f,t,w in self.edges
        ])
        self.G = G
            

class unweightedNetwork(weightedNetwork):
    
    def __init__(self, nodes, edges, weight):
        return super().__init__(
            nodes,
            [[f,t,weight] for f,t in edges]
        )

    @classmethod
    def from_csv(cls, fn, weight):

        from csv import DictReader
        
        with open(fn, 'r', encoding='utf8') as inf:
            rs = list(DictReader(inf))

        ks = sorted(rs[0].keys())
        assert len(ks) == 2 # must have 2 columns

        edges = [ [int(r[ks[0]]),int(r[ks[1]])] for r in rs ]
        nodes = sorted(set( [r[0] for r in edges] ).union(set( [r[1] for r in edges] )))

        return unweightedNetwork( nodes, edges, weight )


class temporalNetwork:

    def range(self, day_start, day_end):

        start = self.mindayT[day_start]
        end = self.maxdayT[day_end-1]
        
        new_timeStart = self.time_start + (day_start)*3600*24
        new_timeStart_diff = (day_start)*3600*24/20
        
        new_edgelist = [
            (t - new_timeStart_diff, a, b) for (t, a, b) in self.edgelist if start <= t <= end
        ]
        new_edgelist = np.array(new_edgelist)

        return temporalNetwork(
            edgelist = new_edgelist,
            time_start = new_timeStart,
            node_attr = self.node_attr
        )


    def _trim_days_fail(self, days):
        def transfer_time(when):
            day_offset = 0
            for i in range(len(self.day_breaks)-1):
                if self.day_breaks[i+1] >= when >= self.day_breaks[i]:
                    # add the thing
                    if i in days:
                        return when - day_offset
                    else:
                        return None
                else:
                    if i not in days:
                        day_offset += self.day_breaks[i+1] - self.day_breaks[i]
            
            raise Exception("when is not in any day... confusing%s"%when, when)

        new_edgelist = [
            (transfer_time(t), a, b) for (t, a, b) in self.edgelist
        ]
        new_edgelist = list(filter(lambda x: x[0] is not None, new_edgelist))
        new_edgelist = np.array(new_edgelist)

        new_day_breaks = [transfer_time(x) for x in self.day_breaks[:-1]]
        new_day_breaks = list(filter(lambda x: x is not None, new_day_breaks))
        
        return temporalNetwork(
            edgelist = new_edgelist,
            time_start = self.time_start,
            times = self.times,
            node_attr = self.node_attr,
            day_breaks = new_day_breaks,
            node_ids = self.node_ids
        )
    
    def to_weighted(self, normalize=False):
        
        nodes = range(len(self.G.nodes))
        wts = Counter([(u,v) for t,u,v in self.edgelist])
        
        TOTAL_TIME = (
            len(set(self.edgelist[:,0])) 
            if normalize else 1
            # amount of time for which we have data
        )
        
        edges = [
            (f,t,c / TOTAL_TIME)
            for (f,t),c in wts.items()
        ]
        
        w = weightedNetwork(nodes, edges)
        w.node_attr = self.node_attr
        w.time_start = self.time_start
        w.day_breaks = self.day_breaks
        w.node_ids = self.node_ids
        w.Tmax = self.Tmax
        w.Tmin = self.Tmin
        w.Nnodes = self.Nnodes
        w.times = self.times
        w.G = self.G
        
        return w

    def __init__(self,edgelist,time_start,node_attr,):
        self.edgelist = np.copy(edgelist)
        self._edgelist_arg = edgelist
        self.time_start = time_start
        self.node_attr = node_attr

        self.execute_preprocessing()
        self.weighted = self.to_weighted()

        self.nodes = list(self.G.nodes)
    
    @classmethod
    def load(cls, dataset):
        
        if dataset == 'high school2':

            import pandas as pd
            edgelist = pd.read_csv( str(Path(DATA_DIR, "high school 2/full.txt")), sep="\t", header=None)
            edgelist = np.array(edgelist, dtype=int)
            edgelist = edgelist[:,:3]

            time_start = dt.datetime(2010,1,14, 6,0,0).timestamp()

            
            roles = pd.read_csv(str(Path(DATA_DIR, "high school 2/roles.txt")), sep="\t", header=None)
            role_map = {
                x[0]: x[1]
                for i,x in roles.iterrows()
            }
            
            node_attr = {
                'role': role_map
            }


        elif dataset == 'high school':

            import pandas as pd
            edgelist = pd.read_csv(str(Path(DATA_DIR, "high school/High-School_data_2013.csv")), sep=" ", header=None)

            node_ids = sorted(set(edgelist[1]).union(set(edgelist[2])))

            classes = {
                c: set(edgelist[ edgelist[3] == c ][1]).union(edgelist[ edgelist[4] == c ][2])
                for c in set(edgelist[3])
            }

            #node_ids.index(p)
            class_map = {
                p: c
                for c in classes
                for p in classes[c]
            }
            
            edgelist = np.array(edgelist.loc[:,:2], dtype=int)
            
            # normalize the times
            time_start = dt.datetime(2013,12,2, 7,0,0).timestamp()
            
            node_attr = {
                'class': class_map
            }
        
        elif dataset == 'workplace':
            
            import pandas as pd
            edgelist = pd.read_csv(str(Path(DATA_DIR, "contacts in the workplace/tij_InVS.DAT")), sep=" ", header=None)
            edgelist = np.array(edgelist)

            time_start = dt.datetime(2013,6,24).timestamp()

            depts = pd.read_csv(str(Path(DATA_DIR, "contacts in the workplace/metadata_InVS13.TXT")), sep="\t", header=None)
            
            node_attr = {
                'department': {
                    x[0]: x[1]
                    for i,x in depts.iterrows()
                }
            }

        elif dataset == 'workplace2':
            
            import pandas as pd
            edgelist = pd.read_csv(str(Path(DATA_DIR, "contacts in the workplace 2/tij_InVS15.DAT")), sep=" ", header=None)
            edgelist = np.array(edgelist)

            time_start = dt.datetime(2014,6,24, 9,0,0).timestamp() # I ACTUALLY CANT FIND THIS ONE...

            depts = pd.read_csv(str(Path(DATA_DIR, "contacts in the workplace 2/metadata_InVS15.TXT")), sep="\t", header=None)

            node_ids = set(edgelist[:,1]).union(set(edgelist[:,2]))
            
            node_attr = {
                'department': {
                    x[0]: x[1]
                    for i,x in depts.iterrows()
                    if x[0] in node_ids
                }
            }


        elif dataset == 'primary school':
            import pandas as pd
            edgelist = pd.read_csv(str(Path(DATA_DIR, "primary school/primaryschool.csv")), sep="\t", header=None)
            edgelist = np.array(edgelist)

            time_start = dt.datetime(2009,10,1, 7,0,0).timestamp()

            depts = pd.read_csv(str(Path(DATA_DIR, "primary school/metadata_primaryschool.txt")), sep="\t", header=None)
            
            node_attr = {
                'class': {
                    x[0]: x[1]
                    for i,x in depts.iterrows()
                },
                'gender': {
                    x[0]: x[2]
                    for i,x in depts.iterrows()
                },
            }
        
        elif dataset == 'gillespie_sample':
            import pandas as pd
            edgelist = pd.read_csv(str(Path(DATA_DIR, "temporal_gillespie_sample.txt")), sep="\t", header=None)
            edgelist = np.array(edgelist)

            time_start = dt.datetime(2020,1,1).timestamp()
            
            node_attr = {}
        
        
        else:
            raise Exception('dataset not found...')

        return temporalNetwork(
            edgelist = edgelist,
            time_start = time_start,
            node_attr = node_attr,
        )
        
    def execute_preprocessing(self):
        
        self.node_ids = sorted(set(self.edgelist[:,1]).union(self.edgelist[:,2]))

        #node_transfer = lambda x: self.node_ids.index(x)

        #self.edgelist[:,1] = np.array(list(map(node_transfer, self.edgelist[:,1])))
        #self.edgelist[:,2] = np.array(list(map(node_transfer, self.edgelist[:,2])))

        # normalize the times
        orig_times = sorted(set(self.edgelist[:,0]))
        first_diff = orig_times[1] - orig_times[0]
        
        self.edgelist[:,0] = self.edgelist[:,0] - self.edgelist[0,0]
        
        if first_diff == 20:
            self.edgelist[:,0] = self.edgelist[:,0] / 20
        elif first_diff != 1:
            raise Exception("time difference must be 20 or 1")

        self.times = sorted(set(self.edgelist[:,0]))


        self.days = sorted(set(
            [
                dt.datetime.fromtimestamp(x).date() 
                for x in np.array(self.times)*20 + self.time_start
            ]
        ))
        
        self.Tmin = np.min(self.edgelist[:,0])
        self.Tmax = np.max(self.edgelist[:,0])
        
        self.day_breaks = [
            (dt.datetime.fromordinal(x.toordinal()).timestamp() - self.time_start) / 20
            for x in self.days
        ]
        self.day_breaks += [self.day_breaks[-1] + 3600*24 / 20]

        self.day_breaks = [int(x) for x in self.day_breaks]
        
        # this is necessary to check whether the network is directed intiitally
        edgelist_counter = Counter([
            tuple([t]+sorted([u,v]))
            for t,u,v in self.edgelist
        ])

        if any(x%2!=0 for x in edgelist_counter.values()):
            # ensure the network is undirected!
            edgelist = np.array(self.edgelist, dtype=int)
            edgelist_directed = np.zeros(( len(edgelist)*2, 3 ))
            for idx in range( len(edgelist) ):
                t, u, v = edgelist[idx]
                edgelist_directed[2*idx, :] = t, u, v
                edgelist_directed[2*idx + 1, :] = t, v, u
            edgelist = edgelist_directed.astype(int)
            self.edgelist = edgelist
        
        
        
        
        self.t_edges = defaultdict(list)
        for t, a, b in self.edgelist:
            #if a>b:
            self.t_edges[t].append((a,b))
        
        self.edgelist = self.edgelist.astype(int)
                
        # make time-aggregated network
        G = nx.Graph()
        wts = Counter([(u,v) for t,u,v in self.edgelist])
        G.add_weighted_edges_from([(u,v, wts[(u,v)]) for (u,v) in wts])
        
        self.Nedges = len(G.edges)
        self.Nnodes = len(G.nodes)
        self.G = G
        
        self.eldf = pd.DataFrame.from_records( self.edgelist )


        
        self.mindayT = {
            day: min( [t for (t,_,_) in self.edgelist 
                if self.day_breaks[day] <= t <= self.day_breaks[day+1]] )
            for day in range( len(self.days) )
        }
        
        self.maxdayT = {
            day: max( [t for (t,_,_) in self.edgelist 
                if self.day_breaks[day] <= t <= self.day_breaks[day+1]] )
            for day in range( len(self.days) )
        }