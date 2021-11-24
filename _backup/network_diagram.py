import sys; sys.path.append("..")
from epi_model import *

def network_diagram(a, groups=None, output_file='temp.png'):

    plt.hist([t for t,u,v in a.edgelist], bins=50);

    [(t,u,v) for t,u,v in a.edgelist if t > 2500][:5]

    [(t,u,v) for t,u,v in a.edgelist][:5]

    [(t,u,v) for t,u,v in a.edgelist if t > 3000][:5]

    ccc = a.node_attr['class']

    classes = {
        c: [k for k,v in ccc.items() if v==c]
        for c in set(ccc.values())
    }

    css = a.node_attr['class']

    nds = sorted(a.node_attr['class'])

    class_map = a.node_attr['class']

    C_DIST = 3

    wt = Counter()
    wt += Counter(u for t,u,v in a.edgelist)
    wt += Counter(v for t,u,v in a.edgelist)

    pos = {}

    wts3 = Counter()

    for i,c in enumerate(sorted(classes)):
        gnodes = [ k for k,v in class_map.items() if v==c ]

        G = nx.Graph()

        wts = Counter([(u,v) for t,u,v in a.edgelist if (u in gnodes and v in gnodes)])
        wts2 = Counter([u for t,u,v in a.edgelist if (u in gnodes and v in gnodes)])
        wts2 += Counter([v for t,u,v in a.edgelist if (u in gnodes and v in gnodes)])

        #G.add_nodes_from([k for k in gnodes if wts2[k] > 0])
        G.add_nodes_from({k for k in gnodes if wt[k]})
        G.add_weighted_edges_from([(u,v, wts[(u,v)]) for (u,v) in wts])
        pos_class = nx.spring_layout(G, k=0.5)

        shift_v = np.array( [
            np.cos((i/len(classes)) * 2 * np.pi),
            np.sin((i/len(classes)) * 2 * np.pi),
        ] ) * C_DIST
        #shift_v = np.array([C_DIST,0]) * i

        x_mx, x_mn = max(v[0] for v in pos_class.values()), min(v[0] for v in pos_class.values())
        y_mx, y_mn = max(v[1] for v in pos_class.values()), min(v[1] for v in pos_class.values())

        print(x_mn, x_mx, y_mn, y_mx, c)
        if c == 'MP*2':
            print(pos_class)

        for k,v in pos_class.items():
            v[0] = (v[0] - x_mn) / ((x_mx-x_mn)/2) + -1
            v[1] = (v[1] - y_mn) / ((y_mx-y_mn)/2) + -1

            pos[k] = v*0.5 + shift_v

    len([x for x in a.edgelist if css[x[1]]==css[x[2]]])

    G = nx.Graph()
    wts = Counter([(u,v) for t,u,v in a.edgelist])
    G.add_weighted_edges_from([(u,v, wts[(u,v)]) for (u,v) in wts])
    #G.add_weighted_edges_from([(u,v, wts[(u,v)]) for (u,v) in wts if a.node_attr['class'][u]==a.node_attr['class'][v]])

    sum( k not in G.nodes for k,v in a.node_attr['class'].items() if v == 'MP*2' )

    CUTOFF = 10 # 15 segments * 20s = 5 minutes

    from matplotlib.patches import ConnectionStyle

    #{k for k,v in a.node_attr['class'].items() if v == 'MP*2'}.difference({z for z in G.nodes if a.node_attr['class'][z] == 'MP*2'})

    fig,ax = plt.subplots(figsize=(10,10))

    #edges,weights = zip(*[((e1,e2),w) for (e1,e2),w in nx.get_edge_attributes(G,'weight').items()])
    edges,weights = zip(*[((e1,e2),w) for (e1,e2),w in nx.get_edge_attributes(G,'weight').items() if w >= CUTOFF])
    #edges,weights = zip(*[((e1,e2),w) for (e1,e2),w in nx.get_edge_attributes(G,'weight').items() if a.node_attr['class'][e1]==a.node_attr['class'][e2]])

    weights = np.array(weights)
    weights = weights / np.max(weights)

    weights = np.log(weights)
    weights = weights + -np.min(weights)
    weights = weights / weights.max()

    #node_color='b', , 
    #arcs.set_color([(0,0,0,w/3) for w in weights])
    nx.draw_networkx_nodes(G, pos, node_size=1, node_color='black')

    for i,c in enumerate(sorted(classes)):
        shift_v = np.array( [
            np.cos((i/len(classes)) * 2 * np.pi),
            np.sin((i/len(classes)) * 2 * np.pi),
        ] ) * C_DIST

        # draw the circles
        circ = plt.Circle(shift_v, 1, color='b', alpha=0.05)
        circ2 = plt.Circle(shift_v, 1, fc='none', ec='#DDD', alpha=0.5)
        plt.text( shift_v[0], shift_v[1]+1.1, c, horizontalalignment='center', backgroundcolor='#FFF' )
        ax.add_patch(circ)
        ax.add_patch(circ2)


    arcs = nx.draw_networkx_edges(G, pos, edgelist=edges, node_size=1, edge_color=weights, width=1, edge_cmap=plt.cm.Blues)
    #plt.scatter(x,y)

    plt.savefig(output_file, dpi=300)