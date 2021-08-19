from .common_imports import *
from . import weightedNetwork,temporalNetwork
import os

def doit(tnet, out_name, sim=None,
    window = 15*60, # in seconds
    FPS = 20, groups=None,
    C_DIST = 3, max_frames=None,
         overwrite=False,
    # one hour every two seconds
    step = 60*60 / (20*2)
):
    
    start = tnet.time_start
    days = [x+start for x in tnet.day_breaks]



    # create the folder
    # delete the files
    fold = Path( out_name )
    fold.mkdir(exist_ok=True)
    if overwrite:
        for f in fold.glob("*"):
            if f.is_dir():
                continue
            f.unlink()
    else:
        print("Warning: not overwriting. Use overwrite=True to fix.")

        

    
    if groups is not None:
        # determine positions based on classes...
        classes = sorted(set(groups.values()))

        node_groups = {
            k: [v for v,gg in groups.items() if gg==k]
            for k in classes
        }
        
        pos = {}
        
        for i,c in enumerate(classes):
            # arrange just this class with the spring formation thing
            gnodes = node_groups[c]

            G = nx.Graph()
            if isinstance(tnet, weightedNetwork):
                G.add_weighted_edges_from([(u,v,w) for (u,v,w) in tnet.edges if (u in gnodes and v in gnodes)])
                
            elif isinstance(tnet, temporalNetwork):
                wts = Counter([(u,v) for t,u,v in tnet.edgelist if (u in gnodes and v in gnodes)])
                G.add_nodes_from(gnodes)
                G.add_weighted_edges_from([(u,v, wts[(u,v)]) for (u,v) in wts])
                
            pos_class = nx.spring_layout(G)

            shift_v = np.array( [
                np.cos((i/len(classes)) * 2 * np.pi),
                np.sin((i/len(classes)) * 2 * np.pi),
            ] ) * C_DIST

            for k,v in pos_class.items():
                pos[k] = v + shift_v
                
    else: # no groupings
        G = nx.Graph()
        if isinstance(tnet, weightedNetwork):
            G.add_weighted_edges_from([(u,v,w) for (u,v,w) in tnet.edges])
            
        elif isinstance(tnet, temporalNetwork):
            wts = Counter([(u,v) for t,u,v in tnet.edgelist])
            G.add_weighted_edges_from([(u,v, wts[(u,v)]) for (u,v) in wts])
        pos = nx.spring_layout(G)
        
    # determine X and Y LIM
    G = nx.Graph()
    if isinstance(tnet, weightedNetwork):
        G.add_weighted_edges_from([(u,v,w) for (u,v,w) in tnet.edges])
    elif isinstance(tnet, temporalNetwork):
        wts = Counter([(u,v) for t,u,v in tnet.edgelist])
        G.add_weighted_edges_from([(u,v, wts[(u,v)]) for (u,v) in wts])
    arcs = nx.draw_networkx_edges(G, pos, node_size=1, width=1)
    nx.draw_networkx_nodes(G, pos, node_size=1)
    X_LIM = plt.xlim()
    Y_LIM = plt.ylim()
    
    ALL_NODES = sorted(G.nodes)

    def generate(FRAME_NUMBER):

        old_outf = fold.joinpath('%03.0f.png' % FRAME_NUMBER)
        new_outf = fold.joinpath('%05.0f.png' % FRAME_NUMBER)
        
        if new_outf.exists(): # only will happen when overwrite=False
            return False
        
        outf = new_outf

        if old_outf.exists():
            #print("Renaming %s to %s" % (old_outf,new_outf))
            old_outf.rename(new_outf)
            return False

        day_begin = start + FRAME_NUMBER*step
        day_end = start + FRAME_NUMBER*step + window
        
        # translate these
        day_begini = (day_begin - start) // 20
        day_endi = (day_end - start) // 20

        if not len([x for x in tnet.times if day_begini <= x <= day_endi]): # eldf.loc[ (day_begini <= tnet.eldf[0])&(tnet.eldf[0] <= day_endi) ].shape[0] == 0
            #print("skipping %s" % i)
            return False

        G = nx.Graph()
        G.add_nodes_from( ALL_NODES )
        
        if isinstance(tnet, weightedNetwork):
            G.add_weighted_edges_from([(u,v,w) for (u,v,w) in tnet.edges])
            
        elif isinstance(tnet, temporalNetwork):
            wts = Counter([(u,v) for t,u,v in tnet.edgelist if day_begini <= t <= day_endi])
            G.add_weighted_edges_from([(u,v, wts[(u,v)]) for (u,v) in wts])

        fig, ax = plt.subplots(figsize=(15,15))
        #print(day_begin, day_end)


        edges,weights = zip(*nx.get_edge_attributes(G,'weight').items())

        weights = np.array(weights)
        weights = weights / np.max(weights)

        weights = np.log(weights)
        weights = weights + -np.min(weights)
        weights = weights / weights.max()

        arcs = nx.draw_networkx_edges(G, pos, edgelist=edges, node_size=1, edge_color=weights, width=1, edge_cmap=plt.cm.Blues)
        #node_color='b', , 
        arcs.set_color([(0,0,0,w) for w in weights])

        nx.draw_networkx_nodes(G, pos, node_size=1, node_color='white')

        # draw the circles
        if groups is not None:
            for ii,c in enumerate(classes):
                shift_v = np.array( [
                    np.cos((ii/len(classes)) * 2 * np.pi),
                    np.sin((ii/len(classes)) * 2 * np.pi),
                ] )*C_DIST

                circ = plt.Circle(shift_v, 1, color='b', alpha=0.1)
                plt.text( shift_v[0], shift_v[1]+1.1, c, horizontalalignment='center' )
                ax.add_patch(circ)

            
        if sim is not None:
            nx.draw_networkx_nodes(G, pos, node_size=100, nodelist=[n for n in G.nodes if sim.stateT(n,FRAME_NUMBER*step/20)=='inf'], node_color='yellow', edgecolors='black')
            nx.draw_networkx_nodes(G, pos, node_size=60, nodelist=[n for n in G.nodes if sim.stateT(n,FRAME_NUMBER*step/20)=='rec'], node_color='brown', edgecolors='black')
            

        sdt = dt.datetime.fromtimestamp(day_begin) + dt.timedelta(hours=6)
        edt = dt.datetime.fromtimestamp(day_end) + dt.timedelta(hours=6)
        
        tx = "%s\n%s - %s" % (
            sdt.strftime("%a, %D"),
            sdt.strftime("%I:%M %p"),
            edt.strftime("%I:%M %p")
        )

        if groups is not None: # put the text in the middle of the groups
            txRect = plt.Rectangle([-1.25,-0.25], 2.5, 1, color='white', ec='black', alpha=1)
            ax.add_patch(txRect)

            plt.text( 0,0, tx, horizontalalignment='center', fontsize=20 )
        else:
            plt.title(tx)


        plt.xlim(*X_LIM)
        plt.ylim(*Y_LIM)


        plt.savefig( outf )
        plt.close()
        return True

    num_steps = int((tnet.Tmax*20 - window) // step)

    #for i in range(24):
    frames_rendered = 0
    for i in range(num_steps):
        
        res = generate(i)
        if res:
            frames_rendered += 1
        
            if frames_rendered % 100 == 0:
                print("%s/%s [%s] frames" % (frames_rendered, num_steps, max_frames))
                pass
        
        if max_frames is not None and frames_rendered >= max_frames:
            break
            
    print("Finished rendering %s frames" % frames_rendered)

    # they aint in a row!
    for i,f in enumerate(sorted(fold.glob("*.png"))):
        #print(f.name, "%05d.png"%i)
        #continue
        f.rename(fold.joinpath("%05d.png"%i))

    # make movie
    outf = fold.joinpath("full_video.mp4")
    if outf.exists():
        outf.unlink()

    f_pattern = str(fold.joinpath("%05d.png"))
    #!ffmpeg -r {FPS} -f image2 -s 1080x1080 -i {f_pattern} -vcodec libx264 -t 83.45 -crf 25  -pix_fmt yuv420p {outf}
    os.system("ffmpeg -framerate {FPS} -i {f_pattern} -c:v libx264 -r 30 {outf}".format(**locals()))