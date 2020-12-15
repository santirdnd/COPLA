#!/usr/bin/env python3

import os
import sys
import pickle
import pathlib
import argparse
import numpy as np
import pandas as pd
import graph_tool.all as gt
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    description='Cluster a set of plasmids into PTUs')
parser.add_argument('outdir', type=str,
                    help='output directory')
parser.add_argument('-a', '--adjacency', type=str,
                    default='databases/Copla_RS84/RS84f_adjacency_matrix_ANIp50.tsv.gz',
                    help='adjacency matrix')
parser.add_argument('-m', '--metadata', type=str,
                    default='databases/Copla_RS84/RS84f_plasmid_metadata.tsv',
                    help='plasmid metadata')
parser.add_argument('-j', '--hierarchical', action='store_true',
                    help='Use hierarchical SBM variant')
parser.add_argument('-d', '--deg_corr', action='store_true',
                    help='Use degree corrected SBM variant')
parser.add_argument('-w', '--weight_model', type=str,
                    choices=['None', 'Exponential', 'Normal', 'LogNormal'], default='None',
                    help='Edge weight distribution model')
parser.add_argument('-n', '--nToss', type=int,
                    default=100,
                    help='number of initial SBM algorithm initializations')
parser.add_argument('-l', '--maxLevels', type=int,
                    default=15,
                    help='maximum number of levels for HSBM')
parser.add_argument('-t', '--transient', type=int,
                    default=2000,
                    help='number of iterations to skip transitory phase')
parser.add_argument('-i', '--iters', type=int,
                    default=10000,
                    help='number of iterations to store posterior statistics')
parser.add_argument('-f', '--filter', action='store_true',
                    help='Filter out connected components with less than 4 members')
parser.add_argument('-r', '--refLabel', type=str,
                    choices=['PTU_Curated', 'PTU_sHSBM'], default='PTU_sHSBM',
                    help='column name for using as reference of PTU labeling')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)

# Load network data
nodes = pd.read_csv(args.metadata, sep="\t")
ani = pd.read_csv(args.adjacency, compression='gzip', sep="\t", header=None, names=nodes.AccessionVersion)
nVertices = len(list(nodes.AccessionVersion))
ani_np = ani.to_numpy() # Complete ANI matrix (symetric, with autoloops)
ani_np = np.triu(ani.to_numpy(), 0) # Upper triangular ANI matrix (with autoloops)
ani_np_idx = ani_np.nonzero()
edge_list = np.transpose(ani_np_idx)

# Create graph object
g = gt.Graph(directed=False)
g.add_vertex(nVertices)
g.add_edge_list(edge_list)

# Use AccessionVersion to name the vertices
v_AccVer = g.new_vertex_property('string', nodes.AccessionVersion.to_list())
g.vertex_properties['AccessionVersion'] = v_AccVer

# Add other relevant properties to vertices
v_MOB = g.new_vertex_property('string', nodes.MOB_60.to_list())
g.vertex_properties['MOB'] = v_MOB
v_MPF = g.new_vertex_property('string', nodes.MPF.to_list())
g.vertex_properties['MPF'] = v_MPF
v_PFinder = g.new_vertex_property('string', nodes.PFinder_80.to_list())
g.vertex_properties['PFinder'] = v_PFinder
v_AMR = g.new_vertex_property('string', nodes.CARD_80.to_list())
g.vertex_properties['AMR'] = v_AMR
v_Size = g.new_vertex_property('int', nodes.Size.to_list())
g.vertex_properties['Size'] = v_Size
v_Topology = g.new_vertex_property('string', nodes.Topology.to_list())
g.vertex_properties['Topology'] = v_Topology
v_TaxKingdom = g.new_vertex_property('string', nodes.TaxSuperkingdom.to_list())
g.vertex_properties['TaxKingdom'] = v_TaxKingdom
v_TaxPhylum = g.new_vertex_property('string', nodes.TaxPhylum.to_list())
g.vertex_properties['TaxPhylum'] = v_TaxPhylum
v_TaxClass = g.new_vertex_property('string', nodes.TaxClass.to_list())
g.vertex_properties['TaxClass'] = v_TaxClass
v_TaxOrder = g.new_vertex_property('string', nodes.TaxOrder.to_list())
g.vertex_properties['TaxOrder'] = v_TaxOrder
v_TaxFamily = g.new_vertex_property('string', nodes.TaxFamily.to_list())
g.vertex_properties['TaxFamily'] = v_TaxFamily
v_TaxGenus = g.new_vertex_property('string', nodes.TaxGenus.to_list())
g.vertex_properties['TaxGenus'] = v_TaxGenus
v_TaxSpecies = g.new_vertex_property('string', nodes.TaxSpecies.to_list())
g.vertex_properties['TaxSpecies'] = v_TaxSpecies
if args.refLabel == 'PTU_Curated':
    v_PtuRef = g.new_vertex_property('string', nodes.PTU_Curated.to_list())
else:
    v_PtuRef = g.new_vertex_property('string', nodes.PTU_sHSBM.to_list())
g.vertex_properties['PtuRef'] = v_PtuRef

# Add edge weights
e_ANI = g.new_edge_property('double')
e_ANI.a = ani_np[ani_np_idx] / 100
g.edge_properties['ANI'] = e_ANI

# Transform weights
if (args.weight_model != 'None'):
    y = g.ep.ANI.copy()
    if (args.weight_model == 'LogNormal'):
        y.a = 0 - np.log(y.a)

# Find out connected components with less than 4 members
comp, hist = gt.label_components(g)
g.vertex_properties['CComp'] = comp
v_CC4 = g.new_vertex_property('bool', np.isin(comp.a, np.where(hist >= 4)))
g.vertex_properties['CC4_filter'] = v_CC4

# Filter out connected components with less than 4 members
if args.filter:
    g.set_vertex_filter(g.vp.CC4_filter)
    g.purge_vertices()
    g.set_vertex_filter(None)

# Save graph (graph-tool format)
fname = 'graph'
g.save(os.path.join(args.outdir,fname+'.gt.gz'))
g.save(os.path.join(args.outdir,fname+'.xml.gz'))

def write_classes(filename, graph):
    f = open(filename, 'w')
    header = "AccessionVersion\tCComp\tBlock\tBlockCC\tBlockCC4"
    if 'sHSBM' in g.vp.keys():
        header += "\tsHSBM\tPTU\tHRange"
    f.write(header+"\n")
    for v in graph.vertices():
        AccessionVersion = graph.vp.AccessionVersion[v]
        CComp = str(graph.vp.CComp[v])
        Block = str(graph.vp.Block[v])
        BlockCC = graph.vp.BlockCC[v]
        BlockCC4 = graph.vp.BlockCC4[v]
        if 'sHSBM' in g.vp.keys():
            sHSBM = graph.vp.sHSBM[v]
            Ptu = graph.vp.Ptu[v]
            HRange = graph.vp.HRange[v]
            f.write("\t".join((AccessionVersion, CComp, Block, BlockCC, BlockCC4, sHSBM, Ptu, HRange)) + "\n")
        else:
            f.write("\t".join((AccessionVersion, CComp, Block, BlockCC, BlockCC4)) + "\n")
    f.close()

def block_annotation(graph, state):
    if args.hierarchical:
        levels = state.get_levels()

        # Find the informative hierarchical levels (i.e. the non-redundant levels, those with non-equivalent block assignment)
        def check_level_redundancy(l):
            x = state.project_partition(l, 0).a
            y = state.project_partition(l+1, 0).a
            return gt.partition_overlap(x, y, norm=True) == 1
        L = len(levels) - 1
        redundant_levels = [False]+list(map(check_level_redundancy, reversed(range(L))))
        nr_levels = [L-i for i, x in enumerate(redundant_levels) if not x]

        b = levels[0].get_blocks()
        bcc = graph.new_vertex_property('string')
        bcc4 = graph.new_vertex_property('string')
        for i in np.unique(b.a):
            b_filter = (b.a == i)
            u = gt.GraphView(graph, vfilt=b_filter)
            tmp = []
            r = u.get_vertices()[0]
            for l in range(len(levels)):
                r = levels[l].get_blocks()[r]
                if l in nr_levels:
                    tmp.append(str(r))
            tmp.reverse()
            comp, hist = gt.label_components(u)
            for v in u.vertices():
                tag = '_'.join(tmp + [str(comp[v])])
                bcc[v] = tag
                bcc4[v] = tag if (hist[comp[int(v)]] >= 4) else '-'
    else:
        b = state.get_blocks()
        bcc = graph.new_vertex_property('string')
        bcc4 = graph.new_vertex_property('string')
        for i in np.unique(b.a):
            b_filter = (b.a == i)
            u = gt.GraphView(graph, vfilt=b_filter)
            comp, hist = gt.label_components(u)
            for v in u.vertices():
                tag = '_'.join([str(i), str(comp[v])])
                bcc[v] = tag
                bcc4[v] = tag if (hist[comp[int(v)]] >= 4) else '-'

    return((b, bcc, bcc4))

def ptu_annotation(graph):
    complex = []
    bcc4 = list(graph.vp.BlockCC4)
    bcc4_nr = np.unique(bcc4)

    # The algorithm will join clusters belonging to the same hierarchical branch. However, to save
    # time, we only check up to certain hierarchical level and not transverse the 5 upper levels
    #common_levels = 7
    if bcc4_nr[0] != '-':
        common_levels = len(bcc4_nr[0].split('_')) - 5
    else:
        common_levels = len(bcc4_nr[1].split('_')) - 5

    for n, i in enumerate(bcc4_nr[:-1]):
        if i == '-':
            continue
        i_filter = (np.array(bcc4) == i)
        u = gt.GraphView(graph, vfilt=i_filter)
        u_vertices = u.num_vertices()
        u_edges = u.num_edges()
        u_size = np.median(list(u.vp.Size))
        for j in bcc4_nr[n+1:]:
            if (j == '-'):
                continue
            if ('_'.join(j.split('_')[:common_levels]) != '_'.join(i.split('_')[:common_levels])):
                continue
            j_filter = (np.array(bcc4) == j)
            w = gt.GraphView(graph, vfilt=j_filter)
            w_vertices = w.num_vertices()
            w_edges = w.num_edges()
            w_size = np.median(list(w.vp.Size))
            if (u_size > w_size):
                s_comp = True if (abs(u_size - w_size) < (u_size * 0.5)) else False
            else:
                s_comp = True if (abs(u_size - w_size) < (w_size * 0.5)) else False
            k_filter = np.logical_or(i_filter, j_filter)
            z = gt.GraphView(graph, vfilt=k_filter)
            z_vertices = z.num_vertices()
            z_edges = z.num_edges()
            # z_comp_p is true if the number of intercluster edges is >50% of posible edges between both clusters taking into account their respective densities
            z_comp_p = True if (z_edges - (u_edges + w_edges) > (u_vertices * w_vertices*((u_edges-u_vertices)/(u_vertices*(u_vertices-1)/2))*((w_edges-w_vertices)/(w_vertices*(w_vertices-1)/2))*0.5)) else False
            if z_comp_p and s_comp:
                # Annotate both clusters as belonging to the same PTU
                # Find the first cluster already included in one PTU and insert the othe one
                for m, c in enumerate(complex):
                    if (i in c) or (j in c):
                        if i not in c:
                            complex[m].append(i)
                        if j not in c:
                            complex[m].append(j)
                        break
                else:
                    complex.append([i, j])
                # Take into account if each cluster were already included in different PTUs
                for m, c in enumerate(complex[:-1]):
                    for l, d in enumerate(complex[m+1:]):
                        for e in c:
                            if e in d:
                                complex[m] = complex[m] + list(set(d) - set(c))
                                complex.pop(l+m+1)
                                break
    cmplx = {}
    for i in bcc4_nr:
        for m, c in enumerate(complex):
            if i in c:
                cmplx[i] = c[0]
                break
        else:
            cmplx[i] = i
    v_sHSBM = graph.new_vertex_property('string')
    for v in graph.vertices():
        v_sHSBM[v] = cmplx[graph.vp.BlockCC4[v]]

    #for i in np.unique(list(v_sHSBM)):
    #    if i == '-':
    #        continue
    #    i_filter = (np.array(list(v_sHSBM)) == i)
    #    u = gt.GraphView(graph, vfilt=i_filter)
    #    d_intra = (u.num_edges() - u.num_vertices()) / ((u.num_vertices() * (u.num_vertices() - 1))/2)
    #    if (d_intra < 0.25):
    #        print(i, str(d_intra))
    #    k_filter = np.logical_not(i_filter)
    #    z = gt.GraphView(graph, vfilt=k_filter)
    #    d_inter = (g.num_edges()-(u.num_edges()+z.num_edges())) / (u.num_vertices()*z.num_vertices())
    #    if (d_inter > 0) and (d_intra/d_inter < 500):
    #        print(i, str(d_intra), str(d_inter), str(d_intra/d_inter))

    # Rename as many PTUs as possible aligning the partition labels
    v_PtuRef = g.vp.PtuRef.copy()
    ptuRef_index={}
    for i, p in enumerate(np.unique(list(v_PtuRef))):
        ptuRef_index[p] = i
    sHSBM_index={}
    for i, p in enumerate(np.unique(list(v_sHSBM))):
        sHSBM_index[p] = i
    ptuRef_num = []
    sHSBM_num = []
    for v in g.vertices():
        ptuRef_num.append(ptuRef_index[v_PtuRef[v]])
        sHSBM_num.append(sHSBM_index[v_sHSBM[v]])
    sHSBM_aligned = gt.align_partition_labels(sHSBM_num, ptuRef_num)
    v_Ptu = graph.new_vertex_property('string')
    goodPtuLabels = {}
    for v in g.vertices():
        if v_sHSBM[v] == '-':
            v_Ptu[v] = '-'
        elif sHSBM_aligned[int(v)] == ptuRef_num[int(v)]:
            v_Ptu[v] = v_PtuRef[v]
            if v_sHSBM[v] not in goodPtuLabels:
                goodPtuLabels[v_sHSBM[v]] = v_PtuRef[v]
        else:
            v_Ptu[v] = '?' + v_sHSBM[v]
    for v in g.vertices():
        if v_Ptu[v].startswith('?'):
            if v_sHSBM[v] in goodPtuLabels:
                v_Ptu[int(v)] = goodPtuLabels[v_sHSBM[v]]

    return(v_sHSBM, v_Ptu)

def hrange_annotation(graph):
    ptus = list(graph.vp.Ptu)

    ptu_hrange = {}
    for n, i in enumerate(np.unique(ptus)):
        if i == '-':
            ptu_hrange[i] = '-'
            continue

        i_filter = (np.array(ptus) == i)
        u = gt.GraphView(graph, vfilt=i_filter)
        for l in ['TaxKingdom', 'TaxPhylum', 'TaxClass', 'TaxOrder', 'TaxFamily', 'TaxGenus', 'TaxSpecies']:
            taxa = set()
            for v in u.vertices():
                if u.vertex_properties[l][v] != '-':
                    taxa.add(u.vertex_properties[l][v])
            if len(taxa) > 1:
                if (l == 'TaxKingdom') or (l == 'TaxPhylum') or (l == 'TaxClass'):
                    ptu_hrange[i] = 'VI'
                    break
                elif l == 'TaxOrder':
                    ptu_hrange[i] = 'V'
                    break
                elif l == 'TaxFamily':
                    ptu_hrange[i] = 'IV'
                    break
                elif l == 'TaxGenus':
                    ptu_hrange[i] = 'III'
                    break
                elif l == 'TaxSpecies':
                    ptu_hrange[i] = 'II'
                    break
        else:
            ptu_hrange[i] = 'I'

    v_HRange = graph.new_vertex_property('string')
    for v in graph.vertices():
        v_HRange[v] = ptu_hrange[graph.vp.Ptu[v]]

    return v_HRange

state_list, entropy_list = [], []
for k in range(args.nToss):
    if args.hierarchical:
        # Nested stochastic block model (hierarchical SBM)
        if (args.weight_model == 'None'):
            state = gt.minimize_nested_blockmodel_dl(g, deg_corr=args.deg_corr)
        elif (args.weight_model == 'Exponential'):
            state = gt.minimize_nested_blockmodel_dl(g, deg_corr=args.deg_corr, state_args=dict(recs=[y], rec_types=['real-exponential']))
        else:
            state = gt.minimize_nested_blockmodel_dl(g, deg_corr=args.deg_corr, state_args=dict(recs=[y], rec_types=['real-normal']))
        state_0 = state.get_levels()[0]
        nClass = len(np.unique(state_0.get_blocks().a))
    else:
        # Flat stochastic block model (SBM)
        if (args.weight_model == 'None'):
            state = gt.minimize_blockmodel_dl(g, deg_corr=args.deg_corr)
        elif (args.weight_model == 'Exponential'):
            state = gt.minimize_blockmodel_dl(g, deg_corr=args.deg_corr, state_args=dict(recs=[y], rec_types=['real-exponential']))
        else:
            state = gt.minimize_blockmodel_dl(g, deg_corr=args.deg_corr, state_args=dict(recs=[y], rec_types=['real-normal']))
        nClass = len(np.unique(state.get_blocks().a))
    entropy = state.entropy()

    # Update state
    state_list.append(state)
    entropy_list.append(entropy)
    print("Toss %d of %d: %d classes, entropy %f" % (k, args.nToss, nClass, entropy))

    # Save graph
    (v_Block, v_BlockCC, v_BlockCC4) = block_annotation(g, state)
    g.vertex_properties['Block'] = v_Block
    g.vertex_properties['BlockCC'] = v_BlockCC
    g.vertex_properties['BlockCC4'] = v_BlockCC4
    fname = "SBM_i%d_%d_%f" % (k, nClass, entropy)
    write_classes(os.path.join(args.outdir,fname+'.tsv'), g)
    pickle.dump([g, state], open(os.path.join(args.outdir,fname+'.pickle'), 'wb'), -1)
k = np.argmin(entropy_list)
state, entropy = state_list[k], entropy_list[k]
if args.hierarchical:
    bs = state.get_bs()
    bs += [np.zeros(1)] * (args.maxLevels - len(bs))
    state = state.copy(bs=bs, sampling=True)
    state_0 = state.get_levels()[0]
    nClass = len(np.unique(state_0.get_blocks().a))
else:
    nClass = len(np.unique(state.get_blocks().a))
print("Selected toss %d: %d classes, entropy %f" % (k, nClass, entropy))

# Avoid the transient state
gt.mcmc_equilibrate(state, wait=args.transient, nbreaks=2, multiflip=True, mcmc_args=dict(niter=10), verbose=False)
entropy = state.entropy()
if args.hierarchical:
    state_0 = state.get_levels()[0]
    nClass = len(np.unique(state_0.get_blocks().a))
else:
    nClass = len(np.unique(state.get_blocks().a))
print("%d classes, entropy %f" % (nClass, entropy))

# Save graph
(v_Block, v_BlockCC, v_BlockCC4) = block_annotation(g, state)
g.vertex_properties['Block'] = v_Block
g.vertex_properties['BlockCC'] = v_BlockCC
g.vertex_properties['BlockCC4'] = v_BlockCC4
fname = "SBM_transient_%d_%f" % (nClass, entropy)
write_classes(os.path.join(args.outdir,fname+'.tsv'), g)
pickle.dump([g, state], open(os.path.join(args.outdir,fname+'.pickle'), 'wb'), -1)

# Callback to collect the vertex marginal probabilities
dls = [] # Description length history
if args.hierarchical:
    pv = [None] * len(state.get_levels()) # Vertex marginals
else:
    pv = None # Vertex marginals
pe = None # Edge marginals
def collect_marginals(s):
    global pv, pe
    if args.hierarchical:
        levels = s.get_levels()
        pv = [sl.collect_vertex_marginals(pv[l], b=gt.perfect_prop_hash([sl.b])[0]) for l, sl in enumerate(levels)]
        pe = levels[0].collect_edge_marginals(pe)
    else:
        b = gt.perfect_prop_hash([s.b])[0]
        pv = s.collect_vertex_marginals(pv, b=b)
        pe = s.collect_edge_marginals(pe)
    dls.append(s.entropy())

# Apply MCMC
gt.mcmc_equilibrate(state, force_niter=args.iters, mcmc_args=dict(niter=10), callback=collect_marginals)
entropy = state.entropy()
if args.hierarchical:
    S_mf = [gt.mf_entropy(sl.g, pv[l]) for l, sl in enumerate(state.get_levels())]
    S_bethe = gt.bethe_entropy(g, pe)[0]
    L = -np.mean(dls)
    state_0 = state.get_levels()[0]
    nClass = len(np.unique(state_0.get_blocks().a))
    print("%d classes, entropy %f, mean_field %f, bethe %f" % (nClass, entropy, L+sum(S_mf), L+S_bethe+sum(S_mf[1:])))
else:
    S_mf = gt.mf_entropy(g, pv)
    S_bethe = gt.bethe_entropy(g, pe)[0]
    L = -np.mean(dls)
    nClass = len(np.unique(state.get_blocks().a))
    print("%d classes, entropy %f, mean_field %f, bethe %f" % (nClass, entropy, L+S_mf, L+S_bethe))

# Save final graph
(v_Block, v_BlockCC, v_BlockCC4) = block_annotation(g, state)
g.vertex_properties['Block'] = v_Block
g.vertex_properties['BlockCC'] = v_BlockCC
g.vertex_properties['BlockCC4'] = v_BlockCC4
(v_sHSBM, v_Ptu) = ptu_annotation(g)
g.vertex_properties['sHSBM'] = v_sHSBM
g.vertex_properties['Ptu'] = v_Ptu
v_HRange = hrange_annotation(g)
g.vertex_properties['HRange'] = v_HRange
fname = "SBM_mcmc_%d_%f" % (nClass, entropy)
write_classes(os.path.join(args.outdir,fname+'.tsv'), g)
pickle.dump([g, state, dls, pv, pe, S_mf, S_bethe], open(os.path.join(args.outdir,fname+'.pickle'), 'wb'), -1)

# Draw final state
if args.hierarchical:
    state_0.draw(output=os.path.join(args.outdir,fname+'.png'))
    state_0.draw(output=os.path.join(args.outdir,fname+'.svg'))
    e = state_0.get_matrix()
    plt.matshow(e.todense())
    plt.savefig(os.path.join(args.outdir,fname+'.blocks.png'))
    plt.savefig(os.path.join(args.outdir,fname+'.blocks.svg'))
else:
    state.draw(output=os.path.join(args.outdir,fname+'.png'))
    state.draw(output=os.path.join(args.outdir,fname+'.svg'))
    e = state.get_matrix()
    plt.matshow(e.todense())
    plt.savefig(os.path.join(args.outdir,fname+'.blocks.png'))
    plt.savefig(os.path.join(args.outdir,fname+'.blocks.svg'))

# Final version has to be edited by hand to rename new PTUs (those starting with ?)
pickle.dump([g, state], open(os.path.join(args.outdir,'sHSBM.pickle'), 'wb'), -1)
