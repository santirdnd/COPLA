#!/usr/bin/env python3

import os
import sys
import shutil
import pickle
import pathlib
import argparse
import subprocess
import numpy as np
import pandas as pd
import graph_tool.all as gt

parser = argparse.ArgumentParser(
    description='Predict a plasmid PTU from its sequence')
parser.add_argument('sequence', type=str,
                    help='nucleotide sequence')
parser.add_argument('refgraph', type=str,
                    help='pickle of graph with the reference sequences to compare')
parser.add_argument('reflist', type=str,
                    help='list of reference sequence filenames to compare')
parser.add_argument('outdir', type=str,
                    help='output directory')
parser.add_argument('-a', '--aminoacid', type=str,
                    help='aminoacid sequence')
parser.add_argument('-t', '--topology', type=str,
                    choices=['circular', 'linear'], default='circular',
                    help='topology of sequence')
parser.add_argument('-k', '--taxKingdom', type=str,
                    default='-',
                    help='taxon at kingdom level')
parser.add_argument('-p', '--taxPhylum', type=str,
                    default='-',
                    help='taxon at phylum level')
parser.add_argument('-c', '--taxClass', type=str,
                    default='-',
                    help='taxon at class level')
parser.add_argument('-o', '--taxOrder', type=str,
                    default='-',
                    help='taxon at order level')
parser.add_argument('-f', '--taxFamily', type=str,
                    default='-',
                    help='taxon at family level')
parser.add_argument('-g', '--taxGenus', type=str,
                    default='-',
                    help='taxon at genus level')
parser.add_argument('-s', '--taxSpecies', type=str,
                    default='-',
                    help='taxon at species level')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

def insert_ani_edges(graph, fname):
    idx = {}
    for v in graph.vertices():
        idx[graph.vp.AccessionVersion[v]] = v

    with open(fname, 'rt') as fh:
        for line in fh:
	    # ANI output format:
            # $ bin/get_ani_identity.sh queries/NC_028464.1.fna CoplaDB.fofn
            # queries/NC_028464.1.fna	NC_028464.1.fna	100.000	0.000	165	165
            # queries/NC_028464.1.fna	NC_010643.1.fna	97.546	3.385	122	165
            # queries/NC_028464.1.fna	NC_009982.1.fna	99.773	0.331	158	165
            # queries/NC_028464.1.fna	NC_010716.1.fna	97.120	2.937	133	165
            items = line.split("\t")
            ani = float(items[2]) / 100
            #if ani <= 0.7:
            #    print(line.strip())
            #    continue
            ref = items[1][0:-4]
            e = graph.add_edge(v_qry, graph.vertex(idx[ref]))
            graph.ep.ANI[e] = ani

def insert_fastani_edges(graph, fname):
    idx = {}
    for v in graph.vertices():
        idx[graph.vp.AccessionVersion[v]] = v

    with open(fname, 'rt') as fh:
        for line in fh:
	    # FastANI output format:
            # $ bin/get_fastani_identity.sh queries/NC_028464.1.fna CoplaDB.fofn
            # queries/NC_028464.1.fna	databases/Copla_RS84/NC_028464.1.fna.gz	100	22	22
            # queries/NC_028464.1.fna	databases/Copla_RS84/NC_009982.1.fna.gz	99.646	22	22
            # queries/NC_028464.1.fna	databases/Copla_RS84/NC_010716.1.fna.gz	97.642	20	22
            # queries/NC_028464.1.fna	databases/Copla_RS84/NC_010643.1.fna.gz	97.4519	18	22
            items = line.split("\t")
            ani = float(items[2]) / 100
            #if ani < 0.8:
            #    print(line.strip())
            #    continue
            ref = items[1][21:-7]
            e = graph.add_edge(v_qry, graph.vertex(idx[ref]))
            graph.ep.ANI[e] = ani

def run_mobscan(fname, outdir):
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

    fname_hmm = os.path.join(outdir, 'hmmscan.log')
    fname_dom = os.path.join(outdir, 'hmmscan.domtblout.log')
    cmd = ['hmmscan', '--cpu', '9', '--incE', '0.01', '--incdomE', '0.01', '-o', fname_hmm, '--domtblout', fname_dom, 'databases/MOBscan_171004/MOBfamDB', fname]
    # Check output of command to handle those plasmids with an empty ORFeome
    try:
        cp = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError:
        mob_label = '-'
        return(mob_label)

    fname_out = os.path.join(outdir, 'results_tab.tsv')
    cmd = ['bin/hmmscan_domtblout_summarize.py', '-e', '0.01', '-i', '0.01', '-c', '0.6', fname_dom, fname_out]
    try:
        cp = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError as msg:
        print("Error:", ' '.join(cmd), 'Failure', msg, file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as msg:
        print("Error:", ' '.join(cmd), 'Failure', msg, file=sys.stderr)
        sys.exit(1)

    mob = []
    if not os.path.isfile(fname_out):
        mob_label = '-'
        return(mob_label)
    with open(fname_out, 'rt') as fh:
        for line in fh:
            # MOBscan output format:
            # $ bin/hmmscan_domtblout_summarize.py -e 0.01 -i 0.01 -c 0.6 queries/NZ_CP014966.1.fna_mobscan/hmm.domtblout.log queries/NZ_CP014966.1.fna_mobscan/results_tab.tsv
            # NZ_CP014966.1_97  MOBH    T4SS_MOBH       1.00    38      241     6.3e-82 1e-81
            # NZ_CP014966.1_191 MOBP    T4SS_MOBP1      0.97    45      293     4.3e-63 4.3e-63
            if line.strip() == '':
                continue
            items = line.strip().split("\t", 2)
            mob.append(items[1])
    if len(mob) == 0:
        mob_label = '-'
    else:
        mob_label = ';'.join(sorted(mob))

    return(mob_label)

def run_conjscan(fname, outdir, type, topology):
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

    # Check why using a list does not work
    #cmd = ['bin/check_conjugation_systems.sh', fname, outdir, type, topology, 'databases/MacSyFinder_190530/Conjugation']
    #cp = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    try:
        cp = subprocess.run('bin/check_conjugation_systems.sh ' + fname + ' ' + outdir + ' ' + type + ' ' + topology + ' databases/MacSyFinder_190530/Conjugation', \
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    except OSError as msg:
        print("Error:", ' '.join(cmd), 'Failure', msg, file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as msg:
        print("Error:", ' '.join(cmd), 'Failure', msg, file=sys.stderr)
        sys.exit(1)

    mpf = []
    fname_out = os.path.join(outdir, 'results_tab.summary.tsv')
    # File does not exists if plasmid has no MPF
    try:
        with open(fname_out, 'rt') as fh:
            for line in fh:
                if line.strip().startswith('#'):
                    continue
                if line.strip() == '':
                    continue
                items = line.strip().split("\t", 3)
                mpf.append(items[2])
    except FileNotFoundError:
        pass
    if len(mpf) == 0:
        mpf_label = '-'
    else:
        mpf_label = ';'.join(sorted(mpf))

    return(mpf_label)

def run_pfinder(fname, outdir):
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

    # To search against last PlasmidFinder database
    #cmd = ['plasmidfinder.py', '-i', fname, '-o', outdir, '-t', '0.80', '-x', '-q']
    # To search against PlasmidFinder database version from 2019/07/31
    cmd = ['plasmidfinder.py', '-i', fname, '-o', outdir, '-t', '0.80', '-x', '-q', '-p', 'databases/PlasmidFinder_190731', '-d', 'enterobacteriaceae,gram_positive']
    try:
        cp = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError as msg:
        print("Error:", ' '.join(cmd), 'Failure', msg, file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as msg:
        print("Error:", ' '.join(cmd), 'Failure', msg, file=sys.stderr)
        sys.exit(1)

    rep = []
    fname_out = os.path.join(outdir, 'results_tab.tsv')
    with open(fname_out, 'rt') as fh:
        next(fh)
        for line in fh:
            if line.strip() == '':
                continue
            items = line.strip().split("\t", 2)
            rep.append(items[1])
    if len(rep) == 0:
        rep_label = '-'
    else:
        rep_label = ';'.join(sorted(rep))

    return(rep_label)

def run_blastn_card(fname, outdir):
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

    fname_out = os.path.join(outdir, 'card.b6')
    cmd = ['blastn', '-task', 'blastn', '-query', fname, '-db', 'databases/CARD_201015/nucleotide_fasta_protein_homolog_model.fasta', '-evalue', '1e-20',  '-perc_identity', '80', '-culling_limit', '1', '-outfmt', '6', '-out', fname_out]
    cp = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    amr = []
    with open(fname_out, 'rt') as fh:
        for line in fh:
            if line.strip() == '':
                continue
            items = line.strip().split("\t", 2)
            amr.append(items[1].split('|')[5])
    if len(amr) == 0:
        amr_label = '-'
    else:
        amr_label = ';'.join(sorted(amr))

    return(amr_label)

def write_classes(graph, filename):
    f = open(filename, 'w')
    header = "AccessionVersion\tPTU_Ref\tCComp\tBlock\tBlockCC\tBlockCC4"
    if 'sHSBM' in g.vp.keys():
        header += "\tsHSBM\tPTU\tHRange"
    f.write(header+"\n")
    for v in graph.vertices():
        AccessionVersion = graph.vp.AccessionVersion[v]
        PtuRef = str(graph.vp.PtuRef[v])
        CComp = str(graph.vp.CComp[v])
        Block = str(graph.vp.Block[v])
        BlockCC = graph.vp.BlockCC[v]
        BlockCC4 = graph.vp.BlockCC4[v]
        if 'sHSBM' in g.vp.keys():
            sHSBM = graph.vp.sHSBM[v]
            Ptu = graph.vp.Ptu[v]
            HRange = graph.vp.HRange[v]
            f.write("\t".join((AccessionVersion, PtuRef, CComp, Block, BlockCC, BlockCC4, sHSBM, Ptu, HRange)) + "\n")
        else:
            f.write("\t".join((AccessionVersion, PtuRef, CComp, Block, BlockCC, BlockCC4)) + "\n")
    f.close()

def block_annotation(graph, state):
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

def plasmid_hrange(graph):
    for l in ['TaxKingdom', 'TaxPhylum', 'TaxClass', 'TaxOrder', 'TaxFamily', 'TaxGenus', 'TaxSpecies']:
        taxa = set()
        for v in graph.vertices():
            if graph.vertex_properties[l][v] != '-':
                taxa.add(graph.vertex_properties[l][v])
        if len(taxa) > 1:
            if (l == 'TaxKingdom') or (l == 'TaxPhylum') or (l == 'TaxClass'):
                hRange = 'VI'
                break
            elif l == 'TaxOrder':
                hRange = 'V'
                break
            elif l == 'TaxFamily':
                hRange = 'IV'
                break
            elif l == 'TaxGenus':
                hRange = 'III'
                break
            elif l == 'TaxSpecies':
                hRange = 'II'
                break
    else:
        hRange = 'I'

    return hRange

def get_related_plasmids(graph, vfilt):
    ref = []
    new = []
    ref_index = {}
    new_index = {}

    u = gt.GraphView(graph, vfilt=vfilt)
    for v in u.vertices():
        if u.vp.AccessionVersion[v] == qry_acc:
            # As the query is not present in the reference network, we get rid of it to calculate the fit of both partitions
            continue
        if u.vp.sHSBMRef2[v] not in ref_index:
            if len(list(ref_index.values())) > 0:
                ref_index[u.vp.sHSBMRef2[v]] = np.max(list(ref_index.values())) + 1
            else:
                ref_index[u.vp.sHSBMRef2[v]] = 0
        ref.append(ref_index[u.vp.sHSBMRef2[v]])
        if u.vp.sHSBM2[v] not in new_index:
            if len(list(new_index.values())) > 0:
                new_index[u.vp.sHSBM2[v]] = np.max(list(new_index.values())) + 1
            else:
                new_index[u.vp.sHSBM2[v]] = 0
        new.append(new_index[u.vp.sHSBM2[v]])

    return(ref, new)

# Copy the query plasmid to the working directory
pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
shutil.copy2(args.sequence, args.outdir)
fname_fna = os.path.join(args.outdir, os.path.basename(args.sequence))

# Get total base pairs of the query. If query is a multifasta file we asume it is not a closed genome
seq_len = 0
multifasta = False
with open(fname_fna, 'rt') as fh:
    next(fh)
    for line in fh:
        if line.startswith('>'):
            multifasta = True
        else:
            seq_len += len(line.strip())

# If not provided use prodigal to get the ORFeome. As per Prokka example we use 100 Kb to use Prodigal autolearning mode
if args.aminoacid:
    shutil.copy2(args.aminoacid, args.outdir)
    fname_faa = os.path.join(args.outdir, os.path.basename(args.aminoacid))
else:
    fname_faa = fname_fna + '.faa'
    if (seq_len >= 100000):
        prodigal_mode = 'single'
    else:
        prodigal_mode = 'meta'
    cmd = ['prodigal', '-p', prodigal_mode, '-i', fname_fna, '-a', fname_faa, '-o', '/dev/null', '-q']
    try:
        cp = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError as msg:
        print("Error:", ' '.join(cmd), 'Failure', msg, file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as msg:
        print("Error:", ' '.join(cmd), 'Failure', msg, file=sys.stderr)
        sys.exit(1)

# Load RefSeq84 plasmid network
fh = open(args.refgraph, 'rb')
[g, state] = pickle.load(fh)
fh.close()

# Insert query into the network
qry_acc = '<Query>'
v_qry = g.add_vertex()
e = g.add_edge(v_qry, v_qry)
g.vp.AccessionVersion[v_qry] = qry_acc
g.ep.ANI[e] = 1.0

# Save taxonomic info of the query if available
g.vp.TaxKingdom[v_qry] = args.taxKingdom
g.vp.TaxPhylum[v_qry] = args.taxPhylum
g.vp.TaxClass[v_qry] = args.taxClass
g.vp.TaxOrder[v_qry] = args.taxOrder
g.vp.TaxFamily[v_qry] = args.taxFamily
g.vp.TaxGenus[v_qry] = args.taxGenus
g.vp.TaxSpecies[v_qry] = args.taxSpecies

# Save as our reference the existing clustering and PTU assignment from the pickle file
qry_null = '<NULL>' # Using a placeholder until PTU is assigned
g.vp.Block[v_qry] = '0' # Block is a numeric value
g.vertex_properties['BlockRef'] = g.vp.Block.copy()
g.vp.BlockCC[v_qry] = qry_null
g.vertex_properties['BlockCCRef'] = g.vp.BlockCC.copy()
g.vp.BlockCC4[v_qry] = qry_null
g.vertex_properties['BlockCC4Ref'] = g.vp.BlockCC4.copy()
g.vp.sHSBM[v_qry] = qry_null
g.vertex_properties['sHSBMRef'] = g.vp.sHSBM.copy()
g.vp.Ptu[v_qry] = qry_null
g.vertex_properties['PtuRef'] = g.vp.Ptu.copy()
g.vp.HRange[v_qry] = qry_null
g.vertex_properties['HRangeRef'] = g.vp.HRange.copy()

# Populate query edges with database plasmids
cmd = ['bin/get_ani_identity.pl', fname_fna, args.reflist]
try:
    cp = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
except OSError as msg:
    print("Error:", ' '.join(cmd), 'Failure', msg, file=sys.stderr)
    sys.exit(1)
except subprocess.CalledProcessError as msg:
    print("Error:", ' '.join(cmd), 'Failure', msg, file=sys.stderr)
    sys.exit(1)
insert_ani_edges(g, fname_fna+'.ani.tsv')

# Additional plasmid characteristics
g.vp.Size[v_qry] = seq_len
g.vp.MOB[v_qry] = run_mobscan(fname_faa, fname_fna+'_mobscan')
genome_type = 'unordered_replicon' if multifasta else 'ordered_replicon'
g.vp.MPF[v_qry] = run_conjscan(fname_faa, fname_fna+'_conjscan', genome_type, args.topology)
g.vp.PFinder[v_qry] = run_pfinder(fname_fna, fname_fna+'_pfinder')
g.vp.AMR[v_qry] = run_blastn_card(fname_fna, fname_fna+'_amr')

with open(fname_fna+'.qry_info.tsv', 'w') as fh:
    fh.write("#Total bp\tMOB\tMPF\tReplicon\tAMR\n")
    fh.write("{}\t{}\t{}\t{}\t{}\n".format(g.vp.Size[v_qry], g.vp.MOB[v_qry], g.vp.MPF[v_qry], g.vp.PFinder[v_qry], g.vp.AMR[v_qry]))

# Check if query belongs to a graph component of less than 4 members
comp, hist = gt.label_components(g)
g.vertex_properties['CComp'] = comp
if hist[comp[v_qry]] <= 4:
    i_filter = (comp.a == comp[v_qry])
    u = gt.GraphView(g, vfilt=i_filter)
    cl_size = u.num_vertices()
    if cl_size < 4:
        ptu_pred = '-'
        print('PTU could not be assigned')
        print('Query is part of a graph component of size {}'.format(cl_size))
        print('However, at least four members are required for PTU assignation')
        print('This plasmid could form part of a new, still unnamed, PTU')
        print('Other info:')
        print("  Size:\t{}".format(g.vp.Size[v_qry]))
        print("  MOB:\t{}".format(g.vp.MOB[v_qry]))
        print("  MPF:\t{}".format(g.vp.MPF[v_qry]))
        print("  Repl:\t{}".format(g.vp.PFinder[v_qry]))
        print("  AMR:\t{}".format(g.vp.AMR[v_qry]))
    else:
        ptu_pred = 'PTU-?'
        print('New (putative) PTU')
        print('Query is part of a graph component of size {}'.format(cl_size))
        print('This plasmid could form part of a new, still unnamed, PTU')
        print('Other info:')
        print("  Size:\t{}".format(g.vp.Size[v_qry]))
        print("  MOB:\t{}".format(g.vp.MOB[v_qry]))
        print("  MPF:\t{}".format(g.vp.MPF[v_qry]))
        print("  Repl:\t{}".format(g.vp.PFinder[v_qry]))
        print("  AMR:\t{}".format(g.vp.AMR[v_qry]))
    with open(fname_fna+'.ptu_prediction.tsv', 'w') as fh:
        fh.write("#Predicted_PTU\tHost_Range\tScore\tNotes\n")
        score = 1.0
        if cl_size < 4:
            hrange = '-'
            notes = 'PTU could not be assigned. Query is part of a graph component of size {}'.format(cl_size)
        else:
            hrange = plasmid_hrange(u)
            notes = 'New (putative) PTU. Query is part of a graph component of size {}'.format(cl_size)
        fh.write("{}\t{}\t{:.4f}\t{}\n".format(ptu_pred, hrange, score, notes))
    with open(fname_fna+'.related_plasmids.tsv', 'w') as fh:
        fh.write("#AccessionVersion\tPTU_Ref\tNote\n")
        for v in u.vertices():
            note = '*'
            fh.write("{}\t{}\t{}\n".format(u.vp.AccessionVersion[v], u.vp.PtuRef[v], note))
    sys.exit()

# Apply SBM to cluster the query
new_state = state.copy()
new_state.multiflip_mcmc_sweep(d=0, psplit=0, pmerge=0, pmergesplit=0, niter=1000)

# Apply sHSBM algorithm to the new graph with the query
(v_Block, v_BlockCC, v_BlockCC4) = block_annotation(g, new_state)
g.vertex_properties['Block'] = v_Block
g.vertex_properties['BlockCC'] = v_BlockCC
g.vertex_properties['BlockCC4'] = v_BlockCC4
(v_sHSBM, v_Ptu) = ptu_annotation(g)
g.vertex_properties['sHSBM'] = v_sHSBM
g.vertex_properties['Ptu'] = v_Ptu
v_HRange = hrange_annotation(g)
g.vertex_properties['HRange'] = v_HRange
write_classes(g, fname_fna+'.sHSBM.tsv')

# TODO: Code is easier if sHSBM column has the data of all clusters and not only those with 4 or more members
v_sHSBM2 = g.vp.sHSBM.copy()
v_sHSBMRef2 = g.vp.sHSBMRef.copy()
for v in g.vertices():
    if v_sHSBM2[v] == '-':
        v_sHSBM2[v] = g.vp.BlockCC[v]
    if v_sHSBMRef2[v] == '-':
        v_sHSBMRef2[v] = g.vp.BlockCCRef[v]
g.vertex_properties['sHSBM2'] = v_sHSBM2
g.vertex_properties['sHSBMRef2'] = v_sHSBMRef2

# Get all plasmids of the PTU assigned to the query
strict_filter = (np.array(list(g.vp.sHSBM2)) == g.vp.sHSBM2[v_qry])
#(ref_strict, new_strict) = get_related_plasmids(g, strict_filter)
#if len(ref_strict) == 0:
#    # Query was assigned to a singleton
#    overlap_strict = 1.0
#else:
#    overlap_strict = gt.partition_overlap(ref_strict, new_strict, norm=True)

# Use the reference clustering of these plasmids to expand the selection to all plasmids of all PTUs with a member in the new PTU assigned to the query
sHSBMRef2_list = set()
u = gt.GraphView(g, vfilt=strict_filter)
for v in u.vertices():
    sHSBMRef2_list.add(u.vp.sHSBMRef2[v])
expand_filter = np.full_like(strict_filter, False)
for i in sHSBMRef2_list:
    i_filter = (np.array(list(g.vp.sHSBMRef2)) == i)
    expand_filter = np.logical_or(expand_filter, i_filter)
(ref_expand, new_expand) = get_related_plasmids(g, expand_filter)
if len(ref_expand) == 0:
    # Query was assigned to a singleton
    overlap_expand = 1.0
else:
    overlap_expand = gt.partition_overlap(ref_expand, new_expand, norm=True)

# Predicted PTU is the most frequent reference label among the new cluster
w = gt.GraphView(g, vfilt=expand_filter)
cl_size = u.num_vertices()
tmp_ptu = np.delete(np.array(list(u.vp.PtuRef)), np.argwhere(np.array(list(u.vp.PtuRef)) == qry_null))
if cl_size < 4:
    ptu_pred = '-'
    ptu_related = '-'
    # Make sure ptu_related is not '-' if there is another option
    # This works because cl_size <= 3 and one of them is the previously deleted qry_null
    for i in tmp_ptu:
        if i != '-':
            ptu_related = i
    print('PTU could not be assigned')
    print('Query is part of a sHSBM cluster of size {}'.format(cl_size))
    print('However, at least four members are required for PTU assignation')
    print('This plasmid could form part of a new, still unnamed, PTU')
    if ptu_related != '-':
        print('Query is related to {} plasmids'.format(ptu_related))
    print('Other info:')
    print("  Size:\t{}".format(g.vp.Size[v_qry]))
    print("  MOB:\t{}".format(g.vp.MOB[v_qry]))
    print("  MPF:\t{}".format(g.vp.MPF[v_qry]))
    print("  Repl:\t{}".format(g.vp.PFinder[v_qry]))
    print("  AMR:\t{}".format(g.vp.AMR[v_qry]))
    with open(fname_fna+'.ptu_prediction.tsv', 'w') as fh:
        fh.write("#Predicted_PTU\tHost_Range\tScore\tNotes\n")
        hrange = '-'
        score = overlap_expand
        notes = 'PTU could not be assigned. Query is part of a sHSBM cluster of size {}'.format(cl_size)
        fh.write("{}\t{}\t{:.4f}\t{}\n".format(ptu_pred, hrange, score, notes))
    with open(fname_fna+'.related_plasmids.tsv', 'w') as fh:
        fh.write("#AccessionVersion\tPTU_Ref\tNote\n")
        for v in w.vertices():
            note = '*' if g.vp.sHSBM[v] == g.vp.sHSBM[v_qry] else ''
            fh.write("{}\t{}\t{}\n".format(w.vp.AccessionVersion[v], w.vp.PtuRef[v], note))
else:
    unique, counts = np.unique(list(tmp_ptu), return_counts=True)
    ptu_pred = unique[counts.argmax()]
    ptu_related = '-'
    if (ptu_pred == '-') or (ptu_pred != u.vp.Ptu[v_qry]):
        if ptu_pred != u.vp.Ptu[v_qry]:
            ptu_related = ptu_pred
        ptu_pred = 'PTU-?'
        print('New (putative) PTU')
        print('Query is part of a sHSBM cluster of size {}'.format(cl_size))
        print('This plasmid could form part of a new, still unnamed, PTU')
        if ptu_related != '-':
            print('Query is related to {} plasmids'.format(ptu_related))
        print('Other info:')
        print("  Size:\t{}".format(g.vp.Size[v_qry]))
        print("  MOB:\t{}".format(g.vp.MOB[v_qry]))
        print("  MPF:\t{}".format(g.vp.MPF[v_qry]))
        print("  Repl:\t{}".format(g.vp.PFinder[v_qry]))
        print("  AMR:\t{}".format(g.vp.AMR[v_qry]))
        with open(fname_fna+'.ptu_prediction.tsv', 'w') as fh:
            fh.write("#Predicted_PTU\tHost_Range\tScore\tNotes\n")
            hrange = plasmid_hrange(u)
            score = overlap_expand
            notes = 'New (putative) PTU. Query is part of a sHSBM cluster of size {}'.format(cl_size)
            fh.write("{}\t{}\t{:.4f}\t{}\n".format(ptu_pred, hrange, score, notes))
        with open(fname_fna+'.related_plasmids.tsv', 'w') as fh:
            fh.write("#AccessionVersion\tPTU_Ref\tNote\n")
            for v in w.vertices():
                note = '*' if g.vp.sHSBM[v] == g.vp.sHSBM[v_qry] else ''
                fh.write("{}\t{}\t{}\n".format(w.vp.AccessionVersion[v], w.vp.PtuRef[v], note))
    else:
        print('Query is a {} plasmid'.format(ptu_pred))
        print('Query is part of a sHSBM cluster of size {}'.format(cl_size))
        print('Other info:')
        print("  Size:\t{}".format(g.vp.Size[v_qry]))
        print("  MOB:\t{}".format(g.vp.MOB[v_qry]))
        print("  MPF:\t{}".format(g.vp.MPF[v_qry]))
        print("  Repl:\t{}".format(g.vp.PFinder[v_qry]))
        print("  AMR:\t{}".format(g.vp.AMR[v_qry]))
        with open(fname_fna+'.ptu_prediction.tsv', 'w') as fh:
            fh.write("#Predicted_PTU\tHost_Range\tScore\tNotes\n")
            p_filter = (np.array(list(g.vp.PtuRef)) == ptu_pred)
            p_filter[int(v_qry)] = True # Include query
            x = gt.GraphView(g, vfilt=p_filter)
            hrange = plasmid_hrange(x)
            score = overlap_expand
            notes = 'Query is a {} plasmid'.format(ptu_pred)
            hr_ref = g.vp.HRangeRef[list(x.vertices())[0]]
            if hr_ref != hrange:
                notes += '. Former PTU host range was {}'.format(hr_ref)
                print('Query inclusion has caused PTU host range to increase from {} to {}'.format(hr_ref, hrange))
            fh.write("{}\t{}\t{:.4f}\t{}\n".format(ptu_pred, hrange, score, notes))
        with open(fname_fna+'.related_plasmids.tsv', 'w') as fh:
            fh.write("#AccessionVersion\tPTU_Ref\tNote\n")
            for v in w.vertices():
                note = '*' if g.vp.sHSBM[v] == g.vp.sHSBM[v_qry] else ''
                fh.write("{}\t{}\t{}\n".format(w.vp.AccessionVersion[v], w.vp.PtuRef[v], note))
