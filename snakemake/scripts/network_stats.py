import os
os.environ["OMP_NUM_THREADS"] = '8' # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = '8' # export OPENBLAS_NUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = '8' # export NUMEXPR_NUM_THREADS=6

import sys,os
import logging, traceback
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )
# Install exception handler
sys.excepthook = handle_exception

from graph_tool.all import *
import pandas as pd
import numpy as np
import scipy as sp
import statsmodels.api as sm

from multipy.fdr import qvalue
from multipy.fdr import lsu

def filterByEdge(g, corr, cutOff, keepOnlyMain):
    # Filtering edges
    corr = g.edge_properties[corr]
    sign = g.new_ep("bool", True)
    sign.a = np.array(np.abs(corr.a) > cutOff)

    tv = GraphView(g, efilt=sign)

    # Keeping largest component
    if keepOnlyMain:
        comp, hist = label_components(tv)
        main_component = tv.new_vp("bool", (comp.a == np.where(hist == max(hist))[0][0]))
        tv.vertex_properties["main_component"] = main_component
        tv.set_vertex_filter(main_component)
    return tv

def filterByFDR(g, level, keepOnlyMain):
    # Filtering edges
    pvals = np.array(g.edge_properties["pvalue"].a)

    fdr_ep = g.new_ep("bool", True)
    fdr_ep.a = lsu(pvals, q=level)

    tv = GraphView(g, efilt=fdr_ep)

    # Keeping largest component
    if keepOnlyMain:
        comp, hist = label_components(tv)
        main_component = tv.new_vp("bool", (comp.a == np.where(hist == max(hist))[0][0]))
        tv.vertex_properties["main_component"] = main_component
        tv.set_vertex_filter(main_component)
    return tv

def getGeneNetworkStats(g):
    genes = g.vertex_properties["genes"]
    corr = g.edge_properties["spearman"]
    block_df = pd.DataFrame(columns=('Gene', "Average_Spearman",
                                     'Degree_thr_0.1', 'WeightedDegree_thr_0.1',
                                     'Degree_thr_0.2', 'WeightedDegree_thr_0.2',
                                     'Degree_thr_0.3', 'WeightedDegree_thr_0.3',
                                     'Degree_thr_0.4', 'WeightedDegree_thr_0.4',
                                     'Degree_thr_0.5', 'WeightedDegree_thr_0.5',
                                     'Degree_fdr_5e-2' , 'WeightedDegree_fdr_5e-2',
                                     'Degree_fdr_1e-2' , 'WeightedDegree_fdr_1e-2',
                                     'Degree_fdr_1e-3' , 'WeightedDegree_fdr_1e-3',
                                     'Degree_fdr_1e-4' , 'WeightedDegree_fdr_1e-4',
                                     'Degree_fdr_1e-5' , 'WeightedDegree_fdr_1e-5',
                                     'Degree_fdr_1e-6' , 'WeightedDegree_fdr_1e-6'))

    tv_01 = filterByEdge(g, "spearman", 0.1, False)
    tv_02 = filterByEdge(g, "spearman", 0.2, False)
    tv_03 = filterByEdge(g, "spearman", 0.3, False)
    tv_04 = filterByEdge(g, "spearman", 0.4, False)
    tv_05 = filterByEdge(g, "spearman", 0.5, False)

    tv_fdr5e2 = filterByFDR(g, 5e-2, False)
    tv_fdr1e2 = filterByFDR(g, 1e-2, False)
    tv_fdr1e3 = filterByFDR(g, 1e-3, False)
    tv_fdr1e4 = filterByFDR(g, 1e-4, False)
    tv_fdr1e5 = filterByFDR(g, 1e-5, False)
    tv_fdr1e6 = filterByFDR(g, 1e-6, False)

    for v in g.vertex_index:
        line = [genes[v]]
        line.append(np.sum(np.abs(g.get_all_edges(v, [corr] )[:,2])))

        line.append(tv_01.get_total_degrees([v])[0])
        line.append(np.sum(np.abs(tv_01.get_all_edges(v, [corr] )[:,2])))

        line.append(tv_02.get_total_degrees([v])[0])
        line.append(np.sum(np.abs(tv_02.get_all_edges(v, [corr] )[:,2])))

        line.append(tv_03.get_total_degrees([v])[0])
        line.append(np.sum(np.abs(tv_03.get_all_edges(v, [corr] )[:,2])))

        line.append(tv_04.get_total_degrees([v])[0])
        line.append(np.sum(np.abs(tv_04.get_all_edges(v, [corr] )[:,2])))
        
        line.append(tv_05.get_total_degrees([v])[0])
        line.append(np.sum(np.abs(tv_05.get_all_edges(v, [corr] )[:,2])))

        line.append(tv_fdr5e2.get_total_degrees([v])[0])
        line.append(np.sum(np.abs(tv_fdr5e2.get_all_edges(v, [corr] )[:,2])))

        line.append(tv_fdr1e2.get_total_degrees([v])[0])
        line.append(np.sum(np.abs(tv_fdr1e2.get_all_edges(v, [corr] )[:,2])))

        line.append(tv_fdr1e3.get_total_degrees([v])[0])
        line.append(np.sum(np.abs(tv_fdr1e3.get_all_edges(v, [corr] )[:,2])))

        line.append(tv_fdr1e4.get_total_degrees([v])[0])
        line.append(np.sum(np.abs(tv_fdr1e4.get_all_edges(v, [corr] )[:,2])))
        
        line.append(tv_fdr1e5.get_total_degrees([v])[0])
        line.append(np.sum(np.abs(tv_fdr1e5.get_all_edges(v, [corr] )[:,2])))

        line.append(tv_fdr1e6.get_total_degrees([v])[0])
        line.append(np.sum(np.abs(tv_fdr1e6.get_all_edges(v, [corr] )[:,2])))

        block_df.loc[v] = line
    return block_df

if __name__ == '__main__':

    logging.info("Reading full graph...")
    g = load_graph(snakemake.input.graph)

    logging.info("Calculating stats...")
    out = getGeneNetworkStats(g)

    logging.info("Writing table")
    out.to_csv(snakemake.output.stats)