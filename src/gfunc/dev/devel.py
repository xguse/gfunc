"""
####################
devel.py
####################
Code under elopment before being placed in a logical location.
"""

import numpy as np
import pandas

from collections import defaultdict

from scipy import stats as sp_stats

from gfunc.maths import weight_d_for_ptci as scale_the_d
from gfunc import maths as m
from gfunc.data_classes import Bunch
from gfunc import xpermutations


# Clever trick to create endlessly dimentional nested dicts on the fly
dict_tree = lambda: defaultdict(dict_tree)

def is_none_or_nan(value):
    """
    Returns True if value is None or nan, False otherwise.
    """
    if value == None:
        return True
    else:
        pass
    
    try:
        return np.all(np.isnan(value))
    except TypeError:
        return False

def edge_correlation(gFunc_edge):
    """
    Returns tuple of the pearson r value and corresponding p-value
    for a given edge's nodes is it is possible, returns None otherwise.
    """
    #try:
        #r_val,p_val = gFunc_edge.data.r_val , gFunc_edge.data.p_val
        #return r_val,p_val
    #except AttributeError:
        #pass
    
    node1,node2 = gFunc_edge.nodes
    try:
        r_val,p_val = sp_stats.pearsonr(node1.data.expression_vector, node2.data.expression_vector)
        if is_none_or_nan(r_val):
            gFunc_edge.data.r_val = None
            gFunc_edge.data.p_val = None
        else:
            gFunc_edge.data.r_val = r_val
            gFunc_edge.data.p_val = p_val
            return r_val,p_val
    
    except AttributeError as err:
        if """'Bunch' object has no attribute""" in err.message:
            # if this executes then one of the nodes did not have an expression_vector which means no r is possible
            # in this case return None
            gFunc_edge.data.r_val = None
            gFunc_edge.data.p_val = None
            return None
        else:
            # if this executes then something ELSE went wrong: thus I will fail.
            raise err



def get_edge_r_and_p_vals(edges,quiet=True):
    """
    set and get r and p vals from list of edges
    """
    # collect all the results using edge_correlation()
    edge_r_and_p_values = [edge_correlation(edge) for edge in edges]
    
    if not quiet:
        print "r_vals before cleaning: %s" % (len(edge_r_and_p_values))

    # get rid of any results that equal None
    edge_r_and_p_values = [x for x in edge_r_and_p_values if not is_none_or_nan(x)]
    
    if not quiet:
        print "Returning %s r_vals." % (len(edge_r_and_p_values))
        
    return edge_r_and_p_values

def get_pairwise_ptci_vals(edges,kind='rpd',quiet=True,w_min=1.0,w_max=1.1):
    """
    Function to calculate and store the 1-to-1 pairwise PTCI values in the graph edges
    """
    if not quiet:
        print "kind: %s" % (kind)
    pairwise_ptci_vals = [get_ptci(edge,kind,w_min,w_max) for edge in edges]
    if not quiet:
        print "ptci_vals before cleaning: %s" % (len(pairwise_ptci_vals))
    # remove any None values
    pairwise_ptci_vals = [ptci for ptci in pairwise_ptci_vals if not is_none_or_nan(ptci)]
    if not quiet:
        print "Returning %s ptci_vals." % (len(pairwise_ptci_vals))
        
    return pairwise_ptci_vals

def get_null_pairwise_ptci_distributions(graphHandler,graphBuilder,n_way_ortho_table,ortho_parser,reps=50,kind='rpd',quiet=True,w_min=1.0,w_max=1.1):
    """
    Function to calculate and store RANDOMIZED 1-to-1 pairwise PTCI values in the graph edges to generate many NULL distributions
    """
    null_paired_ptci_distributions = []

    for rep in range(reps):
        # scramble edges for this rep and set new r&p vals
        reset_random_edges(graphHandler,graphBuilder,n_way_ortho_table,ortho_parser)
        graphHandler.measure_relations()
        
        # do prep
        null_edges = graphHandler.edge_dict.values()
        null_r_and_p_values = get_edge_r_and_p_vals(null_edges,quiet)
        null_r_values = [null_r_and_p_values[i][0] for i in range(len(null_r_and_p_values))]
        null_z_stats = m.get_z_score_stats(null_r_values)
        set_z_vals(graphHandler,null_z_stats,use_center='median')
        
        # calculate null ptci vals
        null_pairwise_ptci_vals = get_pairwise_ptci_vals(null_edges,kind,quiet,w_min,w_max)
        
        # collect null ptci distribution
        null_paired_ptci_distributions.append(null_pairwise_ptci_vals)
        
    
    return null_paired_ptci_distributions


def get_mean_ptcis(graphHandler,n_way_ortho_table,kind='rpd',quiet=True,w_min=1.0,w_max=1.1):
    """
    returns list of n-way averaged PTCI values for N-way ortholog subgraphs if and only if 
    all edges successfully generated non-None value PTCI results.
    
    Stores mean ptci in ortho_edge.data.meanPTCI
    """
    
    # dictionary of nodes/edges indexed by gene names
    node_dict = graphHandler.node_dict
    edge_dict = graphHandler.edge_dict
    
    graph = graphHandler.graph 
    
    mean_ptcis = []
    
    # calculate all pairwise combinations of indexes
    # so that each ortho-edge of n-way orthos are obtained
    
    index_combos = list(xpermutations.xuniqueCombinations(range(len(n_way_ortho_table.columns)),2))
    
    
    for node_list in n_way_ortho_table.itertuples():
        
        node_list = node_list[1:]
        
        ortho_edges = []
        for i in index_combos:
            key = tuple(sorted([node_list[i[0]],node_list[i[1]]]))
            
            try:
                ortho_edges.append(edge_dict[key])
            except KeyError:
                break
                
        ptcis = [get_ptci(edge,kind,w_min,w_max) for edge in ortho_edges]
        
        try:
            mean_ptci = np.mean(ptcis)
            
            #store mean ptci in each member edge
            for ortho_edge in ortho_edges:
                ortho_edge.data.meanPTCI = mean_ptci
                
            mean_ptcis.append(mean_ptci)
        except TypeError:
            pass

    return mean_ptcis



def set_z_vals(graphHandler,z_stats,use_center='median'):
    """
    Function to use z-score stats to calculate and store z-score converted r-values in the gFunc graph
    """
    z_stats = {'mean':z_stats[0],'median':z_stats[1],'stdv':z_stats[2]}
    
    center = z_stats[use_center]
    stdv   = z_stats['stdv']
    
    def z_val(r_val,center,stdv):
        return  (r_val - center) / stdv
    
    edges = graphHandler.edge_dict.values() 
    for edge in edges:
        try:
            edge.data.z_val = z_val(edge.data.r_val,center,stdv)
            
        except (TypeError,AttributeError) as exc:
            if 'TypeError' in str(exc):
                edge.data.z_val = None
            elif 'AttributeError' in str(exc):
                edge_correlation(edge)
                if edge.data.r_val == None:
                    edge.data.z_val = None
                else:
                    edge.data.z_val = z_val(edge.data.r_val,center,stdv)
            


def get_null_mean_ptci_distributions(graphHandler,graphBuilder,n_way_ortho_table,ortho_parser,reps=50,kind='rpd',quiet=True,w_min=1.0,w_max=1.1):
    """
    Function to calculate RANDOMIZED mean PTCI values to generate many NULL distributions
    """
    null_mean_ptci_distributions = []

    for rep in range(reps):
        # scramble edges for this rep and set new r&p vals
        reset_random_edges(graphHandler,graphBuilder,n_way_ortho_table,ortho_parser)
        graphHandler.measure_relations()
        
        
        # do prep
        null_edges = graphHandler.edge_dict.values()
        null_r_and_p_values = get_edge_r_and_p_vals(null_edges,quiet)
        null_r_values = [null_r_and_p_values[i][0] for i in range(len(null_r_and_p_values))]
        null_z_stats = m.get_z_score_stats(null_r_values)
        set_z_vals(graphHandler,null_z_stats,use_center='median')
        
        # calculate null ptci vals
        null_mean_ptci_vals = calc_mean_ptcis(graphHandler,n_way_ortho_table,kind,quiet,w_min,w_max)
        
        # collect null ptci distribution
        null_mean_ptci_distributions.append(null_mean_ptci_vals)
        
    
    return null_mean_ptci_distributions





######################## START CALC  PTCI LAND ######################################

def calc_ptci_rpd(gFunc_edge,w_min=1.0,w_max=1.1):
    """
    calculate the PTCI
    edge_correlation() should already have been run
    """
    try:
        r_val = gFunc_edge.data.r_val
        p_val = gFunc_edge.data.p_val
        d_val,d_min,d_max = gFunc_edge.data.divergence
        
        return r_val * (1-p_val) * scale_the_d(d_val,d_min,d_max,w_min,w_max) 
    
    except (TypeError,AttributeError):
        return None


def calc_ptci_rd(gFunc_edge,w_min=1.0,w_max=1.1):
    """
    calculate the PTCI
    edge_correlation() should already have been run
    """
    try:  
        r_val = gFunc_edge.data.r_val
        d_val,d_min,d_max = gFunc_edge.data.divergence
        
        return r_val * scale_the_d(d_val,d_min,d_max,w_min,w_max) 
    
    except (TypeError,AttributeError):
        return None
        
def calc_ptci_zd(gFunc_edge,w_min=1.0,w_max=1.1):
    """
    calculate the PTCI
    edge_correlation() should already have been run
    """
    try:  
        z_val = gFunc_edge.data.z_val
        d_val,d_min,d_max = gFunc_edge.data.divergence
        
        return z_val * scale_the_d(d_val,d_min,d_max,w_min,w_max) 
    
    except (TypeError,AttributeError):
        return None

def calc_ptci_zpd(gFunc_edge,w_min=1.0,w_max=1.1):
    """
    calculate the PTCI
    edge_correlation() should already have been run
    """
    try:
        z_val = gFunc_edge.data.z_val
        p_val = gFunc_edge.data.p_val
        d_val,d_min,d_max = gFunc_edge.data.divergence
        
        return z_val * (1-p_val) * scale_the_d(d_val,d_min,d_max,w_min,w_max) 
    
    except (TypeError,AttributeError):
        return None

def calc_ptci_z(gFunc_edge,w_min=1.0,w_max=1.1):
    """
    calculate the PTCI
    edge_correlation() should already have been run
    """
    try:
        z_val = gFunc_edge.data.z_val

        return z_val 
    
    except (TypeError,AttributeError):
        return None

def calc_ptci_r(gFunc_edge,w_min=1.0,w_max=1.1):
    """
    calculate the PTCI
    edge_correlation() should already have been run
    """
    try:
        r_val = gFunc_edge.data.r_val

        return r_val 
    
    except (TypeError,AttributeError):
        return None

def calc_ptci_rprd(gFunc_edge,w_min=1.0,w_max=1.1):
    """
    calculate the PTCI
    edge_correlation() should already have been run
    
    NOTE: this adds tfbs similarity to the equation
    """
    try:
        r_val_expn = gFunc_edge.data.r_val
        p_val_expn = gFunc_edge.data.p_val
        r_val_tfbs = gFunc_edge.data.tfbs_vector_similarity.r_val
        d_val,d_min,d_max = gFunc_edge.data.divergence
        
        return (r_val_expn * (1-p_val_expn)) * r_val_tfbs * scale_the_d(d_val,d_min,d_max,w_min,w_max) 
    
    except (TypeError,AttributeError):
        return None

def calc_ptci_rprpd(gFunc_edge,w_min=1.0,w_max=1.1):
    """
    calculate the PTCI
    edge_correlation() should already have been run
    
    NOTE: this adds tfbs similarity AND tfbs p of r to the equation
    """
    try:
        r_val_expn = gFunc_edge.data.r_val
        p_val_expn = gFunc_edge.data.p_val
        r_val_tfbs = gFunc_edge.data.tfbs_vector_similarity.r_val
        p_val_tfbs = gFunc_edge.data.tfbs_vector_similarity.p_val
        d_val,d_min,d_max = gFunc_edge.data.divergence
        
        return (r_val_expn * (1-p_val_expn)) * (r_val_tfbs * (p_val_tfbs)) * scale_the_d(d_val,d_min,d_max,w_min,w_max) 
    
    except (TypeError,AttributeError):
        return None


def calc_ptci_rrd(gFunc_edge,w_min=1.0,w_max=1.1):
    """
    calculate the PTCI
    edge_correlation() should already have been run
    
    NOTE: this adds tfbs similarity to the equation
    """
    try:
        r_val_expn = gFunc_edge.data.r_val
        p_val_expn = gFunc_edge.data.p_val
        r_val_tfbs = gFunc_edge.data.tfbs_vector_similarity.r_val
        d_val,d_min,d_max = gFunc_edge.data.divergence
        
        return (r_val_expn + r_val_tfbs) * scale_the_d(d_val,d_min,d_max,w_min,w_max) 
    
    except (TypeError,AttributeError):
        return None

def calc_ptci_rsrd(gFunc_edge,w_min=1.0,w_max=1.1):
    """
    calculate the PTCI
    edge_correlation() should already have been run
    
    NOTE: this adds tfbs similarity to the equation
    """
    try:
        r_val_expn = gFunc_edge.data.r_val
        r_val_tfbs = gFunc_edge.data.tfbs_vector_similarity.r_val
        d_val,d_min,d_max = gFunc_edge.data.divergence
        
        return (r_val_expn + (r_val_tfbs/2.0) * scale_the_d(d_val,d_min,d_max,w_min,w_max)) 
    
    except (TypeError,AttributeError):
        return None
    
def calc_ptci_srsrd(gFunc_edge,w_min=1.0,w_max=1.1):
    """
    calculate the PTCI
    edge_correlation() should already have been run
    
    NOTE: this adds tfbs similarity to the equation
    """
    try:
        r_val_expn = gFunc_edge.data.r_val
        p_val_expn = gFunc_edge.data.p_val
        r_val_tfbs = gFunc_edge.data.tfbs_vector_similarity.r_val
        d_val,d_min,d_max = gFunc_edge.data.divergence
        
        return (scale_r_val(r_val_expn) + scale_r_val(r_val_tfbs)) * scale_the_d(d_val,d_min,d_max,w_min,w_max) 
    
    except (TypeError,AttributeError):
        return None

######################## END CALC PTCI LAND ######################################

def scale_r_val(r_val):
    """
    returns r_val scaled to between 0 and 1
    """
    return (r_val + 1.0)/3.0


def get_ptci(gFunc_edge,kind='rpd',w_min=1.0,w_max=1.1):
    """
    calculate and store 
    """
    ptci_kind = Bunch({ 'rprpd' : calc_ptci_rprpd,
                        'rprd' : calc_ptci_rprd,
                        'rpd' : calc_ptci_rpd,
                        'zpd' : calc_ptci_zpd,
                        'rd'  : calc_ptci_rd,
                        'rsrd'  : calc_ptci_rsrd,
                        'srsrd'  : calc_ptci_srsrd,
                        'rrd'  : calc_ptci_rrd,
                        'zd'  : calc_ptci_zd,
                        'r'   : calc_ptci_r,
                        'z'   : calc_ptci_z })
    
    ptci = ptci_kind[kind](gFunc_edge,w_min,w_max)
    
    if not is_none_or_nan(ptci):
        # If we get a valid ptci store the value in the gFunc_edge object and also return it
        gFunc_edge.data.PTCI = ptci
        return ptci
    else:
        # If we get an invalid ptci, store it and return it as None
        gFunc_edge.data.PTCI = None
        return None


def get_null_bin_members_counts_dataframe(null_hist_data_list):
    """
    """
    
    n_dict = {}
    for i,nHist_data in enumerate(null_hist_data_list):
        n_dict['nd_%s' % (i)] = nHist_data[0]
        
    null_bin_members_counts = pandas.DataFrame(n_dict)
    null_bin_members_counts = null_bin_members_counts.transpose()
    
    return null_bin_members_counts