"""
####################
devel.py
####################
Code under development before being placed in a logical location.
"""

import numpy as np

from scipy import stats as sp_stats

from gfunc.maths import weight_d_for_ptci as scale_the_d
from gfunc.data_classes import Bunch


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

######################## END CALC PTCI LAND ######################################



def get_ptci(gFunc_edge,kind='rpd',w_min=1.0,w_max=1.1):
    """
    calculate and store 
    """
    ptci_kind = Bunch({ 'rpd' : calc_ptci_rpd,
                        'zpd' : calc_ptci_zpd,
                        'rd'  : calc_ptci_rd,
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
