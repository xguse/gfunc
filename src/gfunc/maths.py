"""
####################
maths.py
####################
Code supporting specialized calculations for gfunc.
"""

import numpy as np
#np.seterr(divide='raise')

def bayesian_score(c,m,n,scores,scale_mod=1):
    """
    | BS = ((c * m) + sum([x for x in scores])) / (n + c)
    | 
    | Where:
    | 
    | ``n``: number of votes for THIS item
    | ``C``: median number of votes for all items that got at least 1 vote (weighting or dampening factor)
    | ``m``: median UNweighted score for all items that got at least 1 vote  
    """
    c = float(c)
    
    bs = ((c * scale_mod * m) + sum([x for x in scores])) / (len(scores) + c)
    if np.isnan(bs):
        raise ValueError
    return bs
    
def weight_d_for_ptci(d_i,d_min,d_max,w_min=1.0,w_max=1.1):
    """
    Scaling function to transform 'd' onto a weight-spectrum to either
    punish or reward the final ptci score based on the phylogenetic distance
    between the two current species as it relates to the range of phylogenetic
    distances in the data set.
    
    Default weight scale is no change for the shortest distance (return 1.0) to
    a 10% reward for the longest distance (return 1.1).
    """
    if (d_i < d_min) or (d_i > d_max):
        raise ValueError("d_i (%s) is outside of range (%s to %s)" % (d_i,d_min,d_max))
    
    return ((d_i-d_min)*(w_max-w_min))/(d_max-d_min) + w_min