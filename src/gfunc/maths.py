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
    BS = ((c * m) + sum([x for x in scores])) / (n + c)
    
    Where:
    
    n: number of votes for THIS item
    C: median number of votes for all items that got at least 1 vote (weighting or dampening factor)
    m: median UNweighted score for all items that got at least 1 vote  
    """
    c = float(c)
    
    bs = ((c * scale_mod * m) + sum([x for x in scores])) / (len(scores) + c)
    if np.isnan(bs):
        raise ValueError
    return bs
    