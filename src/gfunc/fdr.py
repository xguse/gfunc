"""
####################
fdr.py
####################
Code supporting easy generation of emprical FDR for results.
"""

import random 

def shuffle_dict(original_dict,shuffle_count):
    """
    RETURNS:
        - generator that yields randomly shuffled key/value pairs as new dict
          <shuffle_count> times.
    """
    original_keys = tuple(original_dict.keys())
    random_keys = original_dict.keys()
    
    for c in range(shuffle_count):
        new_dict = {}
        random.shuffle(random_keys)
        
        for i,val in enumerate(random_keys):
            new_dict[val] = original_dict[original_keys[i]]
        
        yield new_dict

    
#def get_fdr_full(original_dict,positive_keys,score_threshold,scoring_func,fdr_reps):
    #"""
    #RETURNS:
        #- 
    #"""