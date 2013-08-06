"""
####################
motifs.py
####################
Code supporting the searching, recording, and analysis of sequence motifs for gFunc.
"""

import cPickle

from collections import defaultdict

import numpy as np
from scipy import stats
import pandas


import MOODS

from gfunc.parsers.JASPAR import ParseJasparMatrixOnly
from gfunc.fileIO import tableFile2namedTuple

from rSeq.utils.files import ParseFastA

#def dict_tree():
    #"""
    #GIVEN:
        #- nothing
    #DOES:
        #- creates a lambda fuction to facilitate easy arbitrary depth dict-tree creation
    #RETURNS:
        #- lambda fuction
    #"""
    #dt = lambda: defaultdict(dt)
    #return dt



def motif_profiles_weighted_by_score(processed_moods_result_dict):
    """
    *GIVEN:*
        * processed_moods_result_dict: output from (def ``process_MOODS_results()``)
    *DOES:*
        * iterates through ``processed_moods_result_dict`` and calculates a motif presence score (mps)
          for each motifName:seqName pair by summing the all positive ``site_scores`` for motifName in seqName.
        * mps are recorded in a ``pandas.DataFrame`` with columns = motifName and indexs = seqName.
    *RETURNS:*
        * mps_table: ``pandas.DataFrame`` constructed as above.
    """
    # TODO: test output for correct columns/indexes
    dict_tree = lambda: defaultdict(dict_tree)
    mps_dict = dict_tree()
    
    for seqName,motifNames in processed_moods_result_dict.iteritems():
        for m_name in motifNames.keys():
            scores = processed_moods_result_dict[seqName][m_name].itervalues()
            summed_scores = sum([x for x in scores if x > 0])
            mps_dict[seqName][m_name] = summed_scores
    
    mps_table = pandas.DataFrame.from_dict(mps_dict)
    return mps_table
    
    

def process_MOODS_results(moods_result_dict,motif_names):
    """
    *GIVEN:*
        * ``moods_result_dict``: orderedDict ???is this true??? of moods_results_tuples keyed by seqName/geneName
        * ``motif_names``: correctly ordered motif names (Motifs.motifs.keys())
    *DOES:*
        * converts ``moods_result_dict`` into a three key'd multi-level dict as follows:
          ``processed_moods_result_dict[SeqName][motifName][location] = score``
    *RETURNS:*
        * ``processed_moods_result_dict``
    """
    dict_tree = lambda: defaultdict(dict_tree)
    processed_moods_result_dict = dict_tree()
    
    for seq,motif_data in moods_result_dict.iteritems():
        for i,m_name in enumerate(motif_names):
            for site_loc,site_score in motif_data[i]:
                processed_moods_result_dict[seq][m_name][site_loc] = site_score
    
    return processed_moods_result_dict
    

def save_MOODS_result(moods_hits,out_path):
    """
    *GIVEN:*
        * ``moods_hits``: non-processed result from ``Motif.scan_seq()`` or ``Motif.scan_seqDict()``.
        * ``out_path``: path to store results
    *DOES:*
        * stores ``moods_hits`` as binary pickle to ``out_path``.
    *RETURNS:*
        * ``None``
    """
    
    out_file = open(out_path,'wb')
    
    cPickle.dump(moods_hits,out_file,protocol=-1)
    
def load_MOODS_result(in_path):
    """
    *GIVEN:*
        * ``in_path``: path to store results
    *DOES:*
        * loads ``processed_moods_result_dict`` from binary pickle to ``in_path``.
    *RETURNS:*
        * un-pickled MOODS result object.
    """
    
    in_path = open(in_path,'rb')
    return cPickle.load(in_path)

class Motifs(object):
    def __init__(self):
        """
        Class to manage loading and searching a specific motif set using MOODS.
        """
        self._settings = {}
    
    def __len__(self):
        """
        *RETURNS:*
            * len(``self.motifs``)
        """
        try:
            return len(self.motifs)
        except AttributeError:
            return 0
        
    def load_JASPAR_motifs(self,pwm_file):
        """
        *GIVEN:*
            * ``pwm_file``: path to JASPAR formated PWMs.
        *DOES:*
            * Parses PWMs and stores them in ``self``.
        *RETURNS:*
            * ``None``
        """
        motifs = ParseJasparMatrixOnly(pwm_file)
        self.motifs = motifs.to_dict()
        
    def settings(self):
        """
        *RETURNS:*
            * dict of current settings for ``Motifs`` instance.
        """
        return self._settings.copy()
        
    def set_threshold(self,thresh=0.001):
        """
        *GIVEN:*
            * ``thresh``: threshold cut-off for MOODS to report a 'hit' location.
        *DOES:*
            * Sets and stores the current scanning sensitivity.
        *RETURNS:*
            * ``None``
        """
        self._settings['thresh'] = thresh
        
    def set_both_strands(self,both=False):
        """
        *GIVEN:*
            * ``both``: True or False.
        *DOES:*
            * Sets and stores whether to search both strands of each sequence.
        *RETURNS:*
            * ``None``
        """
        self._settings['both_strands'] = both    

    def scan_seq(self,sequence):
        """
        *GIVEN:*
            * ``sequence``: a single sequence to be scanned.
        *DOES:*
            * Scans seq for each motif in ``self``.
            * Reports the location AND score for each 'hit' in the sequence.
        *RETURNS:*
            * ``moods_hits_for_seq``: nested tuples
        """
        # set up recorded options
        motifs = self.motifs.values()
        thresh = self._settings['thresh']
        both_strands = self._settings['both_strands']       
        
        # search
        moods_hits_for_seq = MOODS.search(sequence=sequence, matrices=motifs, thresholds=thresh, both_strands=both_strands)
        
        # when both_strands=True, MOODS seems to output one empty list for every motif. Remove these
        if both_strands == True:
            moods_hits_for_seq = moods_hits_for_seq[:len(motifs)]
            
        
        return tuple([tuple(x) for x in moods_hits_for_seq])
        
    def scan_seqDict(self,seq_dict):
        """
        *GIVEN:*
            * ``seq_dict``: dict of sequences (key=rec_name,value=rec_sequence).
        *DOES:*
            * Scans each seq in ``seq_dict`` for each motif in ``self``.
            * Reports the location AND score for each 'hit' in each sequence.
        *RETURNS:*
            * ``moods_hits_for_seqDic``t: dict of nested tuples
        """       

        moods_hits_for_seqDict = defaultdict(tuple)

        for name,seq in seq_dict.iteritems():
            moods_hits_for_seq = self.scan_seq(sequence=seq)
            moods_hits_for_seqDict[name] = moods_hits_for_seq
            
        return moods_hits_for_seqDict