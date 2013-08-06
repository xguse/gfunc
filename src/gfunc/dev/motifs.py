import cPickle

from collections import defaultdict

import numpy as np
from scipy import stats
import pandas


import MOODS

from gfunc.parsers.JASPAR import ParseJasparMatrixOnly
from gfunc.fileIO import tableFile2namedTuple
from gfunc.motifs import process_MOODS_results

from spartan.utils.files import ParseFastA


def get_score_dataframe(processed_MOODS_results):
    """
    :returns: pandas.DataFrame obj with seq_names as row headings and
    motif_names as column headings; values are sum of the positive log-odds
    scores for a particular motif_name in that seq_name.
    """
    motif_names = sorted(processed_MOODS_results.values()[0].keys())
    
    dict_of_series_objects = {}  
    
    for seq_name in processed_MOODS_results.keys():
        # intialize pandas.Series(data, index) object and store in dict_of_series_objects
        data = []
        for motif_name in motif_names:
            # for each motif sum all positive scores (log-odds)
            # negative log-odds mean that the site has less than even odds of being the motif given the background
            pos_log_odds = [ x for x in processed_MOODS_results[seq_name][motif_name].values() if x > 0 ]
            data.append(np.sum(pos_log_odds))
            
        dict_of_series_objects[seq_name] = pandas.Series(data=data, index=motif_names)
        
    return pandas.DataFrame.from_dict(dict_of_series_objects,orient='index')

def get_standardized_score_dataframe(score_dataframe,center='mean'):
    """
    :returns: standardized version (by column) of input score_dataframe.
    """
    df = score_dataframe
    
    if center == 'mean':
        return (df - df.mean())/df.std()
    elif center == 'median':
        return (df - df.median())/df.std()
    else:
        raise ValueError('"center" must be either ["median","mean"].  You provided: %s.' % (center))

def get_score_dataframe(processed_MOODS_results):
    """
    """
    motif_names = sorted(processed_MOODS_results.values()[0].keys())
    
    dict_of_series_objects = {}  
    
    for seq_name in sorted(processed_MOODS_results.keys()):
        # intialize pandas.Series(data, index) object and store in dict_of_series_objects
        data = []
        for motif_name in motif_names:
            # for each motif sum all positive scores (log-odds)
            # negative log-odds mean that the site has less than even odds of being the motif given the background
            pos_log_odds = [ x for x in processed_MOODS_results[seq_name][motif_name].values() if x > 0 ]
            data.append(np.sum(pos_log_odds))
            
        dict_of_series_objects[seq_name] = pandas.Series(data=data, index=motif_names)
        
    return pandas.DataFrame.from_dict(dict_of_series_objects,orient='index')

def save_dataframe_to_motif_table(dataframe, species, out_path):
    """
    Headers:
    seq_name,species,motif_name1,motif_name2,...,motif_nameN
    """
    out = open(out_path,'w')
    
    motif_names = sorted(dataframe.columns)
    seq_names = sorted(dataframe.axes[0])
    
    headers = 'seq_name\tspecies\t%s\n' % ('\t'.join(motif_names))
    out.write( headers )
    
    for seq in seq_names:
        out.write( '%s\t%s\t%s\n' % (seq,species,'\t'.join([ str(x) for x in dataframe.ix[seq] ])) )
        
        
def get_and_save_motif_hits(motif_obj,fasta_path,species,out_path):
    """
    """
    
    p = ParseFastA(fasta_path)
    fasta_dict = p.to_dict()    

    # motif_obj has already be initialized with motifs, the p-val threshold, and both_starnds settings
    hit_dict = motif_obj.scan_seqDict(fasta_dict)
    
    processed_MOODS_results = process_MOODS_results(hit_dict,motif_obj.motifs.keys())
    df = get_score_dataframe(processed_MOODS_results)
    
    # standardize the log-odds scores
    df_norm_median = (df - df.median())/df.std()
    
    # write results to out_path as TSV
    save_dataframe_to_motif_table(dataframe=df_norm_median, species=species, out_path=out_path)