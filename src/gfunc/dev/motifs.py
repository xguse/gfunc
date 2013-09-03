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


def get_motif_names_from_meme_file(meme_file_path):
    """
    """
    motif_names = []
    
    for line in open(meme_file_path,'rU'):
        if line.startswith('MOTIF'):
            motif_names.append(line.strip().split()[1])
    
    return motif_names
            

def get_seq_names_from_steme(steme_pwm_scan_seqs):
    """
    """
    pos_name_dict = {}
    index = 0
    for line in open(steme_pwm_scan_seqs,'rU'):
        line = line.strip().split(',')
        pos_name_dict[int(index)] = line[1]
        index += 1
    
    return pos_name_dict

def process_STEME_SCAN_results(steme_pwm_scan_out,steme_pwm_scan_seqs,motif_file_meme):
    """
    *GIVEN:*
        * ``steme_pwm_scan_out``: csv file from steme-pwm-scan (steme-pwm-scan.out)
        * ``steme_pwm_scan_seqs``: csv file from steme-pwm-scan (steme-pwm-scan.seqs)
        * ``motif_file_meme``: meme formated PWM file used to generate ``steme-pwm-scan.out``
    *DOES:*
        * converts ``steme_pwm_scan_out`` into a three key'd multi-level dict as follows:
          ``process_STEME_SCAN_results_dict[SeqName][motifName][location] = score``
    *RETURNS:*
        * ``process_STEME_SCAN_results_dict``
    """
    dict_tree = lambda: defaultdict(dict_tree)
    processed_steme_result_dict = dict_tree()
    
    
    motif_names = get_motif_names_from_meme_file(motif_file_meme)
    seq_pos_name = get_seq_names_from_steme(steme_pwm_scan_seqs)
    
    # initialize process_STEME_SCAN_results_dict with seq_names and motif_names
    for seq_index,seq_name in seq_pos_name.iteritems():
        for motif in motif_names:
            processed_steme_result_dict[seq_name][motif]
    
    for line in open(steme_pwm_scan_out,'rU'):
        try:
            fields = line.strip().split(',')
            m_name = fields[0]
            s_name = seq_pos_name[int(fields[2])]
            location = int(fields[3])
            score  = float(fields[5])
            if score > 0.5:
                processed_steme_result_dict[s_name][m_name][location] = score
        except KeyError:
            print fields
            raise
    return processed_steme_result_dict



def get_score_dataframe(processed_MOTIF_SEARCH_results,search_kind='MOODS'):
    """
    :returns: pandas.DataFrame obj with seq_names as row headings and
    motif_names as column headings; values are sum of the filtered_scores
    for a particular motif_name in that seq_name.
    """
    score_threshold_dict = {'MOODS':0,
                            'STEME_SCAN':0.5}
    motif_names = sorted(processed_MOTIF_SEARCH_results.values()[0].keys())
    
    dict_of_series_objects = {}  
    
    for seq_name in sorted(processed_MOTIF_SEARCH_results.keys()):
        # intialize pandas.Series(data, index) object and store in dict_of_series_objects
        motif_data = []
        for motif_name in motif_names:
            # for each motif sum all scores over threshold (MOODS = log-odds;0  STEME_SCAN = Z;0.5)
            # negative log-odds mean that the site has less than even odds of being the motif given the background
            filtered_scores = [ x for x in processed_MOTIF_SEARCH_results[seq_name][motif_name].values() if x > score_threshold_dict[search_kind] ]
            motif_data.append(np.sum(filtered_scores))
            
        dict_of_series_objects[seq_name] = pandas.Series(data=motif_data, index=motif_names)
        
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
    
def norm_and_save_STEME_SCAN(steme_pwm_scan_out,steme_pwm_scan_seqs,motif_file_meme,species,out_path):
    """
    """
    processed_STEME_SCAN_results = process_STEME_SCAN_results(steme_pwm_scan_out,steme_pwm_scan_seqs,motif_file_meme)
    df = get_score_dataframe(processed_STEME_SCAN_results,search_kind='STEME_SCAN')
    
    # standardize the scores
    df_norm_median = (df - df.median())/df.std()
    
    # write results to out_path as TSV
    save_dataframe_to_motif_table(dataframe=df_norm_median.fillna(0), species=species, out_path=out_path)    
    
    
#for sn,ss in cq_seqs.iteritems():
    #try:
        #start = ss.find(seqs.revcomp('CCCAAAATAG'))
        #if start != -1:
            #print '%s: %s' % (sn,start) 
    #except IndexError:
        #print '%s: length %s' % (sn,len(ss))