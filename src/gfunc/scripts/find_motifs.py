"""
####################
find_motifs.py
####################
Code supporting a VERY simple motif analysis of multiple promoters.  This script depends on the MOODS
library.
"""
import argparse
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

import numpy as np
from scipy import stats

import MOODS

from gfunc.parsers.JASPAR import ParseJasparMatrixOnly
from gfunc.fileIO import tableFile2namedTuple

from rSeq.utils.files import ParseFastA

def write_normalized_table(headers,norm_path,norm_dict):
    """
    """
    axis0 = len(norm_dict)
    axis1 = len(norm_dict.values()[0])
    data = np.zeros([axis0,axis1])
    for i,counts in enumerate(norm_dict.itervalues()):
        data[i,:] = counts
    
    n_data = upper_qrtl_norm(data)
    
    data = n_data
    
    norm_file = open(norm_path,'w')
    norm_file.write(headers)
    labels = norm_dict.iterkeys()
    for row in data:
        norm_file.write('%s\t%s\n' % (labels.next(),'\t'.join([str(round(n,4)) for n in list(row)])))
    norm_file.close()
    

    
def upper_qrtl_norm(data):
    """
    """
    ma = np.ma
    # get array of 25th percentile limits
    lower_quartiles = stats.scoreatpercentile(data,25)
    
    # mask any gene in each column of counts that has less counts than the lower_quartile
    m_data = ma.MaskedArray(data)
    m_data = ma.masked_less(m_data,lower_quartiles)
    
    # get mean and stdDev for each column with lower qrtl masked out
    means = np.array([m_data[:,i].mean() for i in range(len(data[0,:]))])
    stdDevs = np.array([m_data[:,i].std() for i in range(len(data[0,:]))])
    
    return (data-means)/stdDevs
    
 
def load_to_normalize(table_path):
    """
    """
    table = tableFile2namedTuple(table_path,sep='\t',headers=None)
    headers = '%s\n' % ('\t'.join(table[0]._fields))
    norm_dict = OrderedDict()
    for row in table:
        norm_dict['%s\t%s' % (row.seq_name,row.species)] = np.array([float(count) for count in row[2:]])
    return norm_dict,headers
    
    
def main():
    """
    The main loop.  Lets ROCK!
    """
    
    desc = """... ask me later! I'm on a deadline! ..."""
    
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument('--seqs', type=str,
                        help="""Path to a fasta file containing the 'promoter'
                        regions of a single species you wish to scan with motifs.""")
    
    parser.add_argument('--species', type=str,
                        help="""A quoted string of the species name: 'Anophles gambiae'.""")
    
    parser.add_argument('--motifs', type=str,
                        help="""Path to a file containing the motifs you wish to use.  
                        The file must be in JASPAR's 'matrix_only.txt' format.""")
    
    parser.add_argument('--thresh', type=float, required=False, default=0.001,
                        help="""A p-val cut-off above which hits will be ignored. (default = %(default)s)""")
    
    parser.add_argument('--out', type=str, required=False, default='compare_motifs.out',
                        help="""Path to outfile. (default = %(default)s)""")
    
    parser.add_argument('--norm', type=str, required=False, default=False,
                        help="""Optional path to outfile for data normalized by upper quartiles w.r.t. each motif. (default = %(default)s)""")

    parser.add_argument('--to-norm', type=str, required=False, default=False,
                            help="""Optional path to outfile of previous run that needs to be normalized. (default = %(default)s)""")
    

    
    args = parser.parse_args()
    
    
    if not args.to_norm:
        # create parsers
        motifs = ParseJasparMatrixOnly(args.motifs)
        seqs   = ParseFastA(args.seqs)
    
        # Load all motifs at once
        # We will be loading one seq at a time.
        motifs = motifs.to_dict()
    
        # set up output and headers
        headers  = 'seq_name\tspecies\t%s\n' % ('\t'.join(motifs.keys()))
        out_file = open(args.out,'w')
        out_file.write(headers)
    
    # lets start the major looping
    if args.norm and not args.to_norm:
        norm_dict = OrderedDict()
    elif args.to_norm:
        norm_dict,headers = load_to_normalize(args.to_norm)
        write_normalized_table(headers,args.norm,norm_dict)
        exit(0)
        
    for name,seq in seqs:
        hits = MOODS.search(seq,motifs.values(),args.thresh)
        counts = [len(x) for x in hits]
        out_file.write('%s\t%s\t%s\n' % (name,args.species,'\t'.join([str(x) for x in counts])))
        if args.norm:
            norm_dict['%s\t%s' % (name,args.species)] = np.array(counts)
    out_file.close()
    
    if args.norm:
        write_normalized_table(headers,args.norm,norm_dict)
    
    


if __name__ == "__main__":
    main()