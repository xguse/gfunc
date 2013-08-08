"""
####################
gfunc_sort.py
####################
Script supporting the querying of the graph database and gene-set sorting based on target profiles.
"""

import argparse
import cPickle
from collections import defaultdict

import yaml

import numpy as np
np.seterr(all='raise')

from gfunc.data_classes import Bunch,bunchify
from gfunc.analysis_classes import RelationsHandler
from gfunc.analysis_classes import VoteHandler
from gfunc.analysis_classes import ExpressionSimilarity,TFBSSimilarity,BranchLength
from gfunc.graphTools import GraphBuilder
from gfunc.graphTools import GraphHandler
from gfunc.dev.devel import is_none_or_nan
##from gfunc.parsers.Cufflinks import CDiffFpkmTrackerParser
##from gfunc.parsers.ETE import PhyloXMLParser
##from gfunc.parsers.JASPAR import BasicTFBSParser
from gfunc.maths import bayesian_score

def gather_metric_stats(node_list):
    """
    """
    metric_scores = defaultdict(list)
    metric_votes  = defaultdict(list)
    
    for node in node_list:
        for metric in node.poll_results:
            metric_scores[metric].append(node.poll_results[metric])
        for metric in node.voters_per_metric:
            metric_votes[metric].append(len(node.voters_per_metric[metric]))
            
    return metric_scores,metric_votes


    
def calculate_straight_combo_score(node,gHandler):
    """
    """
    sub_scores = node.get_sub_scores(gHandler.target_node,gHandler.graph)
    
    cleaned_sub_scores = []
    

    for sub_score in sub_scores:
        try:
            if not is_none_or_nan(sub_score):
                cleaned_sub_scores.append(sub_score)
            else:
                pass
        except FloatingPointError:
            pass
        
    return np.mean(cleaned_sub_scores)
    
def record_combo_scores(node_list,gHandler):
    """
    """
    comboScores = []
    
    for node in node_list:
        cScore = calculate_straight_combo_score(node,gHandler)
        node.combo_score = cScore
        comboScores.append(cScore)
    gHandler.combo_scores = comboScores

def get_median_votes(node_list):
    """
    for node where there is at least one vote.
    """
    vote_totals = []
    for node in node_list:
        node_votes = node.total_votes()
        if node_votes > 0:
            vote_totals.append(node_votes)
    return np.median(vote_totals)
            

def record_bayesian_combo_score(node_list,median_votes,gHandler):
    """
    """
    scores_with_votes = []
    for node in node_list:
        if node.total_votes() > 0:
            scores_with_votes.append(node.combo_score)
            
    median_score = np.median(scores_with_votes)
    
    target_node = gHandler.target_node
    graph = gHandler.graph
    for node in node_list:
        c=median_votes
        m=median_score
        n=node.total_votes()
        scores=node.get_sub_scores(target_node,graph)
        node.b_score = bayesian_score(c,m,n,scores,scale_mod=1)
    
def write_sorted_table(node_list,table_path):
    """
    """
    node_list.sort(key=lambda node: node.b_score)
    node_list.reverse()
    b_scores = [node.b_score for node in node_list]
    std_dev = np.std(b_scores)
    out_file = open(table_path,'w')
    header = "%s\n" % ('\t'.join(['node','rank','b_score','naive_score','votes','std_dev_from_median_b_score']))
    out_file.write(header)
    for i,node in enumerate(node_list):
        out_file.write('%s\n' % ('\t'.join([node.name,str(i+1),str(node.b_score),str(node.combo_score),str(node.total_votes()),str(node.b_score/std_dev)])))
    out_file.close()
    
def main():
    """
    The main loop.  Lets ROCK!
    """
    
    desc = """TODO: ... ask me later! I'm on a deadline! ..."""
    
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument('config_file', type=str,
                        help="""Path to a yaml formatted config file containing setup options for the graph.""")
    
    args = parser.parse_args()
    
    yopts = bunchify(yaml.load(open(args.config_file,'rU')))
    
    gHandler = cPickle.load(open(yopts.outputs.graph_pickle,'rb'))
    gHandler.clone_node_as_target(yopts.target)
    gHandler.install_target()
    
    node_list = [node for node in gHandler.node_dict.itervalues() if node.species == yopts.genes_to_sort]
    gHandler.take_votes(node_list)
    
    metric_scores,metric_votes = gather_metric_stats(node_list)
    
    record_combo_scores(node_list,gHandler)
    
    median_votes = get_median_votes(node_list)
    record_bayesian_combo_score(node_list,median_votes,gHandler)
    
    cPickle.dump(gHandler,open(yopts.outputs.polled_pickle,'w'))
    write_sorted_table(node_list,table_path=yopts.outputs.sorted_table)

if __name__ == '__main__':
    main()
    print "main() completed."