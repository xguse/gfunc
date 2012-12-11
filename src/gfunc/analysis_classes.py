"""
####################
analysis_classes.py
####################
Code defining controller classes for supported analysis types.
"""
import pdb

import numpy as np
#np.seterr(divide='raise')

from scipy import stats as sp_stats

from gfunc.data_classes import Bunch
from gfunc.maths import bayesian_score
from gfunc.maths import weight_d_for_ptci

######################
# Metrics Handlers
######################
class MetricHandler(object):
    pass

class RelationsHandler(MetricHandler):
    """
    Relations are metrics that characterize how an
    edge's 2 connected nodes relate given some type
    of relationship (branch length, expression profile
    similarity, etc).
    """
    def __init__(self,list_of_metrics):
        """
        TODO: Doc
        """
        self.metrics = list_of_metrics
        self._setup_metric_dict()
        
        
    def _setup_metric_dict(self):
        """
        TODO: Doc
        """
        met_dict = {}
        for met in self.metrics:
            met_dict[met.relation_metric] = met
        self.metrics = met_dict
        
    def get_vote_types(self):
        """
        TODO: Doc
        """
        take_votes_for = []
        for metric,obj in self.metrics.iteritems():
            if obj.poll_me:
                take_votes_for.append(metric)
        return take_votes_for
    
    def measure_relations(self,edge_dict):
        """
        For each gfunc_edge in dict:
           iterates through metrics:
              calculates & stores result in gfunc_edge and metric_handler
        """
        for gfunc_edge in edge_dict.itervalues():
            for metric in self.metrics.itervalues():
                metric.measure_relation(gfunc_edge)
    


class VoteHandler(MetricHandler):
    """
    VoteHandlers determine how well the current GFuncNode's
    neighborhood agrees with it regarding each type of relation
    being used.  The result is a weighted mean, weighted by
    a strength (or trustability metric) between each neighbor
    and the current node (exp: branch length, p-value, etc).
    If no weight relationship is specified, the result is a
    standard mean (weights are equal).
    """
    def __init__(self,graph):
        """
        TODO: Doc
        """
        self._graph = graph
        self.target_node = None
        
    def _weighted_mean(self,scores_and_weights):
        """
        TODO: Doc
        """
        try:
            a = np.array(scores_and_weights)
            scores = a[:,0]
            weights = a[:,1]+1
            #weights = weights/weights.sum() # converts to proportions (weights now sum to 1)
            return (scores*weights).sum()   # calculates weighted mean with proportional weights
        except:
            if len(a) == 0:
                return float('nan')
            else:
                raise
        
    def set_vote_types(self,vote_types,weight_by=None):
        """
        Recieves and stores as LIST 'Metric.relation_metric' string for
        each Metric class that needs a vote taken.
        
        weight_by: one or none of the gfunc_edge.data key strings to use to weight
        the votes of each node's neighbor edge metrics.
        (weight_by=None results in equal weights)
        """
        self.vote_types = vote_types
        self.weight_by = weight_by
    
    def take_votes(self, node_list, poll_func=None):
        """
        TODO: Doc
        """
        if type(node_list) is not type([]):
            raise TypeError("node_list must be type([]).")
        
        graph = self._graph
        
        # Idea is to easily allow custom polling strategies
        if poll_func is not None:
            poll_func(node_list)
        
        else:
            for node in node_list:
                # set up which metrics will be polled for each neighbor edge
                node_polls = Bunch()
                for vote_type in self.vote_types:
                    node_polls[vote_type] = []

                
                # poll neighbors for what they think about the target w.r.t each metric in node_polls
                # but weight it by that neighbors relation to node
                for neighbor in graph[node]:
                    if neighbor.name == 'target':
                        continue
                    neighbor_edge = graph[node][neighbor]['edge']
                    neighbor2target = graph[neighbor][self.target_node]['edge']

                    
                    if self.weight_by is not None:
                        weight = neighbor_edge.data[self.weight_by]
                    else:
                        weight = 1
                    for metric in node_polls:
                        # Check to make sure that the metric is a useful number (not NaN)
                        vote_value = neighbor2target.data[metric]
                        if not np.isnan(vote_value):
                            node_polls[metric].append([vote_value,weight])
                            node.voters_per_metric[metric].append(neighbor)
                        
                
                # convert results lists to weighted mean for each metric type
                for metric in node_polls:
                    node_polls[metric] = self._weighted_mean(node_polls[metric])
                node.poll_results = node_polls
                

######################
# Core Metrics 
######################

class Metric(object):
    """
    TODO: doc
    """
    def __init__(self,poll_me=False):
        """
        TODO: doc
        """
        raise NotImplementedError('You must override this method in subclass.')
    
    def _bootstrap_the_distribution(self):
        """
        Estimate some qualities of the distribution of encountered values.  
        """
        raise NotImplementedError('TODO: You still need to write this.')    
    
    def _calc_metric(self,node1,node2):
        """
        Does this metric's specific calculations.
        """
        raise NotImplementedError('You must override this method in subclass.')
        
    def measure_relation(self,gfunc_edge):
        """
        TODO: doc
        """
        data=self._calc_metric(gfunc_edge)
        data_type=self.relation_metric
        
        gfunc_edge.set_data(data,data_type)
        self.recorded_values.append(data)
    
    def mean(self,greater_than=0):
        """
        Returns the mean of the encountered values greater than the value provided.  
        """
        return np.mean([x for x in self.recorded_values if x > greater_than])

    def median(self,greater_than=0):
        """
        Returns the median of the encountered values greater than the value provided.  
        """
        return np.median([x for x in self.recorded_values if x > greater_than])
    
    
class PhyloExpnCorrelationIndex(Metric):
    """
    TODO: doc
    """
    def __init__(self,poll_me=False):
        """
        TODO: doc
        """
        self.poll_me = poll_me
        self.relation_metric = 'PTCI'
        self.recorded_values = []
        
    def _calc_metric(self,gfunc_edge):
        """
        Does this metric's specific calculations.
        """
        node1,node2 = gfunc_edge.nodes
        
        try:
            r_val,p_val = sp_stats.pearsonr(node1.data.expression_vector, node2.data.expression_vector)
            d_val,d_min,d_max = gfunc_edge.data.divergence
            
            ptci = r_val * (1-p_val) * weight_d_for_ptci(d_val,d_min,d_max)
            
            return ptci
        except AttributeError as err:
            if """'Bunch' object has no attribute""" in err.message:
                return float('nan')
            else:
                raise err

class ExpressionSimilarity(Metric):
    """
    TODO: doc
    """
    def __init__(self,poll_me=False):
        """
        TODO: doc
        """
        self.poll_me = poll_me
        self.relation_metric = 'expression_vector_similarity'
        self.recorded_values = []

    def _calc_metric(self,gfunc_edge):
        """
        Does this metric's specific calculations.
        """
        node1,node2 = gfunc_edge.nodes
        try:
            r_val,p_val = sp_stats.pearsonr(node1.data.expression_vector, node2.data.expression_vector)
            # for now, all r_vals will be used
            # TODO: figure out whether that needs revising...
            ##if p_val > 0.1:
                ##return float('nan')
            
            # return r_val scaled to between 0 and 1
            #scaled_rVal = (r_val+1)/2
            return r_val     
        except AttributeError as err:
            if """'Bunch' object has no attribute""" in err.message:
                # TODO: Should I return these or just leave the value unset?
                # for now its left unset.
                r_val = float('nan')
                p_val = float('nan')
                return r_val
            else:
                raise err
        
        
class TFBSSimilarity(Metric):
    """
    TODO: doc
    """
    def __init__(self,poll_me=False):
        """
        """
        self.poll_me = poll_me
        self.relation_metric = 'tfbs_vector_similarity'
        self.recorded_values = []
        
    def _calc_metric(self,gfunc_edge):
        """
        Does this metric's specific calculations.
        """
        node1,node2 = gfunc_edge.nodes
        try:
            r_val,p_val = sp_stats.pearsonr(node1.data.tfbs_vector, node2.data.tfbs_vector)
            # for now, all r_vals will be used
            # TODO: figure out whether that needs revising...
            ##if p_val > 0.1:
                ##return float('nan')
            
            # return r_val scaled to between 0 and 1
            scaled_rVal = (r_val+1)/2
            return r_val
        except AttributeError as err:
            if """'Bunch' object has no attribute""" in err.message:
                # TODO: Should I return these or just leave the value unset?
                r_val = float('nan')
                p_val = float('nan')
                return r_val
            else:
                raise err
        



class BranchLength(Metric):
    """
    TODO: doc
    """
    def __init__(self,poll_me=False):
        """
        TODO: doc
        """
        self.poll_me = poll_me
        self.relation_metric = 'branch_length'
        self.recorded_values = []

    def measure_relation(self,gfunc_edge):
        """
        TODO: doc
        """
        # I am overriding this method BC this edge data was set when the edge was registered
        # We just need to record the value for the branch_length distribution
        data = gfunc_edge.data.branch_length
        self.recorded_values.append(data)
    


