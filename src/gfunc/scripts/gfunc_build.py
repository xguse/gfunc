"""
####################
gfunc_build.py
####################
Script supporting the construction and saving of a new gfunc graph database from user options.
"""
import pdb
import argparse
import cPickle

import yaml


from gfunc.data_classes import Bunch,bunchify
from gfunc.analysis_classes import RelationsHandler
from gfunc.analysis_classes import VoteHandler
from gfunc.analysis_classes import ExpressionSimilarity,TFBSSimilarity,BranchLength
from gfunc.graphTools import GraphBuilder
from gfunc.graphTools import GraphHandler
from gfunc.parsers.Cufflinks import CDiffFpkmTrackerParser
from gfunc.parsers.ETE import PhyloXMLParser
from gfunc.parsers.JASPAR import BasicTFBSParser

supported_data_types = {'cuffdiff_fpkm_profile':CDiffFpkmTrackerParser,
                        'phyloXML':PhyloXMLParser,
                        'tfbs_profile':BasicTFBSParser}




def confirm_parsers():
    """
    TODO: Make sure provided data types are supported in supported_data_types.
    """
    
def setup_parsers(yopts):
    """
    TODO: make this smarter when you have time.
    """
    y = yopts
    parser_list = []
    
    parser_list.append(CDiffFpkmTrackerParser(cuffdiff_fpkm_path=y.species_info.Anopheles_gambiae.expresion_data,
                                              species=y.species_info.Anopheles_gambiae.name))
    parser_list.append(CDiffFpkmTrackerParser(cuffdiff_fpkm_path=y.species_info.Culex_quinquefasciatus.expresion_data,
                                              species=y.species_info.Culex_quinquefasciatus.name))
    
    parser_list.append(BasicTFBSParser(tfbs_path=y.species_info.Anopheles_gambiae.tfbs_data))
    parser_list.append(BasicTFBSParser(tfbs_path=y.species_info.Culex_quinquefasciatus.tfbs_data))
    
    species = [y.species_info.Anopheles_gambiae.name,y.species_info.Culex_quinquefasciatus.name]
    parser_list.append(PhyloXMLParser(phyloXML_path=y.edge_data.branch_length_data,species=species))
    
    return parser_list

def init_metrics():
    """
    TODO: Make this 'smart' when you have more time. 
    """
    expression = ExpressionSimilarity(poll_me=True)
    tfbs       = TFBSSimilarity(poll_me=True)
    branch     = BranchLength(poll_me=False)
    
    return expression,tfbs,branch
    

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
    
    parser_list = setup_parsers(yopts)
    
    gBuilder = GraphBuilder(parser_list)
    gBuilder.populate_registries()
    gHandler = gBuilder.map_registries_to_graph()

    # init metric objs, RelHandler and Vote Handler
    metrics     = init_metrics()
    relHandler  = RelationsHandler(metrics)
    voteHandler = VoteHandler(gHandler.graph)
    voteHandler.set_vote_types(relHandler.get_vote_types(),weight_by='branch_length')
    
    gHandler.install_metric_handlers(rel_hndler=relHandler,vote_hndlr=voteHandler)
    gHandler.measure_relations()
    
    # dump built construct as pickle
    cPickle.dump(gHandler,open(yopts.outputs.graph_pickle,'w'))

if __name__ == '__main__':
    main()
    print "main() completed."