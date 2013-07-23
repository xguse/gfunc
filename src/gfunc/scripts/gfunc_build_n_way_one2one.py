"""
####################
gfunc_build.py
####################
Script supporting the construction and saving of a new gfunc graph database from user options.
"""
import pdb
import argparse
import cPickle
from collections import defaultdict

from numpy import isnan
import numpy as np

from matplotlib import pylab as plb
import matplotlib.lines as mlines

import pandas

import yaml

from gfunc.stats import cumHypergeoP
from gfunc.stats import benjHochFDR
from gfunc.xpermutations import xuniqueCombinations
from gfunc.fileIO import tableFile2namedTuple
from gfunc.data_classes import Bunch,bunchify
from gfunc.analysis_classes import RelationsHandler
from gfunc.analysis_classes import VoteHandler
from gfunc.analysis_classes import PhyloExpnCorrelationIndex
from gfunc.graphTools import GraphBuilder
from gfunc.graphTools import GraphHandler
from gfunc.parsers.Cufflinks import CDiffFpkmTrackerParser
from gfunc.parsers.ETE import PhyloXMLParser
from gfunc.parsers.JASPAR import BasicTFBSParser
from gfunc.parsers.edge_lists import OneToOneOrthoListParser

def print_LoL_to_file(list_of_lists,path):
    """
    """
    out = open(path,'w')
    for row in list_of_lists:
        out.write("%s\n" % (','.join([str(x) for x in row])))
    out.close()

def count_terms_for_hypergeo(graph,positive_node_set,terms_of_enrichment,term_fuction=None):
    """
    TAKES:
        * graph: gfunc graph obj
        * positive_node_set: list of GFuncNodes that are the "positive" group
        * terms_of_enrichment: set of terms present in at least one of the positive nodes
        * term_fuction: function that takes a GFuncNode as input and returns the SET of terms
          of a specific type associated with the node
        
    *DOES:*
        * uses "positive_node_set" and "term_fuction" to count positives and negatives and
          set up the required inputs for the hypergeometric calculations
        
        
    *YIELDS:*
        * tuple: (term,n,i,m,N)
          term = term_id
          n = # of positives in population
          i = # of positives in sample
          m = # of negatives in population
          N = sample size

    """
    if term_fuction is None:
            raise ValueError("You MUST provide your own value for 'term_fuction'.")    

    for term in terms_of_enrichment:
        term = set([term])
        n = 0
        i = 0
        m = 0
        
        for node in graph.nodes_iter():
            if term.issubset(term_fuction(node)):
                n += 1
                if node in positive_node_set:
                    i += 1
            else:
                m += 1
        
        N = len(positive_node_set)        
        
        yield (term,n,i,m,N)
    
def enrichment_of_terms(graph,positive_node_set,term_fuction=None):
    """
    *GIVEN:*
        * graph: gfunc graph obj
        * positive_node_set: set of GFuncNodes that are the "positive" group
        * term_fuction: function that takes a GFuncNode as input and returns the SET of terms of a specific type associated with the node.
        
    *DOES:*
        * uses "positive_node_set" and "term_fuction" to compile a set of terms present in at least one of the positive nodes.
        * Calculates the hypergeometric enrichment p-value of each term.
        * Calculates the benjamini-hochberg adjusted q-values (FDR).
        
    *RETURNS:*
        * table of results as list of lists
    """
    
    if term_fuction is None:
        raise ValueError("You MUST provide your own value for 'term_fuction'.")
    
    # compile a set of terms present in at least one of the positive nodes
    terms_of_enrichment = set()
    for node in positive_node_set:
        node_terms = term_fuction(node)
        terms_of_enrichment.update(node_terms)
    
    # remove 'blank' terms
    terms_of_enrichment.discard('')
    
    
    results_table = []
    
    # Build results table before multiple hypothesis testing correction
    count_iter = count_terms_for_hypergeo(graph,positive_node_set,terms_of_enrichment,term_fuction=term_fuction)
    for term,n,i,m,N in count_iter:
        p = cumHypergeoP(n,i,m,N)
        results_table.append([term,n,i,m,N,p])
    
    # Perform multiple hypothesis testing correction
    results_table = benjHochFDR(table=results_table,pValColumn=5)
    
    return results_table

def load_prot_domains(graph,tsv_paths):
    """
    *DOES:*
        * reads tsv file generated by:
        
            ``wget -c -t 0 --timeout=60 --waitretry=60 -O 
            genomeName_eg_gene.protDomains.txt 'http://www.biomart.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE Query>
            <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.7" >
            <Dataset name = "genomeName_eg_gene" interface = "default" ><Attribute name = "ensembl_gene_id" />
            <Attribute name = "ensembl_transcript_id" /><Attribute name = "superfamily" />
            <Attribute name = "pfam" /><Attribute name = "transmembrane_domain" /><Attribute name = "signal_domain" />
            <Attribute name = "ncoils" /></Dataset></Query>'``
        
        * builds Bunch tree with
            * tier0_keys = ensembl_gene_id
            * tier1_keys = superfamily, pfam, transmembrane_domain, signal_domain, ncoils
            * tier2_vals = superfamily(set), pfam(set), transmembrane_domain(set), signal_domain(set), ncoils(set)
            
        * adds ref to the tree to each GFuncNode's "data" dict as "prot_domains"
 

    *RETURNS:*
        * prot_domains tree
    """
    # open and convert tsv_paths to namedtuple tables
    data_tables = []
    table_headers = ["ensembl_gene_id","ensembl_transcript_id","superfamily","pfam","transmembrane_domain","signal_domain","ncoils"]
    for tsv_path in tsv_paths:
        data_tables.append(tableFile2namedTuple(tsv_path,headers=table_headers))
    
    # initialize and build prot_domains tree
    prot_domains = defaultdict(lambda: defaultdict(set))
    
    for table in data_tables:
        for row in table:
            for domain_type in table_headers[2:]:
                prot_domains[row.ensembl_gene_id][domain_type].add(row.get(domain_type))
    
    # add prot_domains tree to each node
    for node in graph.nodes_iter():
        node.data["prot_domains"] = prot_domains
    
    # return prot_domains
    return prot_domains

def load_go_slims(graph,tsv_paths):
    """
    *DOES:*
        * reads tsv file generated by:
        
            ``wget -c -t 0 --timeout=60 --waitretry=60 -O 
            genomeName_eg_gene.protDomains.txt 'http://www.biomart.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE Query>
            <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.7" >
            <Dataset name = "genomeName_eg_gene" interface = "default" ><Attribute name = "ensembl_gene_id" />
            <Attribute name = "ensembl_transcript_id" /><Attribute name = "goslim_goa_accession" />
            <Attribute name = "goslim_goa_description" /></Dataset></Query>'``
            
        * builds Bunch tree with
            * tier0_keys = ensembl_gene_id
            * tier1_keys = goslim_goa_accession
            * tier2_vals = goslim_goa_description
            
        * adds ref to the tree to each GFuncNode's "data" dict as "go_slims"


    *RETURNS:*
        * go_slims tree
    """

    # open and convert tsv_paths to namedtuple tables
    data_tables = []
    table_headers = ["ensembl_gene_id","ensembl_transcript_id","goslim_goa_accession","goslim_goa_description"]
    for tsv_path in tsv_paths:
        data_tables.append(tableFile2namedTuple(tsv_path,headers=table_headers))
    
    # initialize and build prot_domains tree
    go_slim_terms = defaultdict(lambda: defaultdict(str))
    
    for table in data_tables:
        for row in table:
            go_slim_terms[row.ensembl_gene_id][row.goslim_goa_accession] = row.goslim_goa_description
    
    # add go_slim_terms tree to each node
    for node in graph.nodes_iter():
        node.data["go_slims"] = go_slim_terms
    
    # return prot_domains
    return go_slim_terms    

def plot_hists(ptci_dicts):
    """
    *DOES:*
        * generates experimental and randomized control histograms
    *RETURNS:*
        * None
    """
    # boilerplate for setting up the figure canvas
    fig = plb.figure()
    ax  = fig.add_subplot(111)
    
    # Calc and plot the histograms
    output = [ax.hist(ptci_dicts[0].values(), bins=max([len(d) for d in ptci_dicts])/50.0,
                      cumulative=0, bottom=None, histtype='step')]
    
    output.extend([ax.hist(ptci_dicts[x].values(), bins=max([len(d) for d in ptci_dicts])/50.0,
                          cumulative=0, bottom=None, histtype='step',
                          color='grey', alpha=.1) for x in range(1,len(ptci_dicts))])
    
    # Label the Fig
    ax.set_xlabel('PTCI')
    ax.set_ylabel('Gene Triplets')
    
    # set up proxy artists to make the legend work
    l1 = mlines.Line2D([0,8], [1,8], lw=5.,alpha=1, c="blue") 
    l2 = mlines.Line2D([0,8], [1,8], lw=5.,alpha=1, c="grey")
    ax.legend([l1,l2],["3-way 1:1 Orthology","Randomized 3-way Triplets"])
    
    plb.show()
    

def calc_deltaFPKM(vector,ref_index=0):
    delta_vector = [x-vector[ref_index] for x in vector]
    return delta_vector

def write_deltaFPKM_vectors(ortho_meanScore_dict,gHandler,outPath):
    gH = gHandler
    names = []
    for ortho_set in ortho_meanScore_dict.keys():
        names.extend(ortho_set)
        
    oFile = open(outPath,'w')
    
    for name in names:
        try:
            deltaFPKM = calc_deltaFPKM(gH.node_dict[name].data.expression_vector)
            oFile.write("%s\t%s\n" % (name,'\t'.join([str(x) for x in deltaFPKM])))
        except:
            print name
        
    oFile.close()
    

def get_top_ortho_sets(ortho_meanScore_dict,percentile=90):
    
    scores = ortho_meanScore_dict.values()
    pcntl_thresh = np.percentile(scores,percentile)
    
    top_orthos = {}
    for orthos,score in ortho_meanScore_dict.iteritems():
        if score >= pcntl_thresh:
            top_orthos[orthos] = score
    
    return top_orthos

#def get_ortho_set_PTCIs(ortho_meanScore_dict,percentile=90):
    
    #scores = ortho_meanScore_dict.values()
    #pcntl_thresh = np.percentile(scores,percentile)
    
    #top_orthos = {}
    #for orthos,score in ortho_meanScore_dict.iteritems():
        #if score >= pcntl_thresh:
            #top_orthos[orthos] = score
    
    #return top_orthos

def get_ortho_set_PTCIs(gHandler,gene_list):
    """
    """
    gH = gHandler
    graph = gH.graph
    
    orthologs = {}
    
    for gene in gene_list:
        ortho_names = [gene]
        ortho_scores = []
        gene_node = gH.node_dict[gene]
        for neighbor in graph[gene_node]:
            ortho_names.append(neighbor.name)
        
        for node_combo in xuniqueCombinations(ortho_names,2):
            n1 = gH.node_dict[node_combo[0]]
            n2 = gH.node_dict[node_combo[1]]
            edge = graph[n1][n2]['edge']
            if not isnan(edge.data.PTCI):
                ortho_scores.append(edge.data.PTCI)
        
        if len(ortho_scores) > 0:
            key = tuple(sorted(ortho_names))
            value = np.mean(ortho_scores)
            orthologs[key] = value
    
    return orthologs
        

def get_starting_nodes(ortho_path):
    """
    """
    rows = tableFile2namedTuple(ortho_path)
    node_names = [row[0] for row in rows]
    return node_names
        


def get_div_info(config_object):
    """
    Use data in yaml config object to create map to divergence
    times with a 2D dict tree structure:
    
    div_map['species a']['species c'] --> divergence 
    
    *RETURNS:* tuple
        div_map,min(div),max(div)
    """
    div_map = defaultdict(lambda:defaultdict(float))
    
    div_vals = []
    for relation in config_object.edge_data.divergence_map:
        spec1,spec2,div_time = relation.split(';')
        div_time = float(div_time)
        div_vals.append(div_time)
        div_map[spec1][spec2] = div_time
        div_map[spec2][spec1] = div_time
        
    return div_map,min(div_vals),max(div_vals)
    
def setup_parsers(yopts):
    """
    TODO: make this smarter when you have time.
    """
    y = yopts
    parser_list = []
    
    parser_list.append(CDiffFpkmTrackerParser(cuffdiff_fpkm_path=y.species_info.Aedes_aegypti.expresion_data,
                                              species=y.species_info.Aedes_aegypti.name,
                                              cuffdiff_exp_path=y.species_info.Aedes_aegypti.sig_data,
                                              name_col='nearest_ref_id',
                                              combine_transcripts=True,tx_2_gene=None))
    
    parser_list.append(CDiffFpkmTrackerParser(cuffdiff_fpkm_path=y.species_info.Anopheles_gambiae.expresion_data,
                                              species=y.species_info.Anopheles_gambiae.name,
                                              cuffdiff_exp_path=y.species_info.Anopheles_gambiae.sig_data,
                                              name_col='nearest_ref_id',
                                              combine_transcripts=True,tx_2_gene=None))
    
    parser_list.append(CDiffFpkmTrackerParser(cuffdiff_fpkm_path=y.species_info.Culex_quinquefasciatus.expresion_data,
                                              species=y.species_info.Culex_quinquefasciatus.name,
                                              cuffdiff_exp_path=y.species_info.Culex_quinquefasciatus.sig_data,
                                              name_col='nearest_ref_id',
                                              combine_transcripts=True,tx_2_gene=None))
    
   

    
    species = [y.species_info.Anopheles_gambiae.name,
               y.species_info.Culex_quinquefasciatus.name,
               y.species_info.Aedes_aegypti.name]
    

    div_info = get_div_info(y)
    
    parser_list.append(OneToOneOrthoListParser(list_path=y.edge_data.one_to_one_ortholog_list,
                                               divergence_info=div_info,
                                               relation_type='one_to_one_ortholog'))
    
    return parser_list

def init_metrics():
    """
    TODO: Make this 'smart' when you have more time. 
    """
    expression = PhyloExpnCorrelationIndex(poll_me=False)
    return (expression,)


def construct_builder_and_handler(parser_list):
    """
    *DOES:*
        * initializes and uses GraphBuilder to parse inputs and
          construct the graph relationships
        * initializes the GraphHandler (gBuilder.map_registries_to_graph())
          and installs the RelationsHandler and VoteHandler classes
        
    *RETURNS:*
        * (gHandler,gBuilder)
    """
    gBuilder = GraphBuilder(parser_list)
    gBuilder.populate_registries()
    gHandler = gBuilder.map_registries_to_graph()

    # init metric objs, RelHandler and Vote Handler
    metrics     = init_metrics()
    relHandler  = RelationsHandler(metrics)
    voteHandler = VoteHandler(gHandler.graph)
    #voteHandler.set_vote_types(relHandler.get_vote_types(),weight_by='branch_length')
    gHandler.install_metric_handlers(rel_hndler=relHandler,vote_hndlr=voteHandler)
    
    return gHandler,gBuilder

def reset_random_edges(gHandler,gBuilder,n_way_ortho_table,ortho_parser):
    """
    * deletes current edge connections/objects in gHandler.graph
    * randomizes n-way 1:1 ortholog sets such that actual orthology
      is lost but each ortho-set still contains one gene from each species.
    * rebuilds edge connections/objects in edge_dict and gHandler.graph
    
    :returns: ``None``
    """
    
    graph = gHandler.graph
    node_dict = gHandler.node_dict
    edge_dict = gHandler.edge_dict
    
    # delete current edge connections/objects in gHandler.graph
    graph.remove_edges_from(graph.edges())
    edge_dict.clear()
    
    # randomize n-way 1:1 ortholog sets
    species = n_way_ortho_table.columns
    for s in species:
        np.random.shuffle(n_way_ortho_table[s])
        
    # rebuild edge connections/objects in edge_dict and gHandler.graph
    for row in n_way_ortho_table.iterrows():
        row = list(row[1])
        ortho_parser._resgister_edge(row,node_dict,edge_dict,graph)
        
    gBuilder.map_registries_to_graph(nodes=False,edges=True)

def to_dataFrame_from_node(graph,data_func,index_func=None,column_func=None):
    """
    *GIVEN:*
        * graph (required)
        * data_func (required)
        * index_func (if None: GFuncNode.name is used)
        * column_func (if None: cols numbered)
        
    *DOES:*
        * Iterates through each GFuncNode
        * Creates a pandas DataFrame using the functions provided as parameters
        
    *RETURNS:*
        * data_frame
    """
    # User value interogation and function definition
    if index_func is None:
        def index_func(gfunc_node):
            """
            *GIVEN:*
                * ``GFuncNode``
            *DOES:*
                * gets ``GFuncNode.name``
            *RETURNS:*
                * ``GFuncNode.name``
            """
            return gfunc_node.name
    
    if column_func is None:
        def column_func(gfunc_node):
            """
            *RETURNS:*
                * ``None`` (pandas.DataFrame will number the cols automatically.)
            """
            return None
    
    
    # Build array to feed DataFrame constructor
    data_array = []
    for node in graph.nodes_iter():
        try:
            data_array.append([index_func(node)] + list(data_func(node)))
        except AttributeError:
            pass
        
    # Enforce that all rows have equal columns and get column Names
    column_lengths = set([len(r[1:]) for r in data_array])

    if not len(column_lengths) == 1:
        raise Exception("Sanity Check failed: not all rows have same number of columns.")
    
    # Build pandas.DataFrame
    data_array = np.array(data_array)
    data    = np.array(data_array[:,1:],dtype=np.float)
    data    = np.ma.masked_invalid(data)
    index   = data_array[:,0]
    columns = column_func(node)
    data_frame = pandas.DataFrame(data=data, index=index, columns=columns) # will use last node from loop above
    
    return data_frame
    
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
    ortho_parser = parser_list[-1]
    
    gHandler,gBuilder = construct_builder_and_handler(parser_list)
    gHandler.measure_relations()
    
    #gene_list = get_starting_nodes(yopts.edge_data.one_to_one_ortholog_list)
    #scored_orthos = get_ortho_set_PTCIs(gHandler,gene_list)    
    
    #ptci_dicts = [scored_orthos]
    
    n_way_ortho_table= pandas.read_table(yopts.edge_data.one_to_one_ortholog_list)
    #for rep in range(100):
        #reset_random_edges(gHandler,gBuilder,n_way_ortho_table,ortho_parser)
        #gHandler.measure_relations()
        
        ##g = gHandler.graph
        ##edges = g.edges()
        
        #scored_random_edges = get_ortho_set_PTCIs(gHandler,gene_list)
        #ptci_dicts.append(scored_random_edges)

        
    
    
    
    
    
    #top_ortho_dict = get_top_ortho_sets(ortho_meanScore_dict=scored_orthos,percentile=90)
    
    #write_deltaFPKM_vectors(top_ortho_dict,gHandler,outPath=yopts.outputs.top_ortho_deltaFPKM)
    
    #def test_enrichment():
        #"""
        #"""
        #ptci = ptci_dicts[0]
        ##pd = load_prot_domains(gHandler.graph,['/home/augustine/Dropbox/common/tmp/test_ensmebl_rest/aaegypti_eg_gene.protDomains.txt',
                                               ##'/home/augustine/Dropbox/common/tmp/test_ensmebl_rest/agambiae_eg_gene.protDomains.txt',
                                               ##'/home/augustine/Dropbox/common/tmp/test_ensmebl_rest/cquinquefasciatus_eg_gene.protDomains.txt'])
        #gs = load_go_slims(gHandler.graph,['/home/augustine/Dropbox/common/tmp/test_ensmebl_rest/aaegypti_eg_gene.goSlim.txt',
                                           #'/home/augustine/Dropbox/common/tmp/test_ensmebl_rest/agambiae_eg_gene.goSlim.txt',
                                           #'/home/augustine/Dropbox/common/tmp/test_ensmebl_rest/cquinquefasciatus_eg_gene.goSlim.txt'])
        #pos_nodes = []
        #dump = [pos_nodes.extend(x) for x in ptci if ptci[x] > 0.]
        #for i,v in enumerate(pos_nodes):
            #pos_nodes[i] = gHandler.node_dict[v]
        
        #pos_nodes = set(pos_nodes)
        #def term_func(node):
            #return frozenset(node.data.go_slims[node.name].keys())
            ##return node.data.prot_domains[node.name]['pfam']
        #table = enrichment_of_terms(graph=gHandler.graph, positive_node_set=pos_nodes, term_fuction=term_func)
        ##runExternalApp('fembot','-r 0 -t "Gus, I need your attention."')
        #return table
        
    #results_table = test_enrichment()
    def data_func(gfunc_node):
        """
        *GIVEN:*
            * ``GFuncNode``
        *RETURNS:*
        *   - ``GFuncNode.data.expression_vector``
        """
        return gfunc_node.data.expression_vector
    
    def index_func(gfunc_node):
        """
        *GIVEN:*
            * ``GFuncNode``
        *RETURNS:*
            * ``GFuncNode.species``
        """
        return gfunc_node.species    
    
    expn_dataFrame = to_dataFrame_from_node(graph=gHandler.graph,data_func=data_func,index_func=index_func,column_func=None)
    
    return gHandler,gBuilder,n_way_ortho_table,ortho_parser
    
    # dump built construct as pickle
    #cPickle.dump(gHandler,open(yopts.outputs.graph_pickle,'w'))

if __name__ == '__main__':
    trap = main()
    #main()
    print "main() completed."