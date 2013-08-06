"""
####################
Cufflinks.py
####################
Code supporting parsing of Cufflinks type raw data files.
"""
from collections import defaultdict

import numpy as np

#from bulbs.model import Node as bNode
#from bulbs.model import Relationship as bRelationship
#from bulbs.property import String as bString
#from bulbs.property import Integer as bInteger

from gfunc.parsers.base import GFuncParserBase
from gfunc.data_classes import GFuncNode
from gfunc.data_classes import Bunch
from gfunc.fileIO import tableFile2namedTuple

class CDiffFpkmTrackerParser(GFuncParserBase):
    """
    Class to accept CuffDiff FKPM data table and init/update the relevant gFuncNode Objects.
    """
    
    def __init__(self, cuffdiff_fpkm_path, species, cuffdiff_exp_path=None, name_col='nearest_ref_id',combine_transcripts=True,tx_2_gene=None):
        """
        Test doc for init'ing CuffDiff parser.
        """
        # everything after rows[9] is FPKM data
        
        self.species = species
        self.data_type = 'expression_vector'
        self._name_col = name_col
        self._expDiff  = cuffdiff_exp_path
        self._tableFile = tableFile2namedTuple(tablePath=cuffdiff_fpkm_path,sep='\t')
        self._combine_tx = combine_transcripts
        if combine_transcripts:
            if tx_2_gene is None:
                self._tx2gene_func = lambda x: x[:-3]
            else:
                self._tx2gene_func = tx_2_gene
            self.gene_vectors = self._sum_transcripts()
            
        if self._expDiff is not None:
            self._expDiff = build_expDiffTable_dict(self._expDiff)
    
    def _setup_expn_vector(self, table_row):
        """
        TODO: Doc
        """
        if not self._combine_tx:
            name = table_row.get(self._name_col)
            xloc = table_row.gene_id
            vector = np.array([float(table_row[x]) for x in range(10,len(table_row),3)])
            return name,xloc,vector
        else:
            gene_rec = table_row # for clairity that this is different process
            name_xloc,vector = gene_rec
            name,xloc = name_xloc     # name_xloc = (name,xloc)
            return name,xloc,vector
        
    def _sum_transcripts(self):
        """
        TODO: doc
        """
        # TODO: comment this
        def zeros():
            return np.zeros(len(self._tableFile[0][10::3]))
        
        gene_vectors = defaultdict(zeros)
        
        for row in self._tableFile:
            name = self._tx2gene_func(row.get(self._name_col)) # Converts Tx symbol into Gene Symbol
            xloc = row.gene_id
            vector = np.array([float(row[x]) for x in range(10,len(row),3)])
            name_xloc = (name,xloc)
            gene_vectors[name_xloc] = gene_vectors[name_xloc] + vector
            
        return gene_vectors
             
    
    def resgister_nodes_and_edges(self,node_dict,edge_dict,graph):
        """
        Parses each row from the CuffDiff FKPM data table and either adds the data to the relevant
        GFuncNode in node_dict or creates one and adds it to that then registers it in node_dict.
        """
        if not self._combine_tx:
            data = self._tableFile
        else:
            data = self.gene_vectors.iteritems()
        for record in data: # may want to close this out
            name,xloc,vector = self._setup_expn_vector(record)
            try:
                node_dict[name].set_data(data=vector,data_type=self.data_type)
            except KeyError:
                node = GFuncNode(name=name, species=self.species, graph=graph, is_target=False, debug=False)
                node.set_data(data=vector,data_type=self.data_type)
                node_dict[name] = node
                
            node_dict[name]._sigDiff = am_i_sigDiff(xloc_number=xloc, expDiffTable_dict=self._expDiff, q_thresh=0.05)

def transfer_nearestRefgeneSymbol_from_isoform_to_gene_tracking(isoform_fpkm_path,gene_fpkm_path):
    """
    TODO: doc
    """
    iso_table = tableFile2namedTuple(isoform_fpkm_path)
    gene_table = tableFile2namedTuple(gene_fpkm_path)
    
    gene_id_2_nearest_ref = defaultdict(set)
    
    for row in iso_table:
        gene_id_2_nearest_ref[row.gene_id].add(row.nearest_ref_id[:-3])
    
    for gene_id,symbol_set in gene_id_2_nearest_ref.iteritems():
        if len(symbol_set) != 1:
            raise ValueError('Gene_id to gene symbol conflict: %s to %s' % (gene_id,symbol_set))
        
    headers = '\t'.join(gene_table[0]._fields) + '\n'
    out_file = open('%s.gene_symbols.csv' % (gene_fpkm_path),'w')
    out_file.write(headers)
    for row in gene_table:
        row.nearest_ref_id = gene_id_2_nearest_ref[row.gene_id].pop()
        out_file.write('%s\n' % ('\t'.join(row)))
    out_file.close()
    
    

# ###################
# Complex Output Model
# ###################

##class IsoformFPKMTrackingTable(Node):

    ##element_type = "isoform_fpkm"

##class GeneFPKMTrackingTable(Node):

    ##element_type = "gene_fpkm"
    
##class IsoformExpDiffTable(Node):

    ##element_type = "isoform_exp_diff"

##class GeneExpDiffTable(Node):

    ##element_type = "gene_exp_diff"


##class Knows(Relationship):
    
    ##label = "knows"
    
    ##created = DateTime(default=current_datetime, nullable=False)

##def create_full_cufflinks_model():
    ##"""
    
    ##"""

def build_expDiffTable_dict(expDiffTable_path):
    """
    Build isoformExpDiffTable_dict:
    
    Keys:
        XLOC_xxxxx
    
    Values:
        namedtuple-ified rows with same XLOC_xxxxx
    """
    
    rows = tableFile2namedTuple(expDiffTable_path)
    expDiffTable_dict = defaultdict(list)
    
    for row in rows:
        expDiffTable_dict[row.gene_id].append(row)
    
    return expDiffTable_dict
    
def am_i_sigDiff(xloc_number, expDiffTable_dict, q_thresh):
    """
    | for lines in cuffdiff_fpkm_table:
    |     Return ``True`` if:
    |         at least one of current line's XLOC_xxxx pair tests in ``expDiffTable_dict`` has *q_val* <= ``q_thresh``
    |     Else Return ``False``
        
    """
    try:
        q_vals = [float(x.q_value) for x in expDiffTable_dict[xloc_number]]
        
        if min(q_vals) <= q_thresh:
            return True
        else:
            return False
    except TypeError:
        return None