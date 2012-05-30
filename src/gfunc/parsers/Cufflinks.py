"""
####################
Cufflinks.py
####################
Code supporting parsing of Cufflinks type raw data files.
"""
from collections import defaultdict

import numpy as np


from gfunc.parsers.base import GFuncParserBase
from gfunc.data_classes import GFuncNode
from gfunc.data_classes import Bunch
from gfunc.fileIO import tableFile2namedTuple

class CDiffFpkmTrackerParser(GFuncParserBase):
    """
    Class to accept CuffDiff FKPM data table and init/update the relevant gFuncNode Objects.
    """
    
    def __init__(self, cuffdiff_path, species, name_col='nearest_ref_id',combine_transcripts=True,tx_2_gene=None):
        """
        Test doc for init'ing CuffDiff parser.
        """
        # everything after rows[9] is FPKM data
        
        self.species = species
        self.data_type = 'expression_vector'
        self._name_col = name_col
        self._tableFile = tableFile2namedTuple(tablePath=cuffdiff_path,sep='\t')
        self._combine_tx = combine_transcripts
        if combine_transcripts:
            if tx_2_gene is None:
                self._tx2gene_func = lambda x: x[:-3]
            else:
                self._tx2gene_func = tx_2_gene
            self.gene_vectors = self._sum_transcripts()
    
    def _setup_expn_vector(self, table_row):
        """
        TODO: Doc
        """
        if not self._combine_tx:
            name = table_row.get(self._name_col)
            vector = np.array([float(table_row[x]) for x in range(10,len(table_row),3)])
            return name,vector
        else:
            gene_rec = table_row # for clairity that this is different process
            name,vector = gene_rec
            return name,vector
        
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
            vector = np.array([float(row[x]) for x in range(10,len(row),3)])
            gene_vectors[name] = gene_vectors[name] + vector
            
        return gene_vectors
             
    
    def resgister_nodes_and_edges(self,node_dict,edge_dict):
        """
        Parses each row from the CuffDiff FKPM data table and either adds the data to the relevant
        GFuncNode in node_dict or creates one and adds it to that then registers it in node_dict.
        """
        if not self._combine_tx:
            data = self._tableFile
        else:
            data = self.gene_vectors.iteritems()
        for record in data: # may want to close this out
            name,vector = self._setup_expn_vector(record)
            try:
                node_dict[name].set_data(data=vector,data_type=self.data_type)
            except KeyError:
                node = GFuncNode(name=name, species=self.species, is_target=False, debug=False)
                node.set_data(data=vector,data_type=self.data_type)
                node_dict[name] = node

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
    
    

        
