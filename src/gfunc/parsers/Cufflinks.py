"""
####################
Cufflinks.py
####################
Code supporting parsing of Cufflinks type raw data files.
"""

from gfunc.data_classes import GFuncNode

from rSeq.utils.files import tableFile2namedTuple

class CDiffFpkmTrackerParser(object):
    """
    Class to accept CuffDiff FKPM data table and init the relevant gFuncNode Objects.
    """
    
    def __init__(self, cuffdiff_path, gFunc_graph, species, name_col='gene_short_name'):
        """
        Test doc for init'ing CuffDiff parser.
        """
        # everything after rows[9] is FPKM data
        
        self.species = species
        
        self._graph = gFunc_graph # may not need
        self._name_col = name_col
        self._tableFile = tableFile2namedTuple(tablePath=cuffdiff_path,sep='\t')
    
    def _setup_expn_vector(self, table_row):
        """
        """
        name = table_row[self._name_col]
        vector = [float(table_row[x]) for x in range(10,len(table_row),3)]
        
        return name,vector
    
    def resgister_nodes(self,node_dict):
        """
        Parses each row from the CuffDiff FKPM data table and either adds the data to the relevant
        GFuncNode in node_dict or creates one and adds it to that then registers it in node_dict.
        """
        for row in self._tableFile: # may want to close this out
            name,vector = self._setup_expn_vector(row)
            try:
                node_dict[name].set_data(data=vector,data_type='expression_vector')
            except KeyError:
                node = GFuncNode(name=name, species=self.species, is_target=False, debug=False)
                node_dict[name].set_data(data=vector,data_type='expression_vector')

        
