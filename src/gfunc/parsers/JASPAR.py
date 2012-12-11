
from collections import OrderedDict

import numpy as np

from gfunc.parsers.base import GFuncParserBase
from gfunc.fileIO import tableFile2namedTuple
from gfunc.data_classes import GFuncNode

class BasicTFBSParser(GFuncParserBase):
    """
    Class to accept TFBS profile data table from 'find_motifs.py' and init/update the relevant gFuncNode Objects.
    """
    
    def __init__(self, tfbs_path):
        """
        Test doc for init'ing TFBS parser.
        """
        
        self.data_type = 'tfbs_vector'
        self._tableFile = tableFile2namedTuple(tablePath=tfbs_path,sep='\t')
    
    def _setup_tfbs_vector(self, table_row):
        """
        """
        name = table_row.seq_name
        vector = np.array([float(x) for x in table_row[2:]])
        species = table_row.species
        
        return name,vector,species
    
    def resgister_nodes_and_edges(self,node_dict,edge_dict,graph):
        """
        Parses each row from the TFBS data table and either adds the data to the relevant
        GFuncNode in node_dict or creates one and adds it to that then registers it in node_dict.
        """
        for row in self._tableFile: # may want to close this out
            name,vector,species = self._setup_tfbs_vector(row)
            try:
                node_dict[name].set_data(data=vector,data_type=self.data_type)
            except KeyError:
                node = GFuncNode(name=name, species=species, graph=graph, is_target=False, debug=False)
                node.set_data(data=vector,data_type=self.data_type)
                node_dict[name] = node


class ParseJasparMatrixOnly(object):
    """
    Returns a record-by-record motif parser for JASPAR matric_only.txt files analogous to file.readline().
    example:
    >MA0001.1 AGL3
    A  [ 0  3 79 40 66 48 65 11 65  0 ]
    C  [94 75  4  3  1  2  5  2  3  3 ]
    G  [ 1  0  3  4  1  0  5  3 28 88 ]
    T  [ 2 19 11 50 29 47 22 81  1  6 ]
    >MA0002.1 RUNX1
    A  [10 12  4  1  2  2  0  0  0  8 13 ]
    C  [ 2  2  7  1  0  8  0  0  1  2  2 ]
    G  [ 3  1  1  0 23  0 26 26  0  0  4 ]
    T  [11 11 14 24  1 16  0  0 25 16  7 ]
    """
    def __init__(self,filePath):
        """Returns a record-by-record motif parser analogous to file.readline().
        Exmpl: parser.next()
        Its ALSO an iterator so "for rec in parser" works too!
        """
        
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
            

        self.bufferLine = None   # stores next headerLine between records.
        self._stop = False
        
    
    def __iter__(self):
        return self
    
    def _recData2array(self,recData):
        """
        Converts lines like this:
        A  [ 0  3 79 40 66 48 65 11 65  0 ]
        C  [94 75  4  3  1  2  5  2  3  3 ]
        G  [ 1  0  3  4  1  0  5  3 28 88 ]
        T  [ 2 19 11 50 29 47 22 81  1  6 ] 
        
        to a 2D list like this:
        [[ 0,  3, 79, 40, 66, 48, 65, 11, 65,  0 ],
        [94, 75,  4,  3,  1,  2,  5,  2,  3,  3 ],
        [.],
        [.]]
         
        Returns 2D list.
        """
        motif_dict = {}
        for line in recData:
            line = line.strip('\n').replace('[','').replace(']','').split()
            motif_dict[line[0]] = [float(x) for x in line[1:]]
            
        motif = []
        for base in ['A','C','G','T']:
            motif.append(motif_dict[base])
            
        return motif
        
    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqName,seqStr)"""
        if not self._stop:
            pass
        else:
            raise StopIteration
        # ++++ Get A Record ++++
        recHead = ''
        recData = []
        # ++++ Check to see if we already have a headerLine ++++
        if self.bufferLine:
            recHead = self.bufferLine
        else:
        # ++++ If not, seek one ++++
            while 1:
                try:
                    line = self._file.next()
                except StopIteration:
                    self._stop = True
                    break
                if line.startswith('>'):
                    recHead = line.lstrip('>').strip('\n').replace(' ','_')
                    break
                elif not line:
                    raise InvalidFileFormatError, "CheckFastaFile: Encountered EOF before any data."
                elif line == '\n':
                    continue
                else:
                    raise InvalidFileFormatError, 'CheckFastaFile: The first line containing text does not start with ">".'
        # ++++ Collect recData ++++
        while 1:
            try:
                line = self._file.next()
            except StopIteration:
                self._stop = True
                break
            if line.startswith('>'):
                self.bufferLine = line.lstrip('>').strip('\n').replace(' ','_').replace('.','_')
                break
            elif not line.startswith('>'):
                recData.append(line.strip('\n'))

        # ++++ Minor Seq Validation ++++
        ## AddHere
        # ++++ Format Rec For Return ++++
        if not recData:
            return (recHead,'')
        else:            
            return (recHead,self._recData2array(recData))   
    
    def to_dict(self):
        """Returns a single OrderedDict populated with the motifRecs
        contained in self._file."""
        motifDict = OrderedDict()
        while 1:
            try:
                motifRec = self.next()
            except StopIteration:
                break
            if motifRec:
                if not motifRec[0] in motifDict:
                    motifDict[motifRec[0]] = motifRec[1]
                else:
                    raise ValueError, "DuplicateFastaRec: %s occurs in your file more than once."
            else:
                break
        return motifDict
    
