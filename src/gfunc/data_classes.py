"""
####################
data_classes.py
####################
Code defining container classes for supported data types.
"""

class GenomeFeature(object):
    """
    Class to store and access data for a single GenomeFeature.
    It should be called and populated by the relevant parser objects.
    """
    def __init__(self, transcript_id, species_id):
        """
        Test comment for initiating TxObj.
        """
        self.transcript_id = str()
        self.species_id = str()
        self.abundance_data = 
        
class TranscriptAbundance(object):
    """
    Class to model a single replicate of mRNA abundance data from a single experiment set.
    """
    def __init__(self, transcript_id, species_id, experiement_id, expression_units, abundance_vector):
        """
        """
        self.transcript_id     = transcript_id
        self.species_id        = species_id        
        self.experiement_id    = experiement_id
        self.expression_units  = expression_units
        self.abundance_vector  = abundance_vector # namedtuple 
        self.condition_headers = abundance_vector._fields
    
    def __repr__(self):
        """
        """
        return "%s;%s" % (self.transcript_id,self.experiement_id)