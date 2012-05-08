"""
####################
parsers.py
####################
Code supporting specific data parsing.
"""

class CuffDiffParser(object):
    """
    Class to accept CuffDiff isoform FKPM data table and init the relevant TranscriptAbundance Objects.
    """
    def __init__(self, cuffdiff_path, options=None):
        """
        Test doc for init'ing CuffDiff parser.
        """