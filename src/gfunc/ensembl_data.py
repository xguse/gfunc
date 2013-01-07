"""
####################
ensembl_data.py
####################
Code supporting automated retrieval/mirroring/processing of ensembl data files and directories.

WARNING:
    To some extent the functionality of this code is dependant on the
    current formating of Ensembl's internal directory structure.
"""
import urllib2

from gfunc import externals
from gfunc.data_classes import Bunch

class DataGrabber(object):
    """
    Class to manage conecting to ensembl-based ftp data dumps and retrieving them.
    """    
    
    
    
    
    def __init__(self,base_url,species,data_types,base_local_path):
        """
        Class to manage conecting to ensembl-based ftp data dumps.
        """

        self._supported_programs = ('curl','wget')
        self._init_map = {'curl':self._init_curl,
                          'wget':self._init_wget}
        
        self._settings = Bunch()
        self._settings.base_url = base_url
        self._settings.species = species
        self._settings.data_types = data_types
        self._settings.base_local_path = base_local_path
        self._avail_data_types = self._get_avail_data_types() # will be tuple
        
            
    
    
    def _get_avail_data_types(self):
        """
        *GIVEN:*
            * self
        *DOES:*
            * [x] queries base_url for the availible directories (data_types)
            * [x] stores a tuple to self._avail_data_types
        *RETURNS:*
            * None
        """
        page = urllib2.urlopen(self._settings.base_url)
        data_types = [x.split()[-1] for x in page]
        
        self._avail_data_types = tuple(data_types)
    
    def _init_curl(self):
        """
        *GIVEN:*
            * self
        *DOES:*
            * [] format the command line string for cURL based on self._settings
        *RETURNS:*
            * command line string
        """
        command_string = ''
        
        
    def _init_wget(self):
        """
        *GIVEN:*
            * self
        *DOES:*
            * [] format the command line string for wget based on self._settings
        *RETURNS:*
            * command line string
        """
        pass
        
    def _execute(self,program):
        """
        *GIVEN:*
            * self
            * supported program
        *DOES:*
            * [x] initiates command line string based on which program is included
            * [x] runs externals.runExternalApp(program,args)
        *RETURNS:*
            * return externals.runExternalApp(program,args)
        """
        command_string = self._init_map[program.split('/')[-1]]()
        
        return externals.runExternalApp(progName=program,argStr=command_string)
        

    
    def _which_transfer_method(self):
        """
        *GIVEN:*
            * self
        *DOES:*
            * [x] queries system for certain external file transfer programs
              (NOTE: for now in order of preference [cURL, wget]. Eventually,
              I may add more and/or include a pure urllib2 solution for
              systems w/o supported external apps.)
            * [x] chooses one that will work on the system
        *RETURNS:*
            * choice of program
        """
        
        accepted_program_path = False
        for prog in self._supported_programs:
            accepted_program_path = externals.whereis(prog)
            if accepted_program_path:
                break
        
        if not bool(accepted_program_path):
            raise Exception("Your system does not seem to have any of the supported download methods availible: %s.  Please install one and try again.") % (self._supported_programs)
        else:
            return accepted_program_path
        
    def settings(self):
        """
        
        *RETURNS:*
            * dict of current settings for ensembl access.
        
        """
        return self._settings.copy()
    
    def transfer_data(self,unzip=False):
        """
        *GIVEN:*
            * sufficently initiated self instance
        *DOES:*
            * [x] decides which method to use to get the data
            * [x] creates local directory if needed based on base_local_path
            * [x] initiates data transfer
            * [] complains if it detects incomplete transfer
            * [] if unzip evaluates to True, recursively unzip any files
              with extentions suggesting they are compressed.
        *RETURNS:*
            * None
        """
        prog = self._which_transfer_method()
        
        externals.mkdirp(self._settings.base_local_path)
        
        execution_result = self._execute(prog)
    