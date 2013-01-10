"""
####################
ensembl_data.py
####################
Code supporting automated retrieval/mirroring/processing of ensembl data files and directories.

WARNING:
    To some extent the functionality of this code is dependant on the
    current formating of Ensembl's internal directory structure and how
    it responts to ``urllib2.urlopen(url)``.
"""
import urllib2
import sys
import os
from collections import defaultdict


from gfunc import externals
from gfunc.data_classes import Bunch


def web_ls(url):
    """
    *GIVEN:*
        * ``url`` = the url of an ensembl-based ftp:// target
    *DOES:*
        * reads the url data and extracts file/directory names
        * stores info in a dict named ``contents``; keyed by 'dirs' or 'files' and pointing to lists of respective urls.
    *RETURNS:*
        * ``contents`` dictionary.
    """
    page = urllib2.urlopen(url)
    results = [(x.split()[0] , x.split()[-1]) for x in page]
    
    contents = defaultdict(list)
    for permissions_str,name in results:
        if permissions_str.startswith('d'):
            contents['dirs'].append(url.rstrip('/') + '/' + name)
        else:
            contents['files'].append(url.rstrip('/') + '/' + name)
    
    return contents

def web_walk(base_url):
    """
    *GIVEN:*
        * ``base_url`` = the url of an ensembl-based ftp:// target directory
    *DOES:*
        * recursively stores directories and file urls (uses ``web_ls``),
          continues to follow directories until all file urls have been
          collected below ``base_url``.
    *RETURNS:*
        * ``file_urls`` = list of file url strings under ``base_url``.
    """
    
    file_urls = []
    
    contents = web_ls(base_url)
    file_urls.extend(contents['files'])
    
    for dir_url in contents['dirs']:
        file_urls.extend(web_walk(dir_url))
    
    return file_urls

class DataGrabber(object):
    """
    Class to manage conecting to ensembl-based ftp data dumps and retrieving them.
    """    

    def __init__(self,base_url,species,data_types,base_local_path):
        """
        Initiate DataGrabber object for conecting to ensembl-based ftp data dumps.
        """

        self._supported_programs = ('aria2c','rsync','curl','wget')
        self._init_map = {'rsync':self._init_rsync,
                          'aria2c':self._init_aria2c,
                          'curl':self._init_curl,
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
            * [x] queries base_url for the availible directories (data types)
            * [x] stores a tuple to self._avail_data_types
        *RETURNS:*
            * None
        """
        contents = web_ls(self._settings.base_url)
        
        self._avail_data_types = tuple([x.split('/')[-1] for x in contents['dir']])
    
    def _get_url_list(self):
        """
        *GIVEN:*
            * sufficently initiated ``self`` instance
        *DOES:*
            * [x] uses:
            
                # ``self._settings.base_url``
                # ``self._settings.species``
                # ``self._settings.data_types``
            
              to construct full urls for every file to be downloaded
              and saves them to self._target_urls.
        *RETURNS:*
            * ``None`` 
        """
        
        target_directories = []
        target_urls = []
        
        base_url = self._settings.base_url
        for species in self._settings.species:
            for data_type in self._settings.data_types:
                target_directories.append('%s/%s/%s' % (base_url,data_type,species))
                
        for target_dir in target_directories:
            target_urls.extend(web_walk(target_dir))
        
        self._target_urls = tuple(target_urls)
    
    def _init_curl(self):
        """
        *GIVEN:*
            * self
        *DOES:*
            * [] format the command line string for cURL based on self._settings
        *RETURNS:*
            * command line string
        """
        raise NotImplementedError()
        
    def _init_wget(self):
        """
        *GIVEN:*
            * self
        *DOES:*
            * [] format the command line string for wget based on self._settings
        *RETURNS:*
            * command line string
        """
        raise NotImplementedError
    
    def _init_aria2c(self):
        """
        *GIVEN:*
            * self
        *DOES:*
            * [] format the command line string for aria2c based on self._settings
        *RETURNS:*
            * command line string
        """
        cmd_string = ' -i %s -j5' % (self._target_url_file_path)
        
        return cmd_string
        
    def _init_rsync(self):
        """
        *GIVEN:*
            * self
        *DOES:*
            * [] format the command line string for rsync based on self._settings
        *RETURNS:*
            * command line string
        """
        raise NotImplementedError    
        
        
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
              (NOTE: for now in order of preference (see ``self._supported_programs``). 
              Eventually, I may add more and/or include a pure urllib2 solution for
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
            * sufficently initiated ``self`` instance
        *DOES:*
            * [x] decides which method to use to get the data
            * [x] creates local directory if needed based on base_local_path
            * [x] initiates data transfer
            * [?] complains if it detects incomplete transfer
            * [?] if unzip evaluates to True, recursively unzip any files
              with extentions suggesting they are compressed.
        *RETURNS:*
            * None
        """
        # choose first program to try 
        prog = self._which_transfer_method()
        
        # create local home for the data
        externals.mkdirp(self._settings.base_local_path)
        
        # write out a log file containing the targeted urls.
        self._get_url_list()
        url_list_file = open('%s/target_urls.txt' % (self._settings.base_local_path), 'w')
        for target_url in self._target_urls:
            url_list_file.write("%s\n" % (target_url))
        url_list_file.close()
        self._target_url_file_path = os.path.abspath(url_list_file.name)
        
        # begin execution of transfer
        execution_result = self._execute(prog)
    