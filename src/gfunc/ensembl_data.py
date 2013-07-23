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
from gfunc.fileIO import walk_dirs_for_fileName
from gfunc.data_classes import Bunch


def to_stdout(txt):
    sys.stdout.write(txt)
    sys.stdout.flush()
    
def to_stderr(txt):
    sys.stderr.write(txt)
    sys.stderr.flush()




def validate_downloads(base_dir):
    """
    *GIVEN:*
        * ``base_dir`` = top-level directory containing files to be validated.
    *DOES:*
        * Decends into ``base_dir`` and records paths to all lower-level files.
        * Uses the CHECKSUMS files it finds to set up external system calls to ``sum``
          for each file listed in the direcory's CHECKSUMS file.
        * A file named VALIDATIONS is created in each directory listing each file in the directory and one of the following outcomes: 
            * PASS = passed checksum match
            * FAIL = failed checksum match
            * NOT_IN_CHECKSUMS = file found in directory but NOT listed in the CHECKSUMS file
            * NOT_IN_DIR = file listed in CHECKSUMS file, but not found in the directory.
    *RETURNS:*
        * ``results`` = a summary of all VALIDATIONS files (list of strings)
    """
    
    results = []
    
    file_paths = walk_dirs_for_fileName(dir_path=base_dir,pattern="*")
    
    checksum_paths = [x for x in file_paths if x.endswith('CHECKSUMS')]
    
    for cksum in checksum_paths:
        cksum_data = [line.strip('\n').split() for line in open(cksum)]
        
        current_dir_path = cksum.replace('/CHECKSUMS','')
        
        scores = check_files(dir_path=current_dir_path , cksum_data=cksum_data)
        results.extend(scores)
        
        vFile = open('%s/VALIDATIONS' % (current_dir_path), 'w')
        for line in scores:
            vFile.write('%s\n' % (line))
        vFile.close()
    
    return results
    
        
    
        
def check_files(dir_path,cksum_data):
    """
    *GIVEN:*
        * ``dir_path`` = path to a directory containing files to be validated
        * ``cksum_data`` = parsed contents of the CHECKSUMS file for this directory
    *DOES:*
        * Compares files in the directory with checksums in the CHECKSUMS file and scores them as PASS/FAIL
        * Documents and classifies discrepancies between files listed in CHECKSUMS file vs files actually
          in the directory: classifies them as NOT_IN_CHECKSUMS or NOT_IN_DIR.
    *RETURNS:*
        * ``results`` = list of strings
    """
    results = []
    
    cksums = {}
    for cksm,blocks,name in cksum_data:
        cksums[name] = (cksm,blocks)
    
    # get the sets of file names from CHECKSUMS file as well as from the actual directory.
    in_cksums  = set([x[-1] for x in cksum_data])
    in_cur_dir = walk_dirs_for_fileName(dir_path=dir_path,pattern='*')
    in_cur_dir = set([x.split('/')[-1] for x in in_cur_dir if x.split('/')[-1]])
    
    in_cur_dir.discard('CHECKSUMS')
    in_cur_dir.discard('VALIDATIONS')
    
    not_in_cksums = list(in_cur_dir - in_cksums) # files in current directory but not listed in CHECKSUMS file
    not_in_dir    = list(in_cksums - in_cur_dir) # files in CHECKSUMS file but not current directory
    in_both       = list(in_cksums & in_cur_dir) # files in both
    
    for f in in_both:
        stdout,stderr = externals.runExternalApp(progName='sum',argStr='%s/%s' % (dir_path,f))
        if tuple(stdout.split()) == cksums[f]:
            results.append('PASS\t%s' % (f))
        else:
            results.append('FAIL\t%s' % (f))
    
    for f in not_in_cksums:
        results.append('NOT_IN_CHECKSUMS\t%s' % (f))
    
    for f in not_in_dir:
        results.append('NOT_IN_DIR\t%s' % (f))
    
    results.sort(key=lambda x: x.split()[-1]) # sort results by file name
    
    return results
        


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
    
    for attempt in range(5):
        try:
            page = urllib2.urlopen(url)
            results = [(x.split()[0] , x.split()[-1]) for x in page]
            break
        except:
            if attempt > 5:
                raise
            else:
                to_stderr("Connection failed: trying again...\n")
            
    
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

    def __init__(self,base_url,species,data_types,base_local_path,verbose=False):
        """
        Initiate DataGrabber object for conecting to ensembl-based ftp data dumps.
        """

        self._verbose = verbose
        self._supported_programs = ('aria2c','rsync','curl','wget')
        self._init_map = {'rsync':self._init_rsync,
                          'aria2c':self._init_aria2c,
                          'curl':self._init_curl,
                          'wget':self._init_wget}
        
        self._settings = Bunch()
        self._settings.base_url = base_url.rstrip('/')
        self._settings.species = species
        self._settings.data_types = data_types
        self._settings.base_local_path = base_local_path
        self._get_avail_data_types() # will be tuple
        
            
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
        if self._verbose:
            to_stdout('Checking server for current data types...')
        contents = web_ls(self._settings.base_url)
        
        self._avail_data_types = tuple([x.split('/')[-1] for x in contents['dirs']])
        if self._verbose:
            to_stdout('\tAvailible types set to: %s\n' % (str(self._avail_data_types)))
    
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
        if self._verbose:
            to_stdout('Collecting file urls for download...')
        for target_dir in target_directories:
            target_urls.extend(web_walk(target_dir))
        
        self._target_urls = tuple(target_urls)
        if self._verbose:
            to_stdout('\tDONE!\n')
    
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
            * format the command line string for aria2c based on self._settings
            * writes a file that encodes the target_urls and their new local paths for aria2c to read.
        *RETURNS:*
            * command line string
        """
        
        # write out a log file containing the targeted urls and their new local paths for aria2c to read  
        file_path = '%s/target_urls_aria2c_input.txt' % (self._settings.base_local_path)
        aria2c_input_file = open(file_path, 'w')
        if self._verbose:
            to_stdout('Writing options file for aria2c...')
        for target_url in self._target_urls:
            aria2c_input_file.write("%s\n  dir=%s\n  out=%s\n" % (target_url,
                                                                  self._settings.base_local_path,
                                                                  target_url.split('://')[1]))
        aria2c_input_file.close()
        if self._verbose:
            to_stdout('\tDONE!\n')        
        
        cmd_string = '-l %s/aria2c.log -c -i %s -j5' % (self._settings.base_local_path,file_path)
        
        
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
        if self._verbose:
            to_stdout('Initializing command string for: %s\n' % (program))
        command_string = self._init_map[program.split('/')[-1]]()
        
        if self._verbose:
            to_stdout('Executing command: %s %s\n' % (program,command_string))
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
        if self._verbose:
            to_stdout('Choosing download program from your system...')
        prog = self._which_transfer_method()
        if self._verbose:
            to_stdout('\t SELECTED: %s\n' % (prog))
        
        # create local home for the data
        externals.mkdirp(self._settings.base_local_path)
        
        # write out a log file containing the targeted urls for reference.
        self._get_url_list()
        url_list_file = open('%s/target_urls.txt' % (self._settings.base_local_path), 'w')
        for target_url in self._target_urls:
            url_list_file.write("%s\n" % (target_url))
        url_list_file.close()
        self._target_url_file_path = os.path.abspath(url_list_file.name)
        
        # begin execution of transfer
        execution_result = self._execute(prog)

if __name__ == '__main__':
    pass