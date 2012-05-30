"""
####################
fileIO.py
####################
Code supporting reading and writing from files not related to specific parsers.
"""

import sys
import os
import fnmatch
import csv
from collections import namedtuple


def tableFile2namedTuple(tablePath,sep='\t',headers=None):
    """
    Returns namedTuple from table file using first row fields as
    col headers or a list supplied by user.
    """

    reader  = csv.reader(open(tablePath,'rU'), delimiter=sep)
    if not headers:
        headers = reader.next()
    headers = [x.replace(' ','_').replace('.','_').replace('(','_').replace(')','_').replace('-','_').replace(':','_') for x in headers]
    Table   = namedtuple('Table', headers)
    # wrap Table.__getattribute__() for less typing
    def get(self,colName):
        return self.__getattribute__(colName)
    Table.get = get
    
    data    = [Table._make(x) for x in reader if x!=[]] # reader kept feeding an empty list at the end that botched everything!  wtf?!
    return data

def walk_dirs_for_fileName(dir_path,pattern="*.xml"):
    """
    Recursively collects file paths in a dir and subdirs.
    """
    file_paths = []
    for root, dirs, files in os.walk(dir_path):
        
        for filename in fnmatch.filter(files, pattern):
            file_paths.append(os.path.join(root, filename))
            
        for sub_dir in dirs:
            file_paths.append(walk_dirs_for_fileName(sub_dir,pattern))
            
    return file_paths