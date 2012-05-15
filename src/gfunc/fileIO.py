"""
####################
fileIO.py
####################
Code supporting reading and writing from files not related to specific parsers.
"""

import sys
import os
import fnmatch


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