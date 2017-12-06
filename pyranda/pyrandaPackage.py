import numpy 
import re
import sys
import time



class pyrandaPackage:
    """
    Case physics package module for adding new physics packages to pyranda
    """
    def __init__(self,name,pyranda):

        self.name = name
        self.pyranda = pyranda


    def get_sMap(self):
        """
        String mappings for this package.  Packages added to the main
        pyranda object will check this map
        """
        sMap = {}
        self.sMap = sMap

        

        
