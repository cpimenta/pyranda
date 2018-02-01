# Run a bunch of test and check answers
import os,sys
import subprocess


def sexe(cmd,ret_output=False,echo = False):
    """ Helper for executing shell commands. """
    if echo:
        print "[exe: %s]" % cmd
    if ret_output:
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        res =p.communicate()[0]
        return p.returncode,res
    else:
        return subprocess.call(cmd,shell=True)


def baseDict(string):
    dbase = {}
    # Make dictionary
    for bb in string.split('\n'):
        if bb:
            name = bb.split('--')[0].strip()
            diff = bb.split('--')[1].strip()
            dbase[name] = diff
    return dbase
            

    

class testObj:

    def __init__(self,name):

        self.name = name
        self.parallel = False
        self.np = 1
        self.script = ''
        self.args = None