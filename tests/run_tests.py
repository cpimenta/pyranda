# Run a bunch of test and check answers
import os,sys
import numpy as npy
import subprocess
from testObj import *


test_dir = os.path.dirname(os.path.abspath(__file__))
bin_dir  = os.path.join( test_dir, '../bin')
root_dir = os.path.join( test_dir, '..')
pyranda_exe = os.path.join( bin_dir, 'pyranda')
pyranda_mpi = os.path.join( bin_dir, 'pympirun')



tests = []    # List of test objects
dbase = {}    # Dictionary of baselines

# Add tests here
execfile('test1DAdvection.py')


# Run tests
for test in tests:

    script = os.path.join(root_dir,test.script)

    exe = pyranda_exe
    if test.parallel:
        exe = pyranda_mpi + ' -n %s %s' % (test.np,exe)
    
    # Args
    sargs = ''
    for arg in test.args:
        sargs += '%s ' % arg

        
    cmd = '%s %s %s' % (exe,script,sargs)

    out = sexe(cmd,ret_output=True,echo=False)
    pout = out[1].split('\n')[-2]

    # Diff against baseline
    try:
        baseline = float(dbase[test.name])
        diff = npy.abs( float(pout)-baseline ) / npy.abs( baseline )
        if diff < 1.0e-4:
            testPass = True
            print 'Pass: (Rel. Error = %s )' % diff
            print '%s -- %s' % (test.name,pout)
        else:
            testPass = False
            print 'Fail: (Rel. Error = %s )' % diff
            print '%s -- %s' % (test.name,pout)

    except:
        testPass = False
        print 'Fail: (No baseline data found )'
        print '%s -- %s' % (test.name,pout)
            

    
        
