###############################################################################
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# All rights reserved.
#
# This file is part of pyranda
#
# 
###############################################################################

"""
 file: install.py

 description: uses spack to install the external third party libs used by a project.

"""

import os
import sys
import subprocess
import shutil
import socket
import platform
import json
import datetime

from optparse import OptionParser
from os import environ as env
from os.path import join as pjoin


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


# Download spack into spack_libraries directory
spack_repo = 'https://github.com/flow-phys/spack.git'   # Britton's fork for now.
cmd = 'git clone %s spack_library' % spack_repo
print "Cloning Spack source for building packages..."
#out = sexe( cmd )
print "Done"


# Build python with float-py support
print "Installing spack libraries ..."
cd_spack = 'cd spack_library/bin;'
spack_bld = cd_spack + './spack install py-floatpy'
#out = sexe( spack_bld )


# Activate python modules
py_mods = ['py-numpy','py-matplotlib','py-floatpy','py-f90wrap','py-mpi4py']
for pym in py_mods:
    sexe(cd_spack + './spack activate %s' % pym)


# Get python version number
py_version = cd_spack + './spack find python | grep python'
out = sexe( py_version , ret_output=True )
Major = out[1].strip().split('@')[1].split('.')[0] 
Minor = out[1].strip().split('@')[1].split('.')[1]
Step  = out[1].strip().split('@')[1].split('.')[2]
pyver = 'python%s.%s' % (Major,Minor)
pyspk = 'python@%s.%s.%s' % (Major,Minor,Step)


# Install pyranda symbolic links into site-packages
print 'Installing pyranda library links to python'
pyloc = sexe( cd_spack + './spack location --install-dir %s' % (pyspk) , ret_output=True)[1].strip()

# Make directory in site-packages
pyPackagesLocation = '%s/lib/%s/site-packages' % (pyloc,pyver)

if not os.path.isdir( '%s/pyranda' % pyPackagesLocation):
    cmd_mkdir_pyranda = 'mkdir %s/pyranda' % pyPackagesLocation
    sexe( cmd_mkdir_pyranda )

cmd_ln_pyrandaLib = 'ln -sf %s/pyranda/* %s/pyranda/' % (os.getcwd() , pyPackagesLocation)
sexe( cmd_ln_pyrandaLib )


# Make link to exec for convenience
print 'Installing pyranda to bin/pyranda'
if not os.path.isdir( 'bin'):
    sexe( 'mkdir bin')

cmd_ln_pyranda = 'ln -sf %s/bin/%s bin/pyranda' % ( pyloc, pyver )
sexe( cmd_ln_pyranda )





