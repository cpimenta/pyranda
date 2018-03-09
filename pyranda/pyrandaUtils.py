import numpy
import re


scalar_key = ":"
vector_key = [':[',']:']
tensor_key = [':{','}:']

def fortran3d(form,sMap):
    """
    Formula translator: Take in a string and return 
    the mapped string for evaluation
    """
    keyS = scalar_key

    # Eliminate useless whitespace
    form = form.replace(' ','').strip()
    
   
    # Fill in form with variable data arrays (scalars)
    varList = findVar(form,'scalar')
    for ivar in varList:
        kvar = '%s%s%s' % (keyS,ivar,keyS)
        if kvar not in sMap:
            sMap[kvar] = 'self.variables["%s"].data'%ivar

    keyT = vector_key[:]
    vecList = findVar(form,'vector')
    for ivar in vecList:
        kvar = '%s%s%s' % (keyT[0],ivar,keyT[1])
        if kvar not in sMap:
            sMap[kvar] = 'self.variables["%s"].data'%ivar

    
    #
    for mvar in sMap:
        if "#arg#" in mvar:
            eMap = sMap[mvar]
            prefix = mvar.split('(')[0]
            try:
                args = re.findall(prefix+'\(.*?\)',form)[0].replace('(','').replace(')','').replace(prefix,'')
                form = form.replace(mvar.replace('#arg#',args) ,eMap.replace("#arg#",args))
            except:
                pass
        else:
            form = form.replace(mvar,sMap[mvar])
    #
    return form

def splitLines(aLongString):
    """
    Take in a long string and split into lines for parsing
    """
    ics = aLongString
    ic_lines = filter(None,ics.split('\n'))
    ic_lines = [el.replace(' ','') for el in ic_lines ]  # Comments work
    ic_lines = filter(None,ic_lines)
    ic_lines = [el for el in ic_lines if el.strip()[0] != '#']  # Comments work
    return ic_lines
    



def findVar(string,kind,unique=True):
    """
    Return a list of the variable names to replace
    """
    if kind == "scalar":
        varStr = findVarScalar(string)
    elif kind == "vector":
        varStr = findVarVecTen(string,vector_key)
    elif kind == "tensor":
        varStr = findVarVecTen(string,tensor_key)
    else:
        varStr = []
    #
    # list(set makes LHS our of order
    if unique:
        varStr = list(set(varStr))
    
    return varStr


def findVarScalar(string):

    keyS = scalar_key

    # Tags to exclude for scalars
    excluded_chars = "[]{}#%*@"

    svars = re.findall('%s.*?%s'%(keyS,keyS),string)
    nvars = []
    for nv in svars:
        skip = False
        for ev in excluded_chars:
            if ev in nv:
                skip = True
        if not skip:
            nvars.append(nv.replace(keyS,''))
    return nvars


def findVarVecTen(string,keyS):


    # Fix for regEq escapes
    exMap = {}
    exMap['['] = '\['
    exMap[']'] = '\]'
    exMap['{'] = '\{'
    exMap['}'] = '\}'

    keyL = keyS[0]
    keyR = keyS[1]
    for em in exMap:
        keyL = keyL.replace(em,exMap[em])
        keyR = keyR.replace(em,exMap[em]) 
    
    svars = re.findall('%s.*?%s'%(keyL,keyR),string)
    for isv in range(len(svars)):
        svars[isv] = svars[isv].replace(keyS[0],'')
        svars[isv] = svars[isv].replace(keyS[1],'')
    return svars
