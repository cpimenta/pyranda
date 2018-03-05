import numpy

            
def fortran3d(form,sMap):
    """
    Formula translator: Take in a string and return 
    the mapped string for evaluation
    """
    keyS = ':'
   
    # Fill in form with variable data arrays
    varList = findVar(form,keyS)
    for ivar in varList:
        kvar = '%s%s%s' % (keyS,ivar,keyS)
        if kvar not in sMap:
            sMap[kvar] = 'self.variables["%s"].data'%ivar
    #
    for mvar in sMap:
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
    



def findVar(string,keyS,unique=True):
    """
    Return a list of the variable names to replace
    """
    # check that there are even numbers
    cnt = string.count(keyS)
    #
    if (cnt % 2) != 0:
        print "Error: unclosed variable definition marker ::"
        raise ValueError(string)
    #
    varStr = []
    s = ''
    inc = False
    for ss in string:
        if ss == keyS:
            inc = not inc
            if not inc:
                varStr.append(s)
                s = ''
            else:
                continue
        if inc:
            s += ss
    #
    # list(set makes LHS our of order
    if unique:
        varStr = list(set(varStr))
    
    return varStr
