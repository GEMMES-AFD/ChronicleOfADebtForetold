import numpy as np
import pandas as pd
from collections import namedtuple

import sys

import cppimport
cppimport.settings['force_rebuild'] = True
# Tell Python where to find the C++ codes
fullwd = wd + '/cppCode'
# sys.path.append(fullwd)
sys.path.insert(0, fullwd)
sys.path.insert(0, fullwd + "/src")
    
def solve(time=None, y0=None, parms=None, samplesExogVar=None, nSamplesVarExogVar=None, nVarExogVar=None, solver="dopri", 
          atol=1e-4, rtol=0, fac=0.85, facMin=0.1, facMax=4, nStepMax=100, hInit = 0.01, hMin=0.0001, hMax=0.5):
    
    #INITIALIZATION
    if time is None:
        nt = solvePy.nt()
        tInit = solvePy.tInit()
        tEnd = solvePy.tEnd()
    else: 
        nt = time.size
        tInit = time[0]
        tEnd=time[time.size-1]

    if hMax>(tEnd-tInit)/(nt-1):
        hMax = (tEnd-tInit)/(nt-1)
    

    if y0 is None:
        y0 = solvePy.yInit()
    if parms is None:
        parms = solvePy.parms()
    if samplesExogVar is None:
        samplesExogVar = solvePy.samplesExogVar()
    if nSamplesVarExogVar is None:
        nSamplesVarExogVar = solvePy.nSamplesVarExogVar()
    if nVarExogVar is None:
        nVarExogVar = solvePy.nVarExogVar()
        
        
    nV = solvePy.nV()
    nIV = solvePy.nIV()
    returnRK4 = solvePy.returnType()
    
    # CALL C++ CODE FOR RK4
    if solver=="dopri":
        outCpp = solvePy.dopri(nt, tInit, tEnd, nV, nIV, y0, parms, samplesExogVar, nSamplesVarExogVar, nVarExogVar,
                               atol, rtol, fac, facMin, facMax, nStepMax, hInit, hMin, hMax)    
    elif solver=="RK4Fixed":
        outCpp = solvePy.RK4Fixed(nt, tInit, tEnd, nV, nIV, y0, parms, samplesExogVar, nSamplesVarExogVar, nVarExogVar)    
    else: 
        outCpp = solvePy.euler(nt, tInit, tEnd, nV, nIV, y0, parms, samplesExogVar, nSamplesVarExogVar, nVarExogVar)    
    
    # FORMAT OUTPTUT
    varNames = solvePy.varNames()
    intermediateVarNames = solvePy.intermediateVarNames()
    derivVarNames = [varName+ "Dot" for varName in varNames]        
    
    
    if(returnRK4 == 3):
        outSize = 2*nV+nIV
        outNames = varNames + derivVarNames + intermediateVarNames
    elif (returnRK4 == 2):
        outSize = nV+nIV
        outNames = varNames + intermediateVarNames
    elif (returnRK4 == 1):
        outSize = 2*nV
        outNames = varNames + derivVarNames
    else:
        outSize = nV
        outNames =varNames

    outUnNamed = np.zeros((nt,outSize))
    for i in range(nt):
        outUnNamed[i,:] = outCpp[(i*(outSize)):(i*(outSize)+(outSize))]
    out = pd.DataFrame(outUnNamed, columns=outNames, index=time)
    return out