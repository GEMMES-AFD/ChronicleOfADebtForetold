#define working directory where all files can be found
wd = "/home/rabbimilligan/PhD/Turkey/CALIBRATION/CodesCalibration_New/SourceCode"


# Source utility functions
execfile(wd + "/source.py")
#Note: execfile is Python's equivalent to the source() function in R


# Run RK4 using all default values defined from the .R equations file
# note that there is a minor bug regarding time: python uses a time sequence of the same duration (in years) and with the same time step as the one specified in R, but that will always be initialized at t=0
# I will correct it later.


# Compile and Import C++ library
solvePy = cppimport.imp('functionsForPy')

outEuler = solve(solver="euler")
outRK4Fixed = solve(solver="RK4Fixed")
outDopri = solve(solver="dopri", atol=1e-4, rtol=0, fac=0.85, facMin=0.1, facMax=4, nStepMax=300, hInit=0.025, hMin=0.025/100, hMax=0.2)

# Plot 3 variables (chosen at random)
outEuler.plot(y=["w"], use_index=True)
outRK4Fixed.plot(y=["w"], use_index=True)
outDopri.plot(y=["w"], use_index=True)


###############################
## ACCESSING ELEMENTS IN OUT ##
###############################
## out is a dataFrame (from the library pandas). 
# rows are dates, columns are variables
# The main functions to access elements of a dataFrame are: 
#ACCESS A COLUMN USING ITS NAME: 
outRK4Fixed["w"]
#ACCESS ROWS:
#You can call a row using either its name in single brackets:
outRK4Fixed.loc[0]
outRK4Fixed.loc[1]
# or using its index with double brackets:
outRK4Fixed.loc[[0]]
outRK4Fixed.loc[[1]]
# (note that, here, name=index...)
# Alternative version without .loc (but weird syntax and does not work with double brackets): 
outRK4Fixed[1:2]
#You can also call multiple rows: 
outRK4Fixed.loc[0:1]
outRK4Fixed.loc[[0,1,2,3,4,5]]
outRK4Fixed[0:5]
# get n last/first rows (works with n=1):
n=2
outRK4Fixed.tail(n)
outRK4Fixed.head (n)

#####################################
## SPECIFY TIME SEQUENCE IN PYTHON ##
#####################################
# Using np.arange (Python equivalent to the R seq(from, to, by) function 
#However, contrarily to R, python excludes last point of an interval, 
# so if you want a time sequence from 2020 to 2050, you need to give him the interval 2020:2050.01 or it will stop at 2049.99
time = np.arange(2009, 2021.01, 0.01) 
out = solve(time=time, solver="dopri")
out.plot(y=["w"], use_index=True)
#Here you get proper date on x axis


######################################
## CHANGE PARAMETER VALUE IN PYTHON ##
######################################
# I will make built-in functions to ease the process. But for now, here is how it can be done: 

# Build default parms vector in python importing data from C++
parmsNames = solvePy.parmsNames()
# Declare a structure that stores parameters values and their names (named tuple)
namedParms = namedtuple("namedParms", parmsNames)

# build named parms vector using the parmsNamed structure and loading parms values from C++
newParms = namedParms(*solvePy.parms())

# You now can access parameters values with parms.nameOfParameter:
newParms.alpha


# edit parms value:
newParms = newParms._replace(interest=0.05)

#Call RK4 and pass new parms vector:
out = solve(parms=newParms, solver="dopri")
out.plot(y=["w", "e", "d"])


#####################################
## CHANGE INITIAL VALUES IN PYTHON ##
#####################################

# This is similar to changing parameters value
varNames = solvePy.varNames()
namedVars = namedtuple("namedVars", varNames)
newY0 = namedVars(*solvePy.yInit())
newY0 = newY0._replace(d=1)

out = solve(y0=newY0, solver="dopri")
out.plot(y=["w", "e", "d"], use_index=True)
