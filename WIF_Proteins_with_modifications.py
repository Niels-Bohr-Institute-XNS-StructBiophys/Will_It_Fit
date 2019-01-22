import os
import fnmatch
import sys
import subprocess
import shutil


# compilation of C-code:
# gcc WillItFit.c -O3 -fopenmp -g -Wall -o WillItFit.com -lm -lgsl -lgslcblas
#
# run python script:
# python WIF_Proteins_with_modifications.py
#
# example C-executable call generated by python script
# ./C/WillItFit.com -c=Data/01/R3T3T3R3-crysol.card -d=Data/01/R3T3T3R3_w.pdb -s=Data/01/Protein.dat -p=Data/01/Protein_free.par -e=N/A -t=0 -r=1 -i=2 -j=0 -k=0 -o=1 -h=0.0 -n=0.0 -x=1.0 -z=1 -l=2


InitialDirectory = os.path.abspath(os.path.dirname(sys.argv[0]))
os.chdir(InitialDirectory)

print 'InitialDirectory = ' + InitialDirectory

print ''



# Q-range in 1/\AA applied for all files
MinQStr = '0.0'
MaxQStr = '1.0'

# default arguments for FittingRoutine
#
# FittingRoutineArgument1 = 50
# FittingRoutineArgument2 = 10
# FittingRoutineArgument3 = 32
#
# 0 ComputeModel	-
#
# 1 LevenbergMarquardt	FittingRoutineArgument1
#			(MaxIterations in LM)
#
# 2 GridsearchLM	FittingRoutineArgument1		FittingRoutineArgument2
#			(MaxIterations in LM)		(NumberOfCycles in Gridsearch)
#
# 3 BFGS		FittingRoutineArgument1
#			(MaxIterations in BFGS)
#
# 4 GridsearchBFGS	FittingRoutineArgument1		FittingRoutineArgument2
#			(MaxIterations in BFGS)		(NumberOfCycles in Gridsearch)
#
# 5 Swarm		FittingRoutineArgument1		FittingRoutineArgument2		FittingRoutineArgument3
#			???				???				???
#
# 6 Genetic		-				FittingRoutineArgument2		FittingRoutineArgument3
#							???				???
#
FittingRoutine = 1
FittingRoutineArgument1 = 15
FittingRoutineArgument2 = 3	# only in use for FittingRoutine == 2, 4, 5, 6 
FittingRoutineArgument3 = 0	# only in use for FittingRoutine == 5, 6

ResolutionFile = 'N/A'
NumberOfSmearingFolds = 0

PrintCovarianceMatrix = 1 # boolean flag 0 or 1
ChiSquareFractile = 0.0

CMD = 1 # do not change


# logging / reporting
#
# WriteLog  < 0 output printed to terminal
# WriteLog == 0 nothing logged
# WriteLog  > 0 output printed to logfile.log in results -folder
#
# abs(WriteLog) == 0 -> no output at all
#
# abs(WriteLog) == 1 -> most important output written, includes:
#			-parameter file
#			-excerpt of spectra, ProteinStructure
#			-ChiSquare and parameters at the fitting iteration steps
#
# abs(WriteLog) == 2 -> full output written, includes in addition to 1:
#			full output of Spectra, ProteinStructure
#			output from fitting algorithm (currently Levenberg-Marquardt and Compute Model supported)
#
# default is -1
WriteLog = 1


ExecutableFile = os.path.join( InitialDirectory, "C/WillItFit.com")
print 'ExecutableFile = ' + ExecutableFile

for subdir in os.listdir(InitialDirectory+'/Data/'):


	print ''
	print ''

	subdir = 'Data/' + subdir + '/'
	print 'subdir = ' + subdir

	DataFile = fnmatch.filter(os.listdir( os.path.join( InitialDirectory, subdir) ), '*.card')
	if not DataFile:
		print('No data *.card file found in directory ' + subdir + '. Continue.')
		continue

	DataFile = subdir + DataFile[0]
	print DataFile

	# print os.path.join( InitialDirectory, subdir) + '*.card-results'
	OldResultsFolder = os.path.join( InitialDirectory, DataFile) + '-results'
	if os.path.isdir( OldResultsFolder ) :
		shutil.rmtree( OldResultsFolder )

	# # rm *LinesNotRead.pdb from potential previous runs
	# PDBFile = fnmatch.filter(os.listdir( os.path.join( InitialDirectory, subdir) ), '*_LinesNotRead.pdb')
	# if PDBFile:
	# 	os.remove( os.path.join( InitialDirectory, subdir) + PDBFile[0])

	PDBFile = fnmatch.filter(os.listdir( os.path.join( InitialDirectory, subdir) ), '*.pdb')
	if not PDBFile:
		print('No protein *.pdb file found in directory ' + subdir + '. Continue.')
		continue

	PDBFile = subdir + PDBFile[0]
	print PDBFile

	SampleFile = fnmatch.filter(os.listdir( os.path.join( InitialDirectory, subdir) ), '*.dat')
	if not DataFile:
		print('No sample *.dat file found in directory ' + subdir + '. Continue.')
		continue

	SampleFile = subdir + SampleFile[0]
	print SampleFile

	ParameterFile = fnmatch.filter(os.listdir( os.path.join( InitialDirectory, subdir) ), '*.par')
	if not ParameterFile:
		print('No parameter *.par file found in directory ' + subdir + '. Continue.')
		continue

	ParameterFile = subdir + ParameterFile[0]
	print ParameterFile


	ProcessToCall = []
	ProcessToCall.append(ExecutableFile)

	ProcessToCall.append('-c=%s' % DataFile)
	ProcessToCall.append('-s=%s' % SampleFile)
	ProcessToCall.append('-p=%s' % ParameterFile)
	ProcessToCall.append('-d=%s' % PDBFile)

	ProcessToCall.append('-n=%s' % MinQStr)
	ProcessToCall.append('-x=%s' % MaxQStr)

	ProcessToCall.append('-r=%d' % FittingRoutine)
	ProcessToCall.append('-i=%d' % FittingRoutineArgument1)
	ProcessToCall.append('-j=%d' % FittingRoutineArgument2)
	ProcessToCall.append('-k=%d' % FittingRoutineArgument3)

	ProcessToCall.append('-e=%s' % ResolutionFile)
	ProcessToCall.append('-t=%d' % NumberOfSmearingFolds)

	ProcessToCall.append('-o=%d' % PrintCovarianceMatrix)
	ProcessToCall.append('-h=%s' % ChiSquareFractile)

	ProcessToCall.append('-z=%d' % CMD)
	ProcessToCall.append('-l=%d' % WriteLog)

	print ""
	print ProcessToCall
	print ""

	Process = subprocess.Popen(ProcessToCall)
	out, err = Process.communicate()
	print out, err
