#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, shutil, tempfile
import subprocess
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Process, freeze_support
import multiprocessing

#nMarkersList = [100, 200, 400, 800, 1600, 3200, 4600]
nMarkersList       = [5, 10, 25, 50, 75, 100]
nMarkersDirs       = ['%.2d' % x for x in nMarkersList]
inputsDir          = 'inputs'
runsDir            = 'runs'
paramsTemplateFile = 'template_nBestMarkers.yml'

#--------------------------------------------------------------------
#--------------------------------------------------------------------
def main ():
	global runsDir
	splitDatasets ()

	createDir (runsDir)

	createRunningDirs (runsDir, nMarkersDirs)

	runMultipleGSSingle (runsDir, nMarkersDirs)

#--------------------------------------------------------------------
#--------------------------------------------------------------------
def splitDatasets ():
	runScript ('run-split.R')

#--------------------------------------------------------------------
#--------------------------------------------------------------------
def createRunningDirs (runsDir, nMarkersDirs):
	for nDir in nMarkersDirs:
		ymlLines = open (paramsTemplateFile).readlines ()
		for i,line in enumerate (ymlLines):
			if  line.strip().startswith ("#"):
				ymlLines [i] = ''
				continue
			if 'nMarkers' in line:
				ymlLines [i] = "nMarkers			: %s\n" % int (nDir)
			if '.csv' in line:
				item	     = line.split (':')[0]
				filename     = line.split (':')[1].strip()
				dataDir      = os.path.join (os.getcwd(), 'datasets', filename)
				ymlLines [i] = f"{item}\t: {dataDir}\n"
				if 'markersFile' in line:
					shutil.copy (os.path.join ('inputs', filename), 'datasets')


		runningDir = os.path.join (runsDir, nDir)
		ymlFile    = os.path.join (runningDir, f'params-nMarkers.yml')
		os.makedirs (runningDir)
		with open (ymlFile, 'w') as fp:
			fp.writelines (ymlLines)

		shutil.copy ("script_gs.R", runningDir)

#--------------------------------------------------------------------
def runMultipleGSSingle (inputDir, nMarkersDirs):
	"""Runs an R script with an argument and captures output."""
	nMarkersPaths = [os.path.join (inputDir, x) for x in nMarkersDirs]
	for path in nMarkersPaths:
		run_r_script (path)

#--------------------------------------------------------------------
def runMultipleGSParallel (inputDir, nMarkersDirs):
	"""Runs an R script with an argument and captures output."""
	nMarkersPaths = [os.path.join (inputDir, x) for x in nMarkersDirs]

	num_workers = os.cpu_count() - 1  # Get the total logical processors	
	with ProcessPoolExecutor (max_workers=1) as executor:
		results = list (executor.map (run_r_script, nMarkersPaths))

	#print(results)

#--------------------------------------------------------------------
#--------------------------------------------------------------------
def run_r_script (runningPath):
	"""Runs an R script from a specific directory with an argument."""
	result = subprocess.run (
		["Rscript", 'script_gs.R'],
		stdout=subprocess.PIPE, 
		stderr=subprocess.PIPE, 
		text=True,
		cwd=runningPath	# Change this to your script's directory
	)
	print("STDOUT:", result.stdout)
	print("STDERR:", result.stderr)  # Check this for errors

	return result.stdout.strip()

#------------------------------------------------------------------
#------------------------------------------------------------------
def runScript (script, runningPath=None):
	"""Runs an R script from a specific directory with an argument."""
	result = subprocess.run (
		["Rscript", script],
		capture_output=True,
		text=True,
		cwd=runningPath	# Change this to your script's directory
	)
	return result.stdout.strip()



#------------------------------------------------------------------
# Create a dir and rename the old one if it exists
#------------------------------------------------------------------
def createDir (dir):
	def checkExistingDir (dir):
		if os.path.lexists (dir):
			headDir, tailDir = os.path.split (dir)
			oldDir = os.path.join (headDir, "old-" + tailDir)
			if os.path.lexists (oldDir):
				checkExistingDir (oldDir)

			os.rename (dir, oldDir)
	checkExistingDir (dir)
	os.system ("mkdir %s" % dir)
#--------------------------------------------------------------------
#--------------------------------------------------------------------
if __name__ == '__main__':
	#multiprocessing.set_start_method("spawn", force=True)  # ✅ Ensure Windows uses spawn
	#multiprocessing.freeze_support()  # ✅ Required for Windows multiprocessing
   
	main ()
