#!/usr/bin/env python3

import os

traitNames = open ('inputs/traitnames-BASE.csv').readlines ()
filename   = 'outputs/out-GS-crossval-varying-nmarkers.csv'
spFilename = filename.split (".")[0] + "-SPANISH.csv"
os.system ("cp %s %s" % (filename, spFilename))

for line in traitNames:
	values = line.split (",")
	spanish, english, description = values [0].strip(), values [1].strip(), values [2].strip()
	print (f"{spanish}, {english}, {description}")

	cmm = f'sed -i "s/{spanish}/{english}/g" {filename} '
	print (cmm, "\n")

	os.system (cmm)
