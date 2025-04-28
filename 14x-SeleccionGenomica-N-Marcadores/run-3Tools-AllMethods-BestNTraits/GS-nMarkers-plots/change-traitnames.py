#!/usr/bin/env python3

import os

traitNames = open ('traitnames-BASE.csv').readlines ()
filename   = 'out-GS-crossval-varying-nmarkers.csv'

for line in traitNames:
	values = line.split (",")
	spanish, english, description = values [0].strip(), values [1].strip(), values [2].strip()
	print (f"{spanish}, {english}, {description}")

	cmm = f'sed -i "s/{spanish}/{english}/g" {filename} '
	print (cmm, "\n")

	os.system (cmm)
