#!/usr/bin/env python

import getopt, os, ast, sys
from mtVariantCaller import VCFoutput

def usage():
	print """Produces VCF file output from VCF_dict_tmp file.
		Version 1.1 - Written by Domenico Simone and Claudia Calabrese - 2013-2014
		
		Options:
		-r		Reference sequence [RSRS|RCRS]
		-i 		input directory
		"""

reference_sequence="RSRS"
in_dir="."

try:
	opts, args = getopt.getopt(sys.argv[1:], "h:r:i:")
except getopt.GetoptError, err:
	print str(err)
	usage()
	sys.exit()

for o,a in opts:
	if o == "-r":
		if a in ('RCRS', 'RSRS'):
			reference_sequence = a
		else:
			print "Reference sequence must be RSRS or RCRS."
			sys.exit()
	if o == "-i":
		in_dir = a

try:
	VCF_dict = ast.literal_eval(open(os.path.join(in_dir, 'VCF_dict_tmp'), 'r').read())
except IOError, err:
	print "VCF_dict_tmp file doesn't exist!"
	sys.exit()
VCFoutput(VCF_dict, reference=reference_sequence)