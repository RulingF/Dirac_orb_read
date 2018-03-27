#!/usr/bin/env python
#
#
# Written by Rulin Feng, Washington State University, Kirk Peterson's group,
#         Washington State University, 2018
#
# Contributors:  Rulin Feng
#
#
# Purpose: Given one single dirac output file with mulliken orbitals printed out,
#          to offer an easier readability by reorganising mulliken outputs,
#          only parse openshell electrons in mulliken analysis
#          No more no less.
#
##############################################################################################################################

# set up the flags of start and end lines of in openshell electrons in dirac mulliken analysis
global startofmo, endofmo
startofmo=['********************** Mulliken population analysis **********************']
endofmo=['*** Total gross population ***']

class dirac_mulliken_orbitals:
    """
    class for mulliken analysis occupations from Dirac output
    """

    def __init__(self,*input_lt):
        '''Initialize the set of variables
        '''

        # set up the features of data that are being stored, can be better understood by comparing to a DIRAC Mulliken print
        self.num=[] #list of the number of orbitals
        self.sym=[] #list of the symmetries of orbitals
        self.energy=[] #list of the energies of orbitals
        self.occ_no=[] #list of the occupation number
        self.character=[] #list of the orbital characters
        self.occ=[] #list of the orbital occupation
        # group evergything
        self.group=[self.num,self.sym,self.energy,self.occ_no,self.character,self.occ]
        # assign values
        self.group=list(input_lt)

def extract_and_process_raw_MOs(txt):
    """
    based on the start and end line of the MOs output region, to exact the unprocessed MOs, read every(if recognisable) set of MOs,
    standardize the the raw MO data and create an instance of dirac_mulliken_orbitals class,
    return the instance
    """
 
# end of the class and function definitions
##############################################################################################################################
import sys

# file reading and threshhold control
if len(sys.argv) > 1:
    txt=sys.argv[1]
else:
    txt= raw_input("inpfile> ")
txt=open(txt)
lines=txt.readlines()

processed_MOs_lt=extract_and_process_raw_MOs(lines)

if len(sys.argv) > 2:
    thresh=sys.argv[2]
else:
    thresh = float(raw_input("threshhold(recommended value:0.05)> "))

# Reprocessing from raw data to organised standard output

print "There are " + str(len(processed_MOs_lt)) + " MO orbitals loaded..."
choose= raw_input(" MOs to print\n Example --- second and fifth : 2 , 5\n         --- all of the orbitals : all or leave blank\n> ")
tmplist=[str(k+1) for k in xrange(len(processed_MOs_lt))]
if len(choose.strip())==0 or 'all' in choose:
    for item in processed_MOs_lt:
        if processed_MOs_lt.index(item)!=0:
            print '\n\n'
            print '   Next set of MOs'
        item.print_MOs(thresh)
else:        
    for k in xrange(len(processed_MOs_lt)):
        if input_error(choose,tmplist):
            print "Check your choice of printed orbitals, there might be an error ...\nExit ..."
            sys.exit()
        elif str(k+1) in choose.split(','):
            processed_MOs_lt[k].print_MOs(thresh)

