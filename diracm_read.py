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
        self.character=[] #list of the orbital characters, with the orbital characters being a list of individual orbital components
        self.occ=[] #list of the orbital occupationm, with the orbital occupation being a list of individual orbital components
        # group evergything
        self.group=[self.num,self.sym,self.energy,self.occ_no,self.character,self.occ]
        # assign values
        self.group=list(input_lt)
    
    def print_options(self):
        '''Print all of the variables in an order, num -> sym -> energy -> occ_no -> character -> occ
           This function needs to have more functionality.
        '''
        for a1,a2,a3,a4,a5,a6 in zip(self.num,self.sym,self.energy,self.occ_no,self.character,self.occ):
            print a1,a2,a3,a4,a5,a6
def find_MOs(lines):
    """
    to find the start and end lines of a Mulliken orbital print in a molpro output, 
    return the start and end lines' indices, 
    return empty list if no MO has been found
    """
    indices_lt=[]
    for i in range(0,len(lines)):
        if lines[i].strip() in startofmo:
            indices_lt.append(i)
        if lines[i].strip() in endofmo:
            indices_lt.append(i)
            break
    return indices_lt

def extract_and_process_raw_MOs(txt):
    """
    based on the start and end line of the MOs output region, to exact the unprocessed MOs, read every(if recognisable) set of MOs,
    standardize the the raw MO data and create an instance of dirac_mulliken_orbitals class,
    return the instance
    """
    indices_lt=find_MOs(txt)
    instance=[]
    raw_MO=txt[indices_lt[0]-1:indices_lt[1]]
    MO_instance=process_one_raw_MOs(raw_MO)
    return MO_instance

def process_one_raw_MOs(txt):
    """
    to process one unprocessed raw Mulliken MOs, create the variable for the dirac_mulliken_orbitals class and return dirac_mulliken_orbitals class instance
    """
    lt_of_num = []
    lt_of_sym = []
    lt_of_energy = []
    lt_of_occ_no = []
    lt_of_charater = []
    lt_of_occ = []
    flag = 0 "This is used for charater and occ read being contineous."
    for line in txt:
        line = line.strip()
        if "* Electronic eigenvalue" in line:
            flag = flag + 1
            line_lt = line.split()
            lt_of_num.append(line_lt[line_lt.index('no.')+1].strip(':'))
            lt_of_sym.append(line_lt[line_lt.index('sym=')+1].strip())
            lt_of_energy.append(line_lt[line_lt.index('no.')+2].strip(':'))
            lt_of_occ_no.append(line_lt[line_lt.index('f')+2].strip(':'))
        if "Gross" and "|" in line:
            line_lt = filter(lambda x: x!='',line.split('|')[1].split('  '))
        if flag == 2:
            flag = 0
            lt_of_charater.append()
            lt_of_occ.append()
    return dirac_mulliken_orbitals(lt_of_num,lt_of_sym,lt_of_energy,lt_of_occ_no,lt_of_charater,lt_of_occ)

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

if len(sys.argv) > 2:
    thresh=sys.argv[2]
else:
    thresh = float(raw_input("threshold(recommended value:0.05)> "))  # This threshold is for orbital occupation, when say 5f-2 occupation is greater than the threshold, it will be printed.

# Reprocessing from raw data to organised standard output
instance_MOs_lt=extract_and_process_raw_MOs(lines)    
instance_MOs_lt.print_options()
