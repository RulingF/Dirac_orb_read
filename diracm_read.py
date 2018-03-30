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
        self.num=input_lt[0] #list of the number of orbitals
        self.sym=input_lt[1] #list of the symmetries of orbitals
        self.energy=input_lt[2] #list of the energies of orbitals
        self.occ_no=input_lt[3] #list of the occupation number
        self.character=input_lt[4] #list of the orbital characters, with the orbital characters being a list of individual orbital components
        self.occ_alpha=input_lt[5] #list of the orbital occupationm, with the orbital occupation being a list of individual orbital components
        self.occ_beta=input_lt[6]
        # group evergything
        self.group=[self.num,self.sym,self.energy,self.occ_no,self.character,self.occ_alpha,self.occ_beta]
   
    def print_all_orbitals(self):
        '''Print all of the variables in an order, num -> sym -> energy -> occ_no -> character -> occ
           This function needs to have more functionality.
        ''' 
        for num,sym,energy,occ_no,character_lt,occ_alpha_lt,occ_beta_lt in zip(self.num,self.sym,self.energy,self.occ_no,self.character,self.occ_alpha,self.occ_beta):
            print num, sym, energy, occ_no
            print character_lt
            print occ_alpha_lt
            print occ_beta_lt
            print '\n\n'
    def print_open_shells(self,thresh):
        '''Print all of the variables of openshell orbitals in an order, num -> sym -> energy -> occ_no -> character -> occ
        '''
        for num,sym,energy,occ_no,character_lt,occ_alpha_lt,occ_beta_lt in zip(self.num,self.sym,self.energy,self.occ_no,self.character,self.occ_alpha,self.occ_beta):
            if float(occ_no) < 1 and float(occ_no) > 0:
                print "This is orbital no.%s with orbital energy: %s, and the symmetry is %s with occupation of %s." %(num,energy,sym,occ_no)
                pt_lt = []
                print "%-10s :" %('Character'),
                for i in xrange(len(occ_alpha_lt)):
                    if float(occ_alpha_lt[i]) > thresh:
                        pt_lt.append(i)
                for i in xrange(len(occ_beta_lt)):
                    if float(occ_beta_lt[i]) > thresh:
                        pt_lt.append(i)
                for i in xrange(len(character_lt)):
                    if i in pt_lt:
                        print "%-15s  " %(character_lt[i]),
                print ''
                print "%-10s :" %('Alpha Occ'),
                for i in xrange(len(occ_alpha_lt)):
                    if i in pt_lt:
                        print "%-15s  " %(occ_alpha_lt[i]),
                print ''
                print "%-10s :" %('Beta Occ'),
                for i in xrange(len(occ_alpha_lt)):
                    if i in pt_lt:
                        print "%-15s  " %(occ_beta_lt[i]),
                print '\n'
    def print_virtual_shells(self,thresh,no_of_virtual = 0):
        '''Print all of the variables of virtual(unoccupied) orbitals in an order, num -> sym -> energy -> occ_no -> character -> occ
        '''
        count = 0
        for num,sym,energy,occ_no,character_lt,occ_alpha_lt,occ_beta_lt in zip(self.num,self.sym,self.energy,self.occ_no,self.character,self.occ_alpha,self.occ_beta):
            if float(occ_no) == 0 and not (count > no_of_virtual):
                count = count + 1
                print "This is orbital no.%s with orbital energy: %s, and the symmetry is %s with occupation of %s." %(num,energy,sym,occ_no)
                pt_lt = []
                print "%-10s :" %('Character'),
                for i in xrange(len(occ_alpha_lt)):
                    if float(occ_alpha_lt[i]) > thresh:
                        pt_lt.append(i) 
                for i in xrange(len(occ_beta_lt)):
                    if float(occ_beta_lt[i]) > thresh:
                        pt_lt.append(i)
                for i in xrange(len(character_lt)):
                    if i in pt_lt:
                        print "%-15s  " %(character_lt[i]),
                print ''
                print "%-10s :" %('Alpha Occ'),
                for i in xrange(len(occ_alpha_lt)):
                    if i in pt_lt:
                        print "%-15s  " %(occ_alpha_lt[i]),
                print ''
                print "%-10s :" %('Beta Occ'),
                for i in xrange(len(occ_alpha_lt)):
                    if i in pt_lt:
                        print "%-15s  " %(occ_beta_lt[i]),
                print '\n'

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
    raw_MO = txt[indices_lt[0]-1:indices_lt[1]+1]
    MO_instance = process_one_raw_MOs(raw_MO)
    return MO_instance

def process_one_raw_MOs(txt):
    """
    to process one unprocessed raw Mulliken MOs, create the variable for the dirac_mulliken_orbitals class and return dirac_mulliken_orbitals class instance
    """
    lt_of_num = []
    lt_of_sym = []
    lt_of_energy = []
    lt_of_occ_no = []
    lt_of_character = []
    lt_of_occ_alpha = []
    lt_of_occ_beta = []
    flag = 0 #Read line control, when flag equals 2, append the charater and occ data.
    for line in txt:
        line = line.replace('no.','no. ')
        line = line.strip()
        if "* Electronic eigenvalue" in line:
            flag = flag + 1
            line_lt = line.split()
            lt_of_num.append(line_lt[line_lt.index('no.')+1].strip(':'))
            try:
                lt_of_sym.append(line_lt[line_lt.index('sym=')+1].strip())
            except ValueError:
                lt_of_sym.append(line_lt[line_lt.index('m_j=')+1].strip())
            lt_of_energy.append(line_lt[line_lt.index('no.')+2].strip(':'))
            lt_of_occ_no.append(line_lt[line_lt.index('f')+2].strip(')'))
        if "Gross" in line and "|" in line:
            line_lt = filter(lambda x: x!='',line.split('|')[1].split('  '))
            try:
                tmp_character_lt = tmp_character_lt + line_lt
            except UnboundLocalError:
                tmp_character_lt = []
        if "alpha" in line and "|" in line:
            line_lt = filter(lambda x: x!='',line.split('|')[1].split('  '))
            try:
                tmp_occ_alpha_lt = tmp_occ_alpha_lt + line_lt
            except UnboundLocalError:
                tmp_occ_alpha_lt = []
        if "beta" in line and "|" in line:
            line_lt = filter(lambda x: x!='',line.split('|')[1].split('  '))
            try:
                tmp_occ_beta_lt = tmp_occ_beta_lt + line_lt
            except UnboundLocalError:
                tmp_occ_beta_lt = []
        if flag == 2 or 'Total gross population' in line:
            flag = flag - 1
            lt_of_character.append(tmp_character_lt)
            lt_of_occ_alpha.append(tmp_occ_alpha_lt)
            lt_of_occ_beta.append(tmp_occ_beta_lt)
            tmp_character_lt = []
            tmp_occ_alpha_lt = []
            tmp_occ_beta_lt = []
    return dirac_mulliken_orbitals(lt_of_num,lt_of_sym,lt_of_energy,lt_of_occ_no,lt_of_character,lt_of_occ_alpha,lt_of_occ_beta)

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

if len(sys.argv) > 3:
    no_of_virtuals = sys.argv[3]
else:
    no_of_virtuals = float(raw_input("number of virtual orbitals to print(default value is 0)> "))  #This is for virtual orbitals print option, if it is not given , then 0 is the default value 

print "This is only for openshell orbitals in Mulliken analysis, if you wish to print other orbitals(say,occupied ones or virtuals), pls contact me, frel.feng@wsu.edu"
# Reprocessing from raw data to organised standard output
instance_MOs_lt=extract_and_process_raw_MOs(lines)    
#instance_MOs_lt.print_all_orbitals()
instance_MOs_lt.print_open_shells(thresh)
instance_MOs_lt.print_virtual_shells(thresh,no_of_virtuals)
