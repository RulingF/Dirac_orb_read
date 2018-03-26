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
startofmo=['* Electronic eigenvalue no.']
endofmo=['* Electronic eigenvalue no.']

class dirac_mulliken:
    """
    class for mulliken analysis occupations from Dirac output
    """

    def __init__(self,input_lt):
        '''Initialize the set of variables
        '''

        # set up the features of data that are being stored, can be better understood by comparing to a molpro MOs print
        self.orb=[]
        self.occ=[]
        self.orbital_energy=[]
        self.configuration_composition=[] #this is list of dicts of orbital name being the key and coefficients being the value
        for item in input_lt[0]:
            self.orb.append(item[0])
            self.occ.append(item[1])
            self.orbital_energy.append(item[2])
            self.configuration_composition.append(item[3])
        # other information about the set of MOs
        self.other_information=input_lt[1]
        # group evergything
        self.group=[self.orb,self.occ,self.orbital_energy,self.configuration_composition]
    
    def print_MOs(self,threshhold,isym=-1,jsym=-1):
        '''Print out the orbitals,occupations,orbital energies and compostitions,
        the threshhold only controls the minimal electronic orbital compostion one wants to see,
        isym is the irrep one wants to print, jsym is orbital within that irrep, the default is to print all the irreps and orbitals,
        this is for organized data only.'''
        for item in self.other_information[0].split():
            print "%s\t"%item,
        for i in xrange(len(self.orb)):
            print
            for j in xrange(len(self.orb[i])):
                self.print_orbs(i,j)
                self.print_occs(i,j)
                self.print_orbital_energies(i,j)
                self.print_configuration_compositions(threshhold,i,j)
    def print_orbs(self,isym=-1,jsym=-1):
        '''Print out the symmetry adapted orbitals numbers, 
        isym is the irrep one wants to print, jsym is orbital within that irrep, the default is to print all the irreps and orbitals,
        this is for organized data only.'''
        self.print_options(self.orb,isym,jsym)

    def print_occs(self,isym=-1,jsym=-1):
        '''Print out the symmetry adapted occupations of the orbital, 
        isym is the irrep one wants to print, jsym is orbital within that irrep, the default is to print all the irreps and orbitals,
        this is for organized data only.'''
        self.print_options(self.occ,isym,jsym)
    
    def print_orbital_energies(self,isym=-1,jsym=-1):
        '''Print out the orbital energies, 
        isym is the irrep one wants to print, jsym is orbital within that irrep, the default is to print all the irreps and orbitals,
        only print out certain number of digits, does not need to give the precision control,
        this is for organized data only.'''
        self.print_options(self.orbital_energy,isym,jsym)

    def print_configuration_compositions(self,threshhold,isym=-1,jsym=-1):
        '''Print out the electronic orbital composition coeffients within the threshhold, 
        isym is the irrep one wants to print, jsym is orbital within that irrep, the default is to print all the irreps and orbitals,
        this is for organized data only.'''
        self.print_options(self.configuration_composition,isym,jsym,threshhold)

    def print_options(self,lt,isym=-1,jsym=-1,threshhold=0):
        '''Print jsym component of isym irrep'''
        if -1!=isym and -1!=jsym:
            if type(lt[isym][jsym]) is dict:
                for item in sorted(lt[isym][jsym],key=lt[isym][jsym].get,reverse=True):
                    if abs(lt[isym][jsym][item]) > threshhold:
                        print "|| %-6s:%6.2f " %(item,lt[isym][jsym][item]),
                print
            else:
                if lt[isym][jsym].count('.')==2:
                    print "%s   \t"%lt[isym][jsym],
                else:
                    print "%-.5s\t"%lt[isym][jsym],
        elif -1!=isym and -1==jsym:
            if type(lt[isym][jsym]) is dict:
                for item in lt[isym]:
                    for jtem in sorted(item,key=item.get,reverse=True):
                        if abs(item[jtem]) > threshhold:
                            print "|| %-6s:%6.3f " %(jtem,item[jtem]),
                    print
                print
            else:
                for item in lt[isym]:
                    print item,
                print
        elif -1==isym and -1!=jsym:
            if type(lt[isym][jsym]) is dict:
                for item in lt:
                    for jtem in sorted(item,key=item.get,reverse=True):
                        if abs(item[jsym][jtem]) > threshhold:
                            print "|| %-6s:%6.3f " %(jtem,item[jsym][jtem]),
                    print
                    print
            else:
                for item in lt:
                    print item[jsym],
                print
        else:
            if type(lt[isym][jsym]) is dict:
                for item in lt:
                    for jtem in item:
                        for ktem in sorted(jtem,key=jtem.get,reverse=True):
                            if abs(jtem[ktem]) > threshhold:
                                print "|| %-6s:%6.3f " %(ktem,jtem[ktem]), 
                        print
                    print
            else:
                for item in lt:
                    for jtem in item:
                        print jtem,
                    print
 
def find_MOs(lines):
    """
    to find the start and end lines of a Mulliken openshell orbitals printed in a DIRAC output, 
    return the start and end lines' indices, 
    return empty list if no MO has been foun
    """
    indices_lt=[]
    startandend_lt=[]
    flag=0
    
    for i in range(0,len(lines)):
        if lines[i].strip() in startofmo and flag==0:
            startandend_lt.append(i)
            flag=1
        if flag==1 and lines[i].strip() and lines[i].strip().split()[0] in endofmo:
            startandend_lt.append(i)
            indices_lt.append(startandend_lt)
            startandend_lt=[]
            flag=0
    return indices_lt

def extract_and_process_raw_MOs(file_lines):
    """
    based on the start and end line of the MOs output region, to exact the unprocessed MOs, read every(if recognisable) set of MOs,
    standardize the the raw MO data and create instances of molecular_orbital class,
    return the list of instances
    """
    indices_lt=find_MOs(file_lines)
    class_instance_lt=[]
    for item in indices_lt:
        one_raw_MOs=file_lines[item[0]-1:item[1]]
        one_MOs=molecular_orbital(process_one_raw_MOs(one_raw_MOs))
        class_instance_lt.append(one_MOs)
    return class_instance_lt

def process_one_raw_MOs(raw_MOs):
    """
    to process one unprocessed raw MOs, and return a list of lists
    """
    block_lt=[]
    tmpline=''
    for line in raw_MOs:
        if line.strip():
            tmpline=tmpline+'  '+line.strip()
        else:
            block_lt.append(tmpline)
            tmpline=''
    block_lt=filter(None,block_lt)
    while True:
        for i,j in zip(range(0,len(block_lt)-1),range(1,len(block_lt))):
            if not multiple_replace(block_lt[i],['.','-',' ','+']).isdigit() and multiple_replace(block_lt[i],['.','-','+',' ','s','p','x','y','z','d','f','g','h','i']).isdigit() and not multiple_replace(block_lt[j],['.','-',' ','+']).isdigit() and multiple_replace(block_lt[j],['.','-','+',' ','s','p','x','y','z','d','f','g','h','i']).isdigit():
                block_lt[i]=block_lt[i]+'  '+block_lt[j]
                del block_lt[j]
        else:
            break
    tmplt=[]
    for i in xrange(len(block_lt)):
        if multiple_replace(block_lt[i],['.','-','+',' ','s','p','x','y','z','d','f','g','h','i']).isdigit():
            tmplt.append(i)
    block_lt_MOs_part=[]     #This is the list that stores only data that represents the MOs
    block_lt_not_MOs_part=[] #This is the list that stores data other than the MOs, such as HOMO LUMO data or title, etc.
    for item in tmplt:
        for i in xrange(len(block_lt)):
            if i==item:
                block_lt_MOs_part.append(block_lt[i])
    for item in block_lt:
        if item not in block_lt_MOs_part:
            block_lt_not_MOs_part.append(item) 
    tmplt=[]
    for i in xrange(len(block_lt_MOs_part)):
        if not multiple_replace(block_lt_MOs_part[i],['.','-',' ','+']).isdigit():
            tmplt.append(i)
        else:
            block_lt_MOs_part[i]=block_lt_MOs_part[i].replace('-',' -')
    symmetry_block_lt=[]   #This is the list that puts MOs from different symmetries in sub-lists.
    for i in xrange(len(tmplt)):
        try:
            symmetry_block_lt.append(block_lt_MOs_part[tmplt[i]:tmplt[i+1]])
        except IndexError:
            symmetry_block_lt.append(block_lt_MOs_part[tmplt[i]:])
    return [extract_info_in_symmetry_block_lt(symmetry_block_lt),block_lt_not_MOs_part]

def multiple_replace(string,lt,new=''):
    """
    read a string and replace every item in list lt with new, the default value of new is to delete
    """
    for item in lt:
        string=string.replace(item,new)
    return string

def extract_info_in_symmetry_block_lt(lt):
    """
    this function is merely for readability purposes, to seperate functionality this out from other function,
    this function reads the symmetry blocklized MOs and seperate informations such as occupation, energies, etc.
    returns a list of organized informations, that can be used as molecular_orbital class initialization input.
    """
    #for item in lt:
    tmplt1=[] #Complete list of titles of compositions
    tmplt2=[] #List of species of titles of compositions, generated by contracting tmplt1 
    tmplt3=[] #List of lists of the coefficients of compositions, each list has the same size as the complete list (tmplt1)
    tmplt4=[] #List of orb,occ,and energies
    one_irrep_lt=[[],[],[],[]] #List of all the MOs information in one irrep, the first element is the list of orb, second the list of occ, third the list of energies
                               #fourth the list of orbital characters(composition)
    all_irrep_lt=[] #List of all the irreps, list of one_irrep_lt's
    for element in lt:
        for i in xrange(len(element)):
            if 0==i:
                tmplt=element[i].split('  ') #use two spaces because one space would separate (2 1s), which is not what one might want.
                tmplt = filter(None, tmplt)
                tmplt1=tmplt
                for item in tmplt:
                    if item.strip() not in tmplt2:
                        tmplt2.append(item.strip())
            else:
                tmplt=element[i].split(' ')
                tmplt = filter(None, tmplt)
                tmplt3.append(tmplt[-len(tmplt1):])
                tmplt4.append(tmplt[0:-len(tmplt1)])
        for item in tmplt4:
            one_irrep_lt[0].append(item[0])
            one_irrep_lt[1].append(item[1])
            one_irrep_lt[2].append(reduce(lambda x,y: x+' '+y,item[2:]))
        for one_MO in tmplt3:
            one_MO_coe_dict={}
            for item in tmplt2:
                count=-1
                for i in xrange(len(tmplt1)):
                    if item==tmplt1[i].strip():
                        count+=1
                        one_MO_coe=float(one_MO[i])
                        one_MO_coe_dict[principal_AG_scan(count,item)]=one_MO_coe
            one_irrep_lt[3].append(one_MO_coe_dict)
        all_irrep_lt.append(one_irrep_lt)
        #initialize everything that's used for next irrep processing
        tmplt1=[]
        tmplt2=[]
        tmplt3=[]
        tmplt4=[]
        one_irrep_lt=[[],[],[],[]]
    return all_irrep_lt        

def principal_AG_scan(num,ag):
    """
    num is the number of same angular momentum function count, ag is the string of the angular momentum function.
    this function takes inputs as indicated, and recognises the principal quantum number 
    and returns the angular momentum function with the correct principal quantum number
    the atom number in the front does not change.
    """
    ag_str={'s':1,'p':2,'d':3,'f':4,'g':5,'h':6,'i':7}
    ag_lt=list(ag)
    ag_lt[2]=str(ag_str[ag[3]]+num)
    return "".join(ag_lt)

def input_error(a,b):
    """
    an easy check if every item indicated in a is included in b, if not, return True for a confirmed error.
    """
    for item in a.split(','):
        if item not in b:
            return True
 
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

