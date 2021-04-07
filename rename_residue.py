from prody import *
import numpy as np
import os
from colorama import init, Fore, Back, Style

init()

def prYellow(skk): 
    print("\033[93m {}\033[00m" .format(skk))

PDBFile = input(Fore.YELLOW + "Please enter input pdb file:\n" + Style.RESET_ALL)
PDBOutFile = input(Fore.YELLOW + "Please enter output pdb file:\n" + Style.RESET_ALL)
last_file = input(Fore.RED + "Is this the last file ? [Y/N]" + Style.RESET_ALL)

chid = PDBOutFile[-5]
#PDBFile = '0.15_80_10_pH7.4_protein_noh_chain_A.result.pdb'
#PDBOutFile = 'check.pdb'
structure = parsePDB(PDBFile)

tleap_file = open('tleap.in','a')

# definining the functions

def Diff(li1, li2):
    return (list(list(set(li1)-set(li2)) + list(set(li2)-set(li1))))
 

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

# Renaming the Histidine residues

try:
    HIS = structure.select('resname HIS')
    HIS_resids = np.unique(HIS.getResnums())
    HD1 = list(np.unique(structure.select('resname HIS name HD1').getResnums()))
except:
    print('No histidines with delta substitution')

try:   

    HE2 = list(np.unique(structure.select('resname HIS name HE2').getResnums()))
    HIP = intersection(HD1,HE2)
except:
    print('No histidnes with epsilon substitution')



try:
     for i in HIP:
          rename = structure.select('resnum {0:d}'.format(i))
          rename.setResnames('HIP')
except:
     print('No HIP')

try:
     HD1 = list(np.unique(structure.select('resname HIS name HD1').getResnums()))
     for i in HD1:
         rename = structure.select('resnum {0:d}'.format(i))
         rename.setResnames('HID')
except:
     print('No HID')

try:
     (structure.select('resname HIS')).setResnames('HIE')
except:
     print('Residues HIE not exist')


try:
     ASH = list(np.unique(structure.select('resname ASP name HD2').getResnums()))
     for i in ASH:
         rename = structure.select('resnum {0:d}'.format(i))
         rename.setResnames('ASH')

except:
     print('Residues ASH not exist')


try:
     GLH = list(np.unique(structure.select('resname GLU name HE2').getResnums()))
     for i in GLH:
         rename = structure.select('resnum {0:d}'.format(i))
         rename.setResnames('GLH')
except:
     print('Residues GLH not exist')


try:
     LYS = list(np.unique(structure.select('resname LYS').getResnums()))
     HZ1 = list(np.unique(structure.select('resname LYS name HZ1').getResnums()))
     LYN = Diff(LYS,HZ1)
     print(LYN)
     for i in LYN:
         rename = structure.select('resnum {0:d}'.format(i))
         rename.setResnames('LYN')
except:
     print('Residues LYN not exist')

try:
 #   CYS = list(np.unique(structure.select('resname CYS').getResnums()))
 #   HG = list(np.unique(structure.select('resname CYS name HG').getResnums()))
 #   CYX = Diff(CYS,HG)
 #   print(CYX)
 #   for i in CYX:
 #       rename = structure.select('resnum {0:d}'.format(i))
 #       rename.setResnames('CYX')

     CYS = list(np.unique(structure.select('resname CYS').getResnums()))
     #CYS = list(np.unique(structure.select('resname CYS name SG').getSerials()))
 #   print('CYS',CYS)
    #atom = structure
     for i in range(len(CYS)):
         for j in range(i+1,len(CYS)):
             atom1 = structure.select('resnum {0:d} name SG'.format(CYS[i]))
             atom2 = structure.select('resnum {0:d} name SG'.format(CYS[j]))
             dist = calcDistance(atom1,atom2)
             if (dist <= 2.90):
                 print('dist',dist)
                 rename = structure.select('resnum {0:d} {1:d}'.format(CYS[i],CYS[j]))
                 rename.setResnames('CYX')
                 print(CYS[i],CYS[j])

     SSBOND = list(structure.select('resname CYX name SG').getSerials())
     print('SSBOND',SSBOND)
     for i in range(0,len(SSBOND)):
         for j in range(i+1,len(SSBOND)):
             atom1 = structure.select('serial {0:d}'.format(SSBOND[i]))
             atom2 = structure.select('serial {0:d}'.format(SSBOND[j]))
             dist = calcDistance(atom1,atom2)
             if (dist <= 2.90):
                 R1 = structure.select('serial {0:d}'.format(SSBOND[i])).getResnums()
                 R2 = structure.select('serial {0:d}'.format(SSBOND[j])).getResnums()
                 print(R1,R2)
                 string = 'bond receptor_{0}.{1}.SG receptor_{0}.{2}.SG'.format(chid,int(R1),int(R2))
                 prYellow(string)
                 tleap_file.write('{0} \n'.format(string))

except:
     print('Residues CYX not exist')

if last_file == 'Y':
    f = open('string.dat','r')
    line = f.readlines()
    string_1 = line[0].rstrip()

    tleap_file.write('\n \ncomplex = combine {{ {0} }} \n'.format(string_1))
    tleap_file.write('solvatebox complex TIP3PBOX 10 \n')
    tleap_file.write('addions complex Cl- 0 \n')
    tleap_file.write('addions complex Na+ 0 \n')
    tleap_file.write('savepdb complex complex_solvated.pdb \n')
    tleap_file.write('saveamberparm complex complex_solvated.prmtop complex_solvated.inpcrd \n')
    tleap_file.write('quit \n')
   # os.system("tleap -s -f tleap.in")

tleap_file.close()

protein = structure.select('protein')

writePDB(PDBOutFile, structure)
with open(PDBOutFile, "a") as myfile:
    myfile.write("END")


# SSBOND information from original pdb id

dirname = os.getcwd()
pdbid = dirname[-4:]
print(pdbid)

with open('{1}/{0}.pdb'.format(pdbid,dirname),'r') as inFile:
    for line in inFile.readlines():
        if 'SSBOND' in line:
            print(line.rstrip())



