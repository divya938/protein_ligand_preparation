from prody import *
import numpy as np
import os
from colorama import init, Fore, Back, Style

init()

def prYellow(skk): 
    print("\033[93m {}\033[00m" .format(skk))

path = '/home/divya/Post-doc-projects/BindingMoaD/bmoadClassifiedProteins/'
index = int(input(Fore.YELLOW + "Enter the index:\n" + Style.RESET_ALL))
pdbid = input(Fore.YELLOW + "Enter the pdbid: \n" + Style.RESET_ALL)
file = str(index)+ '.' + pdbid
folder = os.path.join(path,file)

os.chdir(folder)

def inter_ssbond(sel_chain):
    with open('{1}/{0}.pdb'.format(pdbid,folder),'r') as inFile:
        ssbond_chains = [ ]
        for line in inFile.readlines():
            if 'SSBOND' in line:
                ssbond_chains.extend(line.split()[3],line.split()[6])
            return ssbond_chains

tleap_file = open('tleap.in','a')
#PDBFile = input(Fore.YELLOW + "Please enter input pdb file:\n" + Style.RESET_ALL)
#PDBOutFile = input(Fore.YELLOW + "Please enter output pdb file:\n" + Style.RESET_ALL)
#last_file = input(Fore.RED + "Is this the last file ? [Y/N]" + Style.RESET_ALL)

def Diff(li1, li2):
    return (list(list(set(li1)-set(li2)) + list(set(li2)-set(li1))))
 

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

files = [filename for filename in os.listdir(folder) if filename.startswith('0.15')]
chain_id = [ ]
for i in range(len(files)):
    chain_id.append(files[i].split(sep='.')[2][-1])



def rename_str(f):
    name='0.15_80_10_pH7.4_protein_noh_chain_{0}.result.pdb'.format(f)
    structure = parsePDB(name)

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
                print("Divya!!!")
                atom1 = structure.select('serial {0:d}'.format(SSBOND[i]))
                atom2 = structure.select('serial {0:d}'.format(SSBOND[j]))
                dist = calcDistance(atom1,atom2)
                if (dist <= 2.90):
                    R1 = structure.select('serial {0:d}'.format(SSBOND[i])).getResnums()
                    R2 = structure.select('serial {0:d}'.format(SSBOND[j])).getResnums()
                    print(R1,R2)
                    if name == 'peptide.pdb':
                        string = 'bond ligand_{0}.{1}.SG ligand_{0}.{2}.SG'.format(chid,int(R1),int(R2))
                    else:
                        string = 'bond receptor_{0}.{1}.SG receptor_{0}.{2}.SG'.format(chid,int(R1),int(R2))
                    prYellow(string)
                    tleap_file.write('{0} \n'.format(string))

    except:
        print('Residues CYX not exist')
    return(structure)

def inter_cys_rename(c1,c2):
    str_1 = parsePDB('protein_noh_chain_{0}.pdb'.format(c1))
    str_2 = parsePDB('protein_noh_chain_{0}.pdb'.format(c2))
    try:
        CYS_1 = list(np.unique(str_1.select('resname CYS').getResnums()))
        CYS_2 = list(np.unique(str_2.select('resname CYS').getResnums()))

        for f1 in range(len(CYS_1)):
            for f2 in range(len(CYS_2)):
                atom1 = str_1.select('resnum {0:d} name SG'.format(CYS_1[f1]))
                atom2 = str_2.select('resnum {0:d} name SG'.format(CYS_2[f2]))
                dist = calcDistance(atom1,atom2)
                if (dist <= 2.90):
                    print('dist',dist)
                    rename = str_1.select('resnum {0:d}'.format(CYS[f1]))
                    rename.setResnames('CYX')
                    rename = str_2.select('resnum {0:d}'.format(CYS[f2]))
                    rename.setResnames('CYX')
                    print(c1,':',CYS[f1],c2,':',CYS[f2])
                    R1 = atom1.getResnums()
                    R2 = atom2.getResnums()
                    string = 'bond receptor_{0}.{1}.SG receptor_{2}.{3}.SG'.format(c1,R1,c2,R2)
                    prYellow(string)
                    tleap_file.write('{0} \n'.format(string))
        PDBOutFile = 'protein_noh_chain_{0}.pdb'.format(c1)
        writePDB(PDBOutFile,str_1)
        with open(PDBOutFile, "a") as myfile:
            myfile.write("END")
        PDBOutFile = 'protein_noh_chain_{0}.pdb'.format(c2)
        writePDB(PDBOutFile,str_2)   
        with open(PDBOutFile, "a") as myfile:
            myfile.write("END")     
    except:
        print('Residues CYS not exist')
 

for chid in chain_id:
    structure = rename_str(chid)
    PDBOutFile='protein_noh_chain_{0}.pdb'.format(chid)
    writePDB(PDBOutFile, structure)
    with open(PDBOutFile, "a") as myfile:
        myfile.write("END")

for chid in chain_id:
    ssbond_chains = inter_ssbond(chid)
    if not ssbond_chains:
        prYellow("No inter disulphate bond with chain {0}".format(chid))
    else:
        inter_cys_rename(ssbond_chains[0],ssbond_chains[1]) 


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
tleap_file.close()
os.system("tleap -s -f tleap.in")
os.system("python /home/divya/Post-doc-projects/file_generation/Important_scripts/occupancy.py")


protein = structure.select('protein')




# SSBOND information from original pdb id

dirname = os.getcwd()
pdbid = dirname[-4:]
print(pdbid)

exist = 'false'
with open('{1}/{0}.pdb'.format(pdbid,dirname),'r') as inFile:
    for line in inFile.readlines():
        if 'SSBOND' in line:
            print(line.rstrip())
            exist = 'True'
            
if exist == 'True':
    with open('tleap.in','r') as tfile:
        for line in tfile.readlines():
            if 'bond' in line: 
                prYellow(line.rstrip())




