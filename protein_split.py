from prody import *
import os
import re
from colorama import init, Fore, Back, Style

init()


def prRed(skk): 
	print("\033[91m {}\033[00m" .format(skk))

def prGreen(skk): 
	print("\033[92m {}\033[00m" .format(skk))

def prYellow(skk): 
	print("\033[93m {}\033[00m" .format(skk))

def unique(list1):
	unique_list = [ ]
	for i in list1:
		if i not in unique_list:
			unique_list.append(i)
	return unique_list

def cofactor_status(category):
	if category == '4' or category == '1':
		prRed("Your cofactor is doesn't exist")
		cofactor_code = ' '
	elif category == '2' or category == '3':
		prRed("Your cofactor is a small organic molecule")
		cofactor_code = input("please enter the resname of the cofactor: \n")
	else:
		sys.exit('ligand is not either a peptide or small organic molecule: so this code will not work')
	return cofactor_code

def ligand_status(category):
	if category == '3' or category == '4':
		prRed('Your ligand is a poly-peptide')
		ligand_code = input("Please enter the chain_ID: \n")
	elif category == '1' or category == '2':
		prRed('Your ligand is organic molecule')
		ligand_code = input("Please enter the resname of the ligand:\n")
	else:
		sys.exit('ligand is not either a peptide or small organic molecule: so this code will not work')
	return ligand_code

def ligand_write():
	tleap_file.write('# loading parameters for ligand \n \n')
	tleap_file.write('loadamberprep ligand_clean_h.prepi \n')
	tleap_file.write('loadamberparams ligand.frcmod \n \n')
	tleap_file.write('ligand = loadmol2 ligand_clean_h.mol2 \n')

def cofactor_write():
	tleap_file.write('# loading parameters for cofactor \n \n')
	tleap_file.write('loadamberprep cofactor.prepi \n')
	tleap_file.write('loadamberparams cofactor.frcmod \n')
	tleap_file.write('cofactor = loadmol2 cofactor.mol2 \n \n')

def tleap_source():
	tleap_file.write('source /home/divya/softwares/anaconda3/dat/leap/cmd/leaprc.protein.ff14SB \n')
	tleap_file.write('source /home/divya/softwares/anaconda3/dat/leap/cmd/leaprc.gaff2 \n')
	tleap_file.write('source /home/divya/softwares/anaconda3/dat/leap/cmd/leaprc.water.tip3p \n \n \n')

def ismetal(sub,mcode):
	met = sub.select("resname {0}".format(mcode))
	if met.numAtoms() == 0:
		return False
	else:
		return True

def chainsplit(chain_id):
	substructure = structure.select("chain {0}".format(chain_id))
	protein = substructure.select("protein and noh".format(chain_id))
	water = substructure.select("water chain {0}".format(chain_id))
	protein_file = 'protein_noh_chain_{0}.pdb'.format(chain_id)
	if water.numAtoms() != 0:
		water_file = 'water_{0}.pdb'.format(chain_id)
		writePDB(water_file,water)
		tleap_file.write('water_{0} = loadpdb {1} \n'.format(chain_id,water_file))
		water_list.append('water_{0}'.format(chain_id))	
	writePDB(protein_file,protein)
	tleap_file.write('receptor_{0} = loadpdb {1} \n'.format(chain_id,protein_file))
	receptor_list.append('receptor_{0}'.format(chain_id))

	print(len(metal_resname))
	if len(metal_resname) != 0:
		for ele in metal_resname:
			if ismetal(substructure,ele):
				print(ele)
				print('This protein has {0} metal'.format(ele))
				metal = substructure.select("resname {0}".format(ele))
				metal_file = 'metal_{0}_{1}.pdb'.format(ele,chain_id)
				writePDB(metal_file,metal)
				metal_list.append(metal_file)
				tleap_file.write('metal_{0}_{1} = loadpdb metal_{0}_{1}.pdb'.format(ele,chain_id))


def missing(name):
	sub = parsePDB(name)
	subs = sub.select('calpha')
	res = subs.getResnums()

	cy = 0

	for i in range(len(res)):
		try:
			if res[i+1]-res[i] != 1:
				cy = cy + 1
				if cy == 1:
					prRed('{0} has missing residues in it, please model them'.format(name[-23:]))
					prGreen('The residues adjacent to missing residues are follows:')
				prGreen('{0} {1}'.format(res[i],res[i+1]))
		except:
			print('')

	if cy == 0:
		prGreen('{0} doesnot have missing residues'.format(name[-23:]))


path = '/home/divya/Post-doc-projects/BindingMoaD/bmoadClassifiedProteins/'

index = int(input(Fore.YELLOW + "Enter the index:\n" + Style.RESET_ALL))

pdbid = input(Fore.YELLOW + "Enter the pdbid: \n" + Style.RESET_ALL)

file = str(index)+ '.' + pdbid

tleap_folder = os.path.join(path,file,'tleap.in')

folder = os.path.join(path,file)

cofactor_list = [ ]
receptor_list = [ ]
water_list = [ ]
metal_list = [ ]

try:
	os.mkdir(folder)
except:
	prYellow("Requested file is already exist!!!")

os.chdir(folder)
pathPDBFolder(folder)

tleap_file = open(tleap_folder,'w')

tleap_source()

structure = parsePDB(pdbid,compressed=False)

# Checking for the category of the ligand-cofactor
prGreen("Please select the category of the ligand: \n 1. ligand-organic molecule:cofactor-doesnot exist \
	\n 2. ligand-organic molecule:cofactor-organic_molecule\
	 \n 3. ligand-peptide : cofactor-organic_molecule \
	 \n 4. ligand-peptide : cofactor-doesnot exist \
	 \n 5. Others")

option = input(Fore.YELLOW + "Please type the number corresponding to the ligand-cofactor category: \n" + Style.RESET_ALL)

metal_status = input(Fore.YELLOW + "Is your protein has metal? [Y/N] \n" + Style.RESET_ALL)

if metal_status == 'Y':
	n_types = int(input(Fore.YELLOW + "Please enter how different types of metals you had in the pdb \n" + Style.RESET_ALL))
	prYellow("Please enter the resnames of the metals")
else:
	n_types = 0

metal_resname = [ ]

for i in range(n_types):
	ele = input()
	metal_resname.append(ele)




ligand_file = 'ligand.pdb'
cofactor_file = 'cofactor.pdb'


# This piece of code for simple ligand case and cofactor:

if option == '1' or option == '2':
	ligand_code = ligand_status(option)

	ligand = structure.select("resname {0}".format(ligand_code))
	chain = unique(ligand.getChids())
	chain.sort()

	nchains = unique(structure.select('protein').getChids())
	nchains.sort()

	if chain == nchains:
		substructure = structure.select("chain {0}".format(chain[0]))
		ligand_struct = substructure.select("resname {0}".format(ligand_code))
		ligand_struct.setResnames("LIG")
		writePDB(ligand_file,ligand_struct)

		# Loading amber parameters into tleap program
		ligand_write()

		if option == '2':
			cofactor_code = cofactor_status(option)
			cofact = substructure.select("resname {0}".format(cofactor_code))
			writePDB(cofactor_file,cofact)
			
			# Loading amber parameters into tleap program
			cofactor_write()
			cofactor_list.append('cofactor')

		chainsplit(chain[0])	

	else:
		Chains_IDS = structure.select("protein within 8 of resname {0} chain {1}" \
			.format(ligand_code,chain[0]))
		Chains_to_model = unique(Chains_IDS.getChids())
		Chains_to_model.sort()
		ligand_struct = structure.select("resname {0} chain {1}".format(ligand_code,chain[0]))
		ligand_struct.setResnames('LIG')
		writePDB(ligand_file,ligand_struct)
       
		ligand_write()
	

		if option == '2':
			cofactor_code = cofactor_status(option)
			cofact = structure.select("resname {0} chain {1}".format(cofactor_code,chain[0]))
			writePDB(cofactor_file,cofact)

			# Loading amber parameters into tleap program
			cofactor_write()
			cofactor_list.append('cofactor')			

		for chid in Chains_to_model:
			chainsplit(chid)


elif option == '3' or option == '4': # ligand is peptide and cofactor is organic molecule
	cofactor_code = cofactor_status(option)
	ligand_code = ligand_status(option)
	ligand = structure.select("protein chain {0}".format(ligand_code))
	writePDB('peptide.pdb',ligand)
	tleap_file.write('ligand = loadpdb peptide.pdb \n')


	water = structure.select("water chain {0}".format(ligand_code))
	writePDB('water_{0}.pdb'.format(ligand_code),water)
	tleap_file.write('water_{0} = loadpdb water_{0}.pdb \n'.format(ligand_code))
	water_list.append('water_{0}'.format(ligand_code))

	chain = unique(ligand.getChids())

	nchains = unique(structure.select("protein within 8 of chain {0}".format(chain)).getChids)
	nchains.sort()

	if option == '3':
		cofact = structure.select("resname {0} within 8 of chain {1}".format(cofactor_code,chain))
		writePDB(cofactor_file,cofact)
		# Loading amber parameters into tleap program
		cofactor_write()	
		cofactor_list.append('cofactor')			

	for chid in nchains:
		chainsplit(chid)

else: # ligand and cofactor falls under different category
	cofactor_status(option)

string = ''

for i in receptor_list:
	string = string + i + ' '
for i in cofactor_list:
	string = string + i + ' '

string = string + 'ligand' + ' '

if len(metal_list) != 0:
	for i in metal_list:
		string = string + i + ' '

for i in water_list:
	string = string + i + ' '

with open('{1}/{0}.pdb'.format(pdbid,folder),'r') as inFile:
	text = inFile.read()
	if re.search("SSBOND", text):
		tleap_file.write('\n\n')
		tleap_file.write("# Protein has disulphate bonds \n")
		prRed('Protein has disulphate bonds please include the information in tleap.in program')
	else:
		tleap_file.write("# Protein doesn't have disulphate bonds \n")

f = open('string.dat','w')
f.write(string)


# Identifying the missing residues in protein chain

for prefix in receptor_list:
	f = '{0}/{1}.{2}/protein_noh_chain_{3}.pdb'.format(path,index,pdbid,prefix[-1])
	missing(f)
	
	

