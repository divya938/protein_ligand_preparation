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

def Diff(li1, li2):
    return (list(list(set(li1)-set(li2)) + list(set(li2)-set(li1))))

def cofactor_status(category):
	if category == '4' or category == '1':
		prRed("Your cofactor is doesn't exist")
		cofactor_code = ' '
	elif category == '2' or category == '3':
		prRed("Your cofactor is a small organic molecule \n")
		n_cofactors = int(input("Please enter the number of cofactors in the structure \n"))
		prGreen('Please enter the resname of the cofactors resnames \n')
		cofactor_code = [ ]
		for i in range(n_cofactors):
			cofactor_code.append(input(' '))
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


def ligand_write(f):
	if f == 'ligand.pdb':
		tleap_file.write('# loading parameters for ligand \n \n')
		tleap_file.write('loadamberprep ligand_clean_h.prepi \n')
		tleap_file.write('loadamberparams ligand.frcmod \n \n')
		tleap_file.write('ligand = loadmol2 ligand_clean_h.mol2 \n')
	else:
		tleap_file.write('ligand = loadpdb aminoacid.pdb \n')

def cofactor_write(i):
	tleap_file.write('# loading parameters for cofactor \n \n')
	tleap_file.write('loadamberprep cofactor_{0}.prepi \n'.format(i))
	tleap_file.write('loadamberparams cofactor_{0}.frcmod \n'.format(i))
	tleap_file.write('cofactor_{0} = loadmol2 cofactor_{0}.mol2 \n \n'.format(i))

def tleap_source():
	tleap_file.write('source /home/divya/softwares/anaconda3/dat/leap/cmd/leaprc.protein.ff14SB \n')
	tleap_file.write('source /home/divya/softwares/anaconda3/dat/leap/cmd/leaprc.gaff2 \n')
	tleap_file.write('source /home/divya/softwares/anaconda3/dat/leap/cmd/leaprc.water.tip3p \n \n \n')

def ismetal(sub,mcode):
	met = sub.select("resname {0}".format(mcode))
	print(met,type(met))
	if met is None:
		return False
	else:
		return True

def chainsplit(chain_id):
	substructure = structure.select("chain {0}".format(chain_id))
	protein = substructure.select("protein and noh".format(chain_id))
	water = substructure.select("water chain {0}".format(chain_id))
	protein_file = 'protein_noh_chain_{0}.pdb'.format(chain_id)

	if water is not None:
		water_file = 'water_{0}.pdb'.format(chain_id)
		writePDB(water_file,water)
		tleap_file.write('water_{0} = loadpdb {1} \n'.format(chain_id,water_file))
		water_list.append('water_{0}'.format(chain_id))	
	writePDB(protein_file,protein)
	tleap_file.write('receptor_{0} = loadpdb {1} \n'.format(chain_id,protein_file))
	receptor_list.append('receptor_{0}'.format(chain_id))

	print('metal_resname',len(metal_resname))
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


def inter_ssbond(sel_chain):
	with open('{1}/{0}.pdb'.format(pdbid,folder),'r') as inFile:
			for line in inFile.readlines():
				ssbond_chains = [ ]
			if 'SSBOND' in line:
				ssbond_chains.extend(line.split()[3],line.split()[6])
				if sel_chain in ssbond_chains:
					if ssbond_chains[0] != ssbond_chains[1]:
						imp_chain = Diff(sel_chain,ssbond_chains)
				return imp_chain
			

def chain_preference():
	chain_pre = input(Fore.RED + "Do you have any peference to the chains? [Y/N] \n" + Style.RESET_ALL)
	if chain_pre == 'Y':
		good_chain = input(Fore.YELLOW + "Please enter the chain ID \n" + Style.RESET_ALL)
	else:
		prRed("There is no chain preference")
		big_num = 0
		for chain_id in range(len(chain)):
			if ligand.select("chain {0}".format(chain[chain_id])).numAtoms() > big_num:
				big_num = ligand.select("chain {0}".format(chain[chain_id])).numAtoms()
				good_chain = chain[chain_id]
	return good_chain
        
	


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

pro_chains = ((structure.select("protein")).getChids()).tolist()
p_chains = unique(pro_chains)

prRed('#.....................................................#')
prGreen('The protein chains are {0}'.format(p_chains))
prRed('#.....................................................#')


het = structure.select('hetatm not water')
if het is not None:
	het1 = het.getResnames()
	het = het1.tolist()
	het1 = unique(het)
else:
	het1 = "Protein doesnot have hetero atoms"

if het is not None:
	prRed('#.....................................................#')
	prGreen('The hetero atoms in the protein are {0}'.format(het1))
	prRed('#.....................................................#')
else:
	prRed('#.....................................................#')
	prGreen('{0}'.format(het1))
	prRed('#.....................................................#')

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

aminoacid_list = ['ALA','ARG','ASN','ASP','ASH','CYS','CYX','GLU','GLN','GLX', \
	'GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']


# This piece of code for simple ligand case and cofactor:

if option == '1' or option == '2':
	ligand_code = ligand_status(option)
   
	if ligand_code in aminoacid_list:
		ligand = structure.select("hetatm resname {0}".format(ligand_code))
		ligand_file = 'aminoacid.pdb'
	else:
		ligand = structure.select("resname {0}".format(ligand_code))
	print(ligand_file)
	chain = unique(ligand.getChids())
	chain.sort()

	nchains = unique(structure.select('protein').getChids())
	nchains.sort()

	if chain == nchains:
		good_chain = chain_preference()
		ssbond_chain = inter_ssbond(good_chain)
		if not ssbond_chain:
			prRed("No inter disulphate bond")
		else:
			chainsplit(ssbond_chain)

		substructure = structure.select('chain {0}'.format(good_chain))
		### get ligand and chain which does not have missing atom
		if ligand_code in aminoacid_list:
			ligand_struct = substructure.select("hetatm resname {0}".format(ligand_code))
			ligand_file = 'aminoacid.pdb'
		else:
			ligand_struct = substructure.select("resname {0}".format(ligand_code))
			ligand_struct.setResnames("LIG")

		writePDB(ligand_file,ligand_struct)

		# Loading amber parameters into tleap program
		ligand_write(ligand_file)

		if option == '2':
			cofactor_code = cofactor_status(option)
			n_cofactors = len(cofactor_code)
			for i in range(n_cofactors):
				cofactor_file = 'cofactor_{0}.pdb'.format(i)
				cofact = substructure.select("resname {0}".format(cofactor_code[i]))
				if cofact is not None:
					prGreen("ligand and cofactors are in the same chain")
				else:
					cofactor_chains = unique((structure.select('resname {0}'.format(cofactor_code[i]))).getChids())
					cofactor_chains.sort()
					chains_around_good_chain = unique((structure.select("protein within 8 of chain {0}".format(good_chain))).getChids())
					for chid in chains_around_good_chain:
						if chid in cofactor_chains:
							cofact = structure.select("resname {0} and chain {1}".format(cofactor_code[i],chid))
							chainsplit(chid)
							break

				writePDB(cofactor_file,cofact)
				cofactor_write(i)
				cofactor_list.append('cofactor_{0}'.format(i))

		chainsplit(good_chain)	

	else:
		good_chain = chain_preference()
		ssbond_chain = inter_ssbond(good_chain)
		if not ssbond_chain:
			prRed("No inter disulphate bond")
		else:
			chainsplit(ssbond_chain)
		
		substructure = structure.select("protein within 8 of resname {0} chain {1}" \
			.format(ligand_code,good_chain))
		Chains_to_model = unique(substructure.getChids())
		Chains_to_model.sort()
		
		if ligand_code in aminoacid_list:
			ligand_struct = structure.select("hetatm resname {0} chain {1}".format(ligand_code,good_chain))
			ligand_file = 'aminoacid.pdb'
		else:
			ligand_struct = structure.select("resname {0} chain {1}".format(ligand_code,good_chain))
			ligand_struct.setResnames('LIG')

		writePDB(ligand_file,ligand_struct)
       
		ligand_write(ligand_file)
	

		if option == '2':
			cofactor_code = cofactor_status(option)
			n_cofactors = len(cofactor_code)
			for i in range(n_cofactors):
				cofact = structure.select("resname {0} chain {1}".format(cofactor_code[i],good_chain))
				cofactor_file = 'cofactor_{0}.pdb'.format(i)
				if cofact is not None:
					prGreen("ligand and cofactors are in the same chain")
				else:
					cofactor_chains = unique(structure.select('resname {0}'.format(cofactor_code[i])).getChids())
					cofactor_chains.sort()
					print('cofactor_chains',cofactor_chains)
					chains_around_good_chain = unique((structure.select("protein within 8 of chain {0}".format(good_chain))).getChids())
					print('chains_around_good_chain',chains_around_good_chain)
					for chid in chains_around_good_chain:
						if chid in cofactor_chains:
							print(chid)
							cofact = structure.select("resname {0} and chain {1}".format(cofactor_code[i],chid))
							chainsplit(chid)
							break

				writePDB(cofactor_file,cofact)
				cofactor_write(i)
				cofactor_list.append('cofactor_{0}'.format(i))
			

		for chid in Chains_to_model:
			chainsplit(chid)


elif option == '3' or option == '4': # ligand is peptide and cofactor is organic molecule

	cofactor_code = cofactor_status(option)
	ligand_code = ligand_status(option)
	ligand = structure.select("protein chain {0}".format(ligand_code))
	writePDB('peptide.pdb',ligand)
	tleap_file.write('ligand = loadpdb peptide.pdb \n')


	water = structure.select("water chain {0}".format(ligand_code))
	if water is not None:
		writePDB('water_{0}.pdb'.format(ligand_code),water)
		tleap_file.write('water_{0} = loadpdb water_{0}.pdb \n'.format(ligand_code))
		water_list.append('water_{0}'.format(ligand_code))

	chain = unique(ligand.getChids())

	nchains = structure.select("protein within 8 of chain {0}".format(ligand_code)).getChids()
	nchains1 = unique(nchains.tolist())

	nchains = Diff(nchains1,list(ligand_code))
	nchains.sort()
	print(nchains)

	if option == '3':
		n_cofactors = len(cofactor_code)
		for i in range(n_cofactors):
			cofact = structure.select("resname {0} within 8 of chain {1}".format(cofactor_code[i],ligand_code))
			cofactor_file = 'cofactor_{0}.pdb'.format(i)
			writePDB(cofactor_file,cofact)
			# Loading amber parameters into tleap program
			cofactor_write(i)	
			cofactor_list.append('cofactor_{0}'.format(i))			

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
	
	

