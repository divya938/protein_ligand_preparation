#!/bin/bash

Red='\e[1;31m'
Yel='\e[1;33m'
RCol='\e[0m'

echo "${Red} Warning: These script should run after adding the hydrogen to the ligand and cofactor ${RCol}"

read -p "$(echo ${Yel} index: `echo '\n> '` $RCol)" index
read -p "$(echo ${Yel} pdbid: `echo '\n> '` $RCol)" pdbid



cd /home/divya/Post-doc-projects/BindingMoaD/bmoadClassifiedProteins/${index}.${pdbid}


# ligand parameter genreation
read -p "$(echo ${Yel} ligand_charge: `echo '\n> '` $RCol)" ligand_charge


antechamber -fi pdb -fo prepi -i ligand.pdb -o ligand_clean_h.prepi -rn LIG -nc $ligand_charge -c bcc -pf y -at gaff2 -ek "maxcyc=0"
antechamber -fi pdb -fo mol2 -i ligand.pdb -o ligand_clean_h.mol2 -rn LIG -c bcc -nc $ligand_charge -pf y -at gaff2 -ek "maxcyc=0"
parmchk2 -f prepi -i ligand_clean_h.prepi -o ligand.frcmod

# cofactor parameter generation

read -p "$(echo ${Yel} cofactor_exist?: [Y/N] `echo '\n> '` $Rcol)" exist

if [ "$exist" = "Y" ] 
then
	read -p "$(echo ${Yel} enter the number of cofactors: `echo '\n>'` $RCol)" n_cofactor
	i=0
	while [ "$i" -lt "$n_cofactor" ]
	do
		read -p "$(echo ${Yel} cofactor_charge: `echo '\n:'` $RCol)" cofactor_charge
		antechamber -fi pdb -fo prepi -i cofactor_${i}.pdb -o cofactor_${i}.prepi -c bcc -pf y -nc $cofactor_charge -at gaff2 -ek "maxcyc=0"
		antechamber -fi pdb -fo mol2 -i cofactor_${i}.pdb -o cofactor_${i}.mol2 -c bcc -pf y -nc $cofactor_charge -at gaff2 -ek "maxcyc=0"
		parmchk2 -f prepi -i cofactor_${i}.prepi -o cofactor_${i}.frcmod
		i=$(($i+1))
	done
else
	echo "${Red} Cofactor doesn't exist ${RCol}"
fi

rm sqm.in sqm.out sqm.pdb

