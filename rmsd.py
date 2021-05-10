import MDAnalysis
import MDAnalysis.analysis.rms
import os
import re
import warnings
warnings.filterwarnings("ignore")

address = "."
resnames = []
ls = sorted(os.listdir(address))
segments = []
for i in ls:
    if re.search("protein", i) and re.search("0.15", i) and re.search(".pdb$", i):
#         print(i,"p")
        protein = MDAnalysis.Universe(address+"/"+str(i))
        protein = protein.select_atoms("not (resname HOH or type H)")
        segments.append(protein)

for i in ls:
    if re.search("peptide", i) and re.search(".pdb$", i):
#         print(i,"pep")
        protein = MDAnalysis.Universe(address+"/"+str(i))
        protein = protein.select_atoms("not (resname HOH or type H)")
        segments.append(protein)
     
for i in ls:
    if re.search("aminoacid", i) and re.search(".pdb$", i):
#         print(i,"pep")
        protein = MDAnalysis.Universe(address+"/"+str(i))
        protein = protein.select_atoms("not (resname HOH or type H)")
        segments.append(protein)
     
for i in ls:       
    if re.search("cofactor",i) and re.search(".pdb$", i):
#         print(i,"c")
        cofactor = MDAnalysis.Universe(address+"/"+str(i))
        resnames.extend(cofactor.atoms[:].residues.resnames)
        cofactor = cofactor.select_atoms("not (resname HOH or type H)")
        segments.append(cofactor)


for i in ls:
    if re.search("ligand",i) and re.search(".pdb$", i):
#         print(i,"l")
        ligand = MDAnalysis.Universe(address+"/"+str(i))
        resnames.extend(ligand.atoms[:].residues.resnames)
        ligand = ligand.select_atoms("not (resname HOH or type H)")
        segments.append(ligand)
                
        
# print(resnames)        

for i in ls:         
    if i == "complex_solvated.pdb":
        # print(i,"cmp")  
        cmplx = MDAnalysis.Universe(address+"/complex_solvated.pdb")
        selection_str ="(protein and not (resname WAT or type H)) " 
        for name in resnames:
            if name.replace(" ",""): 
#                 print(name)
                selection_str += " or (resname " + str(name) + " and not (resname WAT or type H))"
#         print(selection_str)
        cmplx = cmplx.select_atoms(selection_str)

        
merged = MDAnalysis.core.universe.Merge(segments[0],segments[1])
for i in range(2,len(segments)):
    merged = MDAnalysis.core.universe.Merge(merged.atoms,segments[i])
    
if merged.atoms.positions.shape[0] == cmplx.atoms.positions.shape[0] :
    print("RMSD =" ,MDAnalysis.analysis.rms.rmsd(merged.atoms.positions, cmplx.positions, weights=None, superposition=True))
else:
    print("INVALID: NO umber of atoms is different for both groups")