#!/usr/bin/env python
# coding: utf-8


# In[1]:


from Bio import SeqIO
from contextlib import redirect_stdout
import sys
import pdb
import numpy as np
import subprocess

filename = "2cji.pdb"
chain = "A"
NUMMODELS = 2
filename = sys.argv[1]
chain = sys.argv[2]
pdbId = filename.split(".")[0]
modelKind = "LOOP"
#Making gluc 

chainSeq = list(SeqIO.parse(filename, "pdb-seqres"))[ord(chain) - ord("A")]

#Writing to gluc file:
with open("gluc.ali", "w") as outFile:
    outFile.write(">P1;gluc\nsequence:gluc:::::::0.00: 0.00\n{}*".format(str(chainSeq.seq)))



# In[3]:


## ALIGNS THE TWO SEQUENCES


from modeller import *

print("Beginning to align the true sequence to the structure residues.")
print("Alignment log can be found in alignLog.txt and alignment in alignment.ali")
with open("alignLog.txt", "w") as outFile:
    with redirect_stdout(outFile):        
        env = Environ()
        aln = alignment(env)
        mdl = model(env, file=pdbId, model_segment=('FIRST:{}'.format(chain),'LAST:{}'.format(chain)))
        aln.append_model(mdl, align_codes=pdbId, atom_files=filename)
        aln.append(file='gluc.ali', align_codes='gluc')
        aln.align2d()
        aln.write(file='alignment.ali', alignment_format='PIR')
        aln.write(file='alignment.pap', alignment_format='PAP')

print("Alignment done. :)")
# In[4]:


def returnMissing(seq):
    mode = "FINDING"
    start = 0
    listOfSel = []
    inBeginning = False
    for index, char in enumerate(seq):
        if index == 0 and char == '-':
            inBeginning = True
            continue
        if char == '-' and inBeginning:
            continue
        if char != '-' and inBeginning:
            inBeginning = False
        if mode == "FINDING" and char == '-':
            start = index+1
            mode = "COUNTING"
        if mode == "COUNTING" and char != '-':
            if char == "*":
                continue
            listOfSel.append((start, index))
            mode = "FINDING"
    return listOfSel

#Finding Selections:
with open("alignment.ali", "r") as inFile:
    alText = inFile.read()

sequenceText = "".join(alText.split("\n\n")[0].split("\n")[3:])


# In[5]:


# In[6]:




# In[ ]:


from modeller import *
from modeller.automodel import *    # Load the AutoModel class

listOfTuples = returnMissing(sequenceText)

for fro, to in listOfTuples:
    if abs(fro - to) > 15 or abs(fro - to) < 5:
        #Loop models would work only for small loops less than 15 and greater
        # than 5.
        # Source: https://salilab.org/archives/modeller_usage/2007/msg00250.html
        modelKind = "AUTO"

if modelKind=="LOOP":
    class MyModel(LoopModel):
        # def __init__(self, *args, **kwargs):
            # listOfArgs = list(args)
            # normalArgs = listOfArgs[:-2]
            # listOfTuples = listOfArgs[-2]
            # self.chain = listOfArgs[-1]
            # print(normalArgs)
            # print(kwargs)
            # super().__init__(*normalArgs, **kwargs)
            # self.selectionList = self.returnSelection(listOfTuples)
            
            
            
        def select_atoms(self):
            var = []
            for fro, to in listOfTuples:
                var.append(self.residue_range('{}:{}'.format(fro, chain), '{}:{}'.format(to, chain)))
            return Selection(*var)

else:
    class MyModel(AutoModel):
        # def __init__(self, *args, **kwargs):
            # listOfArgs = list(args)
            # normalArgs = listOfArgs[:-2]
            # listOfTuples = listOfArgs[-2]
            # self.chain = listOfArgs[-1]
            # print(normalArgs)
            # print(kwargs)
            # super().__init__(*normalArgs, **kwargs)
            # self.selectionList = self.returnSelection(listOfTuples)
            
            
            
        def select_atoms(self):
            var = []
            for fro, to in listOfTuples:
                var.append(self.residue_range('{}:{}'.format(fro, chain), '{}:{}'.format(to, chain)))
            return Selection(*var)


print("Modelling the following ranges:")
print(listOfTuples)

print("Model kind used : {}".format(modelKind))
print("Number of Models being made:{}".format(str(NUMMODELS)))
print("Starting modelling now.")
outFile = open("modelLog.txt", "w")
sysout = sys.stdout

sys.stdout = outFile

log.verbose()
env = Environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']
a = MyModel(env,alnfile = 'alignment.ali',
            knowns = pdbId, sequence = 'gluc')
a.starting_model= 1
a.ending_model  = NUMMODELS 

if modelKind == "LOOP":
    a.loop.starting_model = 1
    a.loop.ending_model   = 1
    a.loop.md_level       = refine.fast

a.make()

sys.stdout = sysout
outFile.close()
print("Done modelling")

with open("modelLog.txt", "r") as inFile:
    logLines = inFile.read().splitlines()

logLines = logLines[-(NUMMODELS+1):-1]

modelNames, vals = [], []

for line in logLines:
    name, val = line.split()
    modelNames.append(name)
    vals.append(float(val))

bestModel = modelNames[np.argmin(vals)]
print(modelNames)
print(vals)

print("Best Model is {}".format(bestModel))

command = "mv {} protein_noh_chain_{}.pdb".format(bestModel, chain)

subprocess.run(command.split())
print("Renamed best model to: protein_noh_chain_{}.pdb".format(chain))

print("Deleting rest of the models of the form gluc..")

subprocess.run("rm gluc.???*", shell=True)


# In[17]:




