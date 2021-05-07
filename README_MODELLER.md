# How to run Modelling Script:


1. Install Modeller after setting the appropriate value to the KEY_MODELLER 
environment variable.

```bash
KEY_MODELLER = "<YOUR MODELLER LICENCE KEY HERE>"
conda install -c salilab modeller
```

2. Install Biopython for allignment of the sequence and ability to make customized
input files for the modeller.

```bash
conda install -c conda-forge biopython
```

3. [NOT COMPULSORY] You can keep an alias to the source of the modelling script
so that you can run from any of the inputfile directory without having to enter
full length of the modellingScript each time:

```bash
alias model='conda activate amber && python /mnt/e/ProteinLigand/protocol/ModellingScript.py'
```

4. Running the script:

Format for running the script:

```bash
model <FILENAME_WITH_EXTENSION> <CHAIN_ID>
```

---

## Additional Notes:

- The stdout lines of the allignment procedure are stored in the `alignLog.txt` 
file.
- The stdout lines of the modelling procedure are stored in the  `modelLog.txt` 
file. These are the files to be examined in case no output shows up.
- The NUM_MODELS variable is by default set to 5, in case the modelling doesn't
take much time and you want more accurate results you can increase it or 
decrease it.
- Loop modelling is applicable only when the number of unmodelled residues is 
greater than 5 and less than 15. The script takes care of it automatically.
- Sometimes, you might see a large string of aminoacids either at the 
beginning or ending of the modelled protein. The detection of unmodelled residues 
doesn't always work. Ideally, these parts are not needed to be modelled but 
modelling them does not harm the entire process in anyway.

