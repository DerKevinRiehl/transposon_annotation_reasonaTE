# Running LTRPred with udocker

## Step 1) Installation of udocker
As docker is not available or usable on all environments / servers / clusters we recommend using the [conda package](https://anaconda.org/bioconda/udocker) [udocker](https://github.com/indigo-dc/udocker). It can be installed using:
```
conda install -c bioconda/label/cf201901 udocker 
```

## Step 2) Get udocker image from LTRpred and install it
This is based on the [tutorial](https://github.com/HajkD/LTRpred/issues/16) suggested by LTRpred's authors.
```
udocker pull drostlab/ltrpred
udocker create --name=ltrpred drostlab/ltrpred
```

### Step 3) Open udocker console
This command will open a new console in which you can run commands inside the udocker:
```
udocker run ltrpred
```
Please note, all following steps are considered to be written in the udocker console.

## Step 4) Copy files into the udocker
For this tutorial we create a folder to store our input fasta files for LTRpred:
```
mkdir ltrpred_data
```

Copying a file into the udocker depends on your system / environment. You could either upload your fasta file to the internet and use *wget* inside the udocker to download it into the udocker:
```
wget https://website.com/myUploads/sequence.fasta
```

Or you could use *SSH* and *SCP*:
```
sshpass -p PWD scp USER@MACHINE:workspace/testProject/sequence.fa ltrpred_data/chromosomes.fa
```

## Step 5) Run LTR pred
Now you should open R, which will open a R console.
```
R
```
Inside the console you apply LTRpred using:
```
LTRpred:: LTRpred(genome.file = "ltrpred_data/sequence.fa", cores=2)
```
After LTRpred finished, you can close the R console:
```
quit()
```

## Step 6) Copy the result files to your local machine
Copy the results back to your system depends on your system / environment. You could either use *SSH* and *SCP*:
```
sshpass -p PWD scp -r sequence_ltrpred USER@MACHINE: workspace/testProject/ltrPred
```
or upload it to the internet, and then downloading it back from there.

Thats it, now you should be able to use LTRpred without docker using the conda package udocker.

