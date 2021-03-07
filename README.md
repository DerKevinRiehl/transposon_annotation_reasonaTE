# Transposon Annotator "reasonaTE"
Transposon annotation tool for the annotation of transposons, transposon characteristic proteins and structural elements of transposons. *reasonaTE*  is part of TransposonUltimate.

- **Input**: Genome assembly (FASTA file).
- **Output**: Lots of transposon annotations (GFF3 file).

## Installation
Installation as [CondaPackage](https://anaconda.org/DerKevinRiehl/transposon_annotator_reasonate):
```
 conda install -c derkevinriehl transposon_annotator_reasonate 
```
*Note: Otherwise you can find all source codes in this Github repository.*

## How to use ''reasonaTE''
**Step 1) Create a project**
```
mkdir workspace
resonaTE -mode createProject -projectFolder workspace -projectName testProject -inputFasta sequence.fasta
```

**Step 2) Annotate genome with annotation tools**
To annotate the genome with different annotation tools, three possible ways exist. 

*Option 1:* annotate with all tools automatically (this does not include ltrPred).
```
reasonate -mode annotate -projectFolder workspace -projectName testProject -tool all
```

*Option 2:* annotate with one specific tool (good for parallelization or rerunning, recommended).
It is mandatory to run the protein annotation tools *transposonPSI* and *NCBICDD1000* for the next steps.
```
reasonate -mode annotate -projectFolder workspace -projectName testProject -tool helitronScanner
reasonate -mode annotate -projectFolder workspace -projectName testProject -tool ltrHarvest
reasonate -mode annotate -projectFolder workspace -projectName testProject -tool mitefind
reasonate -mode annotate -projectFolder workspace -projectName testProject -tool mitetracker
reasonate -mode annotate -projectFolder workspace -projectName testProject -tool must
reasonate -mode annotate -projectFolder workspace -projectName testProject -tool repeatmodel
reasonate -mode annotate -projectFolder workspace -projectName testProject -tool repMasker
reasonate -mode annotate -projectFolder workspace -projectName testProject -tool sinefind
reasonate -mode annotate -projectFolder workspace -projectName testProject -tool sinescan
reasonate -mode annotate -projectFolder workspace -projectName testProject -tool tirvish
reasonate -mode annotate -projectFolder workspace -projectName testProject -tool transposonPSI
reasonate -mode annotate -projectFolder workspace -projectName testProject -tool NCBICDD1000
```

*Option 3:* run annotation tools with specified parameters (for advanced users)
For this purpose, we provide conda packages of all [transposon_annotation_tools](https://github.com/DerKevinRiehl/transposon_annotation_tools) except for ltrPred.
Please use the fasta file with renamed sequence names of the workspace project folder. (e.g. *workspace/testProject/sequence.fasta*)
Please note, as some tools (HelitronScanner, MiteFinderII, MITE-Tracker, SINE-Finder, TIRvish) do not annotate on both strands, we recommend to run these on the reverse complementary as well (e.g. *workspace/testProejct/sequence_rc.fasta*). Once you annotated the genomes with your own specified parameter settings, please copy the result files into the workspace's project Folder as shown in the example project (e.g. results of HelitronScanner to copy into *workspace/testProject/helitronScanner*) and rename the files accordingly. Please note, it is mandatory to run the protein annotation tools *transposonPSI* and *NCBICDD1000* for the next steps using the commands of *option 2*.

*Running ltrPred:* If you want to include ltrPred annotations into the pipeline as well, install and run [ltrPred](https://github.com/HajkD/LTRpred). Later on, copy the result files into the project folder (*workspace/testProject/ltrPred*) and rename the files accordingly. Please find [our tutorial for manually running LTRpred](https://github.com/DerKevinRiehl/transposon_annotation_resonaTE/blob/main/TutorialRunLTRPred.md)  even without docker using the conda package [udocker](https://github.com/indigo-dc/udocker). Based on our experience, ltrPred contributed valuable annotations including transposons and structure features. However, we were not able to create a conda package for easy and automated use, and it takes manual efforts to run it.

*Check status of annotation tools:* If you are running multiple annotation tools in parallel, or run the manually, copied and renamed the result files into the workspace folder, you can check the status of the annotation files by:
```
resonaTE -mode checkAnnotations -projectFolder workspace -projectName testProject
>Checking helitronScanner        ... completed
>Checking ltrHarvest     ... completed
>Checking ltrPred        ... completed
>Checking mitefind       ... completed
>Checking mitetracker    ... completed
>Checking must   ... completed
>Checking repeatmodel    ... completed
>Checking repMasker      ... completed
>Checking sinefind       ... completed
>Checking sinescan       ... completed
>Checking tirvish        ... completed
>Checking transposonPSI  ... completed
>Checking NCBICDD1000    ... completed
```
All files that are reported as "completed" will be considered by **reasonaTE** in the next steps.

**Step 3) Parse annotations**
Each of the tools will produce different output file formats. **reasonaTE** therefore provides a parser module that will unify different output files to one standardized format (GFF3). The parser module will automatically detect annotations that are available as a result from step 2, and only the available files will be considered in the next steps by the pipeline.
```
resonaTE -mode parseAnnotations -projectFolder workspace -projectName testProject
```
If you are unsure about the status of the parsing, you can run following command:
```
reasonaTE -mode checkParsed -projectFolder workspace -projectName testProject
```

**Step 4) Run the pipeline on the genome annotations**
```
reasonaTE -mode pipeline -projectFolder workspace -projectName testProject
```

**Step 5) Calculate final statistics**
```
reasonaTE -mode statistics -projectFolder workspace -projectName testProject
```

## Documentation of output files


## Citations
Please cite our paper if you find transposition event detector "deTEct" useful:
(in progress)
