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

**Introduction**
The outputs of the pipeline consist of mainly two parts:
- Tool Annotations = merging the annotations by annotation software tools
- Pipeline Annotations = Tool annotations + additional copies found in the genome

**Project folder structure**
Inside a project's folder (e.g. *testProject*) there are multiple output folders, that are presented in the following.
The collapsed folders and marked files (by the + symbol in green) represent the relevant output files:
```diff
+├── finalResults
+│   ├── FinalAnnotations_ProteinFeatures.gff3
+│   ├── FinalAnnotations_StructuralFeatures.gff3
+│   ├── FinalAnnotations_TransposonMask.gff3
+│   ├── FinalAnnotations_TransposonSequences.fasta
+│   ├── FinalAnnotations_Transposons.gff3
+│   ├── PipelineAnnotations_ProteinFeatures.gff3
+│   ├── PipelineAnnotations_TransposonMask.gff3
+│   ├── PipelineAnnotations_TransposonSequences.fasta
+│   ├── PipelineAnnotations_Transposons.gff3
+│   ├── ToolAnnotations_ProteinFeatures.gff3
+│   ├── ToolAnnotations_StructuralFeatures.gff3
+│   ├── ToolAnnotations_TransposonMask.gff3
+│   ├── ToolAnnotations_TransposonSequences.fasta
+│   └── ToolAnnotations_Transposons.gff3
├── helitronScanner
├── helitronScanner_rc
├── ltrHarvest
├── ltrPred
├── mitefind
├── mitefind_rc
├── mitetracker
├── mitetracker_rc
├── must
├── NCBICDD1000
+├── parsedAnnotations
+│   ├── helitronScanner.fasta
+│   ├── helitronScanner.gff3
+│   ├── ltrHarvest.fasta
+│   ├── ltrHarvest.gff3
+│   ├── ltrPred.fasta
+│   ├── ltrPred.gff3
+│   ├── mitefind.fasta
+│   ├── mitefind.gff3
+│   ├── mitetracker.fasta
+│   ├── mitetracker.gff3
+│   ├── must.fasta
+│   ├── must.gff3
+│   ├── NCBICDD1000.gff3
+│   ├── proteinfeatures.gff3
+│   ├── proteinfeatures_masked2.gff3
+│   ├── proteinfeatures_masked3.gff3
+│   ├── proteinfeatures_masked.gff3
+│   ├── repeatmodel.fasta
+│   ├── repeatmodel.gff3
+│   ├── repeatmodel_repeats.gff3
+│   ├── repMasker.fasta
+│   ├── repMasker.gff3
+│   ├── repMasker_repeats.gff3
+│   ├── sinefind.fasta
+│   ├── sinefind.gff3
+│   ├── sinescan.fasta
+│   ├── sinescan.gff3
+│   ├── tirvish.fasta
+│   ├── tirvish.gff3
+│   └── transposonPSI.gff3
├── repeatmodel
├── repMasker
+├── sequence.fasta
+├── sequence_heads.txt
+├── sequence_rc.fasta
├── sinefind
├── sinefind_rc
├── sinescan
+├── Statistics_FinalAnnotations.txt
+├── Statistics_ToolAnnotations.txt
├── tirvish
├── tirvish_rc
├── transposonCandA
├── transposonCandB
├── transposonCandC
├── transposonCandD
├── transposonCandE
├── transposonCandF
└── transposonPSI
```

First of all, the fasta file used for the creation of the project was copied to *sequence.fasta*. The sequences in the fasta file were renamed, a matching can be found in *sequence_heads.txt*. Also, the reverse complement sequence was copied to *sequence_rc.fasta* for all softwares that annotate a single strand only.

Moreover, *Statistics_FinalAnnotations.txt* and *Statistics_ToolAnnotations.txt* contain the statistics produced by the statistics mode for the two outputs of **reasonaTE**.

The folder *parsedAnnotations* includes the parsed transposon annotations, structural feature annotations and transposon characteristic protein annotations by the different software tools in GFF3 format, as well as extracted sequences for each annotation in a FASTA file.

The folder *finalResults* includes all results - including the tool and pipeline annotations. The *ToolAnnotations_* files contain the tool annotations, the *PipelineAnnotations_* files contain the additional copies found and the *FinalAnnotations_* include both of the prior merged into one file. There are files of the annotated transposons, transposon characteristic proteins, structural features, the mask of transposon regions and the extracted and classified sequences as FASTA file. As transposon annotations are not intersection free and can include nested or overlapping transposon annotations, the basepairs annotated in the mask represent all base pairs that are annotated by one or more transposons of the transposon annotations.

## Citations
Please cite our paper if you find transposition event detector "deTEct" useful:
(in progress)
