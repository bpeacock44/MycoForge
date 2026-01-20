# Example Infiles
The file needs to have the commented headers as shown. These should be the path to your long reads (they don't need to be cleaned or trimmed), your short reads (raw - the pipeline will clean and trim them), your RNAseq data (raw - the pipeline will clean and trim them), and your desired output directory. 

The columns should be tab delimited as follows: 
  - LR (long read data file, compressed as .gz)
  - R1 (short read data file 1, compressed as .gz)
  - R2 (short read data file 2, compressed as .gz)
  - RNA1 (RNAseq transcripts file 1, compressed as .gz)
  - RNA2 (RNAseq transcripts file 2, compressed as .gz)
  - Designation (could be strain - this will be appended to all output files so it should distinguish that organism somehow)
  - Binomial Species Name (use an underscore - e.g homo_sapiens or hyalorbilia_sp)
  - Assembler (canu or flye if you have long reads, spades if you only have short)
  - Genome_Size (expected size of genome - this will only be used if you are using canu)
  - Busco_Type (what kind of busco analysis you want - protein, genome, or both)
  - Busco_Lineage (run ```module load busco``` and then ```busco --list-datasets``` to see options)

Since they are tab-delimited, please don't put tabs IN the values themselves! (Or regular spaces - just to be safe.) Notice I've used underscores below.

## Long AND Short Reads
```
#Location of Long Reads=/path/to/long_read_data                         
#Location of Short Reads=/path/to/short_read_data                          
#Location of RNAseq Reads=/path/to/RNAseq_data                         
#Desired Output Location=/path/to/output                        
#LR R1  R2  RNA1  RNA2  Designation Species_Name  Assembler Genome_Size Busco_Type  Busco_Lineage     
poch_long.gz poch_short_1.fastq.gz  poch_short_2.fastq.gz  poch_RNA_1.fastq.gz  poch_RNA_2.fastq.gz  poch Hyalorbilia_sp  canu  40m protein ascomycota_odb12
acst_long.fq.gz acst_short_1.fq.gz  acst_short_2.fq.gz  acst_RNA_1.fq.gz  acst_RNA_2.fq.gz  acst  Acremonium_strictum flye  40m protein ascomycota_odb12
```
## Long Reads Only
```
#Location of Long Reads=/path/to/long_read_data
#Location of Short Reads=/path/to/short_read_data
#Location of RNAseq Reads=/path/to/RNAseq_read_data
#Desired Output Location=/path/to/output
#LR R1 R2 RNA1 RNA2 Designation Species_Name
dufl_long.fq.gz NA NA dufl_RNA_R1.fq.gz dufl_RNA_R2.fq.gz dufl Duddingtonia_flagrans  flye  40m protein ascomycota_odb12
```
## Short Reads Only With No RNAseq
```
#Location of Long Reads=/path/to/long_read_data
#Location of Short Reads=/path/to/short_read_data
#Location of RNAseq Reads=/path/to/RNAseq_read_data
#Desired Output Location=/path/to/output
#LR R1 R2 RNA1 RNA2 Designation Species_Name
NA orol_short_R1.fq.gz orol_short_R2.fq.gz NA NA orol Orbilia_oligospora canu 40m protein ascomycota_odb12
```