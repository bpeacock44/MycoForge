# Run Guide
Part 1 is only needed if you have short reads, as it will clean and trim them. All you need as input is an infile ([see example_infiles.md](LINK))
```sh
#!/bin/bash
#SBATCH -c 32
#SBATCH --mem=128G
#SBATCH --job-name=fungi_assm_part1

# Load required module
module load AAFTF

# Run the command
./part1.sh --input infile.txt
```

Part 2 will assemble the reads. 

```sh
#!/bin/bash
#SBATCH -c 32
#SBATCH --mem=128G
#SBATCH --job-name=fungi_assm_part2

# Load required module
module load AAFTF
module load Flye # only if you plan to use flye!
module unload java # only if you plan to use canu!
module load canu # only if you plan to use canu!

# Run the command
# Short reads only
./part2.sh --input infile.txt --assembler short --genomeSize 40m 
# Long reads with flye assembler (faster, slightly less accurate)
./part2.sh --input infile.txt --assembler flye --genomeSize 40m
# Long reasd with canu assembler (slower, slightly more accurate)
./part2.sh --input infile.txt --assembler canu --genomeSize 40m
```



module load AAFTF
module load samtools
module load bcftools
module load medaka

./part3.sh --input infile.txt --assembler flye
./part3.sh --input infile.txt --assembler canu

module purge 
module load funannotate

./part4.sh --input infile.txt --assembler flye
./part4.sh --input infile.txt --assembler canu

module purge
module load busco

busco --list-datasets

./part5.sh --input infile.txt --assembler flye --dataset ascomycota_odb12 --type both
./part5.sh --input infile.txt --assembler canu --dataset ascomycota_odb12 --type both
./part5.sh --input fusarium.txt --assembler flye --dataset nectriaceae_odb12 --type both
./part5.sh --input fusarium.txt --assembler canu --dataset nectriaceae_odb12 --type both


```sh
srun -p highmem \
    -c 32 \
    --time 1-00:00:00 \
    --mem 900G \
    --pty bash -l

module load orthofinder/2.5.5

DIR=/rhome/bpeacock/bigdata/PN113_long_read_assembly/funannotate_genomes/canu_proteins
DIR=/rhome/bpeacock/bigdata/PN113_long_read_assembly/funannotate_genomes/flye_proteins
orthofinder -f ${DIR} -t 32 -M msa -T fasttree



