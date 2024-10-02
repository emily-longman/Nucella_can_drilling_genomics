#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=rename_scaffolds

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 # on one node
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=8:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=30G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Call packages
seqkit=/gpfs1/home/e/l/elongman/software/seqkit
seqtk=/gpfs1/home/e/l/elongman/software/seqtk/seqtk

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#This is the location where the reference genome.
REFERENCE=$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_5/polished_assembly.fasta

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER_SCRATCH

if [ -d "rename_scaffolds" ]
then echo "Working rename_scaffolds folder exist"; echo "Let's move on."; date
else echo "Working rename_scaffolds folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/rename_scaffolds; date
fi

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/rename_scaffolds

if [ -d "scaffolds" ]
then echo "Working scaffolds folder exist"; echo "Let's move on."; date
else echo "Working scaffolds folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/rename_scaffolds/scaffolds; date
fi

#--------------------------------------------------------------------------------

# Create guide file - i.e., a list of the scaffold names

# Move to rename scaffolds directory
cd $WORKING_FOLDER_SCRATCH/rename_scaffolds

# Get all of the scaffold names 
grep ">" $REFERENCE > scaffold_names.txt

# Get rid of the ">"
sed -i 's/>//g' scaffold_names.txt

#--------------------------------------------------------------------------------


# Cat file of scaffold names and start while loop
cat scaffold_names.txt | \
while read scaffold 
do echo ${scaffold}

# Break up the genome into each scaffold
grep -EA 1 "^>${scaffold}$" ${REFERENCE} > ${scaffold}.fasta

# Modify the scaffold name to remove "pilon"
cat ${scaffold}.fasta | $seqkit replace -p "\_pilon_pilon_pilon_pilon_pilon" -r '' | $seqtk seq \
> $WORKING_FOLDER_SCRATCH/rename_scaffolds/scaffolds/${scaffold}.fasta

# Housekeeping - remove intermediate files
rm ${scaffold}.fasta

done

#--------------------------------------------------------------------------------

scaffolds_dir=$WORKING_FOLDER_SCRATCH/rename_scaffolds/scaffolds

# Cat files together to produce a final genome with the scaffolds renamed
for file in $(find ${scaffolds_dir} -name "*.fasta"); do
cmd="cat ${file};"
eval $cmd
done > $WORKING_FOLDER_SCRATCH/rename_scaffolds/N.canaliculata_assembly.fasta