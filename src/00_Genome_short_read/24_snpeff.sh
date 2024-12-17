# SNPeff code

# Create directory 

# Re-save genome as sequences.fa then gzip
awk '/^>/ {$0=$1} 1' N.canaliculata_assembly.fasta.softmasked > sequences.fa
gzip sequences.fa

# Move gtf file, rename and gzip
scp braker.gtf /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/
mv braker.gtf genes.gtf
gzip genes.gtf

# Use gtf to create gff 
perl /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read/24_gtf_to_gff.pl \
< /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/braker/braker_cDNA/braker/braker.gtf \
-o myfile.gff

# Move gff and rename
scp myfile.gff /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/
mv myfile.gff N.can.gff


# Create protein file
module load singularity
cd /netfiles/nunezlab/Shared_Resources/Software/AGAT
singularity run agat_1.0.0--pl5321hdfd78af_0.sif
agat_sp_extract_sequences.pl --gff /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/N.can.gff \
-f /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/N.canaliculata_assembly.fasta.softmasked -p -o /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/N.can.protein.fa

# gzip protein file
gzip N.can.protein.fa