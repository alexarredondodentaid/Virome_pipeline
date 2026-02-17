This pipeline should allow to identify viral contigs within paired reads obtained with short read shotgun metagenomics. Moreover, taxonomy and viral genes are annotated and mapped with the reads, giving an abundance estimation of each sample. The results of this pipeline are meant to be further analyzed with statististical software.

Requirements:

· Anaconda or miniconda.
· The following programs installed in separate conda environments:
  · VirSorter2 (https://github.com/jiarong/VirSorter2) version 2.2.4
  · Bowtie2 (https://www.metagenomics.wiki/tools/bowtie2/install) version 2.5.4
  · CheckV (https://anaconda.org/channels/bioconda/packages/checkv/overview) version 1.0.3
  · Cutadapt (https://cutadapt.readthedocs.io/en/stable/) version 2.6
  · DeepVirFinder (https://github.com/jessieren/DeepVirFinder) version 1.0
  · Genomad (https://anaconda.org/channels/bioconda/packages/genomad/overview) version 1.11.2
  · Kneaddata (https://github.com/biobakery/kneaddata) version 0.12.3
  · Prodigal (https://github.com/hyattpd/Prodigal) version 2.6.3
  · Seqkit (https://github.com/shenwei356/seqkit) version 2.3.0
  · FeatureCounts (https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html) version 2.1.1
  · SPAdes (https://ablab.github.io/spades/) version 4.2.0

The pipeline has been tested in ubuntu 24.04.2 LTS and is functional with the versions provided. Please, bear in mind that some changes and bugs might occur if other versions are used. 

======================
INSTRUCTIONS
======================

1) Place the script in the folder that contains all the .fastq files.
2) Check the code for the correct path to the folder that will be created (lines 13 to 18)
3) Check the number of threads and memory that you want to allocate to the pipeline (lines 20 and 21)
4) Check that your conda installation is in the right path. If not, modify it (line 30)
5) Check the extension of your files, if it is .fastq.gz everything is ok, if not, change line 44 to let the loop choose the files properly.
6) Check the code as a guide to name the conda environments (lines 65, 81, 121, 141, 167, 185, 197, and 214)
7) The result of the pipeline are two main files in the folder 
