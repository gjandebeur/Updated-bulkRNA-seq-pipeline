# Updated-bulkRNA-seq-pipeline
Pipeline from basecall -> alignment -> modkit for genome skimming of modifications. Workflow designed for RNA but can easily be applied to DNA (switching alignment and dorado model)



first part (somewhat difficult) is to download Dorado through apptainer to get the newest version. OSCER can only handle older versions currently and these have a lower accuracy.
This first chunk runs in its own script to set up the apptainer.

```
##Add this into script to create the apptainer cache in your scratch directory.

mkdir -p /scratch/$USER/apptainer_images
apptainer pull /scratch/$USER/apptainer_images/ubuntu20.sif docker://ubuntu:20.04 
export APPTAINER_CACHEDIR=/scratch/$USER/apptainer_cache
export APPTAINER_TMPDIR=/scratch/$USER/apptainer_tmp 
mkdir -p $APPTAINER_CACHEDIR $APPTAINER_TMPDIR


### Heres the path to my installation of the new version as of ~ August 2025 (unlocked permissions) 
"/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-1.0.1-linux-x64/"
```


The next chunk is imortant modules for dependencies / setting the apptainer up correctly. 
This goes in the basecaller script itself separate from the installation one.
```

module load GCC
module load PyTorch
module load FlexiBLAS/3.3.1-GCC-12.3.0  
module load FFmpeg/4.4.2-GCCcore-11.3.0 
module load HTSlib
module load protobuf

export APPTAINER_CACHEDIR=/scratch/$USER/apptainer_cache
export APPTAINER_TMPDIR=/scratch/$USER/apptainer_tmp
mkdir -p $APPTAINER_CACHEDIR $APPTAINER_TMPDIR

cd /scratch/$USER
apptainer pull ubuntu20.sif docker://ubuntu:20.04

```

And this is the basecaller script itself

```
apptainer exec  --nv --bind /ourdisk:/mnt /scratch/$USER/ubuntu20.sif \        #This wraps your env into the apptainer to run correctly. It is important to fix the GLIBC 2.18 version mismatches.
/mnt/hpc/rnafold/gjandebeur/dont_archive/software/dorado-1.0.1-linux-x64/bin/dorado \
basecaller \
--models-directory /mnt/hpc/rnafold/gjandebeur/dont_archive/software/dorado-1.0.1-linux-x64/rna004_130bps_sup@v5.2.0/ \
sup,inosine_m6A_2OmeA,2OmeG,m5C_2OmeC,pseU_2OmeU \   #change model to specific modifications, check dorado github for info. (You can call multiple modifications at once similar to the format given)
--min-qscore 10 \    #qscore filtering, 10 is roughly top 90%
--emit-moves \       #this makes it run much quicker and saves 100s of GB but some downstream analysis want these still. 
--estimate-poly-a \
-r "/mnt/input/pod5/" \        #Importantly, mnt replaces the first folder you call so you have to change whatever name was originally there with mnt (Weird but something related to apptainer)
> "output.bam"
```


Once the data is basecalled, it needs to be aligned to reference genome.

The first fastq step seems redundant but its crucial for keeping modification tags in the data!

```
echo "linking the Conda environment to $ENV_PATH..."
source /opt/oscer/software/Mamba/23.1.0-4/etc/profile.d/conda.sh

module load SAMtools/1.16.1-GCC-11.3.0 #must be this version or newer for this command.

samtools fastq -@ 5 -T "*" "/ourdisk/hpc/rnafold/gjandebeur/dont_archive/hbec_dorado1/rawdata/hbec32r1a_8mod_basecalled.bam" > "/ourdisk/hpc/rnafold/gjandebeur/dont_archive/hbec_dorado1/rawdata/hbec32r1a_8mod_basecalled.fastq"
#The -@ 5 is how many threads, and the -T "*" is the command keeping the modficiation tags.


"minimap2/installation/path" \   #minimap2 is pretty easy to install, doesn't need apptainer.
-ax splice -y --secondary=no \ #check options, likely need to switch more for DNA
"/reference/GRCh38.fa" \   #path to reference genome
"/basecalled.fastq" > "/aligned.sam"


samtools sort -o "/alignedsorted.bam" "/aligned.sam"        #sort and bring back to bam 
```

Installation of modkit is similar to dorado, I have a repository already dedicated to this that follows similar steps to dorado.

[https://github.com/gjandebeur/Modkit_OSCER_installation](url)

Then run modkit as such
```
apptainer exec --bind /ourdisk:/mnt /scratch/gjandebeur/apptainer_images/ubuntu20.sif \
/mnt/hpc/rnafold/gjandebeur/dont_archive/software/modkit_v0.5/dist_modkit_v0.5.0_5120ef7/modkit \
pileup \
--mod-threshold a:0.99 --mod-threshold 17802:0.99 --mod-threshold 17596:0.99 \           #Thresholds are set for RNA, most calls below 0.99 can be from False Discovery (AI hallucinations) but change for DNA
--mod-threshold 69426:0.99 --mod-threshold 19228:0.99 --mod-threshold m:0.99 \
--mod-threshold 19229:0.99 --mod-threshold 19227:0.99 \
--motif A 0 --motif T 0 --motif C 0 --motif G 0 \     #Keep this, it basically says that a mod call can't be nucleotides away from the nucleotide its on (eg an 5mC call on a Thymine would be obv incorrect)
--ref "/mnt/reference/GRCh38.fa" \ #reference genome
"/mnt/alignedsorted.bam" \ 
/pileup.bed
```

Following this you have to run a R script to extract the stats. I'll attach mine to this repository; it's for RNA but should work great for skeleton code.

