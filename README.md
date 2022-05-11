# w-Wessim2
## Description
w-Wessim2 is an *in silico* whole exome sequencing (WES) method for generating simulated reads from a given genome. It defines regions for sequencing by using a BLAT alignment of probes to the genome being sequenced, thereby simulating the exon capture hybridisation during library preparation, and capturing the effects of copynumber mutations. 

w-Wessim2 is adapted from w-Wessim, described in Tanner et al., 2019 (https://doi.org/10.1093/bioinformatics/bty1063). Improvements to w-Wessim2 include using ReSeq (Schmeing S and Robinson M, 2020) for sequencing error modelling, as well as requiring less computational resources.
 
## Citation
Please cite the following when using w-Wessim2:

w-Wessim: 

Tanner G, Westhead DR, Droop A, Stead LF. Simulation of heterogeneous tumour genomes with HeteroGenesis and in silico whole exome sequencing. Bioinformatics. 2019 Jan 4;35(16):2850–2. 

Wessim:

Kim S, Jeong K, Bafna V. Wessim: a whole-exome sequencing simulator based on in silico exome capture. Bioinformatics. 2013;29:1076–7.

ReSeq:

Schmeing S, Robinson D. ReSeq simulates realistic Illumina high-throughput sequencing data. bioRxiv 2020.07.17.209072.


## Versions
v1.0 Initial release


## Requirements

Programs required:

* Python2
* ReSeq
* BLAT

#### BLAT

BLAT is required to generate an alignment of the probe sequences to the genome the user wishes to sequence. It takes a list of hybridisation probe sequences in FASTA format and a genome sequence, and outputs a .psl file listing locations of all alignments (that meet the given stringency parameters) for each probe. 

We reccomend using the multi-threaded pblat version of BLAT, created by Wang Meng and available from http://icebert.github.io/pblat/. Alternatively, the single threaded original BLAT program can be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/ and selecting your computer type. You may need to make the file an executable with:

```
chmod a+x ./blat
```

 The -minScore and -minIdentity parameter for BLAT can be used to adjust how stringent the alignment is. We find that 95 for both gives the most realistic WES coverage results with w-Wessim2 and is also likely to accurately model the effect of variants on probe hybridisation in WES.

### faToTwoBit 

faToTwoBit is required to generate a 2bit format of the genome for use with (p)BLAT and can also be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/.

You may need to make the file an executable with: 

```
chmod a+x ./fatotwobit
```


#### ReSeq
The seqToIllumina module of ReSeq is used to incorporate errors into the reads. ReSeq can be installed using conda:

```
conda install -c bioconda -c conda-forge reseq
```

or by downloading it from github and following installation instructions: https://github.com/schmeing/ReSeq.


## Inputs

* **A genome sequence** in FASTA format.

* **A ReSeq profile** of a dataset from the sequncing technology and organism for which data is being simulated. Pre-trained profiles for some datasets are available from https://github.com/schmeing/ReSeq-profiles.

* **Probe sequences.** These can be either 1) exon capture kit hybridisation probe sequences, or 2) real WES reads. Using real WES reads as probes produces a more realistic read coverage from w-Wessim2 across target and off target regions, as regional biases (such as those due to GC content or imperfect exon enrichment) are captured, but requires higher computational resources for the BLAT alignment. The real read probes can however be randomly downsampled, depending on the number of simulated reads needing to be generated; Using at least the same number of probes as the number of reads being simulated is ideal. Using fewer probe numbers (eg. ten times fewer probes than the number of reads needing to be generated) reduces resource requirements, but may result in clumping of the off-target reads, although it has less impact on on-target regions where there's a higher density of real WES probes.   

 Probe sequences for the Agilent SureSelect Human All Exon V4+UTRs kit (or any other kit for which probe sequences are avaialable) can be downloaded from https://earray.chem.agilent.com/suredesign/index.htm and converted to FASTA format with the Prep\_Probe2Fa.py script from the orginal Wessim tool (http://sak042.github.io/Wessim/)
 
 Real WES reads (from the NCBI Sequence Read Archive, accession no. SRR2103613, captured with the Agilent SureSelect Human All Exon V5+UTRs kit) that have been quality and adapter trimmed by cutadapt and filtered for a high BWA MEM mapping quality, are provided as sample ERR2752113 from the European Nucleotide Archive (http://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752113/real_wes_reads_probes.fastq.gz). These will need to be converted from fastq to fasta format with:
 
```
paste - - - - < file.fastq | cut -f 1,2 | \
sed 's/^@/>/' | tr "\t" "\n" > file.fa
```
Probes can be downsampled to the required number {NUM} with:
 
```
cat file.fa| paste - - > shuf -n {NUM} | \
sed 's/^.//' > shuf_file.fa
sort -k1 -n shuf_file.fa > sort_file.fa 
sed -e 's/^/>/' sort_file.fa | tr '\t' '\n' > \
downsampled_{NUM}_file.fa
```
A test dataset consisting of 1x10^5 of the ERR2752113 reads can be downloaded using:

```
 wget https://github.com/GeorgetteTanner/data/raw/master/real_wes_reads_probes_subsampled.fa.gz
```


## Resouces

The BLAT alignment takes ~200hours and ~5GB memory on a single thread when used with 1x10^8 probe sequences. This can be shortened by 1) multithreading with pBLAT, 2) splitting the probes up and running each set in parallel in separate runs, followed by combining the resulting .psl files, or 3) downsampling probe numbers as discussed above.

w-Wessim2 takes ~2hours and ~32GB memory when run with a BLAT alignment of 1x10^8 probes and the hg38 human reference genome, to generate 1x10^8 fragments. The memory requirement is dependent on both genome size and probe number, and runtime is mainly dependent on the number of simulated fragments. 

ReSeq takes ~40min and ~6GB memory to process the genome with replanceN and illuminaPE, and ~300min and ~500MB memory on a single thread to convert the simulated fragments from w-Wessim2 to reads. This can be multithreaded to reduce runtime.

## Parameters

|Parameter|Description|Default Value| 
|---|---|---|
|-R|Genome sequence in FASTA format - must be indexed with faidex. |Required
|-B|BLAT allignment output file of probes to the genome, in .psl format.|Required
|-S|Systematic errors file for genome.|Required
|-O|Output file name.|Required
|-N|Number of fragments to simulate.|Required
|-f|Mean fragement size (when using paired-end sequencing).|200
|-d|Standard deviation of fragment size.|50
|-m|Minimum fragment lenngth. |read length + 20 
|-y|Minimum required fraction of probe match to be hybridized.|50

   
## Use


```
#Download w-Wessim2
git clone 
https://github.com/GeorgetteTanner/w-Wessim2.git

#Randomly replace all Ns in the genome being sequenced
reseq replaceN -r {GENOME}.fasta -R {GENOME}_noNs.fasta

#Generate systematic errors for the genome 
reseq illuminaPE -r {GENOME}_noNs.fasta -s {ReSeq_profile}.reseq \
--stopAfterEstimation --writeSysError {GENOME}_noNs_syserrors.fq

#Convert genome to 2bit format
faToTwoBit {GENOME}_noNs.fasta {GENOME}_noNs.fasta.2bit

#Generate BLAT alignment 
pblat {GENOME}_noNs.fasta.2bit {PROBES}.txt {PROBES}_{GENOME}_noNs.psl \
-threads=24 -minScore=95 -minIdentity=95

#Generate fragments with w-Wessim2
python2 w-Wessim2.py -R {GENOME}_noNs.fasta -S \
{GENOME}_noNs_syserrors.fq -B {PROBES}_{GENOME}_noNs.psl -N \
{READNUM} -O {GENOME}_w-wessim2_{READNUM}.fa \
-T {ReSeq_profile}.reseq -m 20 -f 170 -d 35

#Create reads with ReSeq
reseq seqToIllumina -j 24 -s {ReSeq_profile}.reseq -i \
{GENOME}_w-wessim2_{READNUM}.fa -o \
{GENOME}_w-wessim2_{READNUM}.fq

#Split reads into 1st and 2nd reads
paste - - - - < {GENOME}_w-wessim2_{READNUM}.fq | sed -n \
'0~2p' |  tr "\t" "\n" > {GENOME}_w-wessim2_{READNUM}_1.fq
paste - - - - < {GENOME}_w-wessim2_{READNUM}.fq | sed -n \
'1~2p' |  tr "\t" "\n" > {GENOME}_w-wessim2_{READNUM}_2.fq
```
