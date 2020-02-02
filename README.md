# snakemake_RNAseq_PDX
This is a pipeline to analyze RNAseq data from PDX samples. This pipeline uses following program:    
* Adapter trimming - cutadapt  
* Alignment- STAR  
* Removing ambiguous reads - ngs-disambiguate  
* Quantification - featureCounts (gene, exon)  
* Quanlity assessment - fastqc, multiqc
* Optional - novel junction (STAR two-pass alignment), star-fusion, RSEM (gene and isoform), DEXSeq (differential exon usage)

Requirements
-
To use this pipeline, you need to make a conda environment with some prerequisite program (snakemake, python, samtools, bedtools, and etc). Then the pipeline will be run on this specified environment.  
```
conda env create -f default_env.yaml
```

How to run this pipeline
-
1. In config file, you can set different modes through the 'options' parametes.  
		- bam2fastq : 'True' if inital data format is bam and 'False' if initial data format is fastq.gz.  
		- paired : 'True' if data is paired-end seq and 'False' if data is single-end seq.  
		- novel_junction : 'True' if you want to identify novel junctions, otherwise 'False'. This execute 2-pass alignment by STAR.  
		- star_fusion : 'True' if you want to identify fusion genes, otherwise 'False'.  
		- rsem : 'True' if you want to perform RSEM, otherwise 'False'. This is usually performed if you want to obtain isoform-level expression data.  
		- dexseq : 'True' if you want to identify differential exon usage, otherwise 'False'.  
		
	Also, 'samples' should be changed to the path of your sample list file. 'path' should be also defined by the path where your input data is located.  
	If initial data format is fastq file, then they should be deposited in the folder with name of 'raw'.  
	If initial data format is bam file, then they should be deposited in the folder with name of 'bam_downloaded'.  
	
2. In config file, you can also set specific parameters for each rule, e.g) number of threads, alignment options, ...
	
3. After changing the parameters in the config file, you can test the pipeline with '-np' option:  
```
snakemake --configfile config_XXX.yaml --use-conda --cores N -np
```

   If it works, then you can run the pipeline without '-np'
```
snakemake --configfile config_XXX.yaml --use-conda --cores N
```
Then the pipeline will run following program:
![diagram](https://raw.githubusercontent.com/bynjh007/snakemake_RNAseq_PDX/master/dag.pdf)
