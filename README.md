# snakemake_RNAseq_PDX
Before you use this pipeline, conda environment with some packages is required. 
- 


This is a pipeline to analyze RNAseq data from PDX samples. You can set different modes by changing the parameters in the config file.

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
	
2. In config file, you can set different parameters for each rule, e.g) number of threads, alignment options, ...
	
