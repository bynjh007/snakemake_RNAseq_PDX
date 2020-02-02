from os import path

import numpy as np
import pandas as pd

################################################################################
# Globals                                                                      #
################################################################################

samples = pd.read_csv(config['samples'])

# default path for analysis
data_path = config["path"]["default"]

# directory for output files
bam_down_dir = path.join(data_path, 'bam_downloaded') 
raw_dir = path.join(data_path, 'raw')
bam_dir = path.join(data_path, 'bam')
fusion_dir = path.join(data_path, 'star_fusion')

rsem_dir = path.join(data_path, 'rsem')
featureCounts_dir=path.join(data_path, 'featureCounts')
dexseq_dir=path.join(data_path, 'dexseq')

qc_dir=path.join(data_path, 'qc')
log_dir=path.join(data_path, 'log')


################################################################################
# Functions                                                                    #
################################################################################

all_samples = samples.SAMPLE_ID.tolist()

def format_options(options):
    return ' '.join(options or [])

################################################################################
# Rules                                                                        #
################################################################################

def all_input(wildcards):
    quant_out = [path.join(featureCounts_dir, 'merged.gene.txt'), 
                    path.join(featureCounts_dir, 'merged.exon.txt')]
    qc_out = [path.join(qc_dir, 'multiqc_report.html')]
    index = expand(path.join(bam_dir, '{sample}', 'Aligned.sort.bam.bai'), sample = all_samples)
    all_list = quant_out + qc_out + index
    
    if config['options']['rsem']:
        rsem = expand(path.join(rsem_dir, '{sample}.genes.results'), sample = all_samples)
        all_list = all_list + rsem

    if config['options']['dexseq']:
        dexseq = expand(path.join(dexseq_dir, '{sample}.txt'), sample = all_samples)
        all_list = all_list + dexseq

    if config['options']['star_fusion']:
        fusion = expand(path.join(fusion_dir, '{sample}',
        'star-fusion.fusion_predictions.abridged.tsv'),
        sample = all_samples)
        all_list = all_list + fusion

    if config['options']['novel_junction']:
        unload = [path.join(bam_dir, 'star_unload_2pass.done')]
    else:
        unload = [path.join(bam_dir, 'star_unload.done')]
    all_list = all_list + unload

    return all_list


rule all:
    input:
        all_input

include: 'rules/cutadapt.smk'
include: 'rules/qc.smk'
include: 'rules/alignment_pdx.smk'
include: 'rules/quantification.smk'

