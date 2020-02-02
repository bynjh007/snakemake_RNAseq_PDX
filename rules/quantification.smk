################################################################################
# RSEM for isoform analysis
################################################################################

if config['options']['rsem']:

    rule prepare_reference:
        output:
            config['rsem_reference']['prefix'] + '.transcripts.fa'
        params:
            fa = config['rsem_reference']['fa'],
            gtf = config['rsem_reference']['gtf'],
            options = format_options(config['rsem_reference']['options']),
            out_prefix = config['rsem_reference']['prefix']
        threads:
            config['rsem_reference']['threads']
        conda:
            '../envs/rsem.yaml'
        log:
            path.join(log_dir, 'rsem_quant', 'rsem_prep_ref.log')
        shell:
            'rsem-prepare-reference -p {threads} --gtf {params.gtf} {params.fa} {params.out_prefix} 2> {log}'

    rule rsem_quant:
        input:
            path.join(bam_dir, '{sample}','Aligned.toTranscriptome.out.bam'),
            config['rsem_reference']['prefix'] + '.transcripts.fa'
        output:
            path.join(rsem_dir, '{sample}.genes.results')
        params:
            ref = config['rsem_quant']['ref'],
            options = format_options(config['rsem_quant']['options']),
            out_prefix = rsem_dir+'/{sample}'
        threads:
            config['rsem_quant']['threads']
        conda:
            '../envs/rsem.yaml'
        log:
            path.join(log_dir, 'rsem_quant', '{sample}.log')
        shell:
            'rsem-calculate-expression --bam -p {threads} {params.options} '
            '{input[0]} {params.ref} {params.out_prefix} && '
            'rm {input[0]} 2> {log}'


################################################################################
# FeatureCounts for raw read counts
################################################################################

rule featureCounts_gene:
    input:
        path.join(bam_dir, '{sample}', 'Aligned.out.bam')
    output:
        path.join(featureCounts_dir, '{sample}_counts.gene.txt')
    params:
        options = format_options(config['featureCounts']['options']),
        gtf = config['featureCounts']['gtf']
    threads:
        config['featureCounts']['threads']
    conda:
        '../envs/featurecounts.yaml'
    log:
        path.join(log_dir, 'featureCounts', '{sample}_gene.txt')
    shell:
        'featureCounts {params.options} -a {params.gtf} '
        '-t exon -g gene_id -o {output} -T {threads} {input} 2> {log}'

rule merge_gene:
    input:
        expand(path.join(featureCounts_dir, '{sample}_counts.gene.txt'),
            sample = all_samples)
    output:
        path.join(featureCounts_dir, 'merged.gene.txt')
    run:
        # Merge count files.
        frames = (pd.read_csv(fp, sep="\t", skiprows=1, index_col=list(range(6)))
                    for fp in input)
        merged = pd.concat(frames, axis=1)
        # Extract sample names
        merged = merged.rename(
                    columns=lambda c: path.splitext(path.basename(c))[0])
        
        merged.to_csv(output[0], sep="\t", index=True)

rule featureCounts_exon:
    input:
        path.join(bam_dir, '{sample}', 'Aligned.out.bam')
    output:
        path.join(featureCounts_dir, '{sample}_counts.exon.txt')
    params:
        options = format_options(config['featureCounts']['options']),
        gtf = config['featureCounts']['gtf']
    threads:
        config['featureCounts']['threads']
    conda:
        '../envs/featurecounts.yaml'
    log:
        path.join(log_dir, 'featureCounts', '{sample}_exon.txt')
    shell:
        'featureCounts {params.options} -a {params.gtf} '
        '-t exon -g exon_id -o {output} -T {threads} {input} 2> {log}'

rule merge_exon:
    input:
        expand(path.join(featureCounts_dir, '{sample}_counts.exon.txt'),
            sample = all_samples)
    output:
        path.join(featureCounts_dir, 'merged.exon.txt')
    run:
        # Merge count files.
        frames = (pd.read_csv(fp, sep="\t", skiprows=1, index_col=list(range(6)))
                    for fp in input)
        merged = pd.concat(frames, axis=1)
        # Extract sample names
        merged = merged.rename(
                    columns=lambda c: path.splitext(path.basename(c))[0])
        
        merged.to_csv(output[0], sep="\t", index=True)


################################################################################
# DEXSeq for differential exon usage
################################################################################

if config['options']['dexseq']:

    # Exon level counting for exon analysis using DEXSeq
    rule DEXSeq_gtf:
        output:
            config['DEXseq']['prefix'] + '.gff'
        params:
            gtf = config['DEXseq']['gtf']
        conda:
            '../envs/dexseq.yaml'
        log:
            path.join(log_dir, 'DEXseq', 'DEXseq_gtf.txt')
        shell:
            'dexseq_prepare_annotation.py '
            '{params.gtf} {output} 2> {log}'

    rule DEXSeq_counting:
        input:
            bam = path.join(bam_dir, '{sample}','Aligned.out.bam'),
            gff = config['DEXseq']['prefix'] + '.gff'
        output:
            path.join(dexseq_dir, '{sample}.txt')
        params:
            options = format_options(config['DEXseq']['options'])
        conda:
            '../envs/dexseq.yaml'
        log:
            path.join(log_dir, 'DEXseq', '{sample}_DEXseq_counting.txt')
        shell:
            'dexseq_count.py -f bam {params.options} {input.gff} {input.bam} {output} 2>{log}'


# sorting bam files based on coordinate
# because featureCounts requires input for qname-sorted bam
# (don't need to sort bam within STAR)
rule sambamba_sort_coordinate: 
    input:
        gene = path.join(featureCounts_dir, 'merged.gene.txt'),
        exon = path.join(featureCounts_dir, 'merged.exon.txt'),
        bam = path.join(bam_dir, '{sample}', 'Aligned.out.bam')
    output:
        path.join(bam_dir, '{sample}', 'Aligned.sort.bam')
    threads:
        config['sambamba_sort']['threads']
    conda:
        '../envs/sambamba.yaml'
    params:
        rm_dir = path.join(bam_dir, '{sample}_host')
    shell:
        'sambamba sort -t {threads} -o {output} {input.bam} '
        '&& rm -rf {params.rm_dir}'

# indexing bam files
rule sambamba_index: 
    input:
        path.join(bam_dir, '{sample}', 'Aligned.sort.bam')
    output:
        path.join(bam_dir, '{sample}', 'Aligned.sort.bam.bai')
    threads:
        config['sambamba_index']['threads']
    conda:
        '../envs/sambamba.yaml'
    shell:
        'sambamba index -t {threads} {input} {output}'

