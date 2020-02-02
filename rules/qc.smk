
if config['options']['paired']:

    rule fastqc_raw:
        input:
            path.join(raw_dir, '{sample}_{pair}.fastq.gz')
        output:
            html = path.join(qc_dir, 'fastqc_raw', '{sample}_{pair}_fastqc.html'),
            zip = path.join(qc_dir, 'fastqc_raw', '{sample}_{pair}_fastqc.zip')
        params:
            output_dir = path.join(qc_dir, 'fastqc_raw')
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'

    rule fastqc_trim:
        input:
            path.join(raw_dir, '{sample}_{pair}_trimmed.fastq.gz')
        output:
            html = path.join(qc_dir, 'fastqc_trim', '{sample}_{pair}_trimmed_fastqc.html'),
            zip = path.join(qc_dir, 'fastqc_trim', '{sample}_{pair}_trimmed_fastqc.zip')
        params:
            output_dir = path.join(qc_dir, 'fastqc_trim')
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'


else:

    rule fastqc_raw:
        input:
            path.join(raw_dir, '{sample}.fastq.gz')
        output:
            html = path.join(qc_dir, 'fastqc_raw', '{sample}_fastqc.html'),
            zip = path.join(qc_dir, 'fastqc_raw', '{sample}_fastqc.zip')
        params:
            output_dir=path.join(qc_dir, 'fastqc_raw')
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'


    rule fastqc_trim:
        input:
            path.join(raw_dir, '{sample}_trimmed.fastq.gz')
        output:
            html = path.join(qc_dir, 'fastqc_trim', '{sample}_trimmed_fastqc.html'),
            zip = path.join(qc_dir, 'fastqc_trim', '{sample}_trimmed_fastqc.zip')
        params:
            output_dir=path.join(qc_dir, 'fastqc_trim')
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'


def multiqc_prerequisite(wildcards):

# quantification results
    quant_out = [path.join(featureCounts_dir, 'merged.gene.txt'),
                    path.join(featureCounts_dir, 'merged.exon.txt')]
    
    if config['options']['rsem']:
        rsem_out = expand(path.join(rsem_dir, '{sample}.genes.results'), 
                            sample = all_samples)
        quant_out = quant_out + rsem_out

# fastqc for raw file
    fastqc_temp_raw = expand(path.join(qc_dir, 'fastqc_raw', '{sample}{{pair}}_fastqc.zip'), sample = all_samples)
    pairs = ["_R1", "_R2"] if config['options']['paired'] else [""]
    fastqc_raw = expand(fastqc_temp_raw, pair=pairs)

# merge list for prerequisite
    fastqc_temp_trim = expand(path.join(qc_dir, 'fastqc_trim', '{sample}{{pair}}_trimmed_fastqc.zip'), sample = all_samples)
    pairs = ["_R1", "_R2"] if config['options']['paired'] else [""]
    fastqc_trim = expand(fastqc_temp_trim, pair=pairs)
    all_list = quant_out + fastqc_raw + fastqc_trim
    return all_list


rule multiqc:
    input:
        multiqc_prerequisite
    output:
        path.join(qc_dir, 'multiqc_report.html')
    params:
        qc_dir = qc_dir,
        out_dir = qc_dir
    shell:
        'multiqc --outdir {params.out_dir} {params.qc_dir}'


