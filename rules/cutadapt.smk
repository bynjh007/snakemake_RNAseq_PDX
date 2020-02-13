################################################################################ 
# trimming by cutadapt
################################################################################ 
# If raw file is bam foramt, it is transformed into fastq format
# Before that, the reads in the bam file should be sorted according to its read name
# This is for cutadapt which requires R1-R2 read matches for paired data

if config['options']['paired']:


    if config['options']['bam2fastq']:
        # sort the bam file accoding to read names (paired end) for cutadapt
        rule sort_bam:
            input:
                path.join(bam_down_dir, '{sample}.bam')
            output:
                temp(path.join(bam_down_dir, '{sample}_sort.bam'))
            threads:
                config['sort_bam']['threads']
            log:
                path.join(log_dir, 'sort_bam', '{sample}.log')
            shell:
                'samtools sort -n -@ {threads} {input} > {output}'

        rule bam2fastq:
            input:
                path.join(bam_down_dir, '{sample}_sort.bam')
            output:
                R1 = path.join(raw_dir, '{sample}_R1.fastq.gz'),
                R2 = path.join(raw_dir, '{sample}_R2.fastq.gz')
            params:
                options = format_options(config['bam2fastq']['options'])
            threads:
                config['bam2fastq']['threads']
            log:
                path.join(log_dir, 'bam2fastq', '{sample}.log')
            shell:
                'samtools bam2fq -c 6 -@ {threads} -1 {output.R1} -2 {output.R2} '
                '{input} 2> {log}'

    rule cutadapt:
        input:
            R1 = path.join(raw_dir, '{sample}_R1.fastq.gz'),
            R2 = path.join(raw_dir, '{sample}_R2.fastq.gz')
        output:
            R1 = path.join(raw_dir, '{sample}_R1_trimmed.fastq.gz'),
            R2 = path.join(raw_dir, '{sample}_R2_trimmed.fastq.gz')
        params:
            options = format_options(config['cutadapt']['options'])
        conda:
            '../envs/cutadapt.yaml'
        log:
            path.join(log_dir, 'cutadapt', '{sample}.log')
        shell:
            'cutadapt {params.options} -o {output.R1} -p {output.R2} '
            '{input.R1} {input.R2} 2> {log}' 

else:
   
    if config['options']['bam2fastq']:
        # sorting is not needed for single-end seq
        rule bam2fastq:
            input:
                path.join(bam_down_dir, '{sample}.bam')
            output:
                path.join(raw_dir, '{sample}.fastq.gz')
            params:
                options = format_options(config['bam2fastq']['options'])
            threads:
                config['bam2fastq']['threads']
            log:
                path.join(log_dir, 'bam2fastq', '{sample}.log')
            shell:
                'samtools bam2fq -@ {threads} {input} > {params.fq} '
                '&& gzip -cvf {params.fq} > {output} '
                '&& rm {params.fq}'

    rule cutadapt:
        input:
            path.join(raw_dir, '{sample}.fastq.gz')
        output:
            temp(path.join(raw_dir, '{sample}_trimmed.fastq.gz'))
        params:
            options = format_options(config['cutadapt']['options'])
        conda:
            '../envs/cutadapt.yaml'
        log:
            path.join(log_dir, 'cutadapt', '{sample}.log'),
        shell:
            'cutadapt {params.options} -o {output} {input} 2> {log}'



