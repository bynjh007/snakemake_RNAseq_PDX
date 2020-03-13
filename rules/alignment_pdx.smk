#####################################################
# preparing STAR alignment 
#####################################################

# input fastq files
def align_inputs(wildcards):
    fastq_path = path.join(raw_dir, "{sample}{{pair}}_trimmed.fastq.gz".format(
        sample=wildcards.sample))
    pairs = ["_R1", "_R2"] if config['options']['paired'] else [""]
    return expand(fastq_path, pair = pairs)

# Using shared memory with genomeLoad
rule star_preload:
    params:
        index_host = config['star_align']['index_host'],
        index_graft = config['star_align']['index_graft']
    output:
        touch(path.join(bam_dir, 'star_preload.done'))
    conda:
        '../envs/star.yaml'
    log:
        path.join(log_dir, 'star_align', 'genome_preload.log')
    shell:
        'STAR --genomeLoad LoadAndExit --genomeDir {params.index_host} '
        '&& STAR --genomeLoad LoadAndExit --genomeDir {params.index_graft} 2> {log}'


#####################################################
# align reads to graft and host genome separately
# extract 'more likely' graft reads by comparing two alignments
#####################################################

# alignnment to host genome
rule star_align_host:
    input:
        fq = align_inputs,
        tmp = path.join(bam_dir, 'star_preload.done')
    output:
        temp(path.join(bam_dir, '{sample}_host','Aligned.out.bam'))
    params:
        index = config['star_align']['index_host'],
        options = format_options(config['star_align']['options']),
        out_prefix = path.join(bam_dir, '{sample}_host/') 
    conda:
        '../envs/star.yaml'
    threads:
        config['star_align']['threads']
    log:
        path.join(log_dir, 'star_align', '{sample}_host.log')
    shell:
        'STAR {params.options} --genomeDir {params.index} '
        '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
        '--readFilesCommand gunzip -c --outSAMtype BAM Unsorted '
        '--chimOutType Junctions WithinBAM --outSAMunmapped Within '
        '--genomeLoad LoadAndKeep '
        '--readFilesIn {input.fq} 2> {log}'

# alignment to graft genome
rule star_align_graft:
    input:
        fq = align_inputs,
        tmp = path.join(bam_dir, 'star_preload.done')
    output:
        temp(path.join(bam_dir, '{sample}_graft','Aligned.out.bam'))
    params:
        index = config['star_align']['index_graft'],
        options = format_options(config['star_align']['options']),
        out_prefix = path.join(bam_dir, '{sample}_graft/')
    conda:
        '../envs/star.yaml'
    threads:
        config['star_align']['threads']
    log:
        path.join(log_dir, 'star_align', '{sample}_graft.log')
    shell:
        'STAR {params.options} --genomeDir {params.index} '
        '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
        '--readFilesCommand gunzip -c --outSAMtype BAM Unsorted '
        '--chimOutType Junctions WithinBAM --outSAMunmapped Within '
        '--genomeLoad LoadAndKeep '
        '--readFilesIn {input.fq} 2> {log}'

# sorting based on qnames to run disambiguate tool
rule sambamba_sort_qname:
    input:
        path.join(bam_dir, '{sample}_{organism}', 'Aligned.out.bam')
    output:
        temp(path.join(bam_dir, '{sample}_{organism}', 'Aligned.sorted.bam'))
    conda:
        '../envs/sambamba.yaml'
    threads:
        config['sambamba_sort']['threads']
    shell:
        'sambamba sort --natural-sort -t {threads} '
        '-o {output} {input}'

# obtaining non-ambigious bam files
rule disambiguate:
    input:
        host = path.join(bam_dir, '{sample}_host', 'Aligned.sorted.bam'),
        graft = path.join(bam_dir, '{sample}_graft', 'Aligned.sorted.bam')
    output:
        host_amb = temp(path.join(bam_dir, '{sample}_host', 'Aligned.amb.bam')),
        graft_amb = temp(path.join(bam_dir, '{sample}_graft', 'Aligned.amb.bam')),
        host_disamb = temp(path.join(bam_dir, '{sample}_host', 'Aligned.disamb.bam')),
        graft_disamb = temp(path.join(bam_dir, '{sample}_graft', 'Aligned.disamb.bam'))
    params:
        algorithm = config['disambiguate']['algorithm'],
        prefix = '{sample}',
        out_dir = bam_dir,
        options = config['disambiguate']['options']
    conda:
        '../envs/disambiguate.yaml'
    shell:
        'ngs_disambiguate {params.options} -o {params.out_dir} -s '
        '{params.prefix} -a {params.algorithm} {input.host} {input.graft} '
        '&& mv {params.out_dir}/{params.prefix}.ambiguousSpeciesA.bam {output.host_amb} '
        '&& mv {params.out_dir}/{params.prefix}.ambiguousSpeciesB.bam {output.graft_amb} '
        '&& mv {params.out_dir}/{params.prefix}.disambiguatedSpeciesA.bam {output.host_disamb} '
        '&& mv {params.out_dir}/{params.prefix}.disambiguatedSpeciesB.bam {output.graft_disamb}'

# unloading host/graft genome & clearing unnecessary folers
rule genome_unload_pdx:
    input:
        expand(path.join(bam_dir, '{sample}_{organism}','Aligned.out.bam'),
            sample = all_samples, organism = ['host', 'graft'])
    output:
        touch(path.join(bam_dir, 'star_unload_pdx.done'))
    params:
        host = config['star_align']['index_host'],
        graft = config['star_align']['index_graft']
    conda:
        '../envs/star.yaml'
    shell:
        'STAR --genomeLoad Remove --genomeDir {params.host} '
        '&& STAR --genomeLoad Remove --genomeDir {params.graft}'

# Identification of novel junction based on two-pass alignment
if config['options']['novel_junction']:

#####################################################
# prepare STAR 2-pass alignment
#####################################################

    if config['options']['paired']:

        # transforming bam to fastq
        rule bam2fastq_disamb_paired:
            input:
                path.join(bam_dir, '{sample}_graft', 'Aligned.disamb.bam')
            output:
                R1 = path.join(raw_dir, 'disamb', '{sample}_R1.fastq.gz'),
                R2 = path.join(raw_dir, 'disamb', '{sample}_R2.fastq.gz'),
                singleton = temp(path.join(raw_dir, 'disamb', '{sample}_sing.fastq.gz'))
            params:
                options = format_options(config['bam2fastq']['options'])
            threads:
                config['bam2fastq']['threads']
            log:
                path.join(log_dir, 'bam2fastq_disamb', '{sample}.log')
            shell:
                'samtools bam2fq -c 6 -@ {threads} -1 {output.R1} -2 {output.R2} '
                '-s {output.singleton} -n {input} 2> {log}'

    else:
        rule bam2fastq_disamb_single:
            input:
                path.join(bam_dir, '{sample}_graft', 'Aligned.disamb.bam')
            output:
                path.join(raw_dir, 'disamb', '{sample}.fastq.gz'),
            params:
                fq = path.join(raw_dir, '{sample}_disamb_temp.fastq'),
                options = format_options(config['bam2fastq']['options'])
            threads:
                config['bam2fastq']['threads']
            log:
                path.join(log_dir, 'bam2fastq_disamb', '{sample}.log')
            shell:
                'samtools bam2fq -@ {threads} {input} > {params.fq} '
                '&& gzip -cvf {params.fq} > {output} '
                '&& rm {params.fq}'


#####################################################
# STAR 2-pass alignment
#####################################################

    # Using shared memory with genomeLoad
    rule star_preload_1st:
        input:
            path.join(bam_dir, 'star_unload_pdx.done')
        params:
            config['star_align']['index_graft']
        output:
            touch(path.join(bam_dir, 'star_preload_1.done'))
        conda:
            '../envs/star.yaml'
        log:
            path.join(log_dir, 'star_align', 'genome_preload.log')
        shell:
            'STAR --genomeLoad LoadAndExit --genomeDir {params}'

    # input fastq files
    def align_inputs_disamb(wildcards):
        fastq_path = path.join(raw_dir, 'disamb', "{sample}{{pair}}.fastq.gz".format(
            sample=wildcards.sample))
        pairs = ["_R1", "_R2"] if config['options']['paired'] else [""]
        return expand(fastq_path, pair = pairs)

    # STAR-1st (normal alignment)
    rule star_align_1st:
        input:
            fq = align_inputs_disamb,
            tmp = path.join(bam_dir, 'star_preload_1.done')
        output:
            path.join(bam_dir, 'star_1', '{sample}', 'SJ.out.tab')
        params:
            index = config['star_align']['index_graft'],
            out_prefix = path.join(bam_dir, 'star_1', '{sample}/')
        conda:
            '../envs/star.yaml'
        threads:
            config['star_align']['threads']
        log:
            path.join(log_dir, 'star_align_two_pass', '{sample}_1.log')
        shell:
            'STAR --genomeDir {params.index} '
            '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
            '--readFilesCommand gunzip -c --outSAMtype BAM Unsorted '
            '--genomeLoad LoadAndKeep '
            '--readFilesIn {input.fq} 2> {log}'

    # filtering out the low confident junctionis
    # remove non-canonical junctions ($5>0)
    # remove junctions with # of uniquely mapped reads crossing the junction <=2 ($7>2)
    # remove annotated junctions because they will be included later stage ($6==0)
    rule filter_SJ:
        input:
            sj = expand(path.join(bam_dir, 'star_1', '{sample}', 'SJ.out.tab'), sample = all_samples)
        output:
            path.join(bam_dir, 'SJ_db', 'SJ.out.comb.tab')
        params:
            star_host = bam_dir,
            star_1 = path.join(bam_dir, 'star_1')
        shell:
            "cat {input.sj} | awk '($5>0 && $7>2 && $6==0)' | cut -f1-6 | sort | uniq > {output} "
            "&& rm -rf {params.star_host}/*_host {params.star_1}"

    rule genome_unload_1st:
        input:
            path.join(bam_dir, 'SJ_db', 'SJ.out.comb.tab')
        output:
            touch(path.join(bam_dir, 'star_unload_1st.done'))
        params:
            config['star_align']['index_graft']
        conda:
            '../envs/star.yaml'
        shell:
            'STAR --genomeLoad Remove --genomeDir {params}'
        
    # indexing genome with the junctions obtained from all samples
    rule star_SJ_indexing:
        input:
            path.join(bam_dir, 'SJ_db', 'SJ.out.comb.tab')
        output:
            path.join(bam_dir, 'SJ_db', 'sjdbList.out.tab')
        params:
            sjdb = path.join(bam_dir, 'SJ_db'),
            fa = config['star_align']['fasta'],
            gtf = config['star_align']['gtf']
        conda:
            '../envs/star.yaml'
        threads:
            config['star_align']['threads_indexing']
        log:
            path.join(log_dir, 'star_align', 'star_SJ_indexing.log')
        shell:
            'STAR --runMode genomeGenerate --genomeDir {params.sjdb} '
            '--genomeFastaFiles {params.fa} --sjdbGTFfile {params.gtf} '
            '--runThreadN {threads} --sjdbFileChrStartEnd {input} 2>{log}'

    # Using shared memory for newly indexed genome
    rule star_preload_2nd:
        input:
            pre = path.join(bam_dir, 'star_unload_1st.done'),
            sj = path.join(bam_dir, 'SJ_db', 'sjdbList.out.tab')
        output:
            touch(path.join(bam_dir, 'star_preload_2.done'))
        params:
            path.join(bam_dir, 'SJ_db')
        conda:
            '../envs/star.yaml'
        log:
            path.join(log_dir, 'star_align', 'star_preload_2nd.log')
        shell:
            'STAR --genomeLoad LoadAndExit --genomeDir {params} 2> {log}'

    # STAR-2nd (considering downstream analysis for chimeric alignment)
    rule star_align_2nd:
        input:
            fq = align_inputs_disamb,
            tmp = path.join(bam_dir, 'star_preload_2.done'),
            sj = path.join(bam_dir, 'SJ_db', 'sjdbList.out.tab')
        output:
            temp(path.join(bam_dir, '{sample}','Aligned.out.bam'))
        params:
            index = path.join(bam_dir, 'SJ_db'),
            options = format_options(config['star_align']['options']),
            out_prefix = path.join(bam_dir, '{sample}/'),
            qc_prefix = path.join(qc_dir, 'star_align','{sample}_final'),
            qc_dir = path.join(qc_dir, 'star_align')
        conda:
            '../envs/star.yaml'
        threads:
            config['star_align']['threads']
        log:
            path.join(log_dir, 'star_align_two_pass', '{sample}_2.log')
        shell:
            'STAR {params.options} --genomeDir {params.index} '
            '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
            '--readFilesCommand gunzip -c --outSAMtype BAM Unsorted '
            '--chimOutType Junctions WithinBAM --outSAMunmapped Within '
            '--genomeLoad LoadAndKeep '
            '--readFilesIn {input.fq} 2> {log} '
            '&& mkdir -p {params.qc_dir} '
            '&& mv {params.out_prefix}Log.final.out {params.qc_prefix}.Log.final.out'

    # unloading the loaded genomes
    rule genome_unload_2nd:
        input:
            expand(path.join(bam_dir, '{sample}','Aligned.out.bam'),
                sample = all_samples)
        output:
            touch(path.join(bam_dir, 'star_unload_2pass.done'))
        params:
            sj = path.join(bam_dir, 'SJ_db'),
            rm = path.join(bam_dir, 'star_1') 
        conda:
            '../envs/star.yaml'
        shell:
            'STAR --genomeLoad Remove --genomeDir {params.sj} '
            '&& rm -rf {params.rm}'


#####################################################
# STAR-fusion
#####################################################
    if config['options']['star_fusion']:
       
        if config['options']['paired']:
            rule star_fusion:
                input:
                    R1 = path.join(raw_dir, 'disamb', '{sample}_R1.fastq.gz'),
                    R2 = path.join(raw_dir, 'disamb', '{sample}_R2.fastq.gz')
                output:
                    path.join(fusion_dir, '{sample}', 'star-fusion.fusion_predictions.abridged.tsv')
                params:
                    reflib = config['star_fusion']['reflib'],
                    out_dir = path.join(fusion_dir, '{sample}'),
                    options = format_options(config['star_fusion']['options'])
                threads:
                    config['star_fusion']['threads']
                conda:
                    '../envs/star_fusion.yaml'
                log:
                    path.join(log_dir, 'star_fusion', '{sample}.log')
                shell:
                    'STAR-Fusion {params.options} --left_fq {input.R1} --right_fq {input.R2} '
                    '--genome_lib_dir {params.reflib} '
                    '--CPU {threads} --output_dir {params.out_dir} 2> {log} '
                    '&& rm -rf {params.out_dir}/_starF_checkpoints '
                    '&& rm -rf {params.out_dir}/star-fusion.preliminary'
        else:
            rule star_fusion:
                input:
                    path.join(raw_dir, 'disamb', '{sample}.fastq.gz'),
                output:
                    path.join(fusion_dir, '{sample}', 'star-fusion.fusion_predictions.abridged.coding_effect.tsv')
                params:
                    reflib = config['star_fusion']['reflib'],
                    out_dir = path.join(fusion_dir, '{sample}'),
                    options = format_options(config['star_fusion']['options'])
                threads:
                    config['star_fusion']['threads']
                conda:
                    '../envs/star_fusion.yaml'
                log:
                    path.join(log_dir, 'star_fusion', '{sample}.log')
                shell:
                    'STAR-Fusion {params.options} --left_fq {input} '
                    '--genome_lib_dir {params.reflib} '
                    '--CPU {threads} --output_dir {params.out_dir} 2> {log}'


else:
    # file location change
    rule change_path:
        input:
            path.join(bam_dir, '{sample}_graft', 'Aligned.disamb.bam')
        output:
            temp(path.join(bam_dir, '{sample}', 'Aligned.out.bam'))
        params:
            path.join(bam_dir, '{sample}')
        shell:
            'mkdir -p {params} && mv {input} {output}'

