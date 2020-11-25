
rule dedup_bam:
    input: bam='bams/{genome}/{run}_trimmed_bismark_bt2.bam'
    output: dedup_bam='bams/{genome}/{run}_trimmed_bismark_bt2.deduplicated.bam'
    log: "logs/alignment/{run}_{genome}_dedup.log"
    params: out_dir=lambda wildcards: f'bams/{wildcards.genome}'
    threads: 2
    shell: "deduplicate_bismark --bam --output_dir {params.out_dir} {input.bam} 2> {log}"


rule sort_bam:
    input: dedup_bam='bams/{genome}/{run}_trimmed_bismark_bt2.deduplicated.bam'
    output: sorted_bam="bams/{genome}/{run}_trimmed_bismark_bt2.deduplicated.sorted.bam"
    log: "logs/alignment/{run}_{genome}_sort_bam.log"
    params: sort_opt='coordinate', create_idx='true'
    shell: """
           picard -Xmx4g SortSam --CREATE_INDEX {params.create_idx} --SORT_ORDER {params.sort_opt} 
           --INPUT {input.dedup_bam} --OUTPUT {output.sorted_bam} 2> {log}
           """

rule alignment_multiqc:
    input: expand("bams/{genome}/{run}_trimmed_bismark_bt2.deduplicated.bam", run=runs, genome=genomes)
    output: html="qc/alignment.html"
    log: "logs/alignment_multiqc.log"
    params: out_dir='qc', name='alignment', bams='bams'
    shell: "multiqc -f -o {params.out_dir} -n {params.name} {params.bams} 2> {log}"

rule alignment:
    input: fq=rules.trim_adapters.output.fq, index_dir=rules.prepare_indexes.output.index_dir
    output: bam=temp("bams/{genome}/{run}_trimmed_bismark_bt2.bam")
    log: "logs/alignment/{run}_{genome}.log"
    params: out_dir=lambda wildcards: f'bams/{wildcards.genome}',
            index_dir=lambda wildcards, input: re.sub('Bisulfite_Genome', '', input.index_dir),
            max_memory='1800G'
    threads: 2
    shell: "bismark --parallel {threads} --max_memory {params.max_memory} --output_dir {params.out_dir} {params.index_dir} {input.fq} 2> {log}"

rule prepare_indexes:
    input: genome_dir='genomes/{genome}'
    output: index_dir='genomes/{genome}/Bisulfite_Genome'
    log: "logs/prepare_indexes/{genome}.log"
    threads: 2
    shell: "bismark_genome_preparation --verbose --bowtie2 --parallel {threads} {input.genome_dir} 2> {log}"