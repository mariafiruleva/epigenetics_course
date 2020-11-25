rule call_meth:
    input: dedup_bam='bams/{genome}/{run}_trimmed_bismark_bt2.deduplicated.bam'
    output: meth="bams/{genome}/CpG_context_{run}_trimmed_bismark_bt2.deduplicated.txt.gz"
    log: "logs/call_meth/{run}_{genome}.log"
    params: out_dir=lambda wildcards: f'bams/{wildcards.genome}'
    shell: "bismark_methylation_extractor --gzip --comprehensive --ample_memory -o {params.out_dir} {input.dedup_bam} 2> {log}"

rule meth_multiqc:
    input: expand("bams/{genome}/CpG_context_{run}_trimmed_bismark_bt2.deduplicated.txt.gz", run=runs, genome='hg38_16')
    output: html='qc/meth_qc/bams.html'
    log: "logs/meth_multiqc.log"
    params: in_dir='bams', out_dir='qc/meth_qc'
    shell: "bismark2report --alignment_report {params.in_dir} --dir {params.out_dir}"

rule meth_export:
    input: meth=rules.call_meth.output.meth
    output: graph='meth_CpG/{genome}/{run}.bedGraph.gz', cov='meth_CpG/{genome}/{run}.bismark.cov.gz'
    log: "logs/meth_export/{run}_{genome}.log"
    params: out_dir=lambda wildcards: f'meth_CpG/{wildcards.genome}',
            name=lambda wildcards: f'{wildcards.run}.bedGraph.gz'
    shell: "bismark2bedGraph --dir {params.out_dir} --output {params.name} {input.meth} 2> {log}"

rule get_bw:
    input: meth=rules.meth_export.output.cov
    output: sorted_cov='meth_CpG/{genome}/{run}.sorted.cov', sorted_m_lvl='meth_CpG/{genome}/{run}.sorted.mlevel.bedGraph',
            m_lvl_bw='meth_CpG/{genome}/{run}.sorted.mlevel.bw', c_cov='meth_CpG/{genome}/{run}.sorted.ccov.bedGraph',
            c_cov_bw='meth_CpG/{genome}/{run}.sorted.ccov.bw'
    log: "logs/get_bw/{run}_{genome}.log"
    params: idx='genomes/hg38.chrom.sizes'
    shell: """
           zcat {input.meth} | sort -k1,1 -k2,2n > {output.sorted_cov}
           cat {output.sorted_cov} | awk '{{ print $1, $2-1, $3, $4 }}' > {output.sorted_m_lvl}
           bedGraphToBigWig {output.sorted_m_lvl} {params.idx} {output.m_lvl_bw}
           cat {output.sorted_cov} | awk '{{ print $1, $2-1, $3, $5+$6 }}' > {output.c_cov}
           bedGraphToBigWig {output.c_cov} {params.idx} {output.c_cov_bw}
           """

rule cov_cyto:
    input: cov=rules.meth_export.output.cov
    output: cov_cyto='meth_CpG/{genome}/{run}.CpG_report.merged_CpG_evidence.cov.gz'
    log: "logs/cov_cyto/{run}_{genome}.log"
    params: ref_dir=lambda wildcards: f'genomes/{wildcards.genome}', in_dir=lambda wildcards: f'meth_CpG/{wildcards.genome}', name=lambda wildcards:wildcards.run
    shell: "coverage2cytosine --genome_folder {params.ref_dir} --merge_CpG --gzip --dir {params.in_dir} -o {params.name} {input.cov} 2> {log}"

rule get_meth:
    input: cov_cyto=rules.cov_cyto.output.cov_cyto
    output: meth='meth_CpG/{genome}/{run}.CpG_report.merged_CpG.meth'
    log: "logs/get_meth/{run}_{genome}.log"
    params: ref_dir='genomes', in_dir=lambda wildcards:f'meth_CpG/{wildcards.genome}'
    shell: "zcat {input.cov_cyto}  | awk -v OFS='\\t' '{{ print $1, $2, "+", "CpG", $4/100, $5+$6 }}' > {output.meth}"