rule trim_multiqc:
    input: expand("data/reads_trim/{run}_trimmed.fq.gz", run=runs),
    output: html="qc/reads_trim.html"
    log: "logs/trim_multiqc.log"
    params: out_dir='qc', name='reads_trim', f_qc='data/reads_trim'
    shell: "multiqc -f -o {params.out_dir} -n {params.name} {params.f_qc} 2> {log}"

rule trim_adapters:
    input: fq=rules.fastq_dump.output.fq
    output: fq=temp("data/reads_trim/{run}_trimmed.fq.gz")
    params: out_dir='data/reads_trim'
    log: "logs/trim_galore/{run}.log"
    shell: "trim_galore --gzip --fastqc {input.fq} -o {params.out_dir} 2> {log}"