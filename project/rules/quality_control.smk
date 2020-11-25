rule raw_multiqc:
    input: expand("fastqc/{run}_fastqc.html", run=runs),
           expand("fastp/{run}.fastp.html", run=runs)
    output: html="qc/reads_qc.html"
    log: "logs/raw_multiqc.log"
    params: out_dir='qc', name='reads_qc', f_qc='fastqc', f_p='fastp'
    shell: "multiqc -f -o {params.out_dir} -n {params.name} {params.f_qc} {params.f_p} 2> {log}"

rule fastqc:
    input: fq="data/{run}.fastq.gz"
    output: html="fastqc/{run}_fastqc.html"
    log: "logs/fastqc/{run}.log"
    params: run_id=lambda wildcards: wildcards.run, out_dir='fastqc'
    shell: "fastqc -o {params.out_dir}/{params.run_id} {input.fq} 2> {log}"

rule fastp:
    input: fq="data/{run}.fastq.gz"
    output: html="fastp/{run}.fastp.html", json="fastp/{run}.fastp.json"
    params: run_id=lambda wildcards: wildcards.run
    log: "logs/fastp/{run}.log"
    threads: 2
    shell: """
           fastp --overrepresentation_analysis --thread {threads} --in1 {input.fq} \
           --html {output.html} --json {output.json} 2> {log}
           """

rule fastq_dump:
    output: fq=temp("data/{run}.fastq.gz")
    params: run_id=lambda wildcards: wildcards.run, out_dir='data'
    log: "logs/fq_dump/{run}.log"
    shell: "fastq-dump -O {params.out_dir} --gzip {wildcards.run} 2> {log}"