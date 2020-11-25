rule get_hmr:
    input: meth=rules.get_meth.output.meth
    output: hmr='meth_CpG_hmr/{genome}/{run}.hmr'
    log: "logs/get_hmr/{run}_{genome}.log"
    params: ref_dir='genomes', in_dir=lambda wildcards:f'meth_CpG/{wildcards.genome}'
    shell: "/opt/methpipe-4.1.1/bin/hmr -o {output.hmr} {input.meth} 2> {log}"

rule get_hypermr:
    input: meth=rules.get_meth.output.meth
    output: hypermr='meth_CpG_hypermr/{genome}/{run}.hypermr'
    log: "logs/get_hypermr/{run}_{genome}.log"
    params: ref_dir='genomes', in_dir=lambda wildcards:f'meth_CpG/{wildcards.genome}'
    shell: "/opt/methpipe-4.1.1/bin/hypermr -o {output.hypermr} {input.meth} 2> {log}"

rule get_part_hmr:
    input: meth=rules.get_meth.output.meth
    output: part_hmr='meth_CpG_hmr/{genome}/{run}.partial.hypermr'
    log: "logs/get_part_hmr/{run}_{genome}.log"
    params: ref_dir='genomes', in_dir=lambda wildcards:f'meth_CpG/{wildcards.genome}'
    shell: "/opt/methpipe-4.1.1/bin/hmr -partial -o {output.part_hmr} {input.meth} 2> {log}"

rule get_part_dmn:
    input: meth=rules.get_meth.output.meth
    output: part_pmd='meth_CpG_hmr/{genome}/{run}.partial.pmd'
    log: "logs/get_part_dmn/{run}_{genome}.log"
    params: ref_dir='genomes', in_dir=lambda wildcards:f'meth_CpG/{wildcards.genome}'
    shell: "/opt/methpipe-4.1.1/bin/pmd -o {output.part_pmd} {input.meth} 2> {log}"