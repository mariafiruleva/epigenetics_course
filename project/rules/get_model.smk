rule merge_meths:
    input: meths=expand("meth_CpG/{genome}/{run}.CpG_report.merged_CpG.meth", run=runs, genome='hg38_16')
    output: part_hmr='meth_CpG_dmrs/proportion_table_all.txt'
    shell: "/opt/methpipe-4.1.1/bin/merge-methcounts -t {input.meths} > meth_CpG_dmrs/proportion_table_all.txt"

rule filter_cpg:
    input: meths=rules.merge_meths.output.part_hmr
    output: cov_5='meth_CpG_dmrs/proportion_table_all_5_cov.txt'
    shell: "cat {input.meths} | awk '( NR == 1) || ($2 >= 5 && $4 >=5 && $6 >= 5 && $8 >= 5)' > {output.cov_5}"

rule regression:
    input: cov_5=rules.filter_cpg.output.cov_5, design_matrix='design_matrix.txt'
    output: bed='meth_CpG_dmrs/dmc.bed'
    shell: "/opt/methpipe-4.1.1/bin/radmeth regression -v -factor infected -o {output.bed} {input.design_matrix} {input.cov_5}"

rule adjust:
    input: bed=rules.regression.output.bed
    output: adj_bed='meth_CpG_dmrs/dmc_adjusted.bed'
    params: bin='1:200:1'
    shell: "/opt/methpipe-4.1.1/bin/radmeth adjust -bins {params.bin} -o {output.adj_bed} {input.bed}"

rule merge:
    input: adj_bed=rules.adjust.output.adj_bed
    output: dmr='meth_CpG_dmrs/dmrs.bed'
    shell: """
           /opt/methpipe-4.1.1/bin/radmeth merge -p 0.05 {input.adj_bed} | awk -v OFS='\t' '{{ print $1,$2,$3,$4 NR, $5, $6 }}' > {output.dmr}
           """

rule filter_dmrs:
    input: dmr=rules.merge.output.dmr
    output: filtered='meth_CpG_dmrs/dmrs_avgdiff_0.025_ncyto_3.bed'
    params: thr='0.025', n='3', base_name='meth_CpG_dmrs/dmrs_avgdiff_0.025_ncyto_3'
    shell: "scripts/split_dmrs.sh {input.dmr} {params.base_name} {params.thr} {params.n}"