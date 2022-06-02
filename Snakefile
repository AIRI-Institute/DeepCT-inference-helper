configfile: "config/config.yml"

rule all:
    input:
        "output/{sample}.vcf"


rule lift_to_hg19_pt1:
    input:
        "input/{sample}.vcf"
    output:
        main="tmp/{sample}.hg19.vcf",
        unmap="tmp/{sample}.hg19.vcf.unmap"
    shell:
        "CrossMap.py vcf {config[chain_hg38_to_hg19]}  {input}  {config[genome_hg19]}  {output.main}"


rule swap_alleles_pt1:
    input:
        "tmp/{sample}.hg19.vcf.unmap"
    output:
        "tmp/{sample}.hg38.pre_fixed_alleles.noheader.vcf"
    shell:
        "set +o pipefail;"
        "grep 'Fail(REF==ALT)' {input} |"
        "scripts/swap_alleles.pl > {output}"


rule lift_to_hg19_pt2:
    input:
        "tmp/{sample}.hg38.pre_fixed_alleles.noheader.vcf"
    output:
        main="tmp/{sample}.hg19.fixed_alleles.noheader.vcf",
        unmap="tmp/{sample}.hg19.fixed_alleles.noheader.vcf.unmap"
    shell:
        "CrossMap.py vcf {config[chain_hg38_to_hg19]} {input} {config[genome_hg19]} {output.main}"


rule combine_and_filter:
    input:
        part1="tmp/{sample}.hg19.vcf",
        part2="tmp/{sample}.hg19.fixed_alleles.noheader.vcf"
    output:
        main="tmp/{sample}.pre_cache.hg19.vcf",
        filtered="tmp/{sample}.filtered.hg19.vcf"
    shell:
        "set +o pipefail;"
        "cat {input.part1} {input.part2} |"
        "scripts/filter_indels_chrM_N.pl 2> {output.filtered} |"
        "bedtools sort -header -faidx {config[genome_hg19]}.fai -i stdin > {output.main}"


rule decouple_uncached:
    input:
        "tmp/{sample}.pre_cache.hg19.vcf"
    output:
        "tmp/{sample}.not_in_cache.hg19.noheader.vcf"
    priority: 10
    shell:
        "bedtools intersect -a {input} -v -b {config[cache]} > {output}"


rule decouple_cached:
    input:
        "tmp/{sample}.pre_cache.hg19.vcf"
    output:
        "tmp/{sample}.in_cache.hg19.vcf.gz"
    threads: 4
    shell:
        "bedtools intersect -header -a {input} -u -b {config[cache]} |"
        "bgzip --force -@ {threads} > {output};"
        "tabix --force -p vcf {output}"


rule annotate_cached:
    input:
        "tmp/{sample}.in_cache.hg19.vcf.gz"
    output:
        "tmp/{sample}.deepct.in_cache.hg19.vcf"
    shell:
        "if [[ $(gunzip -c {input} | wc -l) != 0 ]];"
        "  then bcftools annotate -a {config[cache]} -c INFO/DEEPCT_CHANGE,INFO/DEEPCT_ORGANS,INFO/DEEPCT_CELLS {input} > {output};"
        "  else ln {input} {output};"
        "fi"


rule convert_uncached:
    input:
        "tmp/{sample}.not_in_cache.hg19.noheader.vcf"
    output:
        main="tmp/{sample}.inference.tsv",
        flag=temp("tmp/{sample}.inference_ready")
    priority: 10
    shell:
        "scripts/vcf2tsv.pl < {input} > {output.main};"
        "touch {output.flag}"
        

rule make_selene_template:
    input:
        "templates/inference.template.yml"
    output:
        temp("tmp/{sample}.inference.yml")
    params:
        tsv="tmp/{sample}.inference.tsv"
    priority: 10
    template_engine: "jinja2"


rule run_selene:
    input:
        main="tmp/{sample}.inference.yml",
        flag="tmp/{sample}.inference_ready"
    output:
        "tmp/{sample}.inference/concatenated_flat_predictions.npy"
    priority: 10
    threads:
        workflow.cores
    shell:
        "cd DeepCT; python -m selene_sdk ../{input.main}"


rule process_predictions:
    input:
        predictions="tmp/{sample}.inference/concatenated_flat_predictions.npy",
        vcf="tmp/{sample}.not_in_cache.hg19.noheader.vcf"
    output:
        "tmp/{sample}.deepct.not_in_cache.hg19.noheader.vcf"
    threads:
        workflow.cores
    shell:
        "scripts/apply_logreg.py --include-wo-predictions -p {input.predictions} -v {input.vcf} -o {output} 2>/dev/null"


rule store_in_new_cache:
    input:
        "tmp/{sample}.deepct.not_in_cache.hg19.noheader.vcf"
    output:
        "cache/{sample}.newcache.vcf.gz"
    threads: 4
    shell:
        "gunzip -c {config[cache]} |"
        "cat /dev/stdin {input} |"
        "scripts/clear4cache.pl |"
        "bedtools sort -header -faidx {config[genome_hg19]}.fai -i /dev/stdin |"
        "bgzip --force -@ {threads} > {output};"
        "tabix --force -p vcf {output}"


rule replace_cache:
    input:
        "cache/{sample}.newcache.vcf.gz"
    output:
        temp("tmp/{sample}.cache_changed")
    shell:
        "mv {input} {config[cache]};"
        "mv {input}.tbi {config[cache]}.tbi;"
        "touch {output}"


rule combine_predictions:
    input:
        from_cache="tmp/{sample}.deepct.in_cache.hg19.vcf",
        from_calc="tmp/{sample}.deepct.not_in_cache.hg19.noheader.vcf",
        cache_changed="tmp/{sample}.cache_changed"
    output:
        "tmp/{sample}.deepct.pre_lift.hg19.vcf"
    shell:
        "cat {input.from_cache} {input.from_calc} |"
        "bedtools sort -header -faidx {config[genome_hg19]}.fai -i stdin > {output}"


rule lift_to_hg38_pt1:
    input:
        "tmp/{sample}.deepct.pre_lift.hg19.vcf"
    output:
        main="tmp/{sample}.deepct.hg38.vcf",
        unmap="tmp/{sample}.deepct.hg38.vcf.unmap"
    shell:
        "CrossMap.py vcf {config[chain_hg19_to_hg38]} {input} {config[genome_hg38]} {output.main}"


rule swap_alleles_pt2:
    input:
        "tmp/{sample}.deepct.hg38.vcf.unmap"
    output:
        "tmp/{sample}.deepct.hg19.pre_fixed_alleles.noheader.vcf"
    shell:
        "set +o pipefail;"
        "grep 'Fail(REF==ALT)' {input} |"
        "scripts/swap_alleles.pl > {output}"


rule lift_to_hg38_pt2:
    input:
        "tmp/{sample}.deepct.hg19.pre_fixed_alleles.noheader.vcf"
    output:
        main="tmp/{sample}.deepct.hg38.fixed_alleles.noheader.vcf",
        unmap="tmp/{sample}.deepct.hg38.fixed_alleles.noheader.vcf.unmap"
    shell:
        "CrossMap.py vcf {config[chain_hg19_to_hg38]} {input} {config[genome_hg38]} {output.main}"


rule combine_results:
    input:
        part1="tmp/{sample}.deepct.hg38.vcf",
        part2="tmp/{sample}.deepct.hg38.fixed_alleles.noheader.vcf"
    output:
        "output/{sample}.vcf"
    shell:
        "cat {input.part1} {input.part2} |"
        "bedtools sort -header -faidx {config[genome_hg38]}.fai -i stdin > {output}"
