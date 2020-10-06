localrules: infer_sex_from_coverage, calculate_m2_ratio


rule infer_sex_from_coverage:
    input: f"samples/{{sample}}/mosdepth/{{movie}}.{ref}.mosdepth.summary.txt"
    output: f"samples/{{sample}}/mosdepth/{{movie}}.{ref}.mosdepth.inferred_sex.txt"
    log: f"samples/{{sample}}/logs/quality_control/infer_sex_from_coverage.{{movie}}.{ref}.log"
    benchmark: f"samples/{{sample}}/benchmarks/quality_control/infer_sex_from_coverage.{{movie}}.{ref}.tsv"
    conda: "envs/pandas.yaml"
    message: "Executing {rule}: Inferring chromosomal sex from mosdepth coverage summary {input}."
    shell: "(python3 workflow/scripts/infer_sex_from_coverage.py {input} > {output}) > {log} 2>&1"


rule calculate_m2_ratio:
    input: f"samples/{{sample}}/mosdepth/{{movie}}.{ref}.mosdepth.summary.txt"
    output: f"samples/{{sample}}/mosdepth/{{movie}}.{ref}.mosdepth.M2_ratio.txt"
    log: f"samples/{{sample}}/logs/quality_control/calculate_M2_ratio.{{movie}}.{ref}.log"
    benchmark: f"samples/{{sample}}/benchmarks/quality_control/calculate_M2_ratio.{{movie}}.{ref}.tsv"
    conda: "envs/pandas.yaml"
    message: "Executing {rule}: Calculating chrM:chr2 ratio from mosdepth coverage summary {input}."
    shell: "(python3 workflow/scripts/calculate_M2_ratio.py {input} > {output}) > {log} 2>&1"
