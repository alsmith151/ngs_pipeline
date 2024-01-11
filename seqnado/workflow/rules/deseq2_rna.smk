
if config["run_deseq2"]:
    rule deseq2_report_rnaseq:
        input:
            counts="seqnado_output/feature_counts/read_counts.tsv",
            qmd=f"deseq2_{config['project_name']}.qmd",
        output:
            deseq2=f"deseq2_{config['project_name']}.html",
        log:
            "seqnado_output/logs/deseq2/deseq2.log",
        container:
            "library://asmith151/seqnado/seqnado_report:latest"
        shell:
            """
            input_file=$(realpath "{input.qmd}")
            base_dir=$(dirname $input_file)
            cd "$base_dir"
            quarto render {input.qmd} --no-cache --output {output.deseq2} --log {log}
            """
    
    localrules: deseq2_report_rnaseq
