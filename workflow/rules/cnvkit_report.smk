rule cnvkit_json:
    input:
        cns="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
        cnr="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr",
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.filter.germline.vcf",
        fai=config["reference"]["fai"],
        amp_bed=config["annotate_cnv"]["cnv_amp_genes"],
        loh_bed=config["annotate_cnv"]["cnv_loh_genes"]
    output:
        json="cnv_sv/cnvkit_report/{sample}_{type}_cnvkit.json"
    conda:
        "../envs/cnvkit_json.yaml"
    script:
        "../scripts/cnvkit_json.py"

rule cnvkit_html_report:
    input:
        json="cnv_sv/cnvkit_report/{sample}_{type}_cnvkit.json",
        template=config.get("cnvkit_html_report", {}).get("template", "")
    output:
        html="cnv_sv/cnvkit_report/{sample}_{type}_cnvkit.html"
    script:
        "../scripts/cnvkit_html_report.py"
