rule cnvkit_json:
    input:
        cns="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
        cnr="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cnr",
        fai=config["reference"]["fai"]
    output:
        json="cnv_sv/cnvkit_report/{sample}_{type}_cnvkit.json"
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
