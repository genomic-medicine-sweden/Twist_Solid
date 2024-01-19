from hydra_genetics.utils.io.hotspot_report import generate_hotspot_report

hotspot_file = snakemake.input.hotspots
vcf_file = snakemake.input.vcf
vcf_file_wo_pick = snakemake.input.vcf_file_wo_pick
gvcf_file = snakemake.input.gvcf

report = snakemake.output.report

levels = snakemake.params.levels
config = snakemake.params.report_config
sample_name = snakemake.params.sample_name
chr_translation_file = snakemake.params.chr_translation_file


generate_hotspot_report(sample=sample_name,
                        output=report,
                        levels=levels,
                        hotspot_file=hotspot_file,
                        vcf_file=vcf_file,
                        vcf_file_wo_pick=vcf_file_wo_pick,
                        gvcf_file=gvcf_file,
                        chr_mapping=chr_translation_file,
                        column_yaml_file=config)
