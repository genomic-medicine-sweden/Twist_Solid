from hydra_genetics.utils.io.hotspot_report import generate_hotspot_report

hotspot_file = snakemake.input.hotspots
vcf_file = snakemake.input.vcf
gvcf_file = snakemake.input.gvcf

report = snakemake.output.report

levels = snakemake.params.levels
config = snakemake.params.report_config
sample_name = snakemake.params.sample_name
chr_translation_file = snakemake.params.chr_translation_file


generate_hotspot_report(sample_name,
                        report,
                        levels,
                        hotspot_file,
                        vcf_file,
                        gvcf_file,
                        chr_translation_file,
                        config)
