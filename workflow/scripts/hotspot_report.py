from hydra_genetics.utils.io.filtered_mutations import generate_filtered_mutations

hotspot_file = snakemake.input.hotspot
vcf_file = snakemake.input.vcf
gvcf_file = snakemake.input.gvcf

report = snakemake.output.report

levels = snakemake.params.levels
config = snakemake.params.report_config
sample_name = snakemake.params.sample_name
chr_translation_file = None #snakemake.params.chr_translation_file



generate_filtered_mutations(sample_name,
                            report,
                            levels,
                            hotspot_file
                            vcf_file,
                            gvcf_file,
                            chr_translation_file,
                            config)
