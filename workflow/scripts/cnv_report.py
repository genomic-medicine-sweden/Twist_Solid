
#cvn_vcf = open("/scratch2/wp1/nobackup/ngs/utveckling/analys/2022/hydra/UU_run10/cnv_sv/svdb_query/UP-VAL-80_T.svdb_query.vcf")
cnv_vcf_files = open("svdb_files.txt")
amp_genes = open("config/cnv_amplification_genes.tsv")
loh_genes = open("config/cnv_loh_genes.tsv")
out_amp = open("cnv_amp_report.tsv", "w")
out_loh = open("cnv_loh_report.tsv", "w")

def read_gene_list(gene_file):
    gene_dict = {}
    for line in gene_file:
        columns = line.strip().split("\t")
        gene = columns[0]
        chrom = columns[1]
        start_pos = int(columns[2])
        end_pos = int(columns[3])
        gene_dict[gene] = {"chrom": chrom, "start_pos": start_pos, "end_pos": end_pos}
    return gene_dict

def cnv_in_gene_list(gene_dict, chrom, start_pos, end_pos):
    found_gene = []
    for gene in gene_dict:
        if gene_dict[gene]["chrom"] == chrom:
            if ((start_pos >= gene_dict[gene]["start_pos"] and start_pos <= gene_dict[gene]["end_pos"]) or
                (end_pos >= gene_dict[gene]["start_pos"] and end_pos <= gene_dict[gene]["end_pos"]) or
                (start_pos <= gene_dict[gene]["start_pos"] and end_pos >= gene_dict[gene]["end_pos"]) or
                (start_pos >= gene_dict[gene]["start_pos"] and end_pos <= gene_dict[gene]["end_pos"])):
                found_gene.append(gene)
    return found_gene

amp_gene_dict = read_gene_list(amp_genes)
loh_gene_dict = read_gene_list(loh_genes)


for file in cnv_vcf_files:
    cvn_vcf = open(file.strip())
    sample = file.split("/")[-1].split(".")[0]

    header = True
    amp_gene_sample_dict = {}
    loh_gene_sample_dict = {}
    amp_list = []
    loh_list = []
    for line in cvn_vcf:
        if header:
            if line[:6] == "#CHROM" :
                header = False
            continue
        columns = line.strip().split("\t")
        chrom = columns[0]
        start_pos = int(columns[1])
        end_pos = int(columns[7].split("END=")[1].split(";")[0])
        amp_found = cnv_in_gene_list(amp_gene_dict, chrom, start_pos, end_pos)
        loh_found = cnv_in_gene_list(loh_gene_dict, chrom, start_pos, end_pos)
        caller = ""
        cn = 2.0
        if line.find("CORRECTED_COPY_NUMBER=") != -1:
            cn = float(columns[7].split("CORRECTED_COPY_NUMBER=")[1].split(";")[0])
            caller = "gatk_cnv"
        else :
            cn = float(columns[7].split("COPY_NUMBER=")[1].split(";")[0])
            caller = "cnvkit"
        twist_af = 0.0
        if line.find("Twist_AF=") != -1:
            twist_af = float(columns[7].split("Twist_AF=")[1].split(";")[0])
        if twist_af > 0.15:
            continue
        sv_len = float(columns[7].split("SVLEN=")[1].split(";")[0])
        probes = int(columns[7].split("PROBES=")[1].split(";")[0])
        baf = columns[7].split("BAF=")[1].split(";")[0]
        if baf == "":
            baf = "NA"
        else:
            baf = float(baf)
        baf_probes = "NA"
        if line.find("BAF_PROBES=") != -1:
            baf_probes = int(columns[7].split("BAF_PROBES=")[1].split(";")[0])
        if amp_found != []:
            genes = ','.join(amp_found)
            for gene in genes.split(","):
                if gene not in amp_gene_sample_dict:
                    amp_gene_sample_dict[gene] = {}
                amp_gene_sample_dict[gene][caller] = ""
            if cn >= 4.0:
                amp_list.append([sample, caller, genes, chrom, start_pos, end_pos, cn, twist_af, sv_len, probes, baf, baf_probes])
        if loh_found != []:
            genes = ','.join(loh_found)
            for gene in genes.split(","):
                if gene not in loh_gene_sample_dict:
                    loh_gene_sample_dict[gene] = {}
                loh_gene_sample_dict[gene][caller] = ""
            if cn <= 1.3 and twist_af < 0.05 and sv_len < 1000000:
                loh_list.append([sample, caller, genes, chrom, start_pos, end_pos, cn, twist_af, sv_len, probes, baf, baf_probes])

    for amp in amp_list:
        found_gene = False
        for gene in amp[2].split(","):
            if gene in amp_gene_sample_dict:
                if ((caller == "cnvkit" and "gatk_cnv" in amp_gene_sample_dict[gene]) or
                    (caller == "gatk_cnv" and "cnvkit" in amp_gene_sample_dict[gene])):
                    found_gene = True
        if found_gene:
            print(amp[0])
            out_amp.write(amp[0])
            for a in amp[1:]:
                out_amp.write("\t" + str(a))
            out_amp.write("\n")
    for loh in loh_list:
        found_gene = False
        for gene in loh[2].split(","):
            if gene in loh_gene_sample_dict:
                if ((caller == "cnvkit" and "gatk_cnv" in loh_gene_sample_dict[gene]) or
                    (caller == "gatk_cnv" and "cnvkit" in loh_gene_sample_dict[gene])):
                    found_gene = True
        if found_gene:
            out_loh.write(loh[0])
            for a in loh[1:]:
                out_loh.write("\t" + str(a))
            out_loh.write("\n")
