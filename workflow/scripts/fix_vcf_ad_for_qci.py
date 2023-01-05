
def fix_AD(vcf_in, vcf_out):
    header = True
    for line in vcf_in:
        if header:
            vcf_out.write(line)
            if line[:6] == "#CHROM":
                header = False
            continue
        columns = line.strip().split("\t")
        FORMAT = column[8].split(":")
        AD_index = 0
        AF_index = 0
        DP_index = 0
        VD_index = -1
        i = 0
        for f in FORMAT:
            if f == "AD":
                AD_index = i
            if f == "AF":
                AF_index = i
            if f == "DP":
                DP_index = i
            if f == "VD":
                VD_index = i
            i += 1
        DATA = columns[9].split(":")
        AD = DATA[AD_index].split(",")
        AF = float(DATA[AF_index])
        AF_AD = int(AD[1]) / (int(AD[0]) + int(AD[1]))
        DP = int(DATA[DP_index])
        VD = 0
        if VD_index != -1:
            VD = int(DATA[VD_index])
        AD2 = [1, 1]
        if VD != 0:
            AD2[0] = DP - VD
            AD2[1] = VD
        else:
            AD2[0] = DP - int(round(AF * DP, 0))
            AD2[1] = int(round(AF * DP, 0))
        AF2_AD = AD2[1] / (AD2[0] + AD2[1])
        if abs(AF2_AD - AF) > 0.01:
            print("AD and AF are to different!", AF, AF2_AD, line)
        A1 = AD2[0]
        A2 = AD2[1]
        AD2 = f"{A1},{A2}"
        DATA_fixed = ""
        i = 0
        for d in DATA:
            if i == AD_index:
                if i == 0:
                    DATA_fixed += f"{AD2}"
                else:
                    DATA_fixed += f":{AD2}"
            else:
                if i == 0:
                    DATA_fixed += f"{d}"
                else:
                    DATA_fixed += f":{d}"
            i += 1
        i = 0
        for column in columns:
            if i == 9:
                vcf_out.write(f"\t{DATA_fixed}")
            elif i == 0:
                vcf_out.write(f"{column}")
            else:
                vcf_out.write(f"\t{column}")
            i += 1
        vcf_out.write("\n")
    vcf_out.close()
    vcf_in.close()


vcf_in = open(snakemake.input.vcf)
vcf_out = open(snakemake.output.vcf, "w")
fix_AD(vcf_in, vcf_out)
