

import statistics
import logging

log = logging.getLogger()


aa_table = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H",
    "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*"
    }


def translate_aa(AA_org):
    if len(AA_org) == 3:
        if AA_org in aa_table:
            return aa_table[AA_org]
        else:
            return AA_org
    elif len(AA_org) < 4:
        return AA_org
    else:
        AA1_org = AA_org[:3]
        AA2_org = AA_org[-3:]
        AA1_new = AA1_org
        AA2_new = AA2_org
        if AA1_org in aa_table:
            AA1_new = aa_table[AA1_org]
        if AA2_org in aa_table:
            AA2_new = aa_table[AA2_org]
        AA1_len = len(AA1_new)
        AA2_len = len(AA2_new)
        if AA1_len == 1 and AA2_len == 1:
            return f"{AA1_new}{AA_org[3:-3]}{AA2_new}"
        elif AA1_len == 1:
            return f"{AA1_new}{AA_org[3:]}"
        elif AA2_len == 1:
            return f"{AA_org[:-3]}{AA2_new}"
        else:
            return AA_org


def translate_file(input_file, output_file):
    header_dict = {}
    header = True
    for line in input_file:
        if header:
            i = 0
            for h in line.strip().split("\t"):
                header_dict[h] = i
                i += 1
            header = False
            output_file.write(line)
            continue
        columns = line.strip().split("\t")
        AA_change_index = header_dict["AA_change"]
        AA_org = columns[AA_change_index]
        AA_new = ""
        if AA_org.find("fs") != -1:
            AA_new = f"{translate_aa(AA_org.split('fs')[0])}fs{translate_aa(AA_org.split('fs')[1])}"
        elif AA_org.find("delins") != -1:
            if len(AA_org.split("_")) > 1:
                delins = AA_org.split('delins')[0].split('_')
                AA_new = f"{translate_aa(delins[0])}_{translate_aa(delins[1])}delins"
            else:
                AA_new = f"{translate_aa(AA_org.split('delins')[0])}delins"
            i = 0
            str_len = len(AA_org.split('delins')[1])
            while i+3 <= str_len:
                AA_new += f"{translate_aa(AA_org.split('delins')[1][i:i+3])}"
                i += 3
        elif AA_org.find("del") != -1:
            if len(AA_org.split("_")) > 1:
                deletion = AA_org.split('del')[0].split('_')
                AA_new = f"{translate_aa(deletion[0])}_{translate_aa(deletion[1])}del"
            else:
                AA_new = f"{translate_aa(AA_org.split('del')[0])}del"
        elif AA_org.find("ins") != -1:
            insertion = AA_org.split('ins')[0].split('_')
            AA_new = f"{translate_aa(insertion[0])}_{translate_aa(insertion[1])}ins"
            i = 0
            str_len = len(AA_org.split('ins')[1])
            while i+3 <= str_len:
                AA_new += f"{translate_aa(AA_org.split('ins')[1][i:i+3])}"
                i += 3
        elif AA_org.find("dup") != -1:
            if len(AA_org.split("_")) > 1:
                duplication = AA_org.split('dup')[0].split('_')
                AA_new = f"{translate_aa(duplication[0])}_{translate_aa(duplication[1])}dup"
            else:
                AA_new = f"{translate_aa(AA_org.split('dup')[0])}dup"
        elif AA_org.find("ext") != -1:
            AA_new = f"{translate_aa(AA_org.split('ext')[0])}ext{translate_aa(AA_org.split('ext')[1])}"
        else:
            AA_new = translate_aa(AA_org)
        i = 0
        line_new = ""
        for column in columns:
            if i == AA_change_index:
                column = AA_new
            if i == 0:
                line_new += column
            else:
                line_new += f"\t{column}"
            i += 1
        line_new += "\n"
        output_file.write(line_new)


if __name__ == "__main__":
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    translate_file(
        open(snakemake.input.report),
        open(snakemake.output.report, "w"),
    )
