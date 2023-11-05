import os
import subprocess
import tempfile
from configure import db_list, exe

minimap2_exe = exe['minimap2']
samtools_exe = exe['samtools']
bcftools_exe = exe['bcftools']
esl_ref = db_list['esl_ref']
lineage_snp = db_list['lineage_snp']
variant_snp = db_list['variant_snp']


def parse_vcf(vcf_file):
    vcf_dict = {}
    with open(vcf_file, 'rt') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                parts = line.rstrip().split('\t')
                pos, ref_base, alt_base = parts[1], parts[3], parts[4]
                vcf_dict[pos] = [ref_base, alt_base]
    return vcf_dict


def esl_snp(query, hc):
    with tempfile.TemporaryDirectory() as temp_dir:
        mp_sam = os.path.join(temp_dir, "mp.sam")
        mp_bam = os.path.join(temp_dir, "mp.bam")
        sorted_bam = os.path.join(temp_dir, "sorted.bam")
        vcf_file = os.path.join(temp_dir, "vcf.file")

        step1 = f"{minimap2_exe} -a {esl_ref} {query} -o {mp_sam}"
        step2 = f"{samtools_exe} view -bS -o {mp_bam} {mp_sam}"
        step3 = f"{samtools_exe} sort {mp_bam} -o {sorted_bam}"
        step4 = f"{samtools_exe} index {sorted_bam}"
        step5 = f"{bcftools_exe} mpileup -Ou -f {esl_ref} {sorted_bam} | {bcftools_exe} call -mv -Ov -o {vcf_file}"

        for step in [step1, step2, step3, step4, step5]:
            subprocess.run(step, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        vcf_dict = parse_vcf(vcf_file)
        # print(vcf_dict)
    potential_lineage = {}
    with open(lineage_snp, 'rt') as lineage_file:
        for line in lineage_file:
            parts = line.rstrip().replace(' ', '\t').split('\t')
            lineage = parts[0].split('.SNP')[0]
            site, homoplasy, snp_pair = parts[2], parts[3], parts[4].split('->')
            if site in vcf_dict.keys():
                # print(snp_pair, site)
                if snp_pair == vcf_dict[site]:
                    # print(snp_pair, site)
                    if lineage.startswith(str(hc)):
                        if lineage in potential_lineage.keys():
                            potential_lineage[lineage] += 1
                        else:
                            potential_lineage[lineage] = 1

    phy = ['New lineage', 'New variant']
    if len(potential_lineage.keys()) < 1:
        phy = ['New lineage', 'New variant']
    else:
        max = 0
        pl, pv = '', ''
        for l in potential_lineage.keys():
            count = potential_lineage[l]
            if max <= count:
                max = count
                pl = l
            # print(l, count)
            if pl.startswith('2.4'):
                pv = pl
                pl = '2.4'
            else:
                potential_variant = {}
                with open(variant_snp, 'rt') as variant_file:
                    for vari in variant_file:
                        parts = vari.rstrip().replace(' ', '\t').split('\t')
                        variant = parts[0].split('.SNP')[0]
                        site, homoplasy, snp_pair = parts[2], parts[3], parts[4].split('->')
                        if site in vcf_dict:
                            if snp_pair == vcf_dict[site]:
                                if variant.startswith(pl):
                                    if variant in potential_variant.keys():
                                        potential_variant[variant] += 1
                                    else:
                                        potential_variant[variant] = 1
                if len(potential_lineage.keys()) < 1:
                    phy = [pl, 'New variant']
                else:
                    max = 0
                    for v in potential_variant.keys():
                        count = potential_variant[v]
                        if max <= count:
                            max = count
                            pv = v
                        # print(v, count)
        phy = [pl, pv]
    return phy