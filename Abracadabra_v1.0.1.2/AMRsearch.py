import numpy as np
import pandas as pd
import pickle
from configure import db_list, sulbac_list, amr_prot_dict

dbs = db_list['amr_dbs']


def get_tops(bsn):
    bsn = bsn[np.argsort(-bsn.T[11] * bsn.T[3] / bsn.T[12])]
    covs = {c: [] for c in np.unique(bsn.T[1])}
    for p in bsn:
        s, e = p[8:10] if p[8] < p[9] else (p[9], p[8])
        x = np.zeros(e - s + 1, dtype=np.uint8)
        for c in covs[p[1]]:
            if c[0] > e or c[1] < s:
                continue
            ovl = [max(c[0], s), min(c[1], e)]
            if ovl[1] - ovl[0] + 1 >= 0.3 * (c[1] - c[0] + 1):
                p[1] = ''
                break
            x[ovl[0] - s:ovl[1] - s + 1] = 1
            if np.sum(x) >= 0.3 * (e - s + 1):
                p[1] = ''
                break
        else:
            covs[p[1]].append([s, e])
    return bsn[bsn.T[1] != '']


def amrsearch(query):
    import subprocess, os
    from uberBlast import uberBlast, readFastq
    try:
        gene_info = ['contig', 'accession_number', 'identity', 'coverage', 'mismatch', 'gap', 'start', 'end', 'evalue',
                     'bit_score', 'gene', 'function', 'function_category', 'description']
        amr_df = pd.DataFrame(columns=gene_info)
        for db, (db_name, site) in dbs.items():
            if db.startswith('nucl'):
                sul_list = []
                with open(amr_prot_dict, 'rb') as apd:
                    amp_dict = pickle.load(apd)
                with open(sulbac_list, 'rt') as sulbac:
                    for sul in sulbac.readlines():
                        sul_list.append(sul)
                bsn = uberBlast(
                    '-r {0} -q {1} -f --blastn --min_id {2} --min_ratio {3} -t {4} -p -s 2 -e 0,3 -m --merge_gap 300'.format(
                        query, db_name, 90 / 100., 70 / 100., 1).split())
                bsn = get_tops(bsn)
                for i in bsn:
                    gene_det = i[0].split('|')
                    query_contig = i[1]
                    ide = i[2] * 100
                    cov_len = i[3]
                    seq_len = i[6]
                    cov = round(cov_len / seq_len, 2) * 100
                    mismatch = i[4]
                    gap = i[5]
                    start = i[7]
                    end = i[8]
                    evalue = i[9]
                    score = i[10]
                    prot_acc = gene_det[1]
                    gene_name = gene_det[5]
                    gene_desc = gene_det[7]
                    drug = 'Sulbactam' if gene_name in sul_list else amp_dict[prot_acc][-2]
                    drug_cat = amp_dict[prot_acc][-1]
                    amr_df = pd.concat([amr_df, pd.Series(
                        [query_contig, prot_acc, ide, cov, mismatch, gap, start, end, evalue, score, gene_name, drug,
                         drug_cat, gene_desc], index=gene_info).T], ignore_index=False)
            elif db.startswith('vir'):
                bsn = uberBlast(
                    '-r {0} -q {1} -f --diamondx --min_id {2} --min_ratio {3} -t {4} -p -s 2 -e 0,3 -m --merge_gap 300'.format(
                        query, db_name, 80 / 100., 70 / 100., 1).split())
                bsn = get_tops(bsn)
                for i in bsn:
                    gene_det = i[0]
                    query_contig = i[1]
                    ide = i[2] * 100
                    cov_len = i[3]
                    seq_len = i[6]
                    cov = round(cov_len / seq_len, 2) * 100
                    mismatch = i[4]
                    gap = i[5]
                    start = i[7]
                    end = i[8]
                    evalue = i[9]
                    score = i[10]
                    prot_acc = gene_det.split(' ')[0].split('|')[1].replace(')', '')
                    # print(gene_det)
                    gene_name = gene_det.split(' ')[1].replace(')', '').replace('(', '') if len(gene_det.split(' '))>1 else gene_det.split(' ')[0].split('|')[1].replace('(', '')
                    gene_desc = gene_det.split(gene_name)[1].replace(') ', '').split('-')[0]
                    viru = 'virulent'
                    drug_cat = gene_det.split(gene_name)[1].replace(') ', '').split('-')[0].split(' (')[0]
                    amr_df = pd.concat([amr_df, pd.Series(
                        [query_contig, prot_acc, ide, cov, mismatch, gap, start, end, evalue, score, gene_name, viru,
                         drug_cat, gene_desc], index=gene_info).T], ignore_index=False)
        return amr_df
    except subprocess.CalledProcessError as e:
        print(f'Error running BLAST: {e}')
        exit(1)
