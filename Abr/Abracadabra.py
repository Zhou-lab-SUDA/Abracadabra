import argparse, os, gzip
import datetime
from cgMLST import cgmlst
from Species import species
from AMRsearch import amrsearch
from SNP import esl_snp

def write_cgMLST(prefix, alleles) :
    fout = '{0}.cgMLST.profile.gz'.format(prefix)
    with gzip.open(fout, 'wt') as fout :
        fout.write('#Query\t{0}\n'.format('\t'.join([g for g in sorted(alleles.keys())])))
        fout.write('{0}\t{1}\n'.format(prefix, '\t'.join([allele.get('value_md5', '-') for g, allele in sorted(alleles.items())])))

def spell():
    def check_value(value):
        if value not in ['s', 'r', 'c', 'sr', 'sc', 'rs', 'rc', 'cs', 'cr', 'src', 'scr', 'rsc', 'rcs', 'csr', 'crs']:
            raise argparse.ArgumentTypeError(
                'Mortal! You can only beseech for predictions of s(pecies), r(esistance), and c(gMLST)!')
        return value

    parser = argparse.ArgumentParser(
        description='Abracadabra, an ancient incantation, unveils the mystic power concealed within Acinetobacter baumannii genomes.',
        add_help=True,
        usage='''In the realm of code, behold Abracadabra, the enchanting script, that holds the key to unlocking the enigmatic potential within Acinetobacter baumannii genomes.'''
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '-q', '--query',
        type=str,
        help='''-q or --query: Your genome offering for analysis.'''
    )
    group.add_argument(
        '-ql', '--query_list',
        type=str,
        help='''-ql or --query_list: A list of genome offerings for analysis.'''
    )
    parser.add_argument(
        '-p', '--prefix',
        required=False,
        default='Abr',
        type=str,
        help='''-p or --prefix: The mystical prefix for your results. Default as Abr.'''
    )
    parser.add_argument(
        '-n', '--num_threads',
        type=int,
        default=8,
        help='''-n or --num_threads: Pixies can accelerate the magica, more you use, faster spell will be. (Default as 8, only affect cgMLST spell.)'''
    )
    parser.add_argument(
        '-t', '--tasks',
        type=check_value,
        default='src',
        help='''-t or --tasks: Choose your desires - s for species identification, r for resistance prediction, and c for cgMLST assignment. Default as src.
        Example: -t crs would kindly guide you to obtain the cgMLST, drug resistance, and species for your offered genome.'''
    )

    args = parser.parse_args()

    if args.prefix:
        work_dir = args.prefix.rsplit('/', 1)[0] if len(args.prefix.rsplit('/', 1)) > 1 else os.path.join(os.getcwd(), args.prefix)
    else:
        work_dir = os.path.join(os.getcwd(), 'abr')

    tasks = args.tasks
    if 's' in tasks:
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ '\t---In the realm of magic, species identification begins---')
        species_judge = species(args.query)
        [sp, ani] = species_judge
        if sp == True:
            print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ '\tThe offering genome is revealed as Acinetobacter baumannii')
        else:
            print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ '\tThe offering genome is not Acinetobacter baumannii')
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ '\tAverage Nucleotide Identity: ' + ani)
    if 'r' in tasks:
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ '\t---In the realm of magic, researchers will have ARGs and virulent prediction they pray---')
        amr_info = amrsearch(args.query)
        amr_info.to_csv(work_dir + '.amr', sep='\t', index=False)
    if 'c' in tasks:
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ '\t---In the realm of magic, cgMLST shows the deep connection inside the bacteria---')
        esl = False
        alleles, hiercc = cgmlst(args.query, n_thread=args.num_threads)
        most_closed, allele_differed, hc = hiercc['reference'], hiercc['distance'], hiercc['HC1030']
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ '\tHC1030: '+str(hc)+'\tMost closed genome: '+most_closed+'\tAlleles Differed: '+str(allele_differed))
        if hc == 1 or 2:
            esl = True
            print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ '\t***Behold, the malevolent genome of the ESL has been unveiled!***\n\t\tLet the ancient art of SNP analysis unveil its lineage and variant!')
            [pl, pv] = esl_snp(args.query, hc)
            print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ '\tESL_lineage: '+pl + '\tESL_variant: '+pv)
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ '\t---Writing cgMLST---')
        write_cgMLST(work_dir, alleles)
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+ '\t---Magic fades away leaving knowledge for Acinetobacter baumannii---')


if __name__ == '__main__':
    spell()
