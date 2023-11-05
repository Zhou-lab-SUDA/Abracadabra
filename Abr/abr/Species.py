
from configure import db_list, exe

fastANI = exe['fastANI']
reference_genome = db_list['species_ref']

def species(query_genome):
    import tempfile
    import os
    from subprocess import run, PIPE
    
    with tempfile.TemporaryDirectory() as temp_dir:
        ani_output = os.path.join(temp_dir, "ani_output.txt")
        ani_command = f"{fastANI} -q {query_genome} -r {reference_genome} -o {ani_output}"
        run(ani_command, shell=True, stdout=PIPE, stderr=PIPE, text=True)
        if os.path.getsize(ani_output)==0:
            return [False, '<75% Average Nucleotide Identity']
        else:
            with open(ani_output, 'r') as ani_file:
                for line in ani_file:
                    part = line.rstrip().split('\t')
                    ani = part[2]
                    if float(ani) >= 95:
                        return [True, ani]
                    else:
                        return [False, ani]
                    
                    
    