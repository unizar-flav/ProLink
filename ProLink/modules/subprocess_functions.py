
import subprocess
from tempfile import NamedTemporaryFile
from multiprocessing import cpu_count

from Bio.SeqRecord import SeqRecord

from .. import ProLink_path


def align(muscle_input:str, muscle_output:str):
    subprocess.call(['muscle', '-super5', muscle_input, '-output', muscle_output])

def weblogo3(weblogo_format:str, muscle_output:str, weblogo_output:str):
    subprocess.call(['weblogo', '-D', 'fasta','-F', weblogo_format, '-A', 'protein', '-s', 'large', '-n', '80', '-f', str(muscle_output),  '-o', str(weblogo_output)])

def tree(tree_type:str, bootstrap_replications:int, muscle_output:str, mega_output:str):
    mega_config_input = f"{ProLink_path}/mega_configs/{tree_type}_{bootstrap_replications}.mao"
    subprocess.call(['megacc', '-a', mega_config_input, '-d', muscle_output, '-o', mega_output])

def blastp(seq_record:SeqRecord, blast_filename:str, threads:int=None, **kwargs):
    '''
    Run locally 'blastp' for a single sequence record against a database

    Parameters
    ----------
    seq_record : SeqRecord
        Sequence to search
    blast_filename : str
        Path of the file to write the BLAST results (XML format)
    threads : int, optional
        Number of threads to use (def: all available)
    '''
    threads = str(threads or cpu_count())
    with NamedTemporaryFile(mode='w') as f:
        f.write(seq_record.format('fasta'))
        f.flush()
        additional_args = []
        for i, j in kwargs.items():
            additional_args.append(f'-{i}')
            additional_args.append(f'{j}')
        blastp = subprocess.run(['blastp', '-query', f.name, '-out', blast_filename, '-outfmt', '5', '-num_threads', threads ] + additional_args)
        if blastp.returncode != 0:
            raise Exception(f"Error running 'blastp' for {seq_record.id}")
