
import subprocess
from tempfile import NamedTemporaryFile
from multiprocessing import cpu_count

from Bio.SeqRecord import SeqRecord

from .. import ProLink_path


def align(muscle_input:str, muscle_output:str) -> None:
    '''
    Run a local alignment with MUSCLE v5

    Parameters
    ----------
    muscle_input : str
        Path of the input MUSCLE file
    muscle_output : str
        Path of the output MUSCLE file
    '''
    subprocess.call(['muscle', '-super5', muscle_input, '-output', muscle_output])


def weblogo3(weblogo_input:str, weblogo_output:str, format:str='png') -> None:
    '''
    Run WebLogo v3 to generate a sequence logo

    Parameters
    ----------
    filename_input : str
        Path of the input file (FASTA format from MUSCLE)
    weblogo_output : str
        Path of the output file
    format : str, optional
        Format of the output image (def: png)
    '''
    subprocess.call(['weblogo', '-D', 'fasta','-F', format, '-A', 'protein', '-s', 'large', '-n', '80', '-f', weblogo_input, '-o', weblogo_output])


def tree(tree_type:str, bootstrap_replications:int, muscle_output:str, mega_output:str) -> None:
    '''
    Run MEGA-CC to generate a phylogenetic tree

    Parameters
    ----------
    tree_type : str
        Type of tree to generate
    bootstrap_replications : int
        Number of bootstrap replications
    muscle_output : str
        Path of the input file (FASTA format from MUSCLE)
    mega_output : str
        Path of the MEGA-CC output file
    '''
    mega_config_input = f"{ProLink_path}/mega_configs/{tree_type}_{bootstrap_replications}.mao"
    subprocess.call(['megacc', '-a', mega_config_input, '-d', muscle_output, '-o', mega_output])


def blastp(seq_record:SeqRecord, blast_filename:str, threads:int=0, **kwargs) -> None:
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
    threads = threads or cpu_count()
    with NamedTemporaryFile(mode='w') as f:
        f.write(seq_record.format('fasta'))
        f.flush()
        additional_args = []
        for i, j in kwargs.items():
            additional_args.append(f'-{i}')
            additional_args.append(f'{j}')
        blastp = subprocess.run(['blastp', '-query', f.name, '-out', blast_filename, '-outfmt', '5', '-num_threads', str(threads)] + additional_args)
        if blastp.returncode != 0:
            raise Exception(f"Error running 'blastp' for {seq_record.id}")
