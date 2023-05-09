
import logging
import subprocess
from multiprocessing import cpu_count
from tempfile import NamedTemporaryFile

from Bio.SeqRecord import SeqRecord

from .. import ProLink_path


logger = logging.getLogger()

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
    logging.info(f"\n-- Aligning sequences with MUSCLE")
    muscle_cmd = ['muscle', '-super5', muscle_input, '-output', muscle_output]
    logging.debug(f"Running MUSCLE alignment: {' '.join(muscle_cmd)}")
    muscle_run = subprocess.run(muscle_cmd)
    if muscle_run.returncode != 0:
        logger.error(f"ERROR: MUSCLE failed")
        raise RuntimeError(f"MUSCLE failed")

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
    logging.info(f"\n-- Generating phylogenetic tree with MEGA-CC")
    mega_cmd = ['megacc', '-a', mega_config_input, '-d', muscle_output, '-o', mega_output]
    logging.debug(f"Running MEGA-CC: {' '.join(mega_cmd)}")
    mega_run = subprocess.run(mega_cmd)
    if mega_run.returncode != 0:
        logger.error(f"ERROR: MEGA-CC failed")
        raise RuntimeError(f"MEGA-CC failed")

def blastp_run(seq_record:SeqRecord, blast_filename:str, threads:int=0, **kwargs) -> None:
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
        blastp_cmd = ['blastp', '-query', f.name, '-out', blast_filename, '-outfmt', '5', '-num_threads', str(threads)] + additional_args
        logger.debug(f"Running 'blastp' for {seq_record.id}:  {' '.join(blastp_cmd)}")
        blastp_run = subprocess.run(blastp_cmd)
        if blastp_run.returncode != 0:
            logger.error(f"ERROR: local 'blastp' for {seq_record.id} failed")
            raise RuntimeError(f"Error running 'blastp' for {seq_record.id}")
