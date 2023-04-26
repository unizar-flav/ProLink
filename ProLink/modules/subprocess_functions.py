
import subprocess

from .. import ProLink_path


def align(muscle_input:str, muscle_output:str):
    subprocess.call(['muscle', '-super5', muscle_input, '-output', muscle_output])

def weblogo3(weblogo_format:str, muscle_output:str, weblogo_output:str):
    subprocess.call(['weblogo', '-D', 'fasta','-F', weblogo_format, '-A', 'protein', '-s', 'large', '-n', '80', '-f', str(muscle_output),  '-o', str(weblogo_output)])

def tree(tree_type:str, bootstrap_replications:int, muscle_output:str, mega_output:str):
    mega_config_input = f"{ProLink_path}/mega_configs/{tree_type}_{bootstrap_replications}.mao"
    subprocess.call(['megacc', '-a', mega_config_input, '-d', muscle_output, '-o', mega_output])
