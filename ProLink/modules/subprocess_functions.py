
import subprocess


def align(muscle_input, muscle_output):
    print("in align")
    print(muscle_input)
    print(muscle_output)
    subprocess.call(['muscle', '-super5', muscle_input, '-output', muscle_output])

def weblogo3(weblogo_format, muscle_output, weblogo_output):
    subprocess.call(['weblogo', '-D', 'fasta','-F', weblogo_format, '-A', 'protein', '-s', 'large', '-n', '80', '-f', str(muscle_output),  '-o', str(weblogo_output)])

def tree(tree_type, bootstrap_replications, muscle_output, mega_output):
    mega_config_input = "/ProLink/ProLink/mega_configs/" + tree_type + "_" + bootstrap_replications + ".mao"
    print(mega_config_input)
    mega_config_input = f"/ProLink/ProLink/mega_configs/{tree_type}_{bootstrap_replications}.mao"
    print(mega_config_input)
    subprocess.call(['megacc', '-a', mega_config_input, '-d', muscle_output, '-o', mega_output])
