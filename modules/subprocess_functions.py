import subprocess

def align(muscle_input, muscle_output):
    subprocess.call(['muscle', '-super5', muscle_input, '-output', muscle_output])

def weblogo3(weblogo_format, muscle_output, weblogo_output):
    subprocess.call(['weblogo', '-D', 'fasta','-F', weblogo_format, '-A', 'protein', '-s', 'large', '-f', str(muscle_output),  '-o', str(weblogo_output)])

def tree(mega_config_input, muscle_output, mega_output):
    subprocess.call(['megacc', '-a', mega_config_input, '-d', muscle_output, '-o', mega_output])
