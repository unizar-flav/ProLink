
import weblogo


def weblogo3(sequences_input:str, output:str, format:str='png') -> None:
    '''
    Generate a sequence logo with WebLogo v3

    Parameters
    ----------
    filename_input : str
        Path of the input file (FASTA format from MUSCLE)
    output : str
        Path of the output file from WebLogo
    format : str, optional
        Format of the output image (def: png)
    '''
    with open(sequences_input) as f:
        seqs = weblogo.read_seq_data(f)
    logo_data = weblogo.LogoData.from_seqs(seqs)
    logo_options = weblogo.LogoOptions(
        stacks_per_line = 80,
        resolution = 300,
        )
    logo_format = weblogo.LogoFormat(logo_data, logo_options)
    graph = weblogo.formatters[format](logo_data, logo_format)
    with open(output, 'wb') as f:
        f.write(graph)
