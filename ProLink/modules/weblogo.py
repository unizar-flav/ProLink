
import logging

import weblogo


logger = logging.getLogger()

def weblogo3(sequences_input:str, output:str, format:str='png', dpi:int=300) -> None:
    '''
    Generate a sequence logo with WebLogo v3

    Command line equivalent:
    weblogo -f sequences_input -o output -F format --stacks-per-line 80 --resolution dpi

    Parameters
    ----------
    filename_input : str
        Path of the input file (FASTA format from MUSCLE)
    output : str
        Path of the output file from WebLogo
    format : str, optional
        Format of the output image (def: png)
    dpi : int, optional
        Resolution of the output image (def: 300)
    '''
    dpi_fallback = 96
    with open(sequences_input) as f:
        seqs = weblogo.read_seq_data(f)
    logo_data = weblogo.LogoData.from_seqs(seqs)
    # try-except block to catch WebLogo failing because GhostScript messed up the high DPI
    try:
        logo_options = weblogo.LogoOptions(
            stacks_per_line = 80,
            resolution = dpi,
            )
        logo_format = weblogo.LogoFormat(logo_data, logo_options)
        graph = weblogo.formatters[format](logo_data, logo_format)
    except RuntimeError:
        logger.error(f"ERROR: WebLogo failed, trying with fallback DPI ({dpi_fallback})")
        logo_options = weblogo.LogoOptions(
            stacks_per_line = 80,
            resolution = dpi_fallback,
            )
        logo_format = weblogo.LogoFormat(logo_data, logo_options)
        graph = weblogo.formatters[format](logo_data, logo_format)
    with open(output, 'wb') as f:
        f.write(graph)
