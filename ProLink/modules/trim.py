
import logging
from textwrap import dedent

from clipkit.api import clipkit


logger = logging.getLogger()

def trim_align(alignment_fastafile:str, alignment_trim_fastafile:str, mode:str='smart-gap') -> None:
    '''
    Trim poorly aligned positions and sequences from multiple alignment

    Parameters
    ----------
    alignment_fastafile : str
        Path of the input alignment file (FASTA format)
    alignment_trim_fastafile : str
        Path of the output trimmed alignment file (FASTA format)
    mode : str, optional
        Trimming mode (def: 'smart-gap')
        Modes: 'smart-gap', 'gappy', 'kpic', 'kpic-smart-gap', 'kpic-gappy', 'kpi', 'kpi-smart-gap', 'kpi-gappy'
    '''
    logger.info(f"\n-- Trimming alignment with ClipKIT\n")
    clipkit_stats = clipkit(
        input_file_path = alignment_fastafile,
        output_file_path = alignment_trim_fastafile,
        mode = mode
        )[1]
    logger.info(dedent(f"""\
        Trimming statistics:
        Alignment file:           {alignment_fastafile}
        Trimmed alignment file:   {alignment_trim_fastafile}
        Trimming mode:            {mode}
        Alignment lenght:         {clipkit_stats.alignment_length}
        Number of sites kept:     {clipkit_stats.output_length}
        Number of sites trimmed:  {clipkit_stats.trimmed_length}
        Percentage of alignment trimmed: {clipkit_stats.trimmed_percentage:.3f}%
        """))
