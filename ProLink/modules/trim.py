
import logging

from clipkit.clipkit import execute as clipkit_execute
from clipkit.helpers import SeqType
from clipkit.modes import TrimmingMode


logger = logging.getLogger()

def trim_align(alignment_fastafile:str, alignment_trim_fastafile:str) -> None:
    '''
    Trim poorly aligned positions and sequences from multiple alignment

    Parameters
    ----------
    alignment_fastafile : str
        Path of the input alignment file (FASTA format)
    alignment_trim_fastafile : str
        Path of the output trimmed alignment file (FASTA format)
    '''
    logger.info(f"\n-- Trimming alignment with ClipKIT\n")
    clipkit_execute(
        input_file = alignment_fastafile,
        input_file_format = 'fasta',
        output_file = alignment_trim_fastafile,
        output_file_format = 'fasta',
        sequence_type = SeqType('aa'),
        gaps = 0.9,
        complement = False,
        mode = TrimmingMode('smart-gap'),
        use_log = True
        )
