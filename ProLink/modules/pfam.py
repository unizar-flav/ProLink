
import urllib

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def search_hmmer_pfam(seq:str) -> dict:
    '''
    Use HMMER to search for the Pfam families of a sequence

    HMMER: Biosequence analysis using profile hidden Markov Models

    Adapted from ProDy: prody.database.pfam @ github.com/ProDy/ProDy

    Parameters
    ----------
    seq : str
        Protein sequence to search for

    Returns
    -------
    dict
        Dictionary with the Pfam families and information of them
    '''
    parameters = {
        'hmmdb': 'pfam',
        'seq': f">Seq\n{seq}"}
    request = urllib.request.Request(
        'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan',
        urllib.parse.urlencode(parameters).encode('utf-8'))
    results_url = urllib.request.urlopen(request).geturl()
    res_params = {'format': 'tsv'}
    modified_res_url = results_url.replace('results', 'download') + '?' + urllib.parse.urlencode(res_params)
    result_request = urllib.request.Request(modified_res_url)
    try:
        tsv = urllib.request.urlopen(result_request).read().decode()
    except:
        raise Exception('No matching Pfam domains were found.')
    lines = tsv.split('\n')
    keys = lines[0].split('\t')
    root = dict()
    for i, line in enumerate(lines[1:-1]):
        root[i] = dict()
        for j, key in enumerate(keys):
            root[i][key] = line.split('\t')[j]
    matches = dict()
    for child in root.values():
        accession = child['Family Accession']
        pfam_id = accession.split('.')[0]
        matches[pfam_id] = dict()
        matches[pfam_id]['accession'] = accession
        matches[pfam_id]['class'] = 'Domain'
        matches[pfam_id]['id'] = child['Family id']
        matches[pfam_id]['locations'] = dict()
        matches[pfam_id]['locations']['ali_end'] = child['Ali. End']
        matches[pfam_id]['locations']['ali_start'] = child['Ali. Start']
        matches[pfam_id]['locations']['bitscore'] = child['Bit Score']
        matches[pfam_id]['locations']['end'] = child['Env. End']
        matches[pfam_id]['locations']['cond_evalue'] = child['Cond. E-value']
        matches[pfam_id]['locations']['ind_evalue'] = child['Ind. E-value']
        matches[pfam_id]['locations']['evidence'] = 'hmmer v3.0'
        matches[pfam_id]['locations']['hmm_end'] = child['Model End']
        matches[pfam_id]['locations']['hmm_start'] = child['Model Start']
        matches[pfam_id]['locations']['start'] = child['Env. Start']
        matches[pfam_id]['type'] = 'Pfam-A'
    return matches


def fasta_to_dfasta(seq_record:SeqRecord, fasta_input:str, fasta_output:str) -> None:
    '''
    Find the Pfam domain of a sequence and compare them with the Pfam domains of the sequences in a fasta file

    Parameters
    ----------
    seq_record : SeqRecord
        Sequence to search
    fasta_input : str
        Input fasta file to read
    fasta_output : str
        Output fasta file to write
    '''
    print("Checking Pfam domains")
    try:
        my_sequence_domains = search_hmmer_pfam(str(seq_record.seq)).keys()
    except:
        raise Exception('No matching Pfam domains were found for the query sequence.')
    d_sequences = []
    for seq_record in SeqIO.parse(fasta_input, "fasta"):
        try:
            subject_sequence_domains = search_hmmer_pfam(str(seq_record.seq)).keys()
        except:
            subject_sequence_domains = "No_domains_found"
        if my_sequence_domains != subject_sequence_domains:
            domains = "DD:" + str(subject_sequence_domains).replace("dict_keys", "").replace("([", "").replace("])", "").replace("'", "")
        else:
            domains = "SD"
        rec_d = SeqRecord(Seq(str(seq_record.seq)),
                          id=seq_record.id + seq_record.description + "|" + domains,
                          description="")
        string = rec_d.id
        new_string = string[:string.find("|_") + 1].replace("ref", "") + string[string.find("[") + 1:string.find("]")] + string[string.find("]|") + 1:]
        rec_d.id = new_string
        rec_d.description = ""
        d_sequences.append(rec_d)
    SeqIO.write(d_sequences, fasta_output, "fasta")
