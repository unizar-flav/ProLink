
import json
import logging
from copy import copy

import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


logger = logging.getLogger()

def search_hmmer_pfam(seq:str) -> list[dict]:
    '''
    Use HMMER to search for the Pfam families of a sequence

    HMMER: Biosequence analysis using profile hidden Markov Models

    Parameters
    ----------
    seq : str
        Protein sequence to search for

    Returns
    -------
    list[dict]
        List of dictionaries with the Pfam families and information of them
        Common keys: 'flags', 'nregions', 'ndom', 'name', 'score', 'bias', 'taxid', 'acc',
                     'domains', 'nincluded', 'evalue', 'desc', 'pvalue', 'nreported', 'hindex'
    '''
    # request to HMMER
    hmmscan_url = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'
    header = {'Expect': '', 'Accept': 'application/json'}
    parameters = {
        'hmmdb': 'pfam',
        'seq': f">Seq\n{seq}"
        }
    response = requests.post(hmmscan_url, headers=header, data=parameters)
    # check response status code
    if response.status_code != 200:
        if response.status_code == 500:
            logger.error("ERROR: Bad response from the server side while searching for Pfam domains. Try later or blame 'ebi.ac.uk'")
        elif response.status_code == 400:
            logger.error(f"ERROR: Bad request while searching for Pfam domains. Check your query sequence: {seq}")
        else:
            logger.error(f"ERROR: Something went wrong while searching for Pfam domains. Error code: {response.status_code}")
        logger.debug(f"HMMER response:\n{response.text}")
        raise Exception(f"Something went wrong while searching for Pfam domains: Error {response.status_code}")
    # parse JSON response
    response_json = json.loads(response.text)
    return response_json['results']['hits']

def fasta_to_dfasta(seq_record:SeqRecord, fasta_input:str, fasta_output:str, pfam_output:str=None) -> None:
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
    pfam_output : str, optional
        Output file to write the Pfam domains of the sequences (def: None)
    '''
    def domain_names(domain_hits:list[dict]) -> list[str]:
        '''Get the names of the domains from the results of HMMER'''
        names = [domain['name'] for domain in domain_hits]
        if len(domain_hits) == 0:
            logger.warning(f"WARNING: No Pfam domains found for '{seq_record.id}'")
            names = ["No_domains_found"]
        elif len(domain_hits) > 1:
            logger.warning(f"WARNING: Multiple Pfam domain found for '{seq_record.id}'")
        return names
    # query sequence
    try:
        my_seq_domain_hits = search_hmmer_pfam(str(seq_record.seq))
    except Exception as e:
        logger.error(f"An error occurred while searching for Pfam domains of query sequence")
        raise e
    my_seq_domain_name = domain_names(my_seq_domain_hits)
    logger.info(f"Pfam domains found for query sequence '{seq_record.id}': {', '.join(my_seq_domain_name)}")
    # found sequences
    sequences_domain_hits = []
    sequences_domain = []
    #TODO: parallelize or request all sequences at once
    sequences = list(SeqIO.parse(fasta_input, "fasta"))
    for seq_record in sequences:
        try:
            seq_domain_hits = search_hmmer_pfam(str(seq_record.seq))
        except:
            logger.warning(f"WARNING: An error occurred while searching for Pfam domains of '{seq_record.id}'")
            seq_domain_hits = []
        sequences_domain_hits.append(seq_domain_hits)
        seq_domain_name = domain_names(seq_domain_hits)
        common_domains = set(my_seq_domain_name).intersection(seq_domain_name)
        if common_domains and len(common_domains) == len(my_seq_domain_name):
            domain = "Same_Domains"
        else:
            domain = "Different_Domains"
        seq = copy(seq_record)
        seq.id = f"{seq.description}---{domain}"
        seq.description = ""
        sequences_domain.append(seq)
    logger.debug(f"Writing sequences and Pfam domains to '{fasta_output}'")
    SeqIO.write(sequences_domain, fasta_output, "fasta")
    # write Pfam domains output
    if pfam_output:
        logger.debug(f"Writing Pfam domains to '{pfam_output}'")
        with open(pfam_output, "w") as f:
            f.write(f"#Query\tDomain (Accession)\n\n")
            name_acc = [f"{hit['name']} ({hit['acc']})" for hit in my_seq_domain_hits]
            f.write(f"{seq_record.id}\t{'   '.join(name_acc)}\n")
            f.write("\n\n")
            for seq_record, seq_domain_hits in zip(sequences, sequences_domain_hits):
                name_acc = [f"{hit['name']} ({hit['acc']})" for hit in seq_domain_hits]
                f.write(f"{seq_record.id}\t{'   '.join(name_acc)}\n")
