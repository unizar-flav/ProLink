"""
  Simple interface to Pfam families
  =================================

"""

import urllib


def search_hmmer_pfam(seq) -> dict:
    """
        Use HMMER to search for the Pfam families of a sequence

        HMMER: Biosequence analysis using profile hidden Markov Models

        Adapted from ProDy: prody.database.pfam @ github.com/ProDy/ProDy

        Parameters
        ----------
        seq : str
            protein sequence to search for

        Returns
        -------
        dict
            dictionary with the Pfam families and information of them
    """
    parameters = {
        'hmmdb': 'pfam',
        'seq': f">Seq\n{seq}"}
    request = urllib.request.Request(
        'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan',
        urllib.parse.urlencode(parameters).encode('utf-8'))
    results_url = urllib.request.urlopen(request).geturl()
    res_params = {'format' : 'tsv' }
    modified_res_url = results_url.replace('results', 'download') + '?' + urllib.parse.urlencode(res_params)
    result_request = urllib.request.Request(modified_res_url)
    try:
        tsv = urllib.request.urlopen(result_request).read().decode()
    except:
        raise ValueError('No matching Pfam domains were found.')
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
