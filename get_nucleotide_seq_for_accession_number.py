from Bio import Entrez

Entrez.email = ''


def get_nucleotide_accession_and_start_and_stop_from_protein_accession(accession_number: str)->(str, str, str):
    handle = Entrez.efetch("ipg", id=accession_number, rettype='ipg', retmode='xml')
    raw_read = Entrez.parse(handle)

    nucleotide_accession = ''
    start = ''
    stop = ''

    for element in raw_read:
        nucleotide_accession = element['ProteinList'][0]['CDSList'][0].attributes['accver']
        start = element['ProteinList'][0]['CDSList'][0].attributes['start']
        stop = element['ProteinList'][0]['CDSList'][0].attributes['stop']

    return nucleotide_accession, start, stop


def get_nucleotide_from_nucleotide_accession_and_start_and_stop(nucleotide_accession: str, start: str, stop: str) -> str:

    handle = Entrez.efetch('nucleotide', id=nucleotide_accession, retmode='xml', seq_start=start, seq_stop=stop)
    raw_read = Entrez.parse(handle)
    nucleotide = None

    for element in raw_read:

        nucleotide = element['GBSeq_sequence']

    return nucleotide


def get_genome_from_nucleotide_accession(nucleotide_accession):

    handle = Entrez.efetch('nucleotide', id=nucleotide_accession, report='fasta', retmode='xml')
    raw_read = Entrez.parse(handle)

    nucleotide = None

    for element in raw_read:
        nucleotide = element['TSeq_sequence']

    return nucleotide


def get_genome_from_protein_accession(protein_accession):
    blah = get_nucleotide_accession_and_start_and_stop_from_protein_accession(protein_accession)
    nucleotide_seq = get_genome_from_nucleotide_accession(blah[0])
    return nucleotide_seq


def nucleotide_sequence_from_protein_accession(protein_accession):
    blah = get_nucleotide_accession_and_start_and_stop_from_protein_accession(protein_accession)
    nucleotide_seq = get_nucleotide_from_nucleotide_accession_and_start_and_stop(blah[0], blah[1], blah[2])
    return nucleotide_seq
