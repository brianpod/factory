# path = input('Where is the FASTA file?: ')
# fasta_file = open(path)
# path = input('Where is the accession numbers?: ')

from pandas import *
from Bio import SeqIO

path = ''
path2 = ''

class Protein:
    def __init__(self, name, path=None):
        self.name = name
        self.path = path
        self.strains = []

    def add_strain(self, strain):
        self.strains.append(strain)

    def add_path(self, path):
        self.path = path


class Strain:
    def __init__(self, name, accession_number, strain):
        self.name = name
        self.accession_number = accession_number
        self.serotype = 'Null'
        self.strain = strain
        self.source = 'Null'

    def update_source(self, source):
        self.source = source

    def update_serotype(self, new_serotype):
        if self.serotype == 'Null':
            self.serotype = new_serotype
        else:
            self.serotype = self.serotype+' '+new_serotype

    def to_string(self):
        return self.accession_number + ' ' \
               + self.name \
               + ' serotype ' + self.serotype \
               + ' source ' + self.source \
               + ' str. ' + self.strain


def process_accession_numbers_from_csv(path):
    """
    proccesses default pandas csv into a list of Proteins of Strains
    :param path: path of the csv file
    :return: List of Proteins of Strains
    """
    data = read_csv(path, index_col=0)
    test = data.to_dict()
    filtered_list = []
    index = 0

    for protein, strain in test.items():
        idx_protein = Protein(protein)
        filtered_list.append(idx_protein)
        stuff = strain.items()
        for tuples in stuff:
            if type(tuples[1]) is str:

                filtered_list[index].add_strain(sanitized_strain_format(tuples[0], tuples[1]))
        index += 1

    return filtered_list


def update_proteins_path(proteins, paths):
    for x in range(len(proteins)):
        proteins[x].add_path(paths[x])
    return proteins


def sanitized_strain_format(name, accession):
    name = name.split()
    name1 = name[0]+' '+name[1]

    serotype = None

    if accession[-2] == ',':
        accession = accession[:-2]

    if 'str.' in name and 'serotype' not in name:
        strain = name[3]
    elif 'str.' in name and 'serotype' in name:
        strain = name[5]
        serotype = name[3]
    else:
        strain = name[2]

    result = Strain(name1, accession, strain)

    if serotype is not None:
        result.update_serotype(serotype)

    return result


def extract_serotype_and_strain_and_source_from_ncbi_pathogens_metadata(path):
    data = read_csv(path, sep='\t', header=0, index_col=False, usecols=['serovar', 'strain', 'isolation_source'])
    return data


def update_strains_with_serotype_and_source(database, strain):
    """
    pandas dataframe of database of serotype and source with strain objects and update them
    :param database:
    :param strain:
    :return:
    """
    row = database.loc[database['strain'] == strain.strain]
    if len(row) == 0:
        return
    source = str(row['isolation_source'].iloc[0])
    serovar = str(row['serovar'].iloc[0])
    if source == 'nan':
        pass
    else:
        strain.update_source(source)
    if serovar == 'nan':
        pass
    else:
        strain.update_serotype(serovar)


def process_proteins(protein):
    new_file = open(protein.name+'.txt', "a")
    for strain in protein.strains:
        for seq_record in SeqIO.parse(protein.path, "fasta"):
            if strain.accession_number == seq_record.id:
                seq_record.description = strain.to_string()
                SeqIO.write(seq_record, new_file, 'fasta')


proteins = process_accession_numbers_from_csv(path)
for z in proteins:
    data = extract_serotype_and_strain_and_source_from_ncbi_pathogens_metadata(path2)
    for strain in z.strains:
        update_strains_with_serotype_and_source(data, strain)
test_paths = ['']
proteins = update_proteins_path(proteins, test_paths)
for x in proteins:
    process_proteins(x)


