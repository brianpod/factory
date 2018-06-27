class Sequence:
    __slots__ = ['sequence', 'thymine', 'guanine', 'adenine', 'cytosine']

    def __init__(self, sequence):
        self.sequence = sequence.lower()
        self.thymine = self.count_thymine()
        self.guanine = self.count_guanine()
        self.adenine = self.count_adenine()
        self.cytosine = self.count_cytosine()
        self.destroy_sequence()

    def count_thymine(self):
        return self.sequence.count('t')

    def count_guanine(self):
        return self.sequence.count('g')

    def count_adenine(self):
        return self.sequence.count('a')  

    def count_cytosine(self):
        return self.sequence.count('c')

    def gc_ratio(self):
        return (self.guanine+self.cytosine) / (self.adenine+self.thymine+self.guanine+self.cytosine)

    def at_ratio(self):
        return (self.adenine + self.thymine) / (self.guanine + self.cytosine+self.adenine+self.thymine)

    def destroy_sequence(self):
        del self.sequence


class Organism:
    __slots__ = ['name', 'genome', 'gene']

    def __init__(self, name, gene, genome):
        self.name = name
        self.gene = gene
        self.genome = genome


from Bio import SeqIO
import get_nucleotide_seq_for_accession_number


def read_fasta_file(path: str) -> []:
    organisms = []
    for record in SeqIO.parse(path, "fasta"):
        print(record.description)
        organisms.append(Organism(
            record.description,
            Sequence(get_nucleotide_seq_for_accession_number.nucleotide_sequence_from_protein_accession(record.name)),
            Sequence(get_nucleotide_seq_for_accession_number.get_genome_from_protein_accession(record.name))))
    return organisms


def main():
    import csv
    data = read_fasta_file("")
    data2 = read_fasta_file('')
    file1 = open('', 'a')
    author = csv.writer(file1)
    author.writerow(['Name', 'Genome GC Ratio', 'Gene GC Ratio'])
    file1.close()
    file1 = open('', 'a')
    author = csv.writer(file1)
    author.writerow(['Name', 'Genome GC Ratio', 'Gene GC Ratio'])
    file1.close()
    lp30 = data[0].genome.gc_ratio()
    for entry in data:
        with open('', 'a') as csvfile:
            author = csv.writer(csvfile)
            author.writerow([entry.name, entry.genome.gc_ratio(), entry.gene.gc_ratio(), ])
    for entry in data2:
        with open('', 'a') as csvfile:
            author = csv.writer(csvfile)
            author.writerow([entry.name, entry.genome.gc_ratio(), entry.gene.gc_ratio()])

if __name__ == '__main__':
    main()
