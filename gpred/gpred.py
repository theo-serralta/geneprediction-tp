import argparse
import sys
import os
import csv
import re
import textwrap
from re import Pattern
from pathlib import Path
from typing import List, Union, Optional


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True, 
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int, 
                        default=50, help="Minimum gene length to consider (default 50).")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int, 
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif (default 16).")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes - shine box not included (default 40).")
    parser.add_argument('-p', dest='predicted_genes_file', type=Path, 
                        default=Path("predict_genes.csv"),
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=Path,
                        default=Path("genes.fna"),
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file: Path) -> str:
    """Extract genome sequence from fasta files.

    :param fasta_file: (Path) Path to the fasta file.
    :return: (str) Sequence from the genome. 
    """
    sequence = []
    with open(fasta_file, 'r') as file:
        next(file)
        for line in file:
            sequence.append(line.strip().upper())  # Enlever les sauts de ligne et mettre en majuscules
    return ''.join(sequence)


def find_start(start_regex: Pattern, sequence: str, start: int, stop: int) -> Union[int, None]:
    """Find next start codon before a end position.

    :param start_regexp: A regex object that identifies a start codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :param stop: (int) Stop position of the research
    :return: (int) If exist, position of the start codon. Otherwise None. 
    """
    match = start_regex.search(sequence, start, stop)
    if match:
        return match.start(0)  # Retourne la position de début du match
    return None

def find_stop(stop_regex: Pattern, sequence: str, start: int) -> Union[int, None]:
    """Find next stop codon that should be in the same reading phase as the start.

    :param stop_regexp: A regex object that identifies a stop codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :return: (int) If exist, position of the stop codon. Otherwise None. 
    """
    for match in stop_regex.finditer(sequence, start):
        stop_pos = match.start(0)
        if (stop_pos - start) % 3 == 0:  # Vérifie le cadre de lecture
            return stop_pos
    return None

def has_shine_dalgarno(shine_regex: Pattern, sequence: str, start: int, max_shine_dalgarno_distance: int) -> bool:
    """Find a shine dalgarno motif before the start codon

    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Position of the start in the genome
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno to the start position
    :return: (boolean) true -> has a shine dalgarno upstream to the gene, false -> no
    """
    # Calculer la plage de recherche
    search_start = start - max_shine_dalgarno_distance
    search_end = start - 6
    
    # Vérifier si la position de début de recherche est valide
    if search_start < 0:
        return False
    
    # Rechercher le motif Shine-Dalgarno dans l'intervalle
    if shine_regex.search(sequence, search_start, search_end):
        return True
    return False

def predict_genes(sequence: str, start_regex: Pattern, stop_regex: Pattern, shine_regex: Pattern, 
                  min_gene_len: int, max_shine_dalgarno_distance: int, min_gap: int) -> List:
    """Predict most probable genes

    :param sequence: (str) Sequence from the genome.
    :param start_regexp: A regex object that identifies a start codon.
    :param stop_regexp: A regex object that identifies a stop codon.
    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param min_gene_len: (int) Minimum gene length.
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno to the start position.
    :param min_gap: (int) Minimum distance between two genes.
    :return: (list) List of [start, stop] position of each predicted genes.
    """
    predicted_genes = []
    position_courante = 0
    sequence_length = len(sequence)

    while sequence_length - position_courante >= min_gap:
        # Cherche le prochain codon d'initiation
        start_codon_pos = find_start(start_regex, sequence, position_courante, sequence_length)
        if start_codon_pos is not None:
            # Trouver le codon stop dans le même cadre de lecture
            stop_codon_pos = find_stop(stop_regex, sequence, start_codon_pos)
            if stop_codon_pos is not None:
                # Vérifier si le gène a la longueur minimale requise
                gene_length = stop_codon_pos - start_codon_pos + 3
                if gene_length >= min_gene_len:
                    # Vérifier la présence de la séquence de Shine-Dalgarno
                    if has_shine_dalgarno(shine_regex, sequence, start_codon_pos, max_shine_dalgarno_distance):
                        # Gène probable identifié, ajouter à la liste
                        predicted_genes.append([start_codon_pos + 1, stop_codon_pos + 3])  # Positions 1-based
                        # Avancer la position courante après ce gène
                        position_courante = stop_codon_pos + 3 + min_gap
                    else:
                        # Si pas de Shine-Dalgarno, avancer d'une position
                        position_courante += 1
                else:
                    # Si le gène est trop court, avancer d'une position
                    position_courante += 1
            else:
                # Si pas de codon stop trouvé, avancer d'une position
                position_courante += 1
        else:
            # Si pas de codon d'initiation trouvé, arrêter la recherche
            break

    return predicted_genes

def write_genes_pos(predicted_genes_file: Path, probable_genes: List[List[int]]) -> None:
    """Write list of gene positions.

    :param predicted_genes_file: (Path) Output file of gene positions.
    :param probable_genes: List of [start, stop] position of each predicted genes.
    """
    try:
        with predicted_genes_file.open("wt") as predict_genes:
            predict_genes_writer = csv.writer(predict_genes, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))


def write_genes(fasta_file: Path, sequence: str, probable_genes: List[List[int]], sequence_rc: str, 
                probable_genes_comp: List[List[int]]):
    """Write gene sequence in fasta format

    :param fasta_file: (Path) Output fasta file.
    :param sequence: (str) Sequence of genome file in 5'->3'.
    :param probable_genes: (list) List of [start, stop] position of each predicted genes in 5'->3'.
    :param sequence_rc: (str) Sequence of genome file in 3' -> 5'.
    :param probable_genes_comp: (list)List of [start, stop] position of each predicted genes in 3' -> 5'.
    """
    try:
        with open(fasta_file, "wt") as fasta:
            for i,gene_pos in enumerate(probable_genes):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                    i+1, os.linesep, 
                    fill(sequence[gene_pos[0]-1:gene_pos[1]])))
            i = i+1
            for j,gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            i+1+j, os.linesep,
                            fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])))
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(sequence: str) -> str:
    """Get the reverse complement

    :param sequence: (str) DNA Sequence.
    :return: (str) Reverse complemented sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in sequence[::-1]])


#==============================================================
# Main program
#==============================================================
def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    # Définitions des motifs regex pour les codons d'initiation, stop et la séquence de Shine-Dalgarno
    start_regex = re.compile(r'AT[TG]|[ATCG]TG')
    stop_regex = re.compile(r'TA[GA]|TGA')
    shine_regex = re.compile(r'A?G?GAGG|GGAG|GG.{1}GG')

    # Récupérer les arguments
    args = get_arguments()

    # Lecture de la séquence du génome
    sequence = read_fasta(args.genome_file)

    # Prédiction des gènes dans le sens 5' -> 3'
    probable_genes = predict_genes(
        sequence=sequence,
        start_regex=start_regex,
        stop_regex=stop_regex,
        shine_regex=shine_regex,
        min_gene_len=args.min_gene_len,
        max_shine_dalgarno_distance=args.max_shine_dalgarno_distance,
        min_gap=args.min_gap
    )

    # Calcul du complément inverse pour le sens 3' -> 5'
    sequence_rc = reverse_complement(sequence)

    # Prédiction des gènes dans le sens 3' -> 5' (complément inverse)
    probable_genes_rc = predict_genes(
        sequence=sequence_rc,
        start_regex=start_regex,
        stop_regex=stop_regex,
        shine_regex=shine_regex,
        min_gene_len=args.min_gene_len,
        max_shine_dalgarno_distance=args.max_shine_dalgarno_distance,
        min_gap=args.min_gap
    )

    # Correction des positions des gènes dans le sens 3' -> 5' pour les rendre compatibles avec le sens 5' -> 3'
    probable_genes_comp = [
        [len(sequence) - end + 1, len(sequence) - start + 1] for start, end in probable_genes_rc
    ]
    probable_genes_comp.sort()  # Tri pour garantir un ordre croissant des positions

    # Écriture des résultats dans le fichier CSV des positions
    write_genes_pos(args.predicted_genes_file, probable_genes + probable_genes_comp)

    # Écriture des séquences des gènes dans le fichier FASTA
    write_genes(args.fasta_file, sequence, probable_genes, sequence_rc, probable_genes_comp)

if __name__ == '__main__':
    main()
