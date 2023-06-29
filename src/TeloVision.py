#!/usr/bin/env python

"""
TeloVision is a Python package which determines the presence of telomeres and visualises scaffolds in genome assemblies.
"""

"""Import Statements"""
import pandas as pd
from Bio import SeqIO
import plotly.express as px

"""Authorship Information"""
__author__ = "Tim Verschuren"
__credits__ = ["Tim Verschuren", "Jérôme Collemare"]

__licence__ = "MIT"
__date__ = "29-06-2023"
__version__ = "1.0"
__maintainer__ = "Tim Verschuren"
__email__ = "t.verschuren@wi.knaw.nl"
__status__ = "Development"


class findTelomeres:
    """Search for repeating sequences and determine presence of telomeres.

    Attributes:
        fasta_file (str): Path to fasta file.
        output (str): Name of output file.
    """
    def __init__(self, fasta_file: str, output: str):
        self.fasta = read_fasta(fasta_file)
        self.output = output

    def identify_repeats(self, seq_slice: str) -> int:
        """Iterates over nucleotide sequence and uses
        k-mers of various sizes to determine the presence
        of a repeating sequence.

        Attributes:
            seq_slice (str): A sub portion of a nucleotide sequence.

        Returns:
            repeat_length (int): The length of the repeating
            sequence, if present.
        """
        kmer_dict = {}
        for k in range(6,10):
            for i in range(0, len(seq_slice) - k + 1):
                if seq_slice[i:i+k] == seq_slice[i+k:i+2*k] or \
                    seq_slice[i:i+k] == seq_slice[i+k+1:i+2*k+1]:
                    j = i
                    repeat_count = 1
                    while True:
                        if seq_slice[j:j+k] == seq_slice[j+k:j+2*k]:
                            repeat_count += 1
                            j += k
                        else:
                            if seq_slice[j:j+k] == seq_slice[j+k+1:j+2*k+1]:
                                repeat_count += 1
                                j += k+1
                            else:
                                if seq_slice[i:i+k] in kmer_dict:
                                    if len(kmer_dict[seq_slice[i:i+k]]) - len(seq_slice[i:j+k]) < 0:
                                        kmer_dict[str(seq_slice[i:i+k])] = str(seq_slice[i:j+k])
                                else:
                                    kmer_dict[str(seq_slice[i:i+k])] = str(seq_slice[i:j+k])
                                break
        
        seq_len = [len(seq) for seq in list(kmer_dict.values())]
        if len(seq_len) == 0:
            repeat_sequence = "NA"
            repeat = "NA"
        else:
            position = seq_len.index(max(seq_len))
            repeat_sequence = list(kmer_dict.values())[position]
            repeat = list(kmer_dict.keys())[position]
        return repeat_sequence, repeat
    
    def telomere_position(self, rep_len=30, seq_size=200) -> pd.DataFrame:
        """Loop over the scaffolds of a fasta file and determine
        whether a telomere is present at the beginning and end of
        each scaffold.

        Attributes:
            rep_len (int): Minimum length of a repetitive sequence
            to qualify as a telomeric repeat.
            seq_size (int): Size of sequence taken from the top
            and bottom of the scaffolds. 

        Returns:
            df (pd.Dataframe): Dataframe containing the name, length
            and presence of telomeres for each scaffold.
        
        Yields:
            A tsv file containing information about the repeats.
        """
        telomeres = {}
        telo_pos = []
        scaffolds = []
        lengths = []

        scaf_data = []
        len_data = []
        gc_data = []
        repeat = []
        repeat_sequence = []
        telo_class = []

        for key, value in self.fasta.items():
            telo_bin = []
            scaffolds.append(key)
            lengths.append(len(value))
            five_prime = self.identify_repeats(value[:seq_size])
            three_prime = self.identify_repeats(value[-seq_size:])

            scaf_data.append(f"{key}_5'")
            scaf_data.append(f"{key}_3'")
            repeat.append(five_prime[1])
            repeat.append(three_prime[1])
            len_data.append(len(five_prime[0]))
            len_data.append(len(three_prime[0]))
            gc_data.append(self.calculate_gc(five_prime[0]))
            gc_data.append(self.calculate_gc(three_prime[0]))
            repeat_sequence.append(five_prime[0])
            repeat_sequence.append(three_prime[0])

            if len(five_prime[0]) >= rep_len and self.calculate_gc(five_prime[0]) > 0.2:
                telo_bin.append(1)
                telomeres[key] = telo_bin
                telo_class.append("Y")
            if len(five_prime[0]) < rep_len or self.calculate_gc(five_prime[0]) < 0.2:
                telo_bin.append(0)
                telomeres[key] = telo_bin
                telo_class.append("N")
            if len(three_prime[0]) >= rep_len and self.calculate_gc(three_prime[0]) > 0.2:
                telo_bin.append(1)
                telomeres[key] = telo_bin
                telo_class.append("Y")
            if len(three_prime[0]) < rep_len or self.calculate_gc(three_prime[0]) < 0.2:
                telo_bin.append(0)
                telomeres[key] = telo_bin
                telo_class.append("N")

        for value in telomeres.values():
            telo_pos.append(value)

        telo_data = pd.DataFrame(data={"Repeat": repeat, 
                                "Length": len_data, 
                                "GC%": gc_data, 
                                "Repetitive Sequence": repeat_sequence,
                                "Telomere": telo_class}, 
                                index=scaf_data)
        telo_data.to_csv(f"{self.output}_info.tsv", sep="\t")

        df = pd.DataFrame(data={'Scaffolds': scaffolds, 'Lengths': lengths, 'Telomeres': telo_pos})
        return df

    def calculate_gc(self, sequence: str) -> int:
        """Calculates the GC content of a given sequence.

        Attributes:
        sequence (str): Nucleotide sequence.

        Returns:
        Integer between 0 and 1.
        """
        return (sequence.count("G") + sequence.count("C"))/len(sequence)


class visualiseGC:
    """Calculation and visualisation of GC content and telomere position.

    Attributes:
        fasta_file (str): Path to fasta file.
        telo_df (pd.Dataframe): Dataframe containing the name, length
        and presence of telomeres for each scaffold.
        output (str): Name of output file.
    """
    def __init__(self, fasta_file: str, telo_df: pd.DataFrame, output: str):
        self.fasta = read_fasta(fasta_file)
        self.scaffolds = []
        self.GC_cont = []
        self.length_list = []
        self.telo_df = telo_df
        self.output = output

    def gc_content(self, k=5000) -> px:
        """Calculates the GC content of each scaffold using 
        a k-mer sliding window. GC content and present telomeres
        are subsequently visualised using a plotly bar plot. 
        Telomeres will be visualised as black bars at either 
        the top or bottom of a scaffold.

        Attributes:
            k (int): Size of k-mer sliding window.

        Yields:
            px.bar: html file of plotly bar plot.
        """
        for key, value in self.fasta.items():
            if list(self.telo_df.loc[self.telo_df["Scaffolds"] == key]["Telomeres"])[0][0] == 1:
                self.GC_cont.append(0)
                self.scaffolds.append(key)
                self.length_list.append(self.telo_df["Lengths"].max()*1e-6/75)
            else:
                self.GC_cont.append(100)
                self.scaffolds.append(key)
                self.length_list.append(self.telo_df["Lengths"].max()*1e-6/75)

            for i in range(0, len(value)-k+1, k):
                self.GC_cont.append((value[i:i+k].count('C') + \
                                     value[i:i+k].count('G'))/k*100)
                self.scaffolds.append(key)
                self.length_list.append(k/1e6)
            
            if list(self.telo_df.loc[self.telo_df["Scaffolds"] == key]["Telomeres"])[0][1] == 1:
                self.GC_cont.append(0)
                self.scaffolds.append(key)
                self.length_list.append(self.telo_df["Lengths"].max()*1e-6/75)
            else:
                self.GC_cont.append(100)
                self.scaffolds.append(key)
                self.length_list.append(self.telo_df["Lengths"].max()*1e-6/75)

        fig_df = pd.DataFrame(data={"Scaffolds": self.scaffolds, 
                                    "Length": self.length_list, 
                                    "GC%": self.GC_cont})
        fig = px.bar(fig_df, 
                     x="Scaffolds", 
                     y = "Length", 
                     color="GC%", 
                     color_continuous_scale=px.colors.sequential.Hot, 
                     range_color=(0,fig_df["GC%"].mean()))
        
        fig.update_traces(marker_line_width = 0,
                  selector=dict(type="bar"))
        fig.update_xaxes(showgrid=False, title="Scaffolds", tickangle=45)
        fig.update_yaxes(showgrid=False, title="Physical position (Megabases)")
        fig.write_html(f"{self.output}.html")


def read_fasta(fasta_file) -> dict:
    """Read content of fasta file and store data
    in a dictonary with the scaffold names as the 
    key and the nucleotide sequence as the values.

    Attributes:
        fasta_file (str): Path to the fasta file.

    Returns:
        fasta_dict (dict): Dictionary containing 
        scaffold names as keys and nucleotide sequences
        as the values.
    """
    fasta_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dict[record.id] = record.seq

    return fasta_dict
