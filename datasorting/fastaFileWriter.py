import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#Separate file used to get FASTA sequences from csv file

#Reads column with list of amino acids
sequenceData = pd.read_csv(r'data/Flu_hemagglutinin_sequence_data.csv', usecols = ['POSITION'])

#Creates sequence of amino acids
sequence = ""
for acid in sequenceData.POSITION:
    sequence = sequence + acid

#Creates SeqRecord from amino acid sequence
HG_Flu_record = SeqRecord(
    Seq(sequence),
    id = 'Influenza Hemaglutin'
)

#Uses SeqIO to write SeqRecord to FASTA file
SeqIO.write(HG_Flu_record,'data/HG_Flu.fasta','fasta')