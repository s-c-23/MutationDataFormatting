from math import nan
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re


#Parses through fasta file to get fasta sequence
for record in SeqIO.parse("P0DP23.fasta","fasta"):
    current_seq = record.seq

#Gets column titles for mutations using a spreadsheet for a simple protein
proteinExample = pd.read_csv(r'1L2Y.csv') 
proteinExample_titles = list(proteinExample.columns.values)
proteinExample_titles.remove('Unnamed: 0')

#Creates a new dataframe with the column titles
fitnessData = pd.DataFrame(columns = proteinExample_titles)

#Dict to convert from single letter protein abbreviations to three letter protein abbreviations
proteinAbbrev = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','B':'ASX','C':'CYS','E':'GLU','Q':'GLN','Z':'GLX','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}

#Goes through fasta sequence, converts it to 3 digit codes and fills in the position and wtAA column
wtaaList = list()
positionList = list()
for index, letter in enumerate(current_seq):
    acidName = proteinAbbrev[letter]
    wtaaList.append(acidName)
    positionList.append((index+1))

fitnessData.wtAA = wtaaList
fitnessData.position = positionList

#Reads the excel sheet with the protein mutational data, specifically the mutant and screenscore columns
proteinData = pd.read_excel(r'Protein Mutation Data.xlsx', sheet_name = 'CALM1_HUMAN_Roth2017', usecols = ['mutant','screenscore']) 

#converts mutant and screenscore columns to list
mutationList = list(proteinData.mutant)
screenScoreList = list(proteinData.screenscore)

#Creates lists for the probability of each amino acid, and a dict of lists correlated to the amino acid abbreviation
prALAList=prARGList=prASNList=prASPList=prASXList=prCYSList=prGLUList=prGLNList=prGLXList=prGLYList=prHISList=prILEList=prLEUList=prLYSList=prMETList=prPHEList=prPROList=prSERList=prTHRList=prTRPList=prTYRList=prVALList = list(range(149))
prLists = {'ALA':prALAList, 'ARG':prARGList, 'ASN':prASNList,'ASP':prASPList,'ASX':prASXList,'CYS':prCYSList,'GLU':prGLUList,'GLN':prGLNList,'GLX':prGLXList,'GLY':prGLYList,'HIS':prHISList,'ILE':prILEList,'LEU':prLEUList,'LYS':prLYSList,'MET':prMETList,'PHE':prPHEList,'PRO':prPROList,'SER':prSERList,'THR':prTHRList,'TRP':prTRPList,'TYR':prTYRList,'VAL':prVALList}

#For loop which unzips mutationList and screenscoreList to fill in fitnessData dataframe
for (mutation,screenscore) in zip(mutationList,screenScoreList):
    num = ""
    mutantData = list()
    res = re.sub("\D","",mutation) # Gets the position of the mutation
    for character in mutation:
        ch = character
        if(ch.isalpha()):
            acidName = proteinAbbrev[ch] #Gets the amino acids(original and changed)
            mutantData.append(acidName)
    mutantData.append(int(res))
    mutatedAcid = mutantData[1]
    mutatedList = prLists.get(mutatedAcid) # Gets the list of probabilities using the amino acid it was mutated
    index = (mutantData[2]-1) #Gets the mutation index
    if(screenscore==screenscore):
        mutatedList.insert(index,screenscore) # Adds the fitness if it exists
    else:
        mutatedList.insert(index,nan) # Adds nan value if the fitness value doesn't exist
    colName = 'pr'+mutatedAcid
    fitnessData[colName] = pd.Series(mutatedList) #Sets the dataframe equal to the list of probabilities

print(fitnessData)

fitnessData.to_csv('./proteincheck.csv')