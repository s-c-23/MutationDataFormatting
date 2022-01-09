from math import nan
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re
import numpy as np

#Function to process through mutation file using variables related to the mutation data
def processMutations(fastaRecord: str, mutationDataSheet: str, fitnessCol: str, sequenceLength: int, resultSheetName: str):
    #Parses through fasta file to get fasta sequence
    for record in SeqIO.parse(fastaRecord,"fasta"):
        current_seq = record.seq

    protein_titles = ['pdb','chain_id','position','wtAA','ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','Unknown']

    #Creates a new dataframe with the column titles
    fitnessData = pd.DataFrame(columns = protein_titles)

    #Dict to convert from single letter protein abbreviations to three letter protein abbreviations
    proteinAbbrev = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL','X':'Unknown'}

    #Goes through fasta sequence, converts it to 3 digit codes and fills in the position and wtAA column
    wtaaList = list()
    positionList = list()
    for index, letter in enumerate(current_seq):
        acidName = proteinAbbrev[letter]
        wtaaList.append(acidName)
        positionList.append((index+1))

    fitnessData.wtAA = wtaaList
    fitnessData.position = positionList
    fitnessData.pdb = '#N/A'
    fitnessData.chain_id = '#N/A'

    #Reads the excel sheet with the protein mutational data, specifically the mutant and score columns
    proteinData = pd.read_excel(r'Protein Mutation Data.xlsx', sheet_name = mutationDataSheet, usecols = ['mutant',fitnessCol]) 

    #converts mutant and score columns to list
    mutationList = list(proteinData.mutant)
    screenScoreList = list(proteinData[fitnessCol])


    #Creates lists for the fitness of each amino acid, and a dict of lists correlated to the amino acid abbreviation
    ALAList, ARGList, ASNList, ASPList, ASXList, CYSList, GLUList, GLNList, GLXList, GLYList, HISList, ILEList, LEUList, LYSList, METList, PHEList, PROList, SERList, THRList, TRPList, TYRList, VALList, UnknownList = [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength, [None]*sequenceLength
    prLists = {'ALA':ALAList, 'ARG':ARGList, 'ASN':ASNList,'ASP':ASPList,'ASX':ASXList,'CYS':CYSList,'GLU':GLUList,'GLN':GLNList,'GLX':GLXList,'GLY':GLYList,'HIS':HISList,'ILE':ILEList,'LEU':LEUList,'LYS':LYSList,'MET':METList,'PHE':PHEList,'PRO':PROList,'SER':SERList,'THR':THRList,'TRP':TRPList,'TYR':TYRList,'VAL':VALList,'Unknown':UnknownList}

    #For loop which unzips mutationList and scoreList to fill in fitnessData dataframe
    for (mutation,score) in zip(mutationList,screenScoreList):
        mutantData = list()
        res = re.sub("\D","",mutation) # Gets the position of the mutation
        for character in mutation:
            ch = character
            if(ch.isalpha()):
                acidName = proteinAbbrev[ch] #Gets the amino acids(original and changed)
                mutantData.append(acidName)
        mutantData.append(int(res))
        # print(mutantData)
        mutatedAcid = mutantData[1]
        # print(mutatedAcid)
        mutatedList = prLists.get(mutatedAcid) # Gets the list of probabilities using the amino acid it was mutated
        index = (mutantData[2]-1) #Gets the mutation index
        if(score==score):
            mutatedList[index]=score # Adds the fitness if it exists
        else:
            mutatedList[index]=np.nan # Adds nan value if the fitness value doesn't exist
        # print(mutatedList)
        colName = mutatedAcid
        fitnessData[colName] = pd.Series(mutatedList) #Sets the dataframe equal to the list of probabilities
        # print(fitnessData)

    # print(ALAList)
    fitnessData.to_csv(resultSheetName,na_rep='#N/A')
