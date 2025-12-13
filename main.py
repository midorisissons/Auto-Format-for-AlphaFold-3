import csv
import requests
import os
import json

from datetime import datetime
from Bio import SeqIO
from io import StringIO

import implementation as imp

userChoice = 0
path = os.getcwd() # change to path for save location


while userChoice != 4:
  userChoice = int(input("Enter a number:\n1. Get protein sequence\n2. Fragment protein sequence\n3. Make new AlphaFold task\n4. End program\nEnter: "))

  if userChoice < 3:
    getSeqMethod = int(input("Get protein sequence. Enter a number:\n1. Manual input\n2. Download from Uniprot\n3. Retrieve from CSV store\nEnter: "))
    if getSeqMethod == 1:
      id, sequence = imp.manualInputSequence(path)

    elif getSeqMethod == 2:
      protID = input("Enter protein ID to download from Uniprot: ")
      id, sequence = imp.downloadSequence(path, protID)
    
    elif getSeqMethod == 3:
      print("IDs contained within your CSV file:")
      with open(path + '/store_sequences.csv', 'r+', newline="") as csvfile:
          csvreader = csv.reader(csvfile)
          for row in csvreader:
            print(row[0])
      sequence = ""
      while sequence == "":
        id = input("Enter ID of protein seq to retrieve: ")
        sequence = imp.getSequence(path, id)
        if sequence != "":
          break
        print("Sequence not found in store, try again.")
    
    print(f'Sequence information: ID {id}, {sequence}')

    if userChoice == 2:
      print("Make fragments from selected sequence")
      imp.fragment(sequence, path)
  
  if userChoice == 3:
    makeJobMethod = int(input("Enter a number:\n1. Full length protein sequence\n2. Fragmented sequence\nEnter: "))
    if makeJobMethod == 1:
      imp.newFullLengthJob(path)
    elif makeJobMethod == 2:
      imp.newFragmentJob(path)

print("End program, have a nice day!")

