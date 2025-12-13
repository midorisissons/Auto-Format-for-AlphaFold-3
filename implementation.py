import csv
import requests
import os
import json

from datetime import datetime
from Bio import SeqIO
from io import StringIO

def manualInputSequence(filepath):
  Seq = input("Enter protein seq: ")
  id = input("Enter id for storing the seq: ")

  open(filepath + '/store_sequences.csv', 'a', newline="")

  with open(filepath + '/store_sequences.csv', 'r+', newline="") as csvfile:
    csvreader = csv.reader(csvfile)
    flag = False

    for row in csvreader:
      if id == row[0]:
        flag = True
        print(f"id {id} already exists in csv")
        break

    if flag == False:
      csvwriter = csv.writer(csvfile)
      csvwriter.writerow([id, Seq])
      csvfile.close()

  return id, Seq


def downloadSequence(filepath, protID):
  # make request
  baseUrl = "http://www.uniprot.org/uniprot/"
  currentUrl = baseUrl+protID+".fasta"
  response = requests.post(currentUrl)

  # get sequence
  cData = ''.join(response.text)
  Seq = StringIO(cData)
  pSeq = list(SeqIO.parse(Seq,'fasta'))
  sequence = pSeq[0].seq

  # store protein ID and sequence to csv file. If protein ID already exists, do not duplicate
  entry = [pSeq[0].id, pSeq[0].seq]
  flag = False

  altID = input(f'Current ID is {entry[0]}, enter new ID name (Enter to skip): ')
  if altID != "":
    entry[0] = altID

  open(filepath + '/store_sequences.csv', 'a', newline="")

  with open(filepath + '/store_sequences.csv', 'r+', newline="") as csvfile:
    csvreader = csv.reader(csvfile)
    flag = False
    for row in csvreader:
      if entry[0] == row[0]:
        flag = True
        print(f"id {entry[0]} already exists in csv")
        break

    if flag == False:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(entry)
        csvfile.close()

  return entry[0], entry[1]


def getSequence(filepath, seqID):
  sequence = ""
  with open(filepath + '/store_sequences.csv', 'r+', newline="") as csvfile:
    csvreader = csv.reader(csvfile)

    for row in csvreader:
      if row[0] == seqID: #search protein ID in csv file
        sequence = row[1] #store protein sequence in the variable 'sequence'
        return sequence
  if sequence =="":
    print("ID not found; no sequence retrieved")
    return sequence


"""
JSON information
https://gitlab.rc.uab.edu/rc-data-science/community-containers/Alphafold3/-/blob/6166cdc03cb35acb23815aab7aeb749c213c6d8a/docs/input.md
"""

def newFullLengthJob(filepath):
  # Make job name
  jobName = input("Enter job name (without date): ")
  d = datetime.today()
  date = f"{d.year}-{d.month}-{d.day}"
  jobName = f"{date}_{jobName}"
  print(f"Job Name: {jobName}")

  # Add protein sequences
  sequences = []
  print("Add protein sequences now...")
  while True:
    seqID = input("Enter seqID (press enter if finished with proteins): ")
    if seqID == "":
      # Move on to ligands
      break
    
    # Retrieve sequence
    sequence = getSequence(filepath, seqID)

    count = int(input("Enter sequence count: "))

    # Add to sequences
    sequences.append({"proteinChain": {"sequence": sequence, "count": count}})

  # Add ligands
  ligands = []
  print("Add ligands now...")
  while True:
    ligID = input("Enter ligand ID (press enter if finished with ligands): ")
    if ligID == "":
      # Move on to ions
      break

    count = int(input("Enter ligand count: "))

    # Add to ligands
    ligands.append({"ligand": {"ligand": ligID, "count": count}})

  # Add ions
  ions = []
  print("Add ions now...")
  while True:
    ionID = input("Enter ion ID (press enter if finished with ions): ")
    if ionID == "":
      # Move on to finalising job
      break

    count = int(input("Enter ion count: "))

    # Add to ligands
    ions.append({"ion": {"ion": ligID, "count": count}})

  jobTemplate = """[{
  "name": " """+jobName+""" ",
  "modelSeeds": [],
  "sequences": """+json.dumps(sequences+ligands+ions)+""",
  "dialect": "alphafoldserver",
  "version": 1
  }]"""

  os.makedirs(filepath+'/jobs/' + jobName, exist_ok=True)
  with open(filepath + '/jobs/' + jobName + '/' + jobName + ".json", 'w') as file:
    file.write(jobTemplate)
    file.close

  print(f'JSON file stored at: {filepath}/jobs/{jobName}/{jobName}.json')
  return jobTemplate


def fragment(sequence, filepath):
  # User-entered fragment length may be adjusted to make lengths more even and fit the protein sequence
  id = input("Enter name for existing or new directory to store fragments: ")
  os.makedirs(f'{filepath}/fragments/{id}', exist_ok=True)
  print(f"Your fragments will be stored at: {filepath}/fragments/{id}" )

  fragSize = int(input("Enter desired length of protein segment (int): "))

  overlapSize = int(input("Enter desired length of overlap between segments (0 for no overlap): "))

  sectionStatus = input("Do you wish to use a specific section (Y/N)?: ")
  if sectionStatus == "Y":
    start_residue = int(input("Enter the number of start residue: "))
    start_residue -= 1
    end_residue = int(input("Enter the number of end residue: "))
    end_residue -= 1
  else:
    if sectionStatus != "N":
      print("Incorrect entry, assume full length protein to be used.")
    start_residue = 0
    end_residue = len(sequence)

  sequence_section = sequence[start_residue:end_residue]

  if overlapSize == 0:
    num_blocks = len(sequence_section)// fragSize
    remainder = len(sequence_section)%fragSize

    fragments = [] # store all segments of equal size
    i=0

    while i < len(sequence_section):
      remainder_per_block = 0
      if remainder != 0:
        remainder_per_block = 1
      remainder -= remainder_per_block

      fragments.append(sequence_section[i:i+fragSize+remainder_per_block])
      i=i+fragSize+remainder_per_block

      num_blocks -= 1
  else:
    # make fragments with overlaps

    num_blocks = (len(sequence_section)-overlapSize) // (fragSize - overlapSize)  # calculate block size with overlap length
    remainder = len(sequence_section)+(overlapSize-fragSize)*num_blocks-overlapSize

    fragments = [] # store all segments of equal size
    i=0

    while i < len(sequence_section)-overlapSize:
      remainder_per_block = 0
      if remainder != 0:
        remainder_per_block = 1
      remainder -= remainder_per_block

      fragments.append(sequence_section[i:i+fragSize+remainder_per_block])
      i=i+fragSize+remainder_per_block-overlapSize

  filename = f"{id}_size{fragSize}_overlap{overlapSize}_residues{start_residue+1}-{end_residue}"
  with open(f"{filepath}/fragments/{id}/{filename}.csv", 'w+', newline="") as csvfile:
    csvwriter = csv.writer(csvfile)

    for i in range(len(fragments)):
      csvwriter.writerow( [fragments[i]] )

    csvfile.close()
    print(f'Fragments CSV file stored at: {filepath}/fragments/{id}/{filename}.csv')



def newFragmentJob(filepath):
  # Make job name
  jobName = input("Enter job name (without date): ")
  d = datetime.today()
  date = f"{d.year}-{d.month}-{d.day}"
  jobName = f"{date}_{jobName}"
  print(f"Job Name: {jobName}")

  # Select fragments to use
  fragments_filepath = input("Enter filepath for fragments CSV (/fragments/ ...): ")
  fragments_filepath = filepath + '/fragments/' + fragments_filepath
  
  # Add protein sequences
  sequences = []
  print("Add protein sequences now...")
  while True:
    seqID = input("Enter seqID (press enter if finished with proteins): ")
    if seqID == "":
      # Move on to ligands
      break
    
    # Retrieve sequence
    sequence = getSequence(filepath, seqID)
    count = int(input("Enter sequence count: "))

    # Add to sequences
    sequences.append({"proteinChain": {"sequence": sequence, "count": count}})

  # Add ligands
  ligands = []
  print("Add ligands now...")
  while True:
    ligID = input("Enter ligand ID (press enter if finished with ligands): ")
    if ligID == "":
      # Move on to ions
      break

    count = int(input("Enter ligand count: "))

    # Add to ligands
    ligands.append({"ligand": {"ligand": ligID, "count": count}})

  # Add ions
  ions = []
  print("Add ions now...")
  while True:
    ionID = input("Enter ion ID (press enter if finished with ions): ")
    if ionID == "":
      # Move on to finalising job
      break

    count = int(input("Enter ion count: "))

    # Add to ligands
    ions.append({"ion": {"ion": ligID, "count": count}})

  i = 1
  file_number = 1

  os.makedirs(filepath+'/jobs/' + jobName, exist_ok=True)
  
  tempjobName = ""
  with open(fragments_filepath, 'r+', newline="") as csvfile:  # open fragments file and read per line
    csvreader = csv.reader(csvfile)
    all_tasks = ""
    for row in csvreader: #each fragment
      tempjobName = f"{jobName}_seg{i}"
      insert_seq = row[0]
      insert_frag = [{"proteinChain": {"sequence": insert_seq, "count": 1}}]

      jobTemplate = """{
      "name": " """+tempjobName+""" ",
      "modelSeeds": [],
      "sequences": """+json.dumps(sequences+insert_frag+ligands+ions)+""",
      "dialect": "alphafoldserver",
      "version": 1
      }"""

      all_tasks += jobTemplate + ",\n"

      if i%100 == 0:  # if i=100,200, etc. then the JSON file is full, begin new file
        with open(filepath + '/jobs/' + jobName + '/' + jobName + '_' + str(file_number) + ".json", 'w') as file: 
          all_tasks = "[" + all_tasks[:-2] + "]"
          file.write(all_tasks)
          file.close
          print(f'JSON file stored at: {filepath}/jobs/{jobName}/{jobName}_{str(file_number)}.json')
        file_number += 1
        all_tasks = "[\n" #reset tasks for next doc
    
      i += 1

  if i%100 != 0:
    with open(filepath + '/jobs/' + jobName + '/' + jobName + '_' + str(file_number) + ".json", 'w') as file:
      all_tasks = all_tasks[:-2] + "]"
      file.write(all_tasks)
      file.close
      print(f'JSON file stored at: {filepath}/jobs/{jobName}/{jobName}_{str(file_number)}.json')

  return jobTemplate