!pip install Bio

"""## Download sequence into csv file by uniprot ID

The protein sequence is retrieved from the uniprot server using its ID. The sequence and ID are then stored in a csv file for storage to avoid having to retrieve a sequence each time.

This section can be skipped if you want to manually enter protein sequences.
"""

def downloadSequence(path, protID):
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
  fields = ['id', 'seq']
  entry = [pSeq[0].id, pSeq[0].seq  ]
  flag = False

  with open(path + 'uniprot_sequences.csv', 'r+') as csvfile:
    csvreader = reader(csvfile)

    for row in csvreader:
      if entry[0] == row[0]:
        flag = True
        print(f"id {entry[0]} already exists in csv")
        break

    if flag == False:
        csvwriter = writer(csvfile)
        csvwriter.writerow(entry)
        csvfile.close()

def getSequence(path, seqID):
  sequence = ""
  with open(path + 'uniprot_sequences.csv', 'r+') as csvfile:
    csvreader = reader(csvfile)

    for row in csvreader:
      if row[0] == seqID: #search protein ID in csv file
        sequence = row[1] #store protein sequence in the variable 'sequence'
        break
  if sequence =="":
    print("ID not found; no sequence retrieved")

  return sequence






"""## Get Json file for full protein sequence (no fragments)

### Explanation on making a JSON file template for AlphaFold:
AlphaFold servers can accept JSON files of a certain format to automate the queuing of up to a 100 tasks at once, speeding up the process greatly.

In this code, I have used template for CK2 beta only:
*   2x CK2 beta
*   1x protein sequence of interest
*   2x Zn2+

and a template for CK2 holoenzyme:
*   2x CK2 alpha
*   2x CK2 beta
*   1x protein sequence of interest
*   2x ADP
*   2x Mg2+
*   2x Zn2+

For each file I have also automated the generation of the task name.

Find more information about making your own template:
https://gitlab.rc.uab.edu/rc-data-science/community-containers/Alphafold3/-/blob/6166cdc03cb35acb23815aab7aeb749c213c6d8a/docs/input.md
"""

def newFullLengthJob(path):
  # Make job name
  jobName = input("Enter job name (without date): ")
  date = datetime.today().strftime('%Y/%-m/%-d')
  jobName = f"{date}_{jobName}"
  print(f"Job Name: {jobName}")

  # Add protein sequences
  sequences = []
  print("Add protein sequences now...")
  while True:
    seqID = input("Enter seqID (press enter if finished with proteins)")
    if seqID == "":
      # Move on to ligands
      break
    
    # Retrieve sequence
    sequence = getSequence(path, seqID)

    count = input("Enter sequence count: ")

    # Add to sequences
    sequences.append({"proteinChain": {"sequence": sequence, "count": count}})

  # Add ligands
  ligands = []
  print("Add ligands now...")
  while True:
    ligID = input("Enter ligand ID (press enter if finished with ligands)")
    if ligID == "":
      # Move on to ions
      break

    count = input("Enter ligand count: ")

    # Add to ligands
    ligands.append({"ligand": {"ligand": ligID, "count": count}})

  # Add ions
  ions = []
  print("Add ions now...")
  while True:
    ionID = input("Enter ion ID (press enter if finished with ions)")
    if ionID == "":
      # Move on to finalising job
      break

    count = input("Enter ion count: ")

    # Add to ligands
    ions.append({"ion": {"ion": ligID, "count": count}})

  jobTemplate = """{
  "name": " """+jobName+""" ",
  "modelSeeds": [],
  "sequences": """+(sequences+ligands+ions)+""",
  "dialect": "alphafoldserver",
  "version": 1
  }"""

  return jobTemplate





date = "2025-02-09"
testprotein_name = "NC_NTDLKR"

#make new folder for protein ID within the path you provided - new folder for each protein


os.mkdir( f"{path}{id}" )

#create Json task files for protein of interest + CK2 holoenzyme, protein of interest + CK2beta

all_tasks = "[\n"

#Holoenzyme

taskname = f"{date}_CK2holo_and_full_{testprotein_name}"



all_tasks += template + ",\n"


#Beta only

taskname = f"{date}_CK2beta_and_full_{testprotein_name}"

template = """{
"name": " """+taskname+""" ",
"modelSeeds": [],
"sequences": [
  {"proteinChain": {
        "sequence": "MSSSEEVSWISWFCGLRGNEFFCEVDEDYIQDKFNLTGLNEQVPHYRQALDMILDLEPDEELEDNPNQSDLIEQAAEMLYGLIHARYILTNRGIAQMLEKYQQGDFGYCPRVYCENQPMLPIGLSDIPGEAMVKLYCPKCMDVYTPKSSRHHHTDGAYFGTGFPHMLFMVHPEYRPKRPANQFVPRLYGFKIHPMAYQLQLQAASNFKSPVKTIR",
        "count": 2
    } },
  {"ion": {"ion": "ZN", "count": 2} },
  {"proteinChain": {
        "sequence": """+insert_seq+""",
        "count": 2
    } }
],
"dialect": "alphafoldserver",
"version": 1
}"""

all_tasks += template + ",\n"


# output

all_tasks = all_tasks[:-3] + "} ]"

with open( f"{path}{id}/{id}_full_jobs.json", 'w') as file:
  file.write(all_tasks)
  file.close

"""## Make csv file of uniform length fragments

This section of code is used to fragment your chosen protein sequence into set lengths and store them into a new CSV file. The CSV file of fragments can then be accessed to create more tasks.
"""

#make CSV file to store fragments of uniform size of protein

segment_size = 20 #choose the rough size of each segment
num_blocks = len(sequence)// segment_size
remainder = len(sequence)%segment_size

avg_length = segment_size + remainder//num_blocks

frag_seq = []

i=0

while i < len(sequence):
  remainder_per_block = remainder // num_blocks
  remainder -= remainder_per_block

  frag_seq.append( sequence[i:i+segment_size+remainder_per_block])
  i=i+segment_size+remainder_per_block

  num_blocks -= 1


#store in csv

filename = f" {id}_size:{avg_length}"

with open(f"{path}{id}/{filename}", 'w+') as csvfile:
  csvwriter = writer(csvfile)

  for i in range(len(frag_seq)):
    csvwriter.writerow( [frag_seq[i]] )

  csvfile.close()

#fragments but with custom start and end point - if there is region of interest in sequence
start_res = 1
end_res = 202
sequence_section = sequence[start_res:end_res]

segment_size = 20 #choose the rough size of each segment
num_blocks = len(sequence_section)// segment_size
remainder = len(sequence_section)%segment_size

avg_length = segment_size + remainder//num_blocks

frag_seq = []

i=0

while i < len(sequence_section):
  remainder_per_block = remainder // num_blocks
  remainder -= remainder_per_block

  frag_seq.append( sequence_section[i:i+segment_size+remainder_per_block])
  i=i+segment_size+remainder_per_block

  num_blocks -= 1


#store in csv

filename = f" {id}_size:{avg_length}_{start_res}:{end_res}"

with open(f"{path}{id}/{filename}", 'w+') as csvfile:
  csvwriter = writer(csvfile)

  for i in range(len(frag_seq)):
    csvwriter.writerow( [frag_seq[i]] )

  csvfile.close()

#fragments but with overlapping sections

#note - didn't get around to coding this!

"""2025-02-04_CK2beta_and_p27KIP1_seg1

## Make json file of fragments jobs
"""

date = "2025-02-09"
filepath = "/content/drive/MyDrive/AlphaFold_project/sp|P25054|APC_HUMAN_size:20/ sp|P25054|APC_HUMAN_size:20"
testprotein_name = "APC"

"""Holoenzyme"""

all_tasks = "[\n"
i=0

with open(filepath, 'r+') as csvfile:
  csvreader = reader(csvfile)
  for row in csvreader: #each fragment
    i+=1

    taskname = f"2025-02-09_CK2holo_and_{testprotein_name}_seg{i}"
    insert_seq = "\"" + row[0] + "\""

    template = """{
    "name": " """+taskname+""" ",
    "modelSeeds": [],
    "sequences": [
      {"proteinChain": {
            "sequence": "MSGPVPSRARVYTDVNTHRPREYWDYESHVVEWGNQDDYQLVRKLGRGKYSEVFEAINITNNEKVVVKILKPVKKKKIKREIKILENLRGGPNIITLADIVKDPVSRTPALVFEHVNNTDFKQLYQTLTDYDIRFYMYEILKALDYCHSMGIMHRDVKPHNVMIDHEHRKLRLIDWGLAEFYHPGQEYNVRVASRYFKGPELLVDYQMYDYSLDMWSLGCMLASMIFRKEPFFHGHDNYDQLVRIAKVLGTEDLYDYIDKYNIELDPRFNDILGRHSRKRWERFVHSENQHLVSPEALDFLDKLLRYDHQSRLTAREAMEHPYFYTVVKDQARMGSSSMPGGSTPVSSANMMSGISSVPTPSPLGPLAGSPVIAAANPLGMPVPAAAGAQQ",
            "count": 2
        } },
      {"proteinChain": {
            "sequence": "MSSSEEVSWISWFCGLRGNEFFCEVDEDYIQDKFNLTGLNEQVPHYRQALDMILDLEPDEELEDNPNQSDLIEQAAEMLYGLIHARYILTNRGIAQMLEKYQQGDFGYCPRVYCENQPMLPIGLSDIPGEAMVKLYCPKCMDVYTPKSSRHHHTDGAYFGTGFPHMLFMVHPEYRPKRPANQFVPRLYGFKIHPMAYQLQLQAASNFKSPVKTIR",
            "count": 2
        } },
      {"ligand": {"ligand": "CCD_ADP", "count": 2} },
      {"ion": {"ion": "MG", "count": 2} },
      {"ion": {"ion": "ZN", "count": 2} },
      {"proteinChain": {
            "sequence": """+insert_seq+""",
            "count": 1
        } }
    ],
    "dialect": "alphafoldserver",
    "version": 1
    }"""

    all_tasks += template + ",\n"

  all_tasks = all_tasks[:-3] + "} ]"

with open(filepath + "_holo_jobs.json", 'w') as file:
  file.write(all_tasks)
  file.close

#max tasks=100, make multiple JSON files if more than 100 (max for server upload)

all_tasks = "[\n"
i+=1
num_docs=0 #keep track of which number JSON file you are on

with open(filepath, 'r+') as csvfile:
  csvreader = reader(csvfile)
  for row in csvreader: #each fragment
    i+=1

    taskname = f"2025-02-09_CK2holo_and_{testprotein_name}_seg{i}"
    insert_seq = "\"" + row[0] + "\""

    template = """{
    "name": " """+taskname+""" ",
    "modelSeeds": [],
    "sequences": [
      {"proteinChain": {
            "sequence": "MSGPVPSRARVYTDVNTHRPREYWDYESHVVEWGNQDDYQLVRKLGRGKYSEVFEAINITNNEKVVVKILKPVKKKKIKREIKILENLRGGPNIITLADIVKDPVSRTPALVFEHVNNTDFKQLYQTLTDYDIRFYMYEILKALDYCHSMGIMHRDVKPHNVMIDHEHRKLRLIDWGLAEFYHPGQEYNVRVASRYFKGPELLVDYQMYDYSLDMWSLGCMLASMIFRKEPFFHGHDNYDQLVRIAKVLGTEDLYDYIDKYNIELDPRFNDILGRHSRKRWERFVHSENQHLVSPEALDFLDKLLRYDHQSRLTAREAMEHPYFYTVVKDQARMGSSSMPGGSTPVSSANMMSGISSVPTPSPLGPLAGSPVIAAANPLGMPVPAAAGAQQ",
            "count": 2
        } },
      {"proteinChain": {
            "sequence": "MSSSEEVSWISWFCGLRGNEFFCEVDEDYIQDKFNLTGLNEQVPHYRQALDMILDLEPDEELEDNPNQSDLIEQAAEMLYGLIHARYILTNRGIAQMLEKYQQGDFGYCPRVYCENQPMLPIGLSDIPGEAMVKLYCPKCMDVYTPKSSRHHHTDGAYFGTGFPHMLFMVHPEYRPKRPANQFVPRLYGFKIHPMAYQLQLQAASNFKSPVKTIR",
            "count": 2
        } },
      {"ligand": {"ligand": "CCD_ADP", "count": 2} },
      {"ion": {"ion": "MG", "count": 2} },
      {"ion": {"ion": "ZN", "count": 2} },
      {"proteinChain": {
            "sequence": """+insert_seq+""",
            "count": 1
        } }
    ],
    "dialect": "alphafoldserver",
    "version": 1
    }"""

    all_tasks += template + ",\n"

  all_tasks = all_tasks[:-3] + "} ]"

  if i%100 == 0: #if i=100, 200, etc. then the JSON file is full
    with open(filepath + "_holo_jobs_" + num_docs + "".json", 'w') as file:
      file.write(all_tasks)
      file.close
    all_tasks = "[\n" #reset tasks for next doc

"""CK2beta only"""

all_tasks = "[\n"
i=0

with open(filepath, 'r+') as csvfile:
  csvreader = reader(csvfile)
  for row in csvreader: #each fragment
    i+=1

    taskname = f"{date}_CK2beta_and_{testprotein_name}_seg{i}"
    insert_seq = "\"" + row[0] + "\""

    template = """{
    "name": " """+taskname+""" ",
    "modelSeeds": [],
    "sequences": [
      {"proteinChain": {
            "sequence": "MSSSEEVSWISWFCGLRGNEFFCEVDEDYIQDKFNLTGLNEQVPHYRQALDMILDLEPDEELEDNPNQSDLIEQAAEMLYGLIHARYILTNRGIAQMLEKYQQGDFGYCPRVYCENQPMLPIGLSDIPGEAMVKLYCPKCMDVYTPKSSRHHHTDGAYFGTGFPHMLFMVHPEYRPKRPANQFVPRLYGFKIHPMAYQLQLQAASNFKSPVKTIR",
            "count": 2
        } },
      {"ion": {"ion": "ZN", "count": 2} },
      {"proteinChain": {
            "sequence": """+insert_seq+""",
            "count": 1
        } }
    ],
    "dialect": "alphafoldserver",
    "version": 1
    }"""

    all_tasks += template + ",\n"

  all_tasks = all_tasks[:-3] + "} ]"

with open(filepath + "_beta_jobs.json", 'w') as file:
  file.write(all_tasks)
  file.close