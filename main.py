import csv
import requests
import os

from datetime import datetime
from Bio import SeqIO
from io import StringIO


path = "/content/drive/MyDrive/AlphaFold_project/" # change to path for save location
protID = 'Q8IUZ5' #replace with the ID to get sequence from


#list out all IDs present in csv file. use these IDs to retrieve sequence later

with open(path + 'uniprot_sequences.csv', 'r+') as csvfile:
  csvreader = reader(csvfile)
  for row in csvreader:
    print(row[0])


#get sequence of protein from csv file using protein ID

seqID = 'sp|Q8IUZ5|AT2L2_HUMAN' #ID of protein seq to search in csv


#ALTERNATIVE: use this block if you do not use a stored ID + sequence from the CSV file
id = "NC_NTDLKR"
sequence = "KTEEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEITPAAAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPLIAADGGYAFKYAAGKYDIKDVGVDNAGAKAGLTFLVDLIKNKHMNADTDYSIAEAAFNKGETAMTINGPWAWSNIDTSAVNYGVTVLPTFKGQPSKPFVGVLSAGINAASPNKELAKEFLENYLLTDEGLEAVNKDKPLGAVALKSYEEELVKDPRVAATMENAQKGEIMPNIPQMSAFWYAVRTAVINAASGRQTVDAALAAAQTNAAA"
