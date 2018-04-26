import requests
import csv
from xml.etree import ElementTree as ET

'''
	Simple script to get the citation IDs for all article IDs specified in the "articleIDs.csv" file

	Script generates the URL for each article ID and then writes the citation IDs from the response onto 
	a .csv file named "citations2"
	
	Created for the 490 project by Prajjwol, Eric and Rawan
'''

output = []
nones = 0
notNones = 0
counter = 0
with open('articleIDs.csv',newline='') as csvfile:
	spamReader = csv.reader(csvfile, delimiter=' ', quotechar='|')
	url = ''
	for row in spamReader:
		rowData = []
		url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&linkname=pubmed_pubmed_citedin&id="+row[0]
		response = requests.get(url)
		root = ET.fromstring(response.text)
		count = 0
		try:
			rowData.append(row)
			for citationID in root[0][2]:
				for a in citationID:
					row.append(a.text)
					count += 1
				if count>10:
					break
			notNones += 1
		except Exception as e:
			row.append("None")
			nones += 1
		output.append(row)
		counter += 1
		if counter%100==0:
			print("progress: ",((counter/5000)*100))

print("Nones: ",nones)
print("notNones: ",notNones)

with open('citations2.csv','w', newline='') as csvfile:
	writer = csv.writer(csvfile,delimiter='\n', quotechar='|',quoting=csv.QUOTE_MINIMAL)
	for a in output:
		writer.writerow([a])