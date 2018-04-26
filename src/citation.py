import pandas as pd
import numpy as np
from Bio import Entrez
import csv


def remove_redundant(df):
	
	needed = ["PMID", "ArticleIdList pubmed", "ArticleTitle", "Title", "AbstractText"]
	df = df[needed]
	df.to_csv("data/retracted.csv")


def pretty(d, indent=0):
   for key, value in d.items():
      print('\t' * indent + str(key))
      if isinstance(value, dict):
         pretty(value, indent+1)
      else:
         print('\t' * (indent+1) + str(value))



def fetch_details(id_list):
	Entrez.email = 'ebolboa@gmail.com'
	handle = Entrez.efetch(
		db='pubmed',
		retmode='xml',
		id=id_list
	)

	results = Entrez.read(handle)
	handle.close()
	return results 



def create_csv(path, res):
	with open(path, "w") as output:
		writer = csv.writer(output, lineterminator='\n')
		for val in res:
			writer.writerow([val])    


'''
Extract all the abstract texts of the articles citing each retracted article.
'''
def extract_abstract(cit):

	# first column represents PMID of retracted articles
	retracted = cit.ix[:,0]

	# drop the PMID from the citation list
	cit.drop(cit.columns[0], axis=1, inplace=True)
	citations = cit.as_matrix()

	# save all IDs in a list but remove NaN values
	cleaned_list = []
	for sub in citations:
		new = [x for x in sub if str(x) != 'nan']
		cleaned_list.append(new)

	citations = cleaned_list
	
	# save all retracted IDs
	retracted_ids = []

	# save all abstract texts
	total_abstract = []
	
	count = 0
	for i in range(len(citations)):
		
		# set of cited articles per retracted article
		abstracted_ids = []
		
		# fetch all cited articles for the current retracted article being referenced
		for j in range(len(citations)):
			try:
				text = fetch_details(citations[i][j])["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
				text = [str(text[0])]
				abstracted_ids += text
			except:
				pass

		# the largest number of articles citing the retracted article is only 4		
		if len(abstracted_ids) == 3 and citations[i][0] is not None and citations[i][0] != 'None':	
			
			# keep the same order for cited and retracted so that they correspond to one another
			total_abstract.append(abstracted_ids)
			retracted_ids.append(retracted[i])
		
		# keep track of where we are
		count += 1
		print(count)

	# lengths should be equal
	print(len(retracted_ids))
	print(len(total_abstract))

	# save the data which will be merged later
	create_csv("../data/retracted_ids.csv", retracted_ids)

	# create column header
	df_abstract = pd.DataFrame(columns=['AbstractText'])
	df_abstract['AbstractText'] = total_abstract
		
	# convert the embeded list into multiple columns
	df_abstract = pd.DataFrame(df_abstract['AbstractText'].values.tolist())
	df_abstract.to_csv("../data/abstract_split.csv")


'''
Merge PMIDs with corresponding article abstracts from articles that cite the retracted articles.
'''
def merged_csv(df1, df2):
	merged = pd.merge(df1, df2, left_index=True, right_index=True)
	merged.to_csv("../data/abstract_ids.csv")


'''
Merge the origincal CSV with the attributes corresponding to the retracted articles.
'''
def merge_retract_cited(retracted, abstract):
	df_merge = pd.merge(retracted, abstract, how='inner', on=['PMID'])
	df_merge.to_csv("../data/retracted_citations.csv")



def main():
		
	cit = pd.read_csv("../data/citations.csv")
	extract_abstract(cit)

	df1 = pd.read_csv("../data/retracted_ids.csv", names=['PMID'])
	df2 = pd.read_csv("data/abstract_split.csv", names=['AbstractCited1', 'AbstractCited2', 'AbstractCited3'])
	merged_csv(df1, df2)

	retracted = pd.read_csv("../data/retracted.csv")
	abstract = pd.read_csv("../data/abstract_ids.csv")

	merge_retract_cited(retracted, abstract)
	

main()