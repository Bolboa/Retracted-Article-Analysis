import numpy as np
from Bio import Entrez
import csv
import sys
import collections

from dicttoxml import dicttoxml
from xml.dom.minidom import parseString



'''
Returns all the IDs of a specified article type.
'''
def search(query):
	Entrez.email = 'ebolboa@gmail.com'
	handle = Entrez.esearch(
		db='pubmed',
		sort='relevance',
		retmax='500',
		term=query
	)
	results = Entrez.read(handle)
	handle.close()
	return results


'''
Takes a list of article IDs and returns a dictionary containing all the features of each article.
'''
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


'''
Recursively print an embdedded dictionary in a clean manner.
'''
def pretty(d, indent=0):
   for key, value in d.items():
      print('\t' * indent + str(key))
      if isinstance(value, dict):
         pretty(value, indent+1)
      else:
         print('\t' * (indent+1) + str(value))


'''
Write key and values to CSV file where key represents the header of a column.
'''
def create_csv(path, key, value):
	with open(path, "w") as csv_file:
		writer = csv.writer(csv_file, delimiter=',')
		writer.writerow(key)
		for item in value:
			line = []
			for k in key:
				if k in item:
					line.append(item[k])
				else:
					line.append("None")
			writer.writerow(line)

            	
'''
Recursively extract all the child nodes in an embedded  dictionary.
'''
def get_terminal(data):
	for k, item in data.items():
		if isinstance(item, dict):
			for i in get_terminal(item):
				yield i
		else:
			yield k, item



def main():

	# extract the IDs of all the retracted articles
	results = search("retracted publication")
	ids = results['IdList']
	
	# list to save all child nodes
	all_children = []
	
	current = 0
	for idx in ids:
		try:
			details = fetch_details(idx)
			main = details['PubmedArticle'][0]
			pretty(details['PubmedArticle'][0])
			terminal_nodes = dict(get_terminal(main))
			if "History" in terminal_nodes:
				for terminal in terminal_nodes["History"]:
					for k, v in terminal.items():
						terminal_nodes[str(terminal.attributes) + k] = v
				del terminal_nodes["History"]
			if "AuthorList" in terminal_nodes:
				all_authors_yes = []
				all_authors_no = []
				for terminal in terminal_nodes["AuthorList"]:
					author_info = []
					for k, v in terminal.items():
						author_info.append((k, v))
					if terminal.attributes == {'ValidYN': 'Y'}:
						all_authors_yes.append(author_info)
					else:
						all_authors_no.append(author_info)
				terminal_nodes["AuthorList " + str("Valid Yes")] = all_authors_yes
				terminal_nodes["AuthorList " + str("Valid No")] = all_authors_no
				del terminal_nodes["AuthorList"]
			if "KeywordList" in terminal_nodes:
				all_words_main = []
				all_words_not_main = []
				for terminal in terminal_nodes["KeywordList"]:
					for sub_terminal in terminal:
						if sub_terminal.attributes == {'MajorTopicYN': 'Y'}:
							all_words_main.append(str(sub_terminal))
						else:
							all_words_not_main.append(str(sub_terminal))
				terminal_nodes["KeywordList " + str("Major Topic: Yes")] = all_words_main
				terminal_nodes["KeywordList " + str("Major Topic: No")] = all_words_not_main
				del terminal_nodes["KeywordList"]
			if "ArticleIdList" in terminal_nodes:
				id_pubmed = []
				id_doi = []
				id_pii = []
				id_pmc = []
				id_pmcpid = []
				id_pmpid = []
				id_mid = []
				id_sici = []
				id_medline = []
				id_pmcid = []
				for terminal in terminal_nodes["ArticleIdList"]:
					if terminal.attributes == {'IdType': 'pubmed'}:
						id_pubmed.append(str(terminal))
					elif terminal.attributes == {'IdType': 'doi'}:
						id_doi.append(str(terminal))
					elif terminal.attributes == {'IdType': 'pii'}:
						id_pii.append(str(terminal))
					elif terminal.attributes == {'IdType': 'pmcpid'}:
						id_pmcpid,append(str(terminal))
					elif terminal.attributes == {'IdType': 'pmpid'}:
						id_pmpid.append(str(terminal))
					elif terminal.attributes == {'IdType': 'mid'}:
						id_mid.append(str(terminal))
					elif terminal.attributes == {'IdType': 'sici'}:
						id_sici.append(str(terminal))
					elif terminal.attributes == {'IdType': 'medline'}:
						id_medline.append(str(terminal))
					elif terminal.attributes == {'IdType': 'pmcid'}:
						id_pmcid.append(str(terminal))
					else:
						id_pmc.append((str(terminal)))

				terminal_nodes["ArticleIdList " + str("pubmed")] = id_pubmed
				terminal_nodes["ArticleIdList " + str("doi")] = id_doi
				terminal_nodes["ArticleIdList " + str("pii")] = id_pii
				terminal_nodes["ArticleIdList " + str("pmc")] = id_pmc
				terminal_nodes["ArticleIdList " + str("pmcpid")] = id_pmcpid
				terminal_nodes["ArticleIdList " + str("pmpid")] = id_pmpid
				terminal_nodes["ArticleIdList " + str("mid")] = id_mid
				terminal_nodes["ArticleIdList " + str("medline")] = id_medline
				terminal_nodes["ArticleIdList " + str("sici")] = id_sici
				terminal_nodes["ArticleIdList " + str("pmcid")] = id_pmcid
				del terminal_nodes["ArticleIdList"]
			if "PublicationTypeList" in terminal_nodes:
				all_types = []
				for terminal in terminal_nodes["PublicationTypeList"]:
					all_types.append(str(terminal))
				terminal_nodes["PublicationTypeList"] = all_types
			if "ArticleDate" in terminal_nodes:
				for date in terminal_nodes["ArticleDate"]:
					date_format = ""
					for k, v in date.items():
						date_format += str(v) + "/"
					terminal_nodes["ArticleDate"] = date_format
			if "MeshHeadingList" in terminal_nodes:
				qualifiers_major = []
				qualifiers_non_major = []
				descriptors_major = []
				descriptors_non_major = []
				for terminal in terminal_nodes["MeshHeadingList"]:
					for k, v in terminal.items():
						if k == "QualifierName":
							for el in v:
								if hasattr(el, "attributes"):
									if str("'MajorTopicYN': 'Y'") in str(el.attributes):
										qualifiers_major.append(str(el))
									else:
										qualifiers_non_major.append(str(el))
						else:
							if hasattr(v, "attributes"):
								if str("'MajorTopicYN': 'Y'") in str(v.attributes):
									descriptors_major.append(str(v))
								else:
									descriptors_non_major.append(str(v))
				terminal_nodes["MeshHeadingList Qualifiers " + "major yes"] = qualifiers_major
				terminal_nodes["MeshHeadingList Qualifiers " + "major no"] = qualifiers_non_major
				terminal_nodes["MeshHeadingList Descriptors" + "major yes"] = descriptors_major
				terminal_nodes["MeshHeadingList Descriptors " + "major no"] = descriptors_non_major
				del terminal_nodes["MeshHeadingList"]
			if "CommentsCorrectionsList" in terminal_nodes:
				retraction_in = []
				retraction_of = []
				cites = []
				erratum_in = []
				erratum_for = []
				comment_in = []
				comment_on = []
				assoc_data = []
				assoc_pub = []
				expr_con_in = []
				expr_con_for = []
				rep_from = []
				rep_in = []
				update_in = []
				update_of = []
				sum_pt_in = []
				original_rpt = []
				reprint_in = []
				reprint_of = []
				for terminal in terminal_nodes["CommentsCorrectionsList"]:
					for k, v in terminal.items():
						if {'RefType': 'Cites'} == terminal.attributes:
							cites.append(str(v))
						elif {'RefType': 'RetractionIn'} == terminal.attributes:
							retraction_in.append(str(v))
						elif {'RefType': 'RetractionOf'} == terminal.attributes:
							retraction_of.append(str(v))
						elif {'RefType': 'ErratumIn'} == terminal.attributes:
							erratum_in.append(str(v))
						elif {'RefType': 'ErratumFor'} == terminal.attributes:
							erratum_for.append(str(v))
						elif {'RefType': 'AssociatedDataset'} == terminal.attributes:
							assoc_data.append(str(v))
						elif {'RefType': 'AssociatedPublication'} == terminal.attributes:
							assoc_pub.append(str(v))
						elif {'RefType': 'CommentOn'} == terminal.attributes:
							comment_on.append(str(v))
						elif {'RefType': 'ExpressionOfConcernIn'} == terminal.attributes:
							expr_con_in.append(str(v))
						elif {'RefType': 'ExpressionOfConcernFor'} == terminal.attributes:
							expr_con_for.append(str(v))
						elif {'RefType': 'RepublishedFrom'} == terminal.attributes:
							rep_from.append(str(v))
						elif {'RefType': 'RepublishedIn'} == terminal.attributes:
							rep_in.append(str(v))
						elif {'RefType': 'UpdateIn'} == terminal.attributes:
							update_in.append(str(v))
						elif {'RefType': 'UpdateOf'} == terminal.attributes:
							update_of.append(str(v))
						elif {'RefType': 'SummaryForPatientsIn'} == terminal.attributes:
							sum_pt_in.append(str(v))
						elif {'RefType': 'OriginalReportIn'} == terminal.attributes:
							original_rpt.append(str(v))
						elif {'RefType': 'ReprintIn'} == terminal.attributes:
							reprint_in.append(str(v))
						elif {'RefType': 'ReprintOf'} == terminal.attributes:
							reprint_of.append(str(v))
						else:
							comment_in.append(str(v))

				terminal_nodes["CommentsCorrectionsList Cites"] = cites
				terminal_nodes["CommentsCorrectionsList Retraction-In"] = retraction_in
				terminal_nodes["CommentsCorrectionsList Retraction-Of"] = retraction_of
				terminal_nodes["CommentsCorrectionsList Erratum-In"] = erratum_in
				terminal_nodes["CommentsCorrectionsList Erratum-For"] = erratum_for
				terminal_nodes["CommentsCorrectionsList Comment-In"] = comment_in
				terminal_nodes["CommentsCorrectionsList Associated-Dataset"] = assoc_data
				terminal_nodes["CommentsCorrectionsList Associated-Publication"] = assoc_pub
				terminal_nodes["CommentsCorrectionsList Comment-On"] = comment_on
				terminal_nodes["CommentsCorrectionsList Expression-Of-Concern-In"] = expr_con_in
				terminal_nodes["CommentsCorrectionsList Expression-Of-Concern-For"] = expr_con_for
				terminal_nodes["CommentsCorrectionsList Republished-From"] = rep_from
				terminal_nodes["CommentsCorrectionsList Republished-In"] = rep_in
				terminal_nodes["CommentsCorrectionsList Update-In"] = update_in
				terminal_nodes["CommentsCorrectionsList Update-Of"] = update_of
				terminal_nodes["CommentsCorrectionsList SummaryForPatientsIn"] = sum_pt_in
				terminal_nodes["CommentsCorrectionsList Original-Report-In"] = original_rpt
				terminal_nodes["CommentsCorrectionsList Reprint-In"] = reprint_in
				terminal_nodes["CommentsCorrectionsList Reprint-Of"] = reprint_of

				del terminal_nodes["CommentsCorrectionsList"]


			all_children.append(terminal_nodes)
			current += 1
			print(current)
		except:
			current += 1
			pass
	
	print(len(all_children))
	
	csv_key = []
	for child in all_children:
		tmp_key = list(child.keys())
		if tmp_key != csv_key:
			first = set(csv_key)
			second = set(tmp_key)
			in_one_not_other = second - first
			csv_key += list(in_one_not_other)

	# write the extrracted data to a CSV file
	create_csv("./data/retracted.csv", csv_key, all_children)

main()
