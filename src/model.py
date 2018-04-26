import pandas as pd
import numpy as np
import collections
import csv
from nltk.corpus import stopwords
import string
from nltk import word_tokenize
from nltk.tokenize import wordpunct_tokenize as tokenize
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.feature_extraction.text import TfidfVectorizer
from scipy.spatial.distance import cosine
from nltk.stem.porter import PorterStemmer
import numpy.linalg as LA
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from scipy.spatial.distance import cosine



'''
Use PCA to reduce dimensionality to graph the TFIDF scores in 2D.
'''
def graph_words(df):

	texts = []

	for i in range(0, df.shape[0]):
		row = df['AbstractText'].iloc[i]
		if row != '[]':
			row = str(row)
			texts.append(row)

	stopWords = set(stopwords.words('english'))
	additional_sw = set(['stringelement', 'nlmcategory', 'also', 'may', 'using', 'one', 'two', '001', 'non', 'three', 'ii'])
	stopWords = stopWords.union(additional_sw)
	pipeline = Pipeline([
		('vect', CountVectorizer(stop_words=stopWords)),
		('tfidf', TfidfTransformer())
	])

	X = pipeline.fit_transform(texts).todense()
	print(X)
	pca = PCA(n_components=2).fit(X)

	data2D = pca.transform(X)

	kmeans = KMeans(n_clusters=2).fit(X)

	plt.scatter(data2D[:,0], data2D[:,1])

	plt.show()



'''
Remove rows with empty values
'''
def preprocess(df):
	df = df[df.AbstractText.str.contains("None") == False]
	df = df[df["AbstractText"] != "[]"]
	df.to_csv("../data/retracted_citations.csv")



'''
For every retracted article the clustering function is called to calculate the
similarity bewteen the retracted article and the articles that cite it.
'''
def word_count_clustering(df):

	# loop through every row
	for i in range(0, df.shape[0]):
		
		# store the texts of the articles that cite the retracted article
		texts = []

		# title of the retracted article
		title = df['ArticleTitle'].iloc[i]
		
		# title and abstract text are merged before calculating the TFIDF matrix
		target = [str(df['AbstractText'].iloc[i]) + str(title)]

		# loop through articles that cite the retracted one
		for column in df:

			# append all the abstract texts
			text = df[column].iloc[i]
			texts.append(str(text))
		
		# cluster the documents
		clustering(texts, target)




'''
Tokenize a documents and remove stop words.
'''
def word_tokenizer(sentences):
	tokens = word_tokenize(sentences)
	stemmer = PorterStemmer()
	tokens = [stemmer.stem(t) for t in tokens if t not in stopwords.words('english')]
	return tokens




'''
Calculates which cluster is closest to the retracted article.
'''
def clustering(cited, target):
	
	tfidf_vect = TfidfVectorizer(
			tokenizer=word_tokenizer,
			max_df=0.9,
			min_df=0.1,
			lowercase=True
		)


	tfidf_matrix = tfidf_vect.fit_transform(cited)
	tfidf_test = tfidf_vect.transform(target)

	kmeans = KMeans(n_clusters=3)
	kmeans.fit(tfidf_matrix)

	clusters = collections.defaultdict(list)

	for i, label in enumerate(kmeans.labels_):
		clusters[label].append(i)
	res = dict(clusters)

	for cluster in range(3):
		print("cluster ", cluster, ":")
		for i, sentence in enumerate(res[cluster]):
			print("\tsentence ", i, ": ", cited[sentence][:30] + '...')

	closest = kmeans.predict(tfidf_test)
	print("Closest cluster is " + str(closest[0]))



		

'''
For every retracted article the cosine similarity function is called to calculate the
similarity bewteen the retracted article and the articles that cite it.
'''
def word_count_cs_similar(df):

	total_flagged = 0
	total = 0

	# loop through every row
	for i in range(0, df.shape[0]):
		
		# store the texts of the articles that cite the retracted article
		texts = []

		# title of the retracted article
		title = df['ArticleTitle'].iloc[i]
		
		# title and abstract text are merged before calculating the TFIDF matrix
		target = [str(df['AbstractText'].iloc[i]) + str(title)]

		# loop through articles that cite the retracted one
		for column in df:

			# append all the abstract texts
			text = df[column].iloc[i]
			texts.append(str(text))
		
		# calculate the cosine similarity
		flagged = vectorize(texts, target)
		total_flagged += flagged[0]
		total += flagged[1]

	# calculate number of flagged articles
	print(total_flagged/total)



'''
Calcualate the cosine similarity bewteen a target document and a group of documents
'''
def vectorize(cited, target):

	# lambda definition for the cosine function, btu it is not being used
	cosine_function = lambda a, b : round(np.inner(a, b)/(LA.norm(a)*LA.norm(b)), 3)

	# stop words dictionary
	stopWords = set(stopwords.words('english'))
	
	# some words were added to the dictionary
	additional_sw = set(['stringelement', 'nlmcategory', 'also', 'may', 'using', 'one', 'two', '001', 'non', 'three', 'ii', 'label'])
	stopWords = stopWords.union(additional_sw)

	# get the frequency of the words being used	
	vectorizer = CountVectorizer(stop_words=stopWords)

	# fit the abstract texts of the articles citing the retracted article
	trainVectorizerArray = vectorizer.fit_transform(cited).toarray()

	# fit the retracted articles
	testVectorizerArray = vectorizer.transform(target).toarray()

	# number of flagged articles
	flagged = 0

	# loop through all the articles citing the retracted one
	for vector in trainVectorizerArray:

		# there is only a single target
		for testV in testVectorizerArray:

			# calculate the cosine similarity
			cs = cosine(vector, testV)
			cs = 1 - cs
			
			# multiply the percentage by 100
			cs *= 100
			print(cs)

			# flag any article with a score greater than or equal to 50
			if cs >= 50:
				flagged += 1

	return [flagged, len(trainVectorizerArray)]

	


def main():
	df = pd.read_csv('../data/retracted_citations.csv')
	preprocess(df)

	word_count_clustering(df)
	word_count_cs_similar(df)
	graph_words(df)

main()
