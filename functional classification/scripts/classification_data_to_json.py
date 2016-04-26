#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys
import json, csv, re

class_files = '../data/functional gene classification/'
gene_lis = []
uniq_gene_lis = []
dictionary = {}
hgnc_dictionary = {}
ensembl_dictionary = {}
classification_dictionary = {}

def hgnc_to_gene_symbol():
	with open('../data/HGNC ids.txt', 'r') as infile:
		hgnc_reader = csv.reader(infile, dialect = 'excel-tab', skipinitialspace = True)
		next(hgnc_reader)
		for row in hgnc_reader:
			ids = row[0].split(':')[1]
			dictionary[ids] = row[1]
	# print(len(dictionary))


def ids_to_gene_symbol():
	with open('../data/HGNC ids.txt', 'r') as infile:
		hgnc_reader = csv.reader(infile, dialect = 'excel-tab', skipinitialspace = True)
		next(hgnc_reader)
		for row in hgnc_reader:
			hgnc_ids = row[0].split(':')[1]
			hgnc_dictionary[hgnc_ids] = row[1]

			ensembl_ids = row[4]
			ensembl_dictionary[ensembl_ids] = row[1]

'''
def hgnc_id_classification(path):
	# count = 0
	for file in os.listdir(path):
		print('For file: {0}'.format(str(file)))
		with open(path+file, 'r') as infile:
			lis = []
			reader = csv.reader(infile, dialect = 'excel-tab', skipinitialspace = True)
			for row in reader:
				if 'HGNC' in row[0]:
					hgnc_id = row[0].split('|')[1].split('=')[1]
					if dictionary[hgnc_id] in gene_lis:
						lis.append(dictionary[hgnc_id])
						# count = count + 1
			count_unique_genes(lis)
		class_name = (str(file)).split('_')[1].split('.txt')[0]
		classification_dictionary[class_name] = lis
	print(len(classification_dictionary))
	# print(count)
	jsonify(classification_dictionary, '../output/functional_classification.json')
'''


def classification(path):
	count = 0
	for file in os.listdir(path):
		print('For file: {0}'.format(str(file)))
		with open(path+file, 'r') as infile:
			lis = []
			reader = csv.reader(infile, dialect = 'excel-tab', skipinitialspace = True)
			for row in reader:
				if 'HGNC' in row[0]:
					hgnc_id = row[0].split('|')[1].split('=')[1]
					if hgnc_dictionary[hgnc_id] in gene_lis:
						lis.append(hgnc_dictionary[hgnc_id])
				elif 'Ensembl' in row[0]:
					# count+=1
					ensembl_id = row[0].split('|')[1].split('=')[1]
					if ensembl_id in ensembl_dictionary.keys() and ensembl_dictionary[ensembl_id] in gene_lis:
						lis.append(ensembl_id[ensembl_id])
			count_unique_genes(lis)
		class_name = (str(file)).split('_')[1].split('.txt')[0]
		# print(len(lis) == len(set(lis)))
		classification_dictionary[class_name] = lis
	print(len(classification_dictionary))
	# print(count)


def count_unique_genes(lis):
	for x in lis:
		if x not in uniq_gene_lis:
			uniq_gene_lis.append(x)
	print(len(uniq_gene_lis))
	# uniq_gene_lis.extend(lis)
	# print(len(set(uniq_gene_lis)))


def jsonify(dictionary, filename, text='None'):
	a = json.dumps(dictionary, sort_keys=True, indent=2, separators=(',', ': '))
	with open(str(filename), 'w') as outfile:
		if text == 'None':
			outfile.write(a)
		else:
			outfile.write(text + ' = ')
			outfile.write(a)


if __name__ == '__main__':
	with open('../data/genes.csv', 'r') as infile:
		reader = csv.reader(infile, dialect='excel-tab', skipinitialspace=True)
		for row in reader:
			gene_lis.append(row[0])

	# hgnc_to_gene_symbol()
	# hgnc_id_classification(class_files)

	ids_to_gene_symbol()
	classification(class_files)