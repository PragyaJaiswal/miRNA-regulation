#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json, csv

dictionary = {}
gene_lis = []

def pathway_to_genes():
	with open('../data/CPDB_pathways_genes.tab', 'r') as infile:
		reader = csv.reader(infile, dialect='excel-tab', skipinitialspace=True)
		next(reader)
		for row in reader:
			lis = []
			for gene in row[3].split(','):
				if gene in gene_lis:
					lis.append(gene)
			dictionary[row[0]] = lis
	jsonify(dictionary, '../output/pathway_to_genes.json')

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
	pathway_to_genes()