#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys
import json, csv

ensembl_dictionary = {}
gene_symbols = []

def convert():
	with open('Homo_sapiens_TF_EnsemblID.txt', 'r') as fl:
		data = csv.reader(fl, dialect = 'excel-tab', skipinitialspace = True)
		for gene in data:
			if gene[0] in ensembl_dictionary.keys():
				gene_symbols.append(ensembl_dictionary[gene[0]])
			else:
				print(gene[0])

	with open('converted.txt', 'w') as outfile:
		print(len(gene_symbols))
		json.dump(gene_symbols, outfile)

if __name__ == '__main__':
	with open('gene_ids.txt', 'r') as infile:
		hgnc_reader = csv.reader(infile, dialect = 'excel-tab', skipinitialspace = True)
		next(hgnc_reader)
		for row in hgnc_reader:
			ensembl_ids = row[4]
			ensembl_dictionary[ensembl_ids] = row[1]

	convert()