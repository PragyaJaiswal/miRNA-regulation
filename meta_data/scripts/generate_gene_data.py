#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys
import json, csv, re
from Bio import SeqIO
import Bio

class restructure_data(object):
	"""docstring for restructure_data"""
	def __init__(self):
		pass
	
	def generate_map(self, mirtar):
		next(mirtar)
		gene_map = {}
		for line in mirtar:
			gene_map.setdefault(line[3], []).append(line[1])
		print(len(gene_map.keys()))

		jsonify(gene_map, '../output_data/gene/gene_map_dict.json')
		restructure_data.append_host_target_mirna(gene_map)

	def append_host_target_mirna(gene_map):
		'''
		gene_data_new contains the host as well as target genes as keys.
		Total number of keys are 15112.
		Hence, it has more than 14894 genes since, 14894 are the target
		genes (from mirTarBase) and the remaining are the host genes by
		coordinate matching of ensembl genes and mirna coordinates.
		'''
		with open('../output_data/mirna/mirna_meta_data_complete.json') as infile:
			miRNA_meta_data = json.loads(infile.read())
			gene_data_new = {}
			for mirna in miRNA_meta_data.keys():
				if 'Host Gene' in miRNA_meta_data[mirna].keys():
					gene = miRNA_meta_data[mirna]['Host Gene']
					
					gene_data_new.setdefault(gene, {})
					gene_data_new[gene].setdefault('Host for', []).append(mirna)

					if gene in gene_map.keys():
						gene_data_new[gene]['Target for'] = gene_map[gene]
					else:
						gene_data_new[gene]['Target for'] = ''

			for key, value in gene_map.items():
				gene_data_new.setdefault(key, {})
				gene_data_new[key]['Target for'] = value

			print(len(gene_data_new.keys()))
			jsonify(gene_data_new, '../output_data/gene/gene_meta_data.json')


def jsonify(dictionary, filename, text='None'):
	a = json.dumps(dictionary, sort_keys=True, indent=2, separators=(',', ': '))
	with open(str(filename), 'w') as outfile:
		if text == 'None':
			outfile.write(a)
		else:
			outfile.write(text + ' = ')
			outfile.write(a)


if __name__ == '__main__':
	instance = restructure_data()

	with open('../../data/hsa_MTI.tsv', 'r') as infile:
		mirtar = csv.reader(infile, dialect = 'excel-tab', skipinitialspace = True)
		gene_map = instance.generate_map(mirtar)