#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json, csv

dictionary = {}
mir_gen_dict = {}
tissue_wise_report = []

def process():
	for each in tissue_wise_report:
		dictionary[each['tissue']] = {}
		dictionary[each['tissue']]['mirnas'] = []
		dictionary[each['tissue']]['genes'] = []
		mirna_scores = each['mirnas']
		gene_scores = each['genes']

		sorted_mirna_scores, sorted_mirna_list = zip(*sorted(zip(mirna_scores, mir_gen_dict['mirnas']), reverse=True))
		for index, mirna in enumerate(sorted_mirna_list[:10]):
			dictionary[each['tissue']]['mirnas'].append((sorted_mirna_scores[index], mirna))

		sorted_gene_scores, sorted_gene_list = zip(*sorted(zip(gene_scores, mir_gen_dict['genes']), reverse=True))
		for index, gene in enumerate(sorted_gene_list[:10]):
			dictionary[each['tissue']]['genes'].append((sorted_gene_scores[index], gene))
		# print(dictionary)
		# input('Enter')
	jsonify(dictionary, 'tissue_wise_top_ranks.json')


def jsonify(dictionary, filename, text='None'):
	a = json.dumps(dictionary, sort_keys=True, indent=2, separators=(',', ': '))
	with open(str(filename), 'w') as outfile:
		if text == 'None':
			outfile.write(a)
		else:
			outfile.write(text + ' = ')
			outfile.write(a)

if __name__ == '__main__':
	with open('mir.gen.lists.json', 'r') as fl:
		mir_gen_dict = json.load(fl)

	with open('full_rpt.2016-04-24T200941.284010.json', 'r') as fl:
		tissue_wise_report = json.load(fl)

	process()