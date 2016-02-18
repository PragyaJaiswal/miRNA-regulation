#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import json
import math
import operator


class stats(object):
	"""docstring for stats"""
	def __init__(self):
		pass
	
	def all_target_interactions(self, miRNA_meta_data):
		targets = 0
		for mirna in miRNA_meta_data.keys():
			if 'Target Gene with Transcript Count' in miRNA_meta_data[mirna].keys():
				targets+= len(miRNA_meta_data[mirna]['Target Gene with Transcript Count'])
				# print(targets)
				# raw_input('Enter')
		print('Total target interactions:\n' + str(targets))

	def all_host_interactions(self, miRNA_meta_data):
		hosts = 0
		for mirna in miRNA_meta_data.keys():
			if 'Host Gene' in miRNA_meta_data[mirna].keys():
				hosts+=1

		print('Total host interactions:\n' + str(hosts))

	def range_of_gene_trans_count_affinity(self, miRNA_meta_data):
		target_transcript_count_dict = {}
		target_affinity_dict = {}
		for mirna in miRNA_meta_data:
			if 'Host Gene' in miRNA_meta_data[mirna].keys():
				gene = miRNA_meta_data[mirna]['Host Gene']
				target_transcript_count_dict[gene] = miRNA_meta_data[mirna]['Host Gene Transcript Count']
			
			if 'Target Gene with Transcript Count' in miRNA_meta_data[mirna].keys():
				for target in miRNA_meta_data[mirna]['Target Gene with Transcript Count']:
					if not target[1] == None:
						target_transcript_count_dict[target[0]] = target[1]
					if not target[2] == None:
						target_affinity_dict[target[0]] = target[2]

		mini_trans_count = min(target_transcript_count_dict.values())
		maxi_trans_count = max(target_transcript_count_dict.values())		
		mini_trans_count_lis = []
		maxi_trans_count_lis = []
		
		for key, value in target_transcript_count_dict.items():
			if value == mini_trans_count:
				mini_trans_count_lis.append(key)
			if value == maxi_trans_count:
				maxi_trans_count_lis.append(key)

		trans_range = (mini_trans_count, maxi_trans_count)
		print('\nRange for transcript counts of genes:\n' + str(trans_range))

		# Negative values. Larger negative value means smaller affinity.
		lis = [x for x in target_affinity_dict.values() if x is not None]

		mini_affinity = max(lis)
		maxi_affinity = min(lis)
		mini_affinity_lis = []
		maxi_affinity_lis = []

		for key, value in target_affinity_dict.items():
			if value == mini_affinity:
				maxi_affinity_lis.append(key)
			if value == maxi_affinity:
				mini_affinity_lis.append(key)

		affinity_range = (mini_affinity, maxi_affinity)

		print('Range for affinity values of genes:\n' + str(affinity_range))


class mirna_stats(object):
	"""docstring for mirna_stats"""
	def __init__(self):
		pass
		
	def mirna_with_max_interactions(self, miRNA_meta_data):
		mirna_interact_dict = {}
		lis = []
		for mirna in miRNA_meta_data.keys():
			if 'Target Gene with Transcript Count' in miRNA_meta_data[mirna].keys():
				mirna_interact_dict[mirna] = len(miRNA_meta_data[mirna]['Target Gene with Transcript Count'])

		mirna_having_max_interactions = max(mirna_interact_dict.items(), key=operator.itemgetter(1))
		number_of_interactions = max(mirna_interact_dict.values())

		for key, value in mirna_interact_dict.items():
			if value == number_of_interactions:
				lis.append((key,value))

		print('\nmiRNA with maximum interactions:\n' + str(lis))

	def mirna_interacting_with_host(self, miRNA_meta_data):
		miRNA_lis = []
		mirna_targets = {}
		for mirna in miRNA_meta_data.keys():
			target_lis = []
			if 'Host Gene' in miRNA_meta_data[mirna].keys():
				host = miRNA_meta_data[mirna]['Host Gene']
				if 'Target Gene with Transcript Count' in miRNA_meta_data[mirna].keys():
					mirna_targets[mirna] = len(miRNA_meta_data[mirna]['Target Gene with Transcript Count'])
					for target in miRNA_meta_data[mirna]['Target Gene with Transcript Count']:
						target_lis.append(target[0])

			if host in target_lis:
				miRNA_lis.append(mirna)
		
		print('Number of miRNAs targeting the host genes:\n' + str(len(miRNA_lis)))
		print('miRNAs targeting the host genes:')

		lis = []
		for mirna in miRNA_lis:
			lis.append((mirna, mirna_targets[mirna], miRNA_meta_data[mirna]['Host Gene']))
		print(lis)
		print(len(lis))

	def genes_targetting_mirna_whose_host_is_known(self, miRNA_meta_data):
		target_lis = []
		for mirna in miRNA_meta_data.keys():
			if 'Host Gene' in miRNA_meta_data[mirna].keys():
				if 'Target Gene with Transcript Count' in miRNA_meta_data[mirna].keys():
					for target in miRNA_meta_data[mirna]['Target Gene with Transcript Count']:
						if not target[0] in target_lis:
							target_lis.append(target[0])

		print('\nNumber of genes targetted by miRNAs whose host is known:\n' + str(len(target_lis)))


class gene_stats(object):
	"""docstring for gene_stats"""
	def __init__(self):
		pass

	def gene_with_maximum_interactions(self, gene_meta_data):
		gene_interact_dict = {}
		lis = []
		for gene in gene_meta_data.keys():
			if 'Target for' in gene_meta_data[gene].keys():
				gene_interact_dict[gene] = len(gene_meta_data[gene]['Target for'])

		gene_having_max_interactions = max(gene_interact_dict.items(), key=operator.itemgetter(1))
		number_of_interactions = max(gene_interact_dict.values())

		for key, value in gene_interact_dict.items():
			if value == number_of_interactions:
				lis.append((key,value))

		print('\nGene being targetted most by miRNAs:\n' + str(lis))

		lis = sorted(gene_interact_dict.items(), key=operator.itemgetter(1), reverse=True)
		print(lis[0:20])

	def gene_having_mirna_as_host_and_target(self, gene_meta_data):
		gene_lis = []
		for gene in gene_meta_data.keys():
			host_lis = []
			target_lis = gene_meta_data[gene]['Target for']
			if 'Host for' in gene_meta_data[gene]:
				host_lis = gene_meta_data[gene]['Host for']
			for mirna in target_lis:
				if mirna in host_lis:
					gene_lis.append(gene)

			for mirna in host_lis:
				if mirna in target_lis:
					if gene not in gene_lis:
						gene_lis.append(gene)
		
		gene_lis = list(set(gene_lis))
		print('Number of host genes being targetted by miRNAs:\n' + str(len(gene_lis)))
		print('Host Genes being targetted by miRNAs:\n' + str(gene_lis))


class network(object):
	def __init__(self):
		pass

	def total_edges():
		pass
		

if __name__ == '__main__':
	stats_instance = stats()
	gene_stats_instance = gene_stats()
	mirna_stats_instance = mirna_stats()

	with open('../meta_data/output_data/mirna/mirna_meta_data_complete.json', 'r') as infile:
		miRNA_meta_data = json.loads(infile.read())
		stats_instance.all_target_interactions(miRNA_meta_data)
		stats_instance.all_host_interactions(miRNA_meta_data)
		stats_instance.range_of_gene_trans_count_affinity(miRNA_meta_data)

		mirna_stats_instance.mirna_with_max_interactions(miRNA_meta_data)
		mirna_stats_instance.mirna_interacting_with_host(miRNA_meta_data)
		mirna_stats_instance.genes_targetting_mirna_whose_host_is_known(miRNA_meta_data)

	with open('../meta_data/output_data/gene/gene_meta_data.json', 'r') as infile:
		gene_meta_data = json.loads(infile.read())
		gene_stats_instance.gene_with_maximum_interactions(gene_meta_data)
		gene_stats_instance.gene_having_mirna_as_host_and_target(gene_meta_data)