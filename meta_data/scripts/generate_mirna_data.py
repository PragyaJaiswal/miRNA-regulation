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
		mirna_map = {}
		for line in mirtar:
			mirna_map.setdefault(line[1], []).append(line[3])
		print(len(mirna_map.keys()))

		with open('../../data/chr_coordinates_of_mirna.csv', 'r') as mirna_file:
			mirna_coordinates = csv.reader(mirna_file, dialect = 'excel', skipinitialspace = True)
			intronic_mirna_map = restructure_data.check_intronic(mirna_map, mirna_coordinates)
			return intronic_mirna_map


	def check_intronic(mirna_map, mirna_coordinates):
		final_intronic_dict = {}
		with open('../../data/intron_coordinates_from_ucsc.tsv') as infile:
			intron_reader = csv.reader(infile, dialect = 'excel-tab', skipinitialspace = True)
			for line in mirna_coordinates:
				if line[3] in mirna_map.keys():
					chro = line[0]
					start = int(line[1])
					end = int(line[2])
					mirna = line[3]
					for each_line in intron_reader:
						if str(each_line[0].split('chr')[1]) == str(chro):
							if int(each_line[1]) <= start and end <= int(each_line[2]):
								final_intronic_dict[mirna] = mirna_map[mirna]
							else:
								pass
					infile.seek(0)
			jsonify(final_intronic_dict, '../output_data/mirna/mirna_map_dict.json')
			print(len(final_intronic_dict.keys()))
			return final_intronic_dict


class meta_data(object):
	"""docstring for meta_data"""
	def __init__(self):
		pass

	def ensembl_coordinates_to_py(self, intronic_mirna_map):
		print('Forming ensembl coordinates dictionary.')
		gene_dict = {}
		with open('../../data/gene_coordinates_from_ensembl.tsv', 'r') as gene_file:
			gene_data = csv.reader(gene_file, dialect = 'excel-tab', skipinitialspace = True)
			next(gene_data, None)
			for each_line in gene_data:
				if each_line[1] in gene_dict.keys():
					if each_line[0] in gene_dict[each_line[1]].keys():
						gene_dict[each_line[1]][each_line[0]]['start'] = int(each_line[2])
						gene_dict[each_line[1]][each_line[0]]['end'] = int(each_line[3])
						gene_dict[each_line[1]][each_line[0]]['gene_name'] = str(each_line[4])
						gene_dict[each_line[1]][each_line[0]]['transcript_count'] = int(each_line[5])
					else:
						gene_dict[each_line[1]][each_line[0]] = {}
						gene_dict[each_line[1]][each_line[0]]['start'] = int(each_line[2])
						gene_dict[each_line[1]][each_line[0]]['end'] = int(each_line[3])
						gene_dict[each_line[1]][each_line[0]]['gene_name'] = str(each_line[4])
						gene_dict[each_line[1]][each_line[0]]['transcript_count'] = int(each_line[5])
				else:
					gene_dict[each_line[1]] = {}
					if each_line[0] in gene_dict[each_line[1]].keys():
						gene_dict[each_line[1]][each_line[0]]['start'] = int(each_line[2])
						gene_dict[each_line[1]][each_line[0]]['end'] = int(each_line[3])
						gene_dict[each_line[1]][each_line[0]]['gene_name'] = str(each_line[4])
						gene_dict[each_line[1]][each_line[0]]['transcript_count'] = int(each_line[5])
					else:
						gene_dict[each_line[1]][each_line[0]] = {}
						gene_dict[each_line[1]][each_line[0]]['start'] = int(each_line[2])
						gene_dict[each_line[1]][each_line[0]]['end'] = int(each_line[3])
						gene_dict[each_line[1]][each_line[0]]['gene_name'] = str(each_line[4])
						gene_dict[each_line[1]][each_line[0]]['transcript_count'] = int(each_line[5])
			meta_data.find_host_gene(gene_dict, intronic_mirna_map)


	def find_host_gene(gene_dict, intronic_mirna_map):
		print('Finding hosts.')
		final_dict = {}
		with open('../../data/chr_coordinates_of_mirna.csv', 'r') as mirna_file:
			mirna_coordinates = csv.reader(mirna_file, dialect = 'excel', skipinitialspace = True)
			
			for line in mirna_coordinates:
				if line[3] in intronic_mirna_map.keys():
					final_dict[line[3]] = {}
					print('For miRNA : ' + str(line[3]))
					chro = str(line[0])
					mirna_start = line[1]
					mirna_end = line[2]
					for each in gene_dict[chro].keys():
						if gene_dict[chro][each]['start'] <= int(mirna_start) and int(mirna_end) <= gene_dict[chro][each]['end']:
							if 'MIR' in str(gene_dict[chro][each]['gene_name']):
								final_dict[line[3]]['miRNA Transcript Count'] = gene_dict[chro][each]['transcript_count']
								final_dict[line[3]]['miRNA Name'] = gene_dict[chro][each]['gene_name']
							else:
								final_dict[line[3]]['Host Gene'] = gene_dict[chro][each]['gene_name']
								final_dict[line[3]]['Host Gene Transcript Count'] = gene_dict[chro][each]['transcript_count']
								# print('Gene : ' + str(gene_dict[chro][each]['gene_name']))
								# print('In dict: ' + final_dict[line[3]]['Host Gene'])
								# print('Transcript count : ' + str(gene_dict[chro][each]['transcript_count']))
			for mirna in intronic_mirna_map.keys():
				if mirna in final_dict:
					pass
				else:
					final_dict[mirna] = {}
			
			meta_data.add_target_transcript_count(gene_dict, final_dict, intronic_mirna_map)

	
	def add_target_transcript_count(gene_dict, final_dict, intronic_mirna_map):
		mirmap_dict = meta_data.form_affinity_map()
		print('Appending targets.')
		count = 0
		for mirna in intronic_mirna_map.keys():
			count+=1
			tot = len(intronic_mirna_map.keys())
			print(('For miRNA: ' + str(mirna) + ', {0}/{1}').format(str(count), str(tot)))
			final_dict[mirna]['Target Gene with Transcript Count'] = []
			for target in intronic_mirna_map[mirna]:
				m = meta_data.target_gene_expression(gene_dict, target)		# Transcript Count
				aff = meta_data.append_affinity(mirmap_dict, mirna, target)	# Affinity
				tup = (target, m, aff)
				final_dict[mirna]['Target Gene with Transcript Count'].append(tup)
			# print(final_dict[mirna])
		jsonify(final_dict ,'../output_data/mirna/mirna_meta_data_test.json')
	
	'''
	def find_target_gene_expression(target_gene):
		with open('../../data/gene_coordinates_from_ensembl.tsv', 'r') as gene_file:
			gene_data = csv.reader(gene_file, dialect = 'excel-tab', skipinitialspace = True)
			next(gene_data, None)
			for each_line in gene_data:
				if str(each_line[4]) == str(target_gene):
					value = int(each_line[5])
					gene_file.seek(0)
					return value
				else:
					# Also returns 0 for the genes whose transcript count is not known to us
					value = int(0)
			gene_file.seek(0)
			return value
	'''

	def target_gene_expression(gene_dict, target_gene):
		value = None
		for chro in gene_dict.keys():
			for ensembl_id in gene_dict[chro].keys():
				if 'gene_name' in gene_dict[chro][ensembl_id].keys():
					if str(target_gene) == str(gene_dict[chro][ensembl_id]['gene_name']):
						value = gene_dict[chro][ensembl_id]['transcript_count']
						break
		return value

	def form_affinity_map():
		print('Forming affinity map from mirmap.')
		mirmap_dict = {}
		with open('../../data/sample.csv') as infile:
			mirmap_reader = csv.reader(infile, dialect = 'excel', skipinitialspace = True)
			next(mirmap_reader)
			for line in mirmap_reader:
				if not line[1] in mirmap_dict:
					mirmap_dict[line[1]] = {}
					mirmap_dict[line[1]].setdefault(line[8], []).append(line[20])
				else:
					mirmap_dict[line[1]].setdefault(line[8], []).append(line[20])
		return mirmap_dict

	def append_affinity(mirmap_dict, mirna, target):
		if mirna in mirmap_dict.keys() and target in mirmap_dict[mirna].keys():
			print('Found affinity value.')
			aff = min(mirmap_dict[mirna][target])
			return aff
			# tup = (target_tup[0], target_tup[1], aff)
		else:
			return None


class collect_meta_data_from_mirbase(object):
	"""docstring for collect_meta_data_from_mirbase"""
	def __init__(self):
		pass

	def extract(self, mirna_meta_data):
		mirbase_data = {}
		with open('../../data/miRNA.dat') as data_file:
			data = SeqIO.parse(data_file, 'embl')
			for record in data:
				mirna = record.name
				if 'hsa' in mirna:
					mirbase_data[mirna] = {}
					abc = {}
					if 'accessions' in record.annotations.keys():
						mirbase_data[mirna]['Accession ID'] = record.annotations['accessions']
					else:
						mirbase_data[mirna]['Accession ID'] = record.id

					mirbase_data[mirna]['Name'] = record.name
					mirbase_data[mirna]['Description'] = record.description
					mirbase_data[mirna]['Database cross-references'] = record.dbxrefs
					
					if 'comment' in record.annotations.keys():
						mirbase_data[mirna]['comment'] = record.annotations['comment']
					if 'references' in record.annotations.keys():
						# print type(record.annotations['references'])
						mirbase_data[mirna]['citations'] = {}
						for i in range(0,len(record.annotations['references'])):
							mirbase_data[mirna]['citations'][i] = {}
							mirbase_data[mirna]['citations'][i]['title'] = record.annotations['references'][i].title
							mirbase_data[mirna]['citations'][i]['authors'] = record.annotations['references'][i].authors
							mirbase_data[mirna]['citations'][i]['journal'] = record.annotations['references'][i].journal
					
					product_dict = {}
					for feature in record.features:
						# feature - contains a undecipherable format of the features
						# feature.qualifiers - a Python dictionary of additional decipherable information about the feature.
						# products - returns a list, but generally contains one product only
						# print feature.qualifiers
						if 'product' in feature.qualifiers.keys():
							products = feature.qualifiers['product']
							for product in products:
								product_dict[product] = {}
								for key, value in feature.qualifiers.items():
									product_dict[product][key] = value
								# print feature.qualifiers['experiment']
					mirbase_data[mirna]['products'] = product_dict
			mirna_meta_data_including_mirbase = collect_meta_data_from_mirbase.extend_meta_data(mirbase_data, mirna_meta_data)
			return mirna_meta_data_including_mirbase


	def extend_meta_data(mirbase_data, mirna_meta_data):
		mirna_meta_data_including_mirbase = {}
		for mirna in mirna_meta_data:
			mirna_meta_data_including_mirbase[mirna] = {}
			for fam_mirna in mirbase_data.keys():
				if mirna in mirbase_data[fam_mirna]['products'].keys():
					for keys, info in mirbase_data[fam_mirna]['products'][mirna].items():
						mirna_meta_data_including_mirbase[mirna][keys] = info

					for key, value in mirbase_data[fam_mirna].items():
						if not key == 'products':
							mirna_meta_data_including_mirbase[mirna][key] = value

					for key, value in mirna_meta_data[mirna].items():
						mirna_meta_data_including_mirbase[mirna][key] = value

					mirna_meta_data_including_mirbase[mirna]['family'] = fam_mirna

		jsonify(mirna_meta_data_including_mirbase, '../output_data/mirna/mirna_meta_data_including_mirbase.json')
		return mirna_meta_data_including_mirbase


class cross_references_from_ncbi(object):
	"""docstring for cross_references_from_ncbi"""
	def __init__(self):
		pass

	def id_dict(self, data, mirna_meta_data_including_mirbase):
		dictionary = {}
		for line in data:
			lis = []
			if 'miRBase' in line[5]:
				accession = line[5].split('miRBase:')[1]
				for ele in line[5].split('|'):
					if 'HGNC' in ele:
						hgnc_id = 'HGNC' + str(ele.split('HGNC')[2])
						lis.append(hgnc_id)
					else:
						lis.append(ele)
				dictionary[accession] = lis
		cross_references_from_ncbi.append_ids(dictionary, mirna_meta_data_including_mirbase)

	def append_ids(dictionary, mirna_meta_data_including_mirbase):
		more_ids = {}
		for mirna in mirna_meta_data_including_mirbase.keys():
			if 'Accession ID' in mirna_meta_data_including_mirbase[mirna].keys():
				accession = mirna_meta_data_including_mirbase[mirna]['Accession ID'][0]
				lis = dictionary[accession]
				# print lis
				for ele in lis:
					if not ele in mirna_meta_data_including_mirbase[mirna]['Database cross-references']:
						mirna_meta_data_including_mirbase[mirna]['Database cross-references'].append(ele)
		jsonify(mirna_meta_data_including_mirbase, '../output_data/mirna/mirna_meta_data_complete.json')
		weights.mmi(mirna_meta_data_including_mirbase)


class weights(object):
	"""docstring for weights"""
	def __init__(self):
		pass

	def weights(self, mirna_meta_data_complete):
		mirna_meta_data_with_weights = {}
		for mirna in mirna_meta_data_complete.keys():
			mirna_meta_data_with_weights[mirna] = {}
			for key, value in mirna_meta_data_complete[mirna].items():
				if not key == 'Target Gene with Transcript Count':
					mirna_meta_data_with_weights[mirna][key] = value
			
			if 'Target Gene with Transcript Count' in mirna_meta_data_complete[mirna].keys():
				for each_target in mirna_meta_data_complete[mirna]['Target Gene with Transcript Count']:
					if each_target[2] == None:
						mmi = None
					else:
						del_g_binding = each_target[2]
						keq = float(math.exp(-1 * del_g_binding/(0.008314 * 298)))
						if 'Host Gene'in miRNA_meta_data[mirna].keys() and not miRNA_meta_data[mirna]['Host Gene'] == '':
							m = each_target[1]
							mi = miRNA_meta_data[mirna]['Host Gene Transcript Count']
							mmi = keq * m * mi
							# print(mmi)
					each_target.append(mmi)
					mirna_meta_data_with_weights[mirna]['Target Gene with Transcript Count'].append(each_target)
		jsonify(mirna_meta_data_with_weights, '../output_data/mirna/mirna_meta_data_with_weights.json')
		

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
	meta_data_instance = meta_data()
	mirbase_meta_data_instance = collect_meta_data_from_mirbase()
	cross_references_instance = cross_references_from_ncbi()
	with open('../../data/hsa_MTI.tsv', 'r') as infile:
		mirtar = csv.reader(infile, dialect = 'excel-tab', skipinitialspace = True)
		intronic_mirna_map = instance.generate_map(mirtar)

		mirna_meta_data = meta_data_instance.ensembl_coordinates_to_py(intronic_mirna_map)
	# with open('../output_data/mirna/mirna_map_dict.json', 'r') as infile:
	# 	intronic_mirna_map = json.loads(infile.read())
	# 	mirna_meta_data = meta_data_instance.ensembl_coordinates_to_py(intronic_mirna_map)

		mirna_meta_data_including_mirbase = mirbase_meta_data_instance.extract(mirna_meta_data)
	# with open('../output_data/mirna/mirna_meta_data_test.json', 'r') as infile:
	# 	mirna_meta_data = json.loads(infile.read())
	# 	mirna_meta_data_including_mirbase = mirbase_meta_data_instance.extract(mirna_meta_data)

		with open('../../data/Homo_sapiens.gene_info', 'r') as infile:
			data = csv.reader(infile, 'excel-tab')
			next(data)
			cross_references_instance.id_dict(data, mirna_meta_data_including_mirbase)