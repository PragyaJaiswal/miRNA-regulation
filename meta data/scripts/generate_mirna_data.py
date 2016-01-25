#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys
import json, csv, re
import math

class restructure_data(object):
	"""docstring for restructure_data"""
	def __init__(self):
		pass

	def generate_map(self, mirtar):
		next(mirtar)
		mirna_map = {}
		# gene_map = {}
		for line in mirtar:
			mirna_map.setdefault(line[1], []).append(line[3])
			# gene_map.setdefault(line[3], []).append(line[1])
		print(len(mirna_map.keys()))
		# print(len(gene_map.keys()))

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
			jsonify(final_intronic_dict, '../output_data/mirna_map_dict.json')
			print(len(final_intronic_dict.keys()))
			return final_intronic_dict

def jsonify(dictionary, filename, text='None'):
	a = json.dumps(dictionary, sort_keys=True, indent=4, separators=(',', ': '))
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
		intronic_mirna_map = instance.generate_map(mirtar)

	# meta_data_intance = meta_data()
	# meta_data_intance.ensembl_coordinates_to_py(intronic_mirna_map)