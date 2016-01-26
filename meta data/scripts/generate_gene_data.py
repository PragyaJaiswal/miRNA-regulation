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

		jsonify(final_intronic_dict, '../output_data/gene/gene_map_dict.json')
