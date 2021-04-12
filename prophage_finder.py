import ast
import json
import os

import pandas as pd

from FeatureConstruction import *


def phages_bact():
	count_bacteria = 0
	for phage in data.index:
		if ast.literal_eval(data.loc[phage, 'Host_ID']):
			count_bacteria += 1
	return count_bacteria


data = pd.read_csv('files/NCBI_Phage_Bacteria_Data.csv', header=0, index_col=0)

with open('C:/Users/Pedro/Downloads/pha_in_bac_2_test.json', encoding='utf-8') as F:
	prophage = json.loads(F.read())

for bact in prophage.keys():
	for phage in prophage[bact]:
		if phage in data.index:
			temp = ast.literal_eval(data.loc[phage, 'Host_ID'])
			if bact + '.1' not in temp:
				temp.append(bact+'.1')
				data.loc[phage, 'Host_ID'] = str(temp)

data.to_csv('files/NCBI_Phage_Bacteria_Data.csv')

fc = FeatureConstruction()
phageTails = fc.phageTails

os.system('cd-hit -i files/tails.fasta -o files/cdhit')

temp_cluster = []
with open('files/cdhit.clstr', 'r') as f:
	for line in f.readlines():
		if '>Cluster' in line:
			for prot in temp_cluster:
				for phage in phageTails:
					if prot in phageTails[phage].keys():
						if phage in data.index:
							temp_ref = ast.literal_eval(data.loc[ref_phage, 'Bacteria ID'])
							temp = ast.literal_eval(data.loc[phage, 'Bacteria ID'])
							for i in temp_ref:
								if i not in temp:
									temp.append(i)
									data.loc[phage, 'Bacteria ID'] = str(temp)
						break
			temp_cluster = []
		elif line[0] == '0':
			pos_i = line.find('>') + 1
			pos_f = line.find('...')
			ref_prot = line[pos_i:pos_f]
			for phage in phageTails:
				if ref_prot in phageTails[phage].keys():
					ref_phage = phage
					break
		else:
			pos_i = line.find('>') + 1
			pos_f = line.find('...')
			temp_cluster.append(line[pos_i:pos_f])


with open('files/bactDNA.json', encoding='utf-8') as F:
	bacProt = json.loads(F.read())

listDone = []
for bact in bacProt:
	if bact in listDone:
		pass
	else:
		listDone.append(bact)
		with open('files/temp_genome.fasta', 'w') as F:
			F.write('>' + bact + '\n' + bacProt[bact] + '\n')
		os.system('phigaro -f files/temp_genome.fasta --not-open -d -o files/temp_phigaro')  # Phigaro
		with open('files/temp_phigaro.html', 'r') as Ph:
			tempPhigaro = Ph.readlines()
		for line in tempPhigaro:
			if '<div class="accordion-body collapse"' in line:
				VOGs = line[line.find('>')+1:].strip('\n').split(', ')
				for vog in VOGs:
					with open('files/VOG_tables/' + vog + '.txt', 'r') as f:
						temp_phages = f.readlines()
					for i in range(len(temp_phages)):
						if i != 0:
							phage = temp_phages[i].split('\t')[2]
							if phage in data.index:
								temp = ast.literal_eval(data.loc[phage, 'Bacteria ID'])
								if bact not in temp:
									temp.append(bact)
									data.loc[phage, 'Bacteria ID'] = str(temp)
		print('Number of phages with associated bacteria strains:', phages_bact(), end="\r")

'''os.system('wget --post-file="files/temp_genome.fasta" "http://phaster.ca/phaster_api?contigs=1" -O files/temp_phaster') # Phaster
with open('files/temp_phaster', encoding='utf-8') as F:
	temp = json.loads(F.read())
os.system('wget "http://phaster.ca/phaster_api?acc=' + temp['job_id'] + 'Z" -O files/temp_phaster') # servidor cheio
os.system('PhiSpy.py files/temp_genome.fasta -o files/temp_phipsy') # Phipsy - não possível com fastas'''
