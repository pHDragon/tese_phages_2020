
class FeatureConstruction:

	def __init__(self):
		"""
		In development. Extract features from proteins.
		"""
		import pandas as pd
		import json
		import ast
		from pathlib import Path
		import os
		from random import randint
		data = pd.read_csv('files/NCBI_Phage_Bacteria_Data.csv', header=0, index_col=0)
		with open('files/phagesProteins.json', encoding='utf-8') as F:
			self.phagesProteins = json.loads(F.read())
		self._filter_phage_domains()
		# with open('files/bactProteins.json', encoding='utf-8') as F:
		# 	self.bactProteins = json.loads(F.read())
		# self._filter_bacteria()
		all_phages = {}
		ecoli = {}
		kpneumoniae = {}
		abaumannii = {}
		my_file = Path("files/FeatureDataset")
		if not my_file.is_file():
			for phage in self.phageTails:
				if phage in data.index and self.phageTails[phage]:
					for bact in ast.literal_eval(data.loc[phage, 'Host_ID']):
						bact = bact[:-2]
						if bact + '.json' in os.listdir('files/bacteria'):
							# if self.externalProts[bact]: # This verification is not necessary for carbohydrates
							all_phages[phage + '--' + bact] = 'Yes'
							name = data.loc[phage, 'Host']
							if 'escherichia' in name.lower() or 'coli' in name.lower():
								ecoli[bact] = 0
							elif 'klebsiella' in name.lower() or 'pneumoniae' in name.lower():
								kpneumoniae[bact] = 0
							elif 'acinetobacter' in name.lower() or 'baumannii' in name.lower():
								abaumannii[bact] = 0
			for phage in self.phageTails:
				if phage in data.index and self.phageTails[phage]:
				# if self.phageTails[phage]:
					name = data.loc[phage, 'Host']
					if 'escherichia' in name.lower() or 'coli' in name.lower():
						i = 0
						while i < 12:
							bact = list(kpneumoniae.keys())[randint(0, len(kpneumoniae.keys()) - 1)]
							all_phages[phage + '--' + bact] = 'No'
							i += 1
					elif 'klebsiella' in name.lower() or 'pneumoniae' in name.lower():
						i = 0
						while i < 12:
							bact = list(ecoli.keys())[randint(0, len(ecoli.keys()) - 1)]
							all_phages[phage + '--' + bact] = 'No'
							i += 1
					elif 'acinetobacter' in name.lower() or 'baumannii' in name.lower():
						i = 0
						while i < 12:
							bact = list(kpneumoniae.keys())[randint(0, len(kpneumoniae.keys()) - 1)]
							all_phages[phage + '--' + bact] = 'No'
							i += 1
			self.features_data = pd.DataFrame({'ID': list(all_phages.keys()), 'Infects': list(all_phages.values())})
			self.features_data = self.features_data.set_index('ID')
		else:
			self.import_feat_data()

	def _filter_phage_domains(self):
		import json
		from pathlib import Path
		'''
		Filters out unwanted proteins. Domains that are unknown or are not associated with fibers, spikes, tails, enzymatic or binding are not considered.
		Still in development.
		:return: phageTails, a dictionary containing only
		'''
		my_file = Path("files/phageTails.json")
		if not my_file.is_file():
			self.phageTails = {}
			for phage in self.phagesProteins:
				self.phageTails[phage] = {}
				for protein in self.phagesProteins[phage]:
					if any(z in self.phagesProteins[phage][protein][0].lower() for z in ['fiber', 'fibre', 'spike', 'hydrolase', 'bind', 'depolymerase', 'peptidase', 'lyase', 'sialidase', 'dextranase', 'lipase', 'adhesin', 'baseplate', 'protein h', 'recognizing'
																						 'protein j', 'protein g', 'gpe', 'duf4035', 'host specifity', 'cor protein', 'specificity', 'baseplate component', 'gp38', 'gp12 tail',  'receptor', 'recognition', 'tail']) \
							and not any(z in self.phagesProteins[phage][protein][0].lower() for z in ['nucle', 'dna', 'rna', 'ligase', 'transferase', 'inhibitor', 'assembly', 'connect', 'nudix', 'atp', 'nad', 'transpos', 'ntp', 'molybdenum', 'hns',
																								 'gtp', 'riib', 'inhibitor', 'replicat', 'codon', 'pyruvate', 'catalyst', 'hinge', 'sheath completion', 'head', 'capsid', 'tape', 'tip', 'strand', 'matur', 'portal'
																									'terminase', 'nucl', 'promot', 'block', 'olfact', 'wedge', 'lysozyme', 'mur', 'sheat']):
						self.phageTails[phage][protein] = self.phagesProteins[phage][protein]
					'''else:
						for i in self.phagesProteins[phage][protein]:
							if type(i) == str:
								if any(z in str(i).lower() for z in ['fiber', 'fibre', 'spike', 'hydrolase', 'bind', 'depolymerase', 'peptidase', 'lyase', 'sialidase', 'dextranase', 'lipase', 'adhesin', 'baseplate', 'protein h', 'recognizing'
																	'protein j', 'protein g', 'gpe', 'duf4035', 'host specifity', 'cor protein', 'specificity', 'baseplate component', 'gp38', 'gp12 tail',  'receptor', 'recognition', 'tail']) \
										and not any(z in str(i).lower() for z in ['nucle', 'dna', 'rna', 'ligase', 'transferase', 'inhibitor', 'assembly', 'connect', 'nudix', 'atp', 'nad', 'transpos', 'ntp', 'molybdenum', 'hns', 'gtp',
																				  'riib', 'inhibitor', 'replicat', 'codon', 'pyruvate', 'catalyst', 'hinge', 'sheath completion', 'head', 'capsid', 'tape', 'tip', 'strand', 'matur', 'portal'
																				'terminase', 'nucl']):
									self.phageTails[phage][protein] = self.phagesProteins[phage][protein]
							else:
								for j in i:
									if any(z in str(j).lower() for z in ['fiber', 'fibre', 'spike', 'hydrolase', 'bind', 'depolymerase', 'peptidase', 'lyase', 'sialidase', 'dextranase', 'lipase', 'adhesin', 'baseplate', 'protein h', 'recognizing'
																		'protein j', 'protein g', 'gpe', 'duf4035', 'host specifity', 'cor protein', 'specificity', 'baseplate component', 'gp38', 'gp12 tail',  'receptor', 'recognition', 'tail']) \
											and not any(z in str(j).lower() for z in ['nucle', 'dna', 'rna', 'ligase', 'transferase', 'inhibitor', 'assembly', 'connect', 'nudix', 'atp', 'nad', 'transpos', 'ntp', 'molybdenum', 'hns', 'gtp',
																					  'riib', 'inhibitor', 'replicat', 'codon', 'pyruvate', 'catalyst', 'hinge', 'sheath completion', 'head', 'capsid', 'tape', 'tip', 'strand', 'matur', 'portal'
																					'terminase', 'nucl']):
										self.phageTails[phage][protein] = self.phagesProteins[phage][protein]'''
			with open('files/phageTails.json', 'w') as f:
				json.dump(self.phageTails, f)
			self.__create_phage_fasta()
		else:
			with open('files/phageTails.json', encoding='utf-8') as F:
				self.phageTails = json.loads(F.read())
		return self.phageTails

	def _filter_bacteria(self):
		import json
		from pathlib import Path
		import pandas as pd
		my_file = Path("files/externalProts.json")
		if not my_file.is_file():
			self.externalProts = {}
			predictions = pd.read_csv('files/results_psort.txt', sep='\t', index_col=False)
			predictions = predictions.set_index('SeqID')
			predictions = predictions.drop_duplicates()
			for bac in self.bactProteins:
				self.externalProts[bac] = {}
				for protein in self.bactProteins[bac]:
					if protein + ' ' in predictions.index:
						maxScore = 0.0
						for loc in ['Cytoplasmic_Score', 'CytoplasmicMembrane_Score', 'Periplasmic_Score', 'OuterMembrane_Score', 'Extracellular_Score']:
							if predictions.loc[protein + ' ', loc] > maxScore:
								maxScore = predictions.loc[protein + ' ', loc]
								location = loc
						if location == 'CytoplasmicMembrane_Score' or location == 'OuterMembrane_Score' or location == 'Extracellular_Score':
							self.externalProts[bac][protein] = self.bactProteins[bac][protein][1]
			if self.externalProts != {}:
				del self.bactProteins
			with open('files/externalProts.json', 'w') as f:
				json.dump(self.externalProts, f)
		else:
			with open('files/externalProts.json', encoding='utf-8') as F:
				self.externalProts = json.loads(F.read())
		return self.externalProts

	def __create_phage_fasta(self):
		"""
		Creates a fasta file containing every protein sequence for every phage.
		:return:
		"""
		with open('files/tails.fasta', 'w') as F:
			for phage in self.phageTails:
				for prot in self.phageTails[phage]:
					F.write('>' + prot + '\n' + self.phageTails[phage][prot][1] + '\n')

	def add_kmers(self):
		from skbio import Sequence
		import json
		groups = '0123456'
		freqs = {}
		for i in groups:
			for j in groups:
				freqs[i+j] = 0.0
		for i in freqs:
			exec('phage_group_{0} = []'.format(i))
			exec('bact_group_{0} = []'.format(i))
		phage = ''
		bact = ''
		for ID in self.features_data.index:
			done_phage = False
			done_bact = False
			if ID[:ID.find('--')] == phage:
				for i in freqs.keys():
					exec('phage_group_{0}.append(phage_group_{0}[-1])'.format(i))
				done_phage = True
			if ID[ID.find('--') + 2:] == bact:
				for i in freqs.keys():
					exec('bact_group_{0}.append(bact_group_{0}[-1])'.format(i))
				done_bact = True
			bact = ID[ID.find('--') + 2:]
			phage = ID[:ID.find('--')]

			if not done_phage:
				totalKmers = freqs.copy()
				count_prots = 0
				for prot in self.list_prot[phage]:
					max_freq = 0.0
					min_freq = 1000000.0
					count_prots += 1
					seq = self.__get_conjoint_triad(self.list_prot[phage][prot][1].replace('X', 'D').replace('B', 'N').replace('Z', 'E').replace('U', 'C').replace('J', 'L'))
					seq = Sequence(seq)
					temp = seq.kmer_frequencies(2, overlap=True, relative=True)
					for i in temp.keys():  # para normalizar
						if temp[i] < min_freq:
							min_freq = temp[i]
						if temp[i] > max_freq:
							max_freq = temp[i]
					for i in temp.keys():
						totalKmers[i] += temp[i] - (min_freq / max_freq)
				if count_prots != 0:
					for i in totalKmers.keys():
						totalKmers[i] = totalKmers[i] / count_prots
						temp_value = totalKmers[i]
						exec('phage_group_{0}.append(temp_value)'.format(i))
				else:
					for i in totalKmers.keys():
						exec('phage_group_{0}.append(0.0)'.format(i))

			if not done_bact:
				totalKmers = freqs.copy()
				count_prots = 0
				with open('files/bacteria/' + bact + '.json', encoding='utf-8') as F:
					bact_prots = json.loads(F.read())
				for prot in bact_prots:
					max_freq = 0.0
					min_freq = 1000000.0
					count_prots += 1
					seq = bact_prots[prot][1]
					seq = seq[:seq.find('"')]
					seq = self.__get_conjoint_triad(seq.replace('X', 'D').replace('B', 'N').replace('Z', 'E').replace('U', 'C').replace('J', 'L'))
					seq = Sequence(seq)
					temp = seq.kmer_frequencies(2, overlap=True, relative=True)
					for i in temp.keys():  # para normalizar
						if temp[i] < min_freq:
							min_freq = temp[i]
						if temp[i] > max_freq:
							max_freq = temp[i]
					for i in temp.keys():
						totalKmers[i] += temp[i] - (min_freq / max_freq)
				if count_prots != 0:
					for i in totalKmers.keys():
						totalKmers[i] = totalKmers[i] / count_prots
						temp_value = totalKmers[i]
						exec('bact_group_{0}.append(temp_value)'.format(i))
				else:
					for i in freqs.keys():
						exec('bact_group_{0}.append(0.0)'.format(i))

		for i in freqs.keys():
			exec('self.features_data["phage_kmer_{0}"] = phage_group_{0}'.format(i))
			exec('self.features_data["bact_kmer_{0}"] = bact_group_{0}'.format(i))

	def get_kmers(self, phage, bacteria):
		from skbio import Sequence
		solution = []
		groups = '0123456'
		freqs = {}
		for i in groups:
			for j in groups:
				freqs[i+j] = 0.0
		for i in freqs:
			exec('phage_group_{0} = 0.0'.format(i))
			exec('bact_group_{0} = 0.0'.format(i))

		totalKmers = freqs.copy()
		count_prots = 0
		for prot in phage:
			max_freq = 0.0
			min_freq = 1000000.0
			count_prots += 1
			seq = self.__get_conjoint_triad(phage[prot][1].replace('X', 'D').replace('B', 'N').replace('Z', 'E').replace('U', 'C'))
			seq = Sequence(seq)
			temp = seq.kmer_frequencies(2, overlap=True, relative=True)
			for i in temp.keys(): #  para normalizar
				if temp[i] < min_freq:
					min_freq = temp[i]
				if temp[i] > max_freq:
					max_freq = temp[i]
			for i in temp.keys():
				totalKmers[i] += temp[i] - (min_freq / max_freq)
		if count_prots != 0:
			for i in totalKmers.keys():
				totalKmers[i] = totalKmers[i] / count_prots
				temp_value = totalKmers[i]
				exec('phage_group_{0} += temp_value'.format(i))

		totalKmers = freqs.copy()
		count_prots = 0
		for prot in bacteria:
			max_freq = 0.0
			min_freq = 1000000.0
			count_prots += 1
			seq = bacteria[prot][1]
			seq = seq[:seq.find('"')]
			seq = self.__get_conjoint_triad(seq.replace('X', 'D').replace('B', 'N').replace('Z', 'E').replace('U', 'C'))
			seq = Sequence(seq)
			temp = seq.kmer_frequencies(2, overlap=True, relative=True)
			for i in temp.keys():  # para normalizar
				if temp[i] < min_freq:
					min_freq = temp[i]
				if temp[i] > max_freq:
					max_freq = temp[i]
			for i in temp.keys():
				totalKmers[i] += temp[i] - (min_freq / max_freq)
		if count_prots != 0:
			for i in totalKmers.keys():
				totalKmers[i] = totalKmers[i] / count_prots
				temp_value = totalKmers[i]
				exec('bact_group_{0} += temp_value'.format(i))

		for i in freqs.keys():
			exec('solution.append(phage_group_{0})'.format(i))
			exec('solution.append(bact_group_{0})'.format(i))
		return solution

	def add_composition(self):
		from skbio import Sequence
		import json
		bact_comp = {}
		phage_comp = {}
		groups = '0123456'
		for i in groups:
			bact_comp['comp_' + i] = []
			phage_comp['comp_' + i] = []
		phage = ''
		bact = ''
		count = -1
		for ID in self.features_data.index:
			done_phage = False
			done_bact = False
			count += 1
			if ID[:ID.find('--')] == phage:
				for i in groups:
					phage_comp['comp_' + i].append(phage_comp['comp_' + i][-1])
				done_phage = True
			if ID[ID.find('--') + 2:] == bact:
				for i in groups:
					bact_comp['comp_' + i].append(bact_comp['comp_' + i][-1])
				done_bact = True
			bact = ID[ID.find('--') + 2:]
			phage = ID[:ID.find('--')]

			if not done_phage:
				count_prots = 0
				for i in groups:
					phage_comp['comp_' + i].append(0)
				for prot in self.list_prot[phage]:
					max_comp = 0.0
					min_comp = 1000000.0
					count_prots += 1
					seq = self.__get_conjoint_triad(self.list_prot[phage][prot][1].replace('X', 'D').replace('B', 'N').replace('Z', 'E').replace('U', 'C'))
					seq = Sequence(seq)
					for i in groups: #  para normalizar
						if seq.count(i) < min_comp:
							min_comp = seq.count(i)
						if seq.count(i) > max_comp:
							max_comp = seq.count(i)
					for i in groups:
						phage_comp['comp_' + i][count] += seq.count(i) - (min_comp / max_comp)
				total = 0
				if count_prots != 0:
					for i in groups:
						phage_comp['comp_' + i][count] = phage_comp['comp_' + i][count] / count_prots
						total += phage_comp['comp_' + i][count]
					for i in groups:
						phage_comp['comp_' + i][count] = phage_comp['comp_' + i][count] / total
				else:
					for i in groups:
						phage_comp['comp_' + i][count] = 0.0

			if not done_bact:
				count_prots = 0
				for i in groups:
					bact_comp['comp_' + i].append(0)
				with open('files/bacteria/' + bact + '.json', encoding='utf-8') as F:
					bact_prots = json.loads(F.read())
				for prot in bact_prots:
					max_comp = 0.0
					min_comp = 1000000.0
					count_prots += 1
					seq = bact_prots[prot][1]
					seq = seq[:seq.find('"')]
					seq = self.__get_conjoint_triad(seq.replace('X', 'D').replace('B', 'N').replace('Z', 'E').replace('U', 'C'))
					seq = Sequence(seq)
					for i in groups:
						if seq.count(i) < min_comp:
							min_comp = seq.count(i)
						if seq.count(i) > max_comp:
							max_comp = seq.count(i)
					for i in groups:
						bact_comp['comp_' + i][count] += seq.count(i) - (min_comp / max_comp)
				total = 0
				if count_prots != 0:
					for i in groups:
						bact_comp['comp_' + i][count] = bact_comp['comp_' + i][count] / count_prots
						total += bact_comp['comp_' + i][count]
				else:
					for i in groups:
						bact_comp['comp_' + i][count] = 0.0
				if total != 0:
					for i in groups:
						bact_comp['comp_' + i][count] = bact_comp['comp_' + i][count] / total
				else:
					for i in groups:
						bact_comp['comp_' + i][count] = 0.0

		for i in groups:
			self.features_data['bact_comp_' + i] = bact_comp['comp_' + i]
			self.features_data['phage_comp_' + i] = phage_comp['comp_' + i]

	def get_composition(self, phage, bacteria):
		from skbio import Sequence
		solution = []
		bact_comp = {}
		phage_comp = {}
		phage_comp_carb = {}
		groups = '0123456'
		for i in groups:
			bact_comp['comp_' + i] = 0
			phage_comp['comp_' + i] = 0
		count_prots = 0
		for prot in phage:
			max_comp = 0.0
			min_comp = 1000000.0
			count_prots += 1
			seq = self.__get_conjoint_triad(phage[prot][1].replace('X', 'D').replace('B', 'N').replace('Z', 'E').replace('U', 'C'))
			seq = Sequence(seq)
			for i in groups: #  para normalizar
				if seq.count(i) < min_comp:
					min_comp = seq.count(i)
				if seq.count(i) > max_comp:
					max_comp = seq.count(i)
			for i in groups:
				phage_comp['comp_' + i] += seq.count(i) - (min_comp / max_comp)
		total = 0
		if count_prots != 0:
			for i in groups:
				phage_comp['comp_' + i] = phage_comp['comp_' + i] / count_prots
				total += phage_comp['comp_' + i]
			for i in groups:
				phage_comp['comp_' + i] = phage_comp['comp_' + i] / total
		else:
			for i in groups:
				phage_comp['comp_' + i] = 0.0

		count_prots = 0
		for prot in bacteria:
			max_comp = 0.0
			min_comp = 1000000.0
			count_prots += 1
			seq = bacteria[prot][1]
			seq = seq[:seq.find('"')]
			seq = self.__get_conjoint_triad(seq.replace('X', 'D').replace('B', 'N').replace('Z', 'E').replace('U', 'C'))
			seq = Sequence(seq)
			for i in groups:
				if seq.count(i) < min_comp:
					min_comp = seq.count(i)
				if seq.count(i) > max_comp:
					max_comp = seq.count(i)
			for i in groups:
				bact_comp['comp_' + i] += seq.count(i) - (min_comp / max_comp)
		total = 0
		if count_prots != 0:
			for i in groups:
				bact_comp['comp_' + i] = bact_comp['comp_' + i] / count_prots
				total += bact_comp['comp_' + i]
			for i in groups:
				bact_comp['comp_' + i] = bact_comp['comp_' + i] / total
		else:
			for i in groups:
				bact_comp['comp_' + i] = 0.0

		for i in groups:
			solution.append(bact_comp['comp_' + i])
			solution.append(phage_comp['comp_' + i])
		return solution

	def add_grouping(self):
		from skbio import Sequence
		import json
		bact_group = {}
		phage_group = {}
		groups = '0123456'
		letters = 'ABCDEFGHIJ'
		for i in groups:
			for j in letters:
				bact_group['group' + j + '_' + i] = []
				phage_group['group' + j + '_' + i] = []
		phage = ''
		bact = ''
		count = -1
		for ID in self.features_data.index:
			done_phage = False
			done_bact = False
			count += 1
			if ID[:ID.find('--')] == phage:
				for i in groups:
					for j in letters:
						phage_group['group' + j + '_' + i].append(phage_group['group' + j + '_' + i][-1])
				done_phage = True
			if ID[ID.find('--') + 2:] == bact:
				for i in groups:
					for j in letters:
						bact_group['group' + j + '_' + i].append(bact_group['group' + j + '_' + i][-1])
				done_bact = True
			bact = ID[ID.find('--') + 2:]
			phage = ID[:ID.find('--')]

			if not done_phage:
				count_prots = 0
				for i in groups:
					for j in letters:
						phage_group['group' + j + '_' + i].append(0)
				for prot in self.list_prot[phage]:
					count_prots += 1
					seq = self.__get_conjoint_triad(self.list_prot[phage][prot][1].replace('X', 'D').replace('B', 'N').replace('Z', 'E').replace('U', 'C'))
					seq = Sequence(seq)
					for j in letters:
						group = self.__get_grouping(seq, j)
						for i in groups:
							phage_group['group' + j + '_' + i][count] += group[i]
				if count_prots != 0:
					for i in groups:
						for j in letters:
							phage_group['group' + j + '_' + i][count] = phage_group['group' + j + '_' + i][count] / count_prots
				else:
					for i in groups:
						for j in letters:
							phage_group['group' + j + '_' + i][count] = 0.0

			if not done_bact:
				count_prots = 0
				for i in groups:
					for j in letters:
						bact_group['group' + j + '_' + i].append(0)
				with open('files/bacteria/' + bact + '.json', encoding='utf-8') as F:
					bact_prots = json.loads(F.read())
				for prot in bact_prots:
					count_prots += 1
					seq = bact_prots[prot][1]
					seq = seq[:seq.find('"')]
					seq = self.__get_conjoint_triad(seq.replace('X', 'D').replace('B', 'N').replace('Z', 'E').replace('U', 'C'))
					seq = Sequence(seq)
					for j in letters:
						group = self.__get_grouping(seq, j)
						for i in groups:
							bact_group['group' + j + '_' + i][count] += group[i]
				if count_prots != 0:
					for i in groups:
						for j in letters:
							bact_group['group' + j + '_' + i][count] = bact_group['group' + j + '_' + i][count] / count_prots
				else:
					for i in groups:
						for j in letters:
							bact_group['group' + j + '_' + i][count] = 0.0

		for i in groups:
			for j in letters:
				self.features_data['bact_group' + j + '_' + i] = bact_group['group' + j + '_' + i]
				self.features_data['phage_group' + j + '_' + i] = phage_group['group' + j + '_' + i]

	def get_grouping(self, phage, bacteria):
		from skbio import Sequence
		bact_group = {}
		phage_group = {}
		groups = '0123456'
		letters = 'ABCDEFGHIJ'
		for i in groups:
			for j in letters:
				bact_group['group' + j + '_' + i] = 0
				phage_group['group' + j + '_' + i] = 0
		solution = []
		count_prots = 0
		for prot in phage:
			count_prots += 1
			seq = self.__get_conjoint_triad(phage[prot][1].replace('X', 'D').replace('B', 'N').replace('Z', 'E').replace('U', 'C'))
			seq = Sequence(seq)
			for j in letters:
				group = self.__get_grouping(seq, j)
				for i in groups:
					phage_group['group' + j + '_' + i] += group[i]
		if count_prots != 0:
			for i in groups:
				for j in letters:
					phage_group['group' + j + '_' + i] = phage_group['group' + j + '_' + i] / count_prots
		else:
			for i in groups:
				for j in letters:
					phage_group['group' + j + '_' + i] = 0.0

		count_prots = 0
		for prot in bacteria:
			count_prots += 1
			seq = bacteria[prot][1]
			seq = seq[:seq.find('"')]
			seq = self.__get_conjoint_triad(seq.replace('X', 'D').replace('B', 'N').replace('Z', 'E').replace('U', 'C'))
			seq = Sequence(seq)
			for j in letters:
				group = self.__get_grouping(seq, j)
				for i in groups:
					bact_group['group' + j + '_' + i] += group[i]
		if count_prots != 0:
			for i in groups:
				for j in letters:
					bact_group['group' + j + '_' + i] = bact_group['group' + j + '_' + i] / count_prots
		else:
			for i in groups:
				for j in letters:
					bact_group['group' + j + '_' + i] = 0.0

		for i in groups:
			for j in letters:
				solution.append(bact_group['group' + j + '_' + i])
				solution.append(phage_group['group' + j + '_' + i])
		return solution

	def __get_conjoint_triad(self, prot):
		ctm = {'A':'0', 'G':'0', 'V':'0','C':'1', 'F':'2', 'I':'2', 'L':'2', 'P':'2', 'M':'3', 'S':'3', 'T':'3', 'Y':'3', 'H':'4', 'N':'4', 'Q':'4', 'W':'4', 'K':'5', 'R':'5', 'D':'6', 'E':'6'}
		for i, j in ctm.items():
			prot = prot.replace(i, j)
		return prot

	def __get_grouping(self, prot, let='A'):
		from skbio import Sequence
		groups = '0123456'
		group = {}
		for i in groups:
			group[i] = 0.0
		if let == 'A':
			seq = Sequence(prot[:int(len(prot) * 0.25)])
			for i in groups:
				group[i] += seq.count(i) / len(seq)
		elif let == 'B':
			seq = Sequence(prot[int(len(prot) * 0.25):int(len(prot) * 0.5)])
			for i in groups:
				group[i] += seq.count(i) / len(seq)
		elif let == 'C':
			seq = Sequence(prot[int(len(prot) * 0.5):int(len(prot) * 0.75)])
			for i in groups:
				group[i] += seq.count(i) / len(seq)
		elif let == 'D':
			seq = Sequence(prot[int(len(prot) * 0.75):])
			for i in groups:
				group[i] += seq.count(i) / len(seq)
		elif let == 'E':
			seq = Sequence(prot[:int(len(prot) * 0.5)])
			for i in groups:
				group[i] += seq.count(i) / len(seq)
		elif let == 'F':
			seq = Sequence(prot[int(len(prot) * 0.5):])
			for i in groups:
				group[i] += seq.count(i) / len(seq)
		elif let == 'G':
			seq = Sequence(prot[int(len(prot) * 0.25):int(len(prot) * 0.75)])
			for i in groups:
				group[i] += seq.count(i) / len(seq)
		elif let == 'H':
			seq = Sequence(prot[:int(len(prot) * 0.75)])
			for i in groups:
				group[i] += seq.count(i) / len(seq)
		elif let == 'I':
			seq = Sequence(prot[int(len(prot) * 0.25):])
			for i in groups:
				group[i] += seq.count(i) / len(seq)
		elif let == 'J':
			seq = Sequence(prot[int(len(prot) * 0.125):int(len(prot) * 0.875)])
			for i in groups:
				group[i] += seq.count(i) / len(seq)
		return group

	def set_output(self):
		import pandas as pd
		output = []
		data = pd.read_csv('files/NCBI_Phage_Bacteria_Data.csv', header=0, index_col=0, names=['Phage Name', 'Bacteria Name', 'Bacteria ID'])
		for phage in self.features_data['ID']:
			phage = phage[:phage.find('--')]
			bact = data.loc[phage, 'Bacteria Name']
			if 'escherichia' in bact.lower():
				output.append('Escherichia coli')
			elif 'klebsiella' in bact.lower():
				output.append('Klebsiella pneumoniae')
			elif 'acinetobacter' in bact.lower():
				output.append('Acinetobacter baumannii')
		self.features_data = self.features_data.set_index('ID')
		self.features_data['Bacteria'] = output

	def save_feat_data(self):
		import pickle
		with open('files/FeatureDataset', 'wb') as f:
			pickle.dump(self.features_data, f)
		return self.features_data

	def import_feat_data(self):
		import pickle
		with open('files/FeatureDataset', 'rb') as f:
			self.features_data = pickle.load(f)
		return self.features_data


if __name__ == '__main__':
	test = FeatureConstruction()
	# test.process_net_surf()
	test.add_grouping()
	test.add_composition()
	test.add_kmers()
	# test.set_output()
	test.save_feat_data()
	'''
	test.process_net_surf()
	test.add_aa_freq()
	test.add_aromaticity()
	test.add_flexibility()
	test.add_molecular_weight()'''
	# test.import_feat_data()
	# test.netSurf()
