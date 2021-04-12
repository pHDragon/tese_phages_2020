

class GalaxyPrediction:

	def __init__(self, phage_input_type='ID', bact_input_type='ID', phage='', bacteria='', ml_model='RandomForests', run_interpro=False):
		import pickle
		import os
		import re
		with open('files/feature_dataset', 'rb') as f:
			dataset = pickle.load(f)
		self.all_phages = []
		self.all_bacteria = []
		for ID in dataset.index:
			temp_phage = ID[:ID.find('--')]
			temp_bacteria = ID[ID.find('--')+2:]
			if temp_phage not in self.all_phages:
				self.all_phages.append(temp_phage)
			if temp_bacteria not in self.all_bacteria:
				self.all_bacteria.append(temp_bacteria)
		if phage_input_type == 'ID':
			phage = re.split('\W', phage.replace(' ', ''))
			len_phage_id = len(phage)
			phage_seqs = self._retrieve_from_phage_id(phage)
		elif phage_input_type == 'seq_file':
			phage_seqs = {}
			phage_seqs['PhageFasta'] = {}
			with open(phage, 'r') as f:
				temp = f.readlines()
			count_prot = 0
			prot = ''
			i=0
			while i < len(temp):
				if '>' in temp[i]:
					if prot:
						phage_seqs['PhageFasta']['Protein' + str(count_prot)] = ['Unknown', prot]
					count_prot += 1
					prot = ''
					i+=1
				else:
					prot += temp[i].strip()
					i+=1
		if bact_input_type == 'ID':
			bacteria = re.split('\W', bacteria.replace(' ', ''))
			if len(bacteria) > 1 and len_phage_id == 1 or len(bacteria) == 1:
				bact_seqs = self._retrieve_from_bact_id(bacteria)
		elif bact_input_type == 'seq_file':
			bact_seqs = {}
			bact_seqs['BacteriaFasta'] = {}
			with open(bacteria, 'r') as f:
				temp = f.readlines()
			count_prot = 0
			prot = ''
			i=0
			while i < len(temp):
				if '>' in temp[i]:
					if prot:
						bact_seqs['BacteriaFasta']['Protein' + str(count_prot)] = ['Unknown', prot]
					count_prot += 1
					prot = ''
					i+=1
				else:
					prot += temp[i].strip()
					i+=1
		phage_seqs = self._find_phage_functions(phage_seqs, run_interpro)
		phage_seqs = self._find_phage_tails(phage_seqs)
		
		list_remove = []
		for org in phage_seqs:
			if not phage_seqs[org]:
				print('Could not find tails for phage ' + org + '. Deleting entry...')
				list_remove.append(org)
		for org in list_remove:
			del phage_seqs[org]
		
		if phage_seqs:
			output = self.run_prediction(phage_seqs, bact_seqs, ml_model)
			self.create_output(output, phage_seqs, bact_seqs)
		else:
			with open(or_location + '/output.tsv', 'w') as f:
				f.write('No phage tails found in query')
		for file in os.listdir('files'):
			if file.startswith('temp'):
				os.remove('files/' + file)

	def _retrieve_from_phage_id(self, phage):
		temp_phage = {}
		for ID in phage:
			temp_phage[ID] = {}
			if ID in self.all_phages:
				import json
				with open('files/phageTails.json', encoding='utf-8') as f:
					phage_tails = json.loads(f.read())
				temp_phage[ID] = phage_tails[ID]
			else:
				from Bio import Entrez
				from Bio import SeqIO
				phage = {}
				Entrez.email = 'insert@email.com'
				try:
					with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ID) as handle:
						genome = SeqIO.read(handle, "gb")
					for feat in genome.features:
						if feat.type == 'CDS':
							try: temp_phage[ID][feat.qualifiers['protein_id'][0]] = [feat.qualifiers['product'][0], feat.qualifiers['translation'][0]]
							except: pass
				except:
					print(ID, 'not found in GenBank')
		return temp_phage

	def _retrieve_from_bact_id(self, bacteria):
		temp_bacteria = {}
		for ID in bacteria:
			temp_bacteria[ID] = {}
			if '.' in ID:
				ID = ID[:ID.find('.')]
			#if ID in self.all_bacteria:
			#	import json
			#	with open('files/bacteria/' + ID + '.json', encoding='utf-8') as f:
			#		temp_bacteria[ID] = json.loads(f.read())
			#else:
			from Bio import Entrez
			from Bio import SeqIO
			bacteria = {}
			Entrez.email = 'insert@email.com'
			try:
				with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ID+'.1') as handle:
					genome = SeqIO.read(handle, "gb")
				for feat in genome.features:
					if feat.type == 'CDS':
						try: temp_bacteria[ID][feat.qualifiers['protein_id'][0]] = [feat.qualifiers['product'][0], feat.qualifiers['translation'][0]]
						except: pass
				if len(genome.features) <= 5:
					with Entrez.efetch(db="nucleotide", rettype="gbwithparts", retmode="text", id=ID) as handle:
						genome = handle.readlines()
					for i in range(len(genome)):
						if ' CDS ' in genome[i]:
							j = i
							protDone = False
							while j < len(genome):
								if protDone:
									break
								if '/product=' in genome[j]:
									product = genome[j].strip()[10:]
									j += 1
								elif '_id=' in genome[j]:
									protKey = genome[j].strip()[13:-1]
									j += 1
								elif '/translation=' in genome[j]:
									protSeq = genome[j].strip()[14:]
									j += 1
									for k in range(j, len(genome)):
										if genome[k].islower():
											j = k
											protDone = True
											break
										else:
											protSeq += genome[k].strip()
								else:
									j += 1
							temp_bacteria[ID][protKey] = [product, protSeq[:protSeq.find('"')]]
			except:
				print(ID, 'not found in GenBank')
		return temp_bacteria

	def _find_phage_functions(self, phage_dict, run_interpro):
		import os
		import json
		with open('files/known_function.json', encoding='utf-8') as F:
			known_function = json.loads(F.read())
		with open('files/temp_database.fasta', 'w') as F:
			for phage in known_function:
				for prot in known_function[phage]:
					F.write('>' + phage + '-' + prot + '\n' + known_function[phage][prot][1] + '\n')
		os.system('makeblastdb -in files/temp_database.fasta -dbtype prot -title PhageProts -parse_seqids -out files/temp_database -logfile files/temp_log')
		for org in phage_dict:
			with open('files/temp.fasta', 'w') as F:
				for prot in phage_dict[org]:
					F.write('>' + prot + '\n' + phage_dict[org][prot][1] + '\n')
			os.system('blastp -db files/temp_database -query files/temp.fasta -out files/temp_blast -num_threads 2 -outfmt 6')
			phage_dict[org] = self.process_blast(phage_dict[org], known_function)
			if run_interpro: phage_dict[org] = self.interpro(phage_dict[org])
		return phage_dict

	def process_blast(self, phage_dict, known_function):
		import pandas as pd
		import re
		blast_domains = pd.read_csv('files/temp_blast', sep='\t', header=None)
		for prot in phage_dict:
			func = phage_dict[prot][0]
			known = False
			if (not any(i in func.lower() for i in ['hypothetical', 'unknown', 'kda', 'uncharacterized', 'hyphothetical']) and len(func) > 3) and not ('gp' in func.lower() and len(func.split(' ')) < 2) and not (len(func.split(' ')) == 1 and len(func) < 5):
				known = True
			if not known:
				evalue = []
				bitscore = []
				pred = blast_domains[blast_domains[0] == prot]
				if pred.shape[0] == 0: break
				for i in pred[10]:
					evalue.append(float(i))
				for i in pred[11]:
					bitscore.append(float(i))
				if min(evalue) < 1.0 and max(bitscore) > 30.0:
					ind = evalue.index(min(evalue))
					if ind != bitscore.index(max(bitscore)):
						ind = bitscore.index(max(bitscore))
					temp = pred.iloc[ind, 1]
					known_phage = temp[:temp.find('-')]
					known_prot = temp[temp.find('-') + 1:]
					if known_function[known_phage][known_prot]:
						new_func = known_function[known_phage][known_prot][0]
					# for j in known_function.keys():
					# 	if pred.iloc[ind, 1] in known_function[j].keys():
					# 		new_func = known_function[j][pred.iloc[ind, 1]][0]
					# 		break
					x = re.findall('(Gp\d{2,}[^,\d -]|Gp\d{1}[^,\d -])', temp)  # se tiver hits, remover
					if not any(z in new_func.lower() for z in ['unknown', 'ucp', 'uncharacterized', 'consensus']) and len(new_func) > 3 and not x:
						phage_dict[prot][0] = new_func
		return phage_dict

	def interpro(self, phage_dict):
		import os
		import pandas as pd
		import re
		os.system('interproscan.sh -b ' + 'files/temp_interpro -i ' + 'files/temp.fasta -f tsv > files/temp_interpro_log')
		domains = pd.read_csv('files/temp_interpro.tsv', sep='\t', index_col=0, header=None, names=list(range(13)))
		domains = domains.fillna('-')
		domains = domains[domains.loc[:, 3] != 'Coils']
		domains = domains[domains.loc[:, 3] != 'MobiDBLite']
		for prot in phage_dict:
			func = phage_dict[prot][0]
			known = False
			if (not any(i in func.lower() for i in ['hypothetical', 'unknown', 'kda', 'uncharacterized', 'hyphothetical']) and len(func) > 3) and not ('gp' in func.lower() and len(func.split(' ')) < 2) and not (len(func.split(' ')) == 1 and len(func) < 5):
				known = True
			if prot in domains.index and not known:
				temp = '-'
				try:
					for i in range(domains.loc[prot, :].shape[0]):
						if '-' not in domains.loc[prot, 12].iloc[i].lower():
							if float(domains.loc[prot, 8].iloc[i]) < 1.0:
								temp = domains.loc[prot, 12].iloc[i]
							break
				except:
					if float(domains.loc[prot, 8]) < 1.0:
						temp = domains.loc[prot, 12]
				x = re.findall('(Gp\d{2,}[^,\d -]|Gp\d{1}[^,\d -])', temp)  # se tiver hits, remover
				if temp != '-' and not any(z in temp.lower() for z in ['unknown', 'ucp', 'uncharacterized', 'consensus']) and len(temp) > 3 and not x:
					phage_dict[prot][0] = temp
				else:
					try:
						for i in range(domains.loc[prot, :].shape[0]):
							if '-' not in domains.loc[prot, 5].iloc[i].lower():
								temp = domains.loc[prot, 5].iloc[i]
								break
					except:
						temp = domains.loc[prot, 5]
					x = re.findall('(Gp\d{2,}[^,\d -]|Gp\d{1}[^,\d -])', temp)
					if temp != '-' and not any(z in temp.lower() for z in ['unknown', 'ucp', 'uncharacterized', 'consensus']) and len(temp) > 3 and not x:
						phage_dict[prot][0] = temp
		return phage_dict

	def _find_phage_tails(self, phage_dict):
		for org in phage_dict:
			list_remove = []
			for protein in phage_dict[org]:
				if any(z in phage_dict[org][protein][0].lower() for z in ['fiber', 'fibre', 'spike', 'hydrolase', 'bind', 'depolymerase', 'peptidase', 'lyase', 'sialidase', 'dextranase', 'lipase', 'adhesin', 'baseplate', 'protein h', 'recognizing', 'protein j', 'protein g', 'gpe', 'duf4035', 'host specifity', 'cor protein', 'specificity', 'baseplate component', 'gp38', 'gp12 tail', 'receptor', 'recognition', 'tail']) \
						and not any(z in phage_dict[org][protein][0].lower() for z in ['nucle', 'dna', 'rna', 'ligase', 'transferase', 'inhibitor', 'assembly', 'connect', 'nudix', 'atp', 'nad', 'transpos', 'ntp', 'molybdenum', 'hns', 'gtp', 'riib', 'inhibitor', 'replicat', 'codon', 'pyruvate', 'catalyst', 'hinge', 'sheath completion', 'head', 'capsid', 'tape', 'tip', 'strand', 'matur', 'portal', 'terminase', 'nucl', 'promot', 'block', 'olfact', 'wedge', 'lysozyme', 'mur', 'sheat']):
					pass
				else:
					list_remove.append(protein)
			for protein in list_remove:
				del phage_dict[org][protein]
		return phage_dict

	def run_prediction(self, phage_dict, bact_dict, ml_model):
		from feature_construction import FeatureConstruction
		import pickle
		from sklearn.preprocessing import LabelEncoder
		from sklearn.preprocessing import StandardScaler
		import numpy as np

		if ml_model == 'RandomForests':
			with open('files/dataset_reduced', 'rb') as f:
				dataset = pickle.load(f)
			columns_remove = [3, 7, 9, 11, 24, 28, 32, 34, 38, 42, 45, 52, 53, 61, 65, 73, 75, 79, 104, 122, 141, 151, 154, 155, 157, 159, 160, 161, 163, 165, 169, 170, 173, 176, 178, 180, 182, 183, 185, 186, 187, 190, 193, 194, 196, 197, 201, 202, 203, 206, 207, 209, 210, 212, 216, 217, 221, 223, 225, 226, 230, 233, 235, 236, 245, 251]
		elif ml_model == 'SVM':
			with open('files/feature_dataset', 'rb') as f:
				dataset = pickle.load(f)
			columns_remove = []

		dataset = dataset.dropna()
		le = LabelEncoder()
		le.fit(['Yes', 'No'])
		output = le.transform(dataset['Infects'].values)
		dataset = dataset.drop('Infects', 1)
		scaler = StandardScaler()
		scaler.fit(dataset)
		data_z = scaler.transform(dataset)

		fc = FeatureConstruction()
		solution = []
		for phage in phage_dict:
			for bacteria in bact_dict:
				temp_solution = np.array([])
				temp_solution = np.append(temp_solution, fc.get_grouping(phage_dict[phage], bact_dict[bacteria]))
				temp_solution = np.append(temp_solution, fc.get_composition(phage_dict[phage], bact_dict[bacteria]))
				temp_solution = np.append(temp_solution, fc.get_kmers(phage_dict[phage], bact_dict[bacteria]))
				temp_solution = temp_solution.reshape(1, -1)
				if columns_remove:
					temp_solution = np.delete(temp_solution, columns_remove, 1)
				if phage in self.all_phages:
					for ID in dataset.index:
						if phage in ID:
							for i in range(len(dataset.loc[ID].index)):
								if 'phage' in dataset.loc[ID].index[i]:
									temp_solution[0][i] = dataset.loc[ID, dataset.loc[ID].index[i]]
							break
				if bacteria in self.all_bacteria:
					for ID in dataset.index:
						if bacteria in ID:
							for i in range(len(dataset.loc[ID].index)):
								if 'bact' in dataset.loc[ID].index[i]:
									temp_solution[0][i] = dataset.loc[ID, dataset.loc[ID].index[i]]
							break
				if type(solution) != np.ndarray:
					solution = temp_solution
				else:
					solution = np.append(solution, temp_solution, axis=0)
		# solution = solution.reshape(1, -1)
		solution = scaler.transform(solution)

		if ml_model == 'RandomForests':
			from sklearn.ensemble import RandomForestClassifier
			clf = RandomForestClassifier(n_estimators=200, bootstrap=False, criterion='gini', min_samples_leaf=2, min_samples_split=4, oob_score=False)
			clf = clf.fit(data_z, output)
		elif ml_model == 'SVM':
			from sklearn.svm import SVC
			clf = SVC(C=10, degree=2, gamma='auto', kernel='rbf')
			clf = clf.fit(data_z, output)
		pred = clf.predict(solution)
		pred = list(le.inverse_transform(pred))
		return pred

	def create_output(self, output, phage_seqs, bact_seqs):
		import pandas as pd
		list_orgs = []
		for phage in phage_seqs:
			for bact in bact_seqs:
				list_orgs.append(phage + ' - ' + bact)
		file = pd.DataFrame({'Phage - Bacteria': list_orgs, 'Infects': output})
		file.to_csv('files/output.tsv', sep='\t', index=False, header=True)
		file.to_csv(or_location + '/output.tsv', sep='\t', index=False, header=True)


if __name__ == '__main__':
	import sys
	import os
	global or_location
	or_location = os.getcwd()
	os.chdir(os.path.dirname(__file__))

	phage_input_type = sys.argv[1]
	Phages = sys.argv[2]
	bact_input_type = sys.argv[3]
	Bacts = sys.argv[4]
	run_interpro = sys.argv[5]
	if run_interpro == 'True':
		run_interpro = True
	else:
		run_interpro = False
	model = sys.argv[6]
	GalaxyPrediction(phage_input_type=phage_input_type, bact_input_type=bact_input_type, phage=Phages, bacteria=Bacts, ml_model=model, run_interpro=run_interpro)
	#rg = GalaxyPrediction(phage_input_type='ID', bact_input_type='ID', phage='NC_050154', bacteria='NC_007414,NZ_MK033499,NZ_CP031214')
	# GalaxyPrediction(phage_input_type='ID', bact_input_type='ID', phage='NC_031087,NC_049833,NC_049838,NC_049444', bacteria='LR133964', ml_model='SVM')
