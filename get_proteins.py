
class PhageBacteriaInformation:

	def __init__(self, dataset=None):
		"""
		Imports a dataset from NCBI Virus, where the columns are Phage ID, Phage Name, Bacteria Name, Bacteria ID.
		If a phage entry does not have a bacteria associated, it is deleted
		:param dataset:
		"""
		import pandas as pd
		import ast
		self.phagesProteins = {}
		# self.phagesDNA = {}
		self.bactProteins = {}
		# self.bactDNA = {}
		self.data = pd.read_csv('files/'+dataset, header=0, index_col=0)
		self.data = self.data.dropna(how='any')
		self.data = self.data[self.data['Host_ID'] != '[]']
		index_remove = []
		for i in range(len(self.data)):
			temp = self.data['Host'][i].split(' ')
			if len(temp) <= 1:
				index_remove.append(i)
		self.data = self.data.drop(self.data.index[index_remove])
		self.all_bact = []
		for i in self.data.index:
			for bact in ast.literal_eval(self.data.loc[i, 'Host_ID']):
				if bact[:-2] not in self.all_bact:
					self.all_bact.append(bact[:-2])
		self.data.to_csv('files/Filtered_Phage_Bacteria.csv')

	def addFeatures(self):
		"""
		For each phage in the data, it saves its DNA sequence and all proteins, as provided by NCBI. It saves them into two variables
		Each bacteria associated with the phage is also searched for its DNA and proteins sequences.
		:return:
		"""
		from Bio import Entrez
		from Bio import SeqIO
		import json
		import ast
		Entrez.email = 'pedro_araujo97@hotmail.com'
		print('Working...')
		for phageID in self.data.index:
			with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=phageID) as handle:
				genomePhage = SeqIO.read(handle, "gb")
			protsPhage = {}
			for feat in genomePhage.features:
				if feat.type == 'CDS':
					try: protsPhage[feat.qualifiers['protein_id'][0]] = [feat.qualifiers['product'][0], feat.qualifiers['translation'][0]]
					except: pass
			self.phagesProteins[phageID] = protsPhage

		for bact in self.all_bact:
			protsBac = {}
			with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=bact) as handle:
				genomeBac = SeqIO.read(handle, "gb")
			for feat in genomeBac.features:
				if feat.type == 'CDS':
					try: protsBac[feat.qualifiers['protein_id'][0]] = [feat.qualifiers['product'][0], feat.qualifiers['translation'][0]]
					except: pass
			self.bactProteins[bact] = protsBac

		with open('files/phagesProteins.json', 'w') as f:
			json.dump(self.phagesProteins, f)
		self.__createFasta(self.phagesProteins, 'phagesProteins')
		with open('files/bactProteins.json', 'w') as f:
			json.dump(self.bactProteins, f)
		self.__createFasta(self.bactProteins, 'bactProteins')
		print('Done')

	def addBacteriaFeatures(self):
		"""
		For each unique bacteria present in the dataset, the DNA and protein sequences are saved in two variables.
		:return:
		"""
		from Bio import Entrez
		from Bio import SeqIO
		import json
		import ast
		Entrez.email = 'pedro_araujo97@hotmail.com'
		print('Working...')
		for bact in self.all_bact:
			if bact not in self.bactProteins.keys():
				protsBac = {}
				try:
					with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=bact) as handle:
						genomeBac = SeqIO.read(handle, "gb")
					for feat in genomeBac.features:
						if feat.type == 'CDS':
							try: protsBac[feat.qualifiers['protein_id'][0]] = [feat.qualifiers['product'][0], feat.qualifiers['translation'][0]]
							except: pass
					if len(genomeBac.features) <= 5:
						with Entrez.efetch(db="nucleotide", rettype="gbwithparts", retmode="text", id=bact) as handle:
							genomeBac = handle.readlines()
						for i in range(len(genomeBac)):
							if ' CDS ' in genomeBac[i]:
								j = i
								protDone = False
								while j < len(genomeBac):
									if protDone:
										break
									if '/product=' in genomeBac[j]:
										product = genomeBac[j].strip()[10:]
										j += 1
									elif '_id=' in genomeBac[j]:
										protKey = genomeBac[j].strip()[13:-1]
										j += 1
									elif '/translation=' in genomeBac[j]:
										protSeq = genomeBac[j].strip()[14:]
										j += 1
										for k in range(j, len(genomeBac)):
											if genomeBac[k].islower():
												j = k
												protDone = True
												break
											else:
												protSeq += genomeBac[k].strip()
									else:
										j += 1
								protsBac[protKey] = [product, protSeq[:protSeq.find('"')]]
					self.bactProteins[bact] = protsBac
				except:
					print(bact + ' failed')
		with open('files/bactProteins.json', 'w') as f:
			json.dump(self.bactProteins, f)
		self.__createFasta(self.bactProteins, 'bactProteins')
		print('Done')

	def add_individual_bacteria(self):
		"""
		For each unique bacteria present in the dataset, the DNA and protein sequences are saved in two variables.
		:return:
		"""
		from Bio import Entrez
		from Bio import SeqIO
		import json
		from pathlib import Path
		Entrez.email = 'pedro_araujo97@hotmail.com'
		print('Working...')
		for bact in self.all_bact:
			my_file = Path('files/bacteria/' + bact + ".json")
			if not my_file.is_file():
				protsBac = {}
				try:
					with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=bact) as handle:
						genomeBac = SeqIO.read(handle, "gb")
					for feat in genomeBac.features:
						if feat.type == 'CDS':
							try: protsBac[feat.qualifiers['protein_id'][0]] = [feat.qualifiers['product'][0], feat.qualifiers['translation'][0]]
							except: pass
					if len(genomeBac.features) <= 5:
						with Entrez.efetch(db="nucleotide", rettype="gbwithparts", retmode="text", id=bact) as handle:
							genomeBac = handle.readlines()
						for i in range(len(genomeBac)):
							if ' CDS ' in genomeBac[i]:
								j = i
								protDone = False
								while j < len(genomeBac):
									if protDone:
										break
									if '/product=' in genomeBac[j]:
										product = genomeBac[j].strip()[10:]
										j += 1
									elif '_id=' in genomeBac[j]:
										protKey = genomeBac[j].strip()[13:-1]
										j += 1
									elif '/translation=' in genomeBac[j]:
										protSeq = genomeBac[j].strip()[14:]
										j += 1
										for k in range(j, len(genomeBac)):
											if genomeBac[k].islower():
												j = k
												protDone = True
												break
											else:
												protSeq += genomeBac[k].strip()
									else:
										j += 1
								protsBac[protKey] = [product, protSeq[:protSeq.find('"')]]
					with open('files/bacteria/' + bact + '.json', 'w') as f:
						json.dump(protsBac, f)
				except:
					print(bact + ' failed')
		# with open('files/bactProteins.json', 'w') as f:
		# 	json.dump(self.bactProteins, f)
		# self.__createFasta(self.bactProteins, 'bactProteins')
		print('Done')

	def importData(self):
		"""
		Imports the previously saved DNA and protein sequences. This needs to be improved so the user can specify which data to import.
		:return:
		"""
		import json
		with open('files/phagesProteins.json', encoding='utf-8') as F:
			self.phagesProteins = json.loads(F.read())
		# with open('files/bactProteins.json', encoding='utf-8') as F:
		# 	self.bactProteins = json.loads(F.read())

	def addPhageProt(self):
		"""
		For each unique phage present in the dataset, the DNA and protein sequences are saved in two variables.
		:return:
		"""
		from Bio import Entrez
		from Bio import SeqIO
		import json
		Entrez.email = 'pedro_araujo97@hotmail.com'
		print('Working...')

		for phageID in self.data.index:
			with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=phageID) as handle:
				genomePhage = SeqIO.read(handle, "gb")
			protsPhage = {}
			for feat in genomePhage.features:
				if feat.type == 'CDS':
					try: protsPhage[feat.qualifiers['protein_id'][0]] = [feat.qualifiers['product'][0], feat.qualifiers['translation'][0]]
					except: pass
			self.phagesProteins[phageID] = protsPhage

		with open('files/phagesProteins.json', 'w') as f:
			json.dump(self.phagesProteins, f)
		self.__createFasta(self.phagesProteins, 'phagesProteins')
		return self.phagesProteins

	def add_missing_phage(self):
		"""
		For each unique phage present in the dataset, the DNA and protein sequences are saved in two variables.
		:return:
		"""
		from Bio import Entrez
		from Bio import SeqIO
		import json
		Entrez.email = 'pedro_araujo97@hotmail.com'
		print('Working...')

		for phageID in self.data.index:
			if phageID not in self.phagesProteins.keys():
				print(phageID)
				with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=phageID) as handle:
					genomePhage = SeqIO.read(handle, "gb")
				protsPhage = {}
				for feat in genomePhage.features:
					if feat.type == 'CDS':
						try: protsPhage[feat.qualifiers['protein_id'][0]] = [feat.qualifiers['product'][0], feat.qualifiers['translation'][0]]
						except: pass
				self.phagesProteins[phageID] = protsPhage
		with open('files/phagesProteins.json', 'w') as f:
			json.dump(self.phagesProteins, f)
		self.__createFasta(self.phagesProteins, 'phagesProteins')
		return self.phagesProteins

	def __createFasta(self, var, name):
		with open('files/' + name + '.fasta', 'w') as F:
			for spec in var:
				try:
					for prot in var[spec]:
						F.write('>' + prot + '-' + spec + '\n' + var[spec][prot][1] + '\n')
				except:
					F.write('>' + spec + '\n' + var[spec] + '\n')


if __name__ == '__main__':
	test = PhageBacteriaInformation('NCBI_Phage_Bacteria_Data.csv')
	test.add_individual_bacteria()
	test.addPhageProt()
	test.add_missing_phage()
	test.importData()
