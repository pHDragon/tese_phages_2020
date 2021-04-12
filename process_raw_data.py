"""
https://www.ncbi.nlm.nih.gov/Sequin/acc.html
https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly
"""


class PhageBacteriaData:

	def __init__(self, dataset=None):
		"""
		Imports a dataset from NCBI Virus, where the columns are Phage ID, Phage Name, Bacteria Name, Bacteria ID.
		If a phage entry does not have a bacteria associated, it is deleted
		:param dataset:
		"""
		import pandas as pd
		import json
		self.listBacID = []
		if dataset is None:
			file = False
			while not file:
				try:
					name = input('File name: ')
					self.data = pd.read_csv('files/' + name, header=0, index_col=0)
					file = True
				except:
					print('Couldn\'t find file')
		else:
			self.data = pd.read_csv('files/' + dataset, header=0, index_col=0)
			self.data = self.data.dropna(how='any')
			self.data = self.data[self.data['Host'] != 'unclassified bacterium']
			index_remove = []
			for i in range(len(self.data)):
				if 'uncultured' in self.data.iloc[i]['Species']:
					index_remove.append(i)
				elif 'virus' not in self.data.iloc[i]['Species'] and 'phage' not in self.data.iloc[i]['Species']:
					index_remove.append(i)
			self.data = self.data.drop(self.data.index[index_remove])
			index_remove = []
			for i in range(len(self.data)):
				temp = self.data['Host'][i].split(' ')
				if len(temp) <= 1:
					index_remove.append(i)
			self.data = self.data.drop(self.data.index[index_remove])
		if 'Host_ID' not in self.data.columns:
			temp = []
			for i in self.data.index:
				temp.append([])
			self.data['Host_ID'] = temp
		try:
			with open('files/searched_accessions', 'r') as f:
				self.searched = json.loads(f.read())
		except:
			self.searched = {}

	def addPhageName(self):
		"""
		Using the entrez service, from NCBI, for each phage the name is added, from its features.
		:return:
		"""
		from Bio import Entrez
		from Bio import SeqIO
		Entrez.email = "pedro_araujo97@hotmail.com"
		listNames = []
		with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=self.data.index) as handle:
			for seq_record in SeqIO.parse(handle, "gb"):
				listNames.append(seq_record.annotations['organism'])
		self.data['Species'] = listNames

	def addBacteriaName(self):
		"""
		Using the entrez service, from NCBI, for each phage the infecting bacteria name is added, from its features.
		:return:
		"""
		from Bio import Entrez
		from Bio import SeqIO
		Entrez.email = "pedro_araujo97@hotmail.com"
		for phage in self.data.index:
			with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=phage) as handle:
				seq_record = SeqIO.read(handle, "gb")
				try:
					if len(seq_record.features[0].qualifiers['host'][0].split(' ')) > 2:
						self.data.loc[phage, 'Host'] = seq_record.features[0].qualifiers['host'][0]
				except:
					if len(seq_record.features[0].qualifiers['lab_host'][0].split(' ')) > 2:
						self.data.loc[phage, 'Host'] = seq_record.features[0].qualifiers['lab_host'][0]
		self.save_data()

	def addBacteriaGenome(self):
		"""
		For each phage, the associated scientific articles in pubmed are looked. From theses pubmed articles, all associated IDs are extracted.
		If the ID corresponds to a bacterial strain, its accession ID is added to the Bacteria ID column.
		:return:
		"""
		from Bio import Entrez
		import ast
		import pickle
		from pathlib import Path
		Entrez.email = "pedro_araujo97@hotmail.com"
		my_file = Path("files/searched_hosts")
		if my_file.is_file():
			with open('files/searched_hosts', 'rb') as f:
				list_done = pickle.load(f)
		else:
			list_done = []
		count = 0
		try:
			for phageID in self.data.index:
				if phageID in list_done: continue
				listBactID = ast.literal_eval(self.data.loc[phageID, 'Host_ID'])
				try:
					with Entrez.elink(dbfrom='nuccore', db='pubmed', id=phageID) as handle:
						pubmed = Entrez.read(handle)
					for i in pubmed[0]['LinkSetDb']:
						if 'weighted' not in i["LinkName"]:
							for link in i["Link"]:
								try:
									with Entrez.elink(dbfrom='pubmed', db="nucleotide", id=link['Id']) as handle:
										genomes = Entrez.read(handle)
									for i in genomes[0]['LinkSetDb']:
										if 'weighted' not in i['LinkName']:
											for id in i['Link']:
												with Entrez.esummary(db='nucleotide', id=id['Id']) as handle:
													bacorg = Entrez.read(handle)
													if bacorg[0]['Caption'] != phageID and 'phage' not in bacorg[0][
														'Title'].lower() \
															and bacorg[0][
														'AccessionVersion'] not in listBactID and 'cds' not in \
															bacorg[0]['Title'].lower() and 'shotgun' not in bacorg[0][
														'Title'].lower():
														if any(z in bacorg[0]['AccessionVersion'][:3] for z in
															   ['NC_', 'AC_', 'NZ_', 'CP', 'AE', 'CY', 'AP']):
															listBactID.append(bacorg[0]['AccessionVersion'])
															self.searched[bacorg[0]['AccessionVersion']] = 'yes'
															count += 1
														elif not any(z in bacorg[0]['AccessionVersion'][:3] for z in
																	 ['MN', 'FM', 'MQ', 'MR', 'MK', 'AB', 'MF', 'KP',
																	  'NM_', 'KC', 'MH', 'AY', 'FN', 'AY']) \
																and 'complete' in bacorg[0]['Title'].lower():
															if bacorg[0]['AccessionVersion'] in self.searched.keys():
																add = self.searched[bacorg[0]['AccessionVersion']]
															else:
																add = input('Check ' + bacorg[0][
																	'AccessionVersion'] + '\nDo you wish to add it? (yes/no)')
															if 'y' in add.lower():
																listBactID.append(bacorg[0]['AccessionVersion'])
																self.searched[bacorg[0]['AccessionVersion']] = 'yes'
																count += 1
															else:
																self.searched[bacorg[0]['AccessionVersion']] = 'no'
								except:
									pass
				except:
					pass
				self.data.loc[phageID, 'Host_ID'] = listBactID
				list_done.append(phageID)
				with open('files/searched_hosts', 'wb') as f:
					pickle.dump(list_done, f)
				self.save_data()
				print(phageID)
		except:
			print('Bacterial host name missing. Searching from phage id')
			pass
		print('For future reference,', count, "new bacterial ID's were added.")
		self.save_data()

	def checkAbstracts(self):
		"""
		For each phage, the associated scientific articles in pubmed are looked. From theses pubmed articles, the abstracted is searched for mentions of bacterial strains.
		If bacterial strains are found, its accession IDs are added to the Bacteria ID column.
		:return:
		"""
		from Bio import Entrez
		import re
		import ast
		Entrez.email = 'pedro_araujo97@hotmail.com'
		count = 0
		for phageID in self.data.index:
			if len(self.data.loc[phageID, 'Host'].split()) < 3:
				with Entrez.elink(dbfrom='nuccore', db='pubmed', id=phageID) as handle:
					pubmed = Entrez.read(handle)
				listBactID = ast.literal_eval(self.data.loc[phageID, 'Host_ID'])
				for i in pubmed[0]['LinkSetDb']:
					if 'weighted' not in i["LinkName"]:
						for link in i["Link"]:
							try:
								with Entrez.efetch(db="pubmed", rettype="medline", retmode="xml",
												   id=link['Id']) as handle:
									article = Entrez.read(handle)
								abstract = \
								article['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
								x = re.findall('\w{0,1}[A-Z]{1,5}[0-9]{1,5}[-,:]{0,1}[A-Z]{0,5}[1-9]{0,5}', abstract)
								for i in range(len(x)):
									x[i] = x[i].strip(',;')
								x = list(set(x))
								for i in x:
									if 'ORF' in i:
										x.remove(i)
								for strain in x:
									with Entrez.esearch(db='nucleotide', term='((((' + self.data.loc[
										phageID, 'Host'] + ' ' + strain + ' AND complete sequence) NOT shotgun[Title]) NOT phage[Title]) NOT cds[Title]) NOT gene[Title]',
														idtype="acc") as handle:
										species = Entrez.read(handle)
										strains = species['IdList']
										for i in strains:
											if any(z in i for z in ['NC_', 'NZ_', 'AC_', 'CP', 'AE', 'CY',
																	'AP']) and i not in listBactID:
												listBactID.append(i)
												self.searched[i] = 'yes'
												count += 1
											elif not any(z in i[:3] for z in
														 ['MN', 'FM', 'MQ', 'MR', 'MK', 'AB', 'MF', 'KP', 'NM_', 'KC',
														  'MH', 'AY', 'FN', 'AY']) and i not in listBactID:
												if i in self.searched.keys():
													add = self.searched[i]
												else:
													add = input('Check ' + i + '\nDo you wish to add it? (yes/no)')
												if 'y' in add.lower():
													listBactID.append(i)
													self.searched[i] = 'yes'
													count += 1
												else:
													self.searched = 'no'
							except:
								pass
				self.data.loc[phageID, 'Host_ID'] = listBactID
		print('For future reference,', count, "new bacterial ID's were added.")
		self.save_data()

	def saveAbstracts(self):
		"""
		For each phage, the associated scientific articles in pubmed are looked.
		From theses pubmed articles, the abstracted is saved in a dictionary structure, where each phage ID is associated with a list of abstracts.
		These abstracts can be used for later processing.
		:return:
		"""
		from Bio import Entrez
		import json
		Entrez.email = 'pedro_araujo97@hotmail.com'
		dicPhageAbstracts = {}

		for phageID in self.data.index:
			print(phageID, end='\n')
			dicPhageAbstracts[phageID] = []
			with Entrez.elink(dbfrom='nuccore', db='pubmed', id=phageID) as handle:
				pubmed = Entrez.read(handle)
			lista = []
			for i in pubmed[0]['LinkSetDb']:
				if 'weighted' not in i["LinkName"]:
					for link in i["Link"]:
						try:
							with Entrez.efetch(db="pubmed", rettype="medline", retmode="xml", id=link['Id']) as handle:
								article = Entrez.read(handle)
							abstract = \
							article['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
							dicPhageAbstracts[phageID].append([link['Id'], abstract])
						except:
							pass
		with open('files/phageAbstracts.json', 'w') as f:
			json.dump(dicPhageAbstracts, f)

	def searchBacName(self):
		from Bio import Entrez
		import ast
		Entrez.email = "pedro_araujo97@hotmail.com"
		count = 0
		for phageID in self.data.index:
			if len(self.data.loc[phageID, 'Host'].split()) > 2:
				listBactID = ast.literal_eval(self.data.loc[phageID, 'Host_ID'])
				with Entrez.esearch(db='nucleotide', term='((((' + self.data.loc[
					phageID, 'Host'] + ' AND complete sequence) NOT shotgun[Title]) NOT phage[Title]) NOT cds[Title]) NOT gene[Title]',
									idtype="acc") as handle:
					species = Entrez.read(handle)
					strains = species['IdList']
					for j in strains:
						if any(z in j[:3] for z in
							   ['NC_', 'NZ_', 'AC_', 'CP', 'AE', 'CY', 'AP']) and j not in listBactID:
							listBactID.append(j)
							self.searched[j] = 'yes'
							count += 1
						elif not any(z in j[:3] for z in
									 ['MN', 'FM', 'MQ', 'MR', 'MK', 'AB', 'MF', 'KP', 'NM_', 'KC', 'MH', 'AY', 'FN',
									  'AY']) and j not in listBactID:
							if j in self.searched.keys():
								add = self.searched[j]
							else:
								add = input('Check ' + j + '\nDo you wish to add it? (yes/no)')
							if 'y' in add.lower():
								listBactID.append(j)
								self.searched[j] = 'yes'
								count += 1
							else:
								self.searched[j] = 'no'
				self.data.loc[phageID, 'Host_ID'] = listBactID
		print('For future reference,', count, "new bacterial ID's were added.")
		self.save_data()

	def createListBacID(self, lower=0, upper=100):
		"""
		More sequential than previous methods. Maybe include every single one...
		:param lower: lower index from the phage list (numeric)
		:param upper: upper index from the phage list (numeric)
		:return:
		"""
		from Bio import Entrez
		Entrez.email = 'pedro_araujo97@hotmail.com'
		for i in range(lower, upper):
			phageID = self.data.index[i]
			BactID = []
			name = test.data.loc[phageID]['Bacteria Name']
			try:
				if name != 'unclassified bacterium' and not name != name:  # Verificação de hosts válidos
					with Entrez.elink(dbfrom='nuccore', db='pubmed', id=phageID) as handle:
						pubmed = Entrez.read(handle)
					for link in pubmed[0]["LinkSetDb"][0]["Link"]:
						try:
							with Entrez.elink(dbfrom='pubmed', db="nucleotide", id=link['Id']) as handle:
								genomes = Entrez.read(handle)
							for id in genomes[0]['LinkSetDb'][0]['Link']:
								with Entrez.esummary(db='nucleotide', id=id['Id']) as handle:
									bacorg = Entrez.read(handle)
									if 'NC_' in bacorg[0]['AccessionVersion'] or 'NZ_' in bacorg[0]['AccessionVersion']:
										if bacorg[0]['Caption'] != phageID:
											BactID.append(bacorg[0]['AccessionVersion'])
						except:
							pass
				else:
					pass
			except:
				pass
			self.listBacID.append(BactID)

	def check_bacteria(self):
		from Bio import Entrez
		from Bio import SeqIO
		import ast
		Entrez.email = "pedro_araujo97@hotmail.com"
		all_bact = []
		for i in self.data.index:
			for bact in ast.literal_eval(self.data.loc[i, 'Host_ID']):
				if bact[:-2] not in all_bact:
					all_bact.append(bact[:-2])
		list_remove = []
		for bact in all_bact:
			if bact not in list_remove:
				try:
					with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=bact) as handle:
						seq_record = SeqIO.read(handle, "gb")
					if not any(i in seq_record.description.lower() for i in ['pneumoniae', 'coli', 'baumannii']) or not any(i in seq_record.description.lower() for i in ['escherichia', 'acinetobacter', 'klebsiella']) \
							or 'phage' in seq_record.description.lower() or 'virus' in seq_record.description.lower():
						list_remove.append(bact)
				except:
					list_remove.append(bact)
		print(list_remove)
		for phage in self.data.index:
			listBactID = ast.literal_eval(self.data.loc[phage, 'Host_ID'])
			for bact in listBactID:
				if bact[:-2] in list_remove:
					listBactID.remove(bact)
			self.data.loc[phage, 'Host_ID'] = listBactID
		self.data.to_csv('files/NCBI_Phage_Bacteria_Data.csv')

	def save_data(self):
		"""
		Saves the data in csv format.
		:return:
		"""
		import json
		self.data.to_csv('files/NCBI_Phage_Bacteria_Data.csv')
		with open('files/searched_accessions', 'w') as f:
			json.dump(self.searched, f)


if __name__ == '__main__':
	test = PhageBacteriaData('NCBI_Phage_Bacteria_Data.csv')  # sequences
	test.addBacteriaName()
	test.addBacteriaGenome()
	test.searchBacName()  # 2266 bacteria added
	test.checkAbstracts()
	test.searchBacName()
	# test.createListBacID(0, 100)
	# test.data = test.data.iloc[:, 0:3]
	test.check_bacteria()
	test.save_data()
# test.extractProtein()
# test.importProtein('Phage')
