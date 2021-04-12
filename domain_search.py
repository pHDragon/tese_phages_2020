
class DomainSearch:

	def __init__(self):
		'''
		This still needs a bit of modifications
		:param phagesProteins: protein function and sequences, as provided in NCBI. Each phage ID has every protein represented with a dicionary with keys as protein IDs
		:param phageDomains: for each phage and each of it's proteins, a list of predicted domains is given. If unavailable, it returns an empty list
		'''
		import json
		import pandas as pd
		with open('files/phagesProteins.json', encoding='utf-8') as F:
			self.phagesProteins = json.loads(F.read())
		# with open('files/bactProteins.json', encoding='utf-8') as F: # For later use, implement the same way as phage, more or less. Include psort
		# 	self.bacProt = json.loads(F.read())
		data = pd.read_csv('files/NCBI_Phage_Bacteria_Data.csv', header=0, index_col=0)
		for phage in data.index:
			if data.loc[phage, 'Host_ID'] == '[]':
				try: del self.phagesProteins[phage]
				except: pass
		self._filter_phage()

	def _create_fasta(self, dic, name):
		'''
		Creates a fasta file containing every protein sequence for a given dictionary.
		:return:
		'''
		with open('files/' + name, 'w') as F:
			for org in dic:
				for prot in dic[org]:
					F.write('>' + org + '-' + prot + '\n' + dic[org][prot][1] + '\n')

	def _filter_phage(self):
		self.known_function = {}
		self.unknown_function = {}
		for phage in self.phagesProteins.keys():
			self.known_function[phage] = {}
			self.unknown_function[phage] = {}
			for prot in self.phagesProteins[phage].keys():
				func = self.phagesProteins[phage][prot][0]
				if (not any(i in func.lower() for i in ['hypothetical', 'unknown', 'kda', 'uncharacterized', 'hyphothetical']) and len(func) > 3) and not ('gp' in func.lower() and len(func.split(' ')) < 2) and not (len(func.split(' ')) == 1 and len(func) < 5):
					self.known_function[phage][prot] = self.phagesProteins[phage][prot]
				else:
					self.unknown_function[phage][prot] = self.phagesProteins[phage][prot]

	def scanInterPro(self, InterPro_path='/home/pedro-linux/Downloads/interproscan-5.46-81.0/', out_path='/home/pedro-linux/OneDrive/UMinho/Cenas_de_tese_idk/test_tese_process/files/'):
		'''
		Creates a fasta file containing every protein and scans it using Interproscan. Creates a tsv file
		:param InterPro_path: path to the interproscan executable
		:param out_path: path to save the tsv output
		:return: domains_output.tsv, a file that contains the domain associated with each protein
		'''
		import os
		self._create_fasta(self.unknown_function, 'unknown_phages.fasta')
		os.system(InterPro_path + 'interproscan.sh -b ' + out_path + '/interpro/domains_output -i ' + out_path + 'unknown_phages.fasta -f tsv')

	def iter_interpro(self, InterPro_path='/home/pedro-linux/Downloads/interproscan-5.46-81.0/', out_path='/home/pedro-linux/OneDrive/UMinho/Cenas_de_tese_idk/test_tese_process/files/interpro/'):
		import os
		from pathlib import Path
		count = 0
		F = open('files/interpro/temp_100.fasta', 'w')
		for phage in self.unknown_function:
			for prot in self.unknown_function[phage]:
				count += 1
				my_file = Path("files/interpro/domains_output" + str(count) + ".tsv")
				if count % 100 == 0 and not my_file.is_file():
					F.write('>' + prot + '\n' + self.unknown_function[phage][prot][1] + '\n')
					F.close()
					os.system(InterPro_path + 'interproscan.sh -b ' + out_path + 'domains_output' + str(count) + ' -i ' + out_path + 'temp_100.fasta -f tsv')
					F = open('files/interpro/temp_100.fasta', 'w')
				else:
					F.write('>' + prot + '\n' + self.unknown_function[phage][prot][1] + '\n')
		if count % 100 != 0:
			F.close()
			os.system(InterPro_path + 'interproscan.sh -b ' + out_path + 'domains_output' + str(count) + ' -i ' + out_path + 'temp_100.fasta -f tsv')

	def processInterPro(self):
		'''
		Processes the tsv file created from scanInterPro. Domains are saved in the protdomains variable.
		:return: phageDomains, a dictionary that, for each protein in a given species, has domains associated
		'''
		import os
		from pathlib import Path
		import pandas as pd
		import re
		my_file = Path("files/interpro/domains_output.tsv")
		if not my_file.is_file():
			with open('files/interpro/domains_output.tsv', 'w') as F:
				for file in os.listdir('files/interpro/'):
					if 'temp_100' not in file:
						with open('files/interpro/' + file, 'r') as f:
							F.write(f.read())
		domains = pd.read_csv('files/interpro/domains_output.tsv', sep='\t', index_col=0, header=None, names=list(range(13)))
		domains = domains.fillna('-')
		domains = domains[domains.loc[:, 3] != 'Coils']
		domains = domains[domains.loc[:, 3] != 'MobiDBLite']
		# domains = domains.groupby(domains.index).last()
		add_domains = {}
		for spec in self.phagesProteins:
			for prot in self.phagesProteins[spec]:
				if prot in domains.index:
					temp = '-'
					try:
						for i in range(domains.loc[prot, :].shape[0]):
							if '-' not in domains.loc[prot, 12].iloc[i].lower():
								if float(domains.loc[id, 8].iloc[i]) < 1.0:
									temp = domains.loc[id, 12].iloc[i]
								break
					except:
						if float(domains.loc[id, 8]) < 1.0:
							temp = domains.loc[id, 12]
					x = re.findall('(Gp\d{2,}[^,\d -]|Gp\d{1}[^,\d -])', temp) # se tiver hits, remover
					if temp != '-' and not any(z in temp.lower() for z in ['unknown', 'ucp', 'uncharacterized', 'consensus']) and len(temp) > 3 and not x:
						if temp not in add_domains.keys():
							add_domains[temp] = input('Add function: ' + temp).lower()
						if 'y' in add_domains[temp]:
							self.phagesProteins[spec][prot][0] = temp
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
							if temp not in add_domains.keys():
								add_domains[temp] = input('Add function: ' + temp).lower()
							if 'y' in add_domains[temp]:
								self.phagesProteins[spec][prot][0] = temp

	def find_domains_interpro(self, dic):
		import os
		import pandas as pd
		import re
		InterPro_path='/home/pedro-linux/Downloads/interproscan-5.46-81.0/'
		out_path='/home/pedro-linux/OneDrive/UMinho/Cenas_de_tese_idk/WholeProcess/files/'
		with open('files/SinglePhageProteins.fasta', 'w') as F:
			for prot in dic.keys():
				F.write('>' + dic[prot][0] + '\n' + dic[prot][1] + '\n')
		os.system(InterPro_path + 'interproscan.sh -b ' + out_path + 'single_phage_domains -i ' + out_path + 'SinglePhageProteins.fasta -f tsv')

		domains = pd.read_csv('files/single_phage_domains.tsv', sep='\t', index_col=0, header=None, names=list(range(13)))
		domains = domains.fillna('-')
		for prot in dic:
			if prot in domains.index:
				temp = '-'
				try:
					for i in range(domains.loc[prot, :].shape[0]):
						if 'coil' not in domains.loc[prot, 12].iloc[i].lower() and '-' not in domains.loc[prot, 12].iloc[i].lower():
							temp = domains.loc[prot, 12].iloc[i]
							break
				except:
					temp = domains.loc[prot, 12]
				x = re.findall('(Gp\d{2,}[^,\d -]|Gp\d{1}[^,\d -])', temp) # se tiver hits, remover
				if temp != '-' and 'unknown' not in temp and 'UCP' not in temp and len(temp)>3 and not x:
					dic[prot][0] = temp
				else:
					try:
						for i in range(domains.loc[prot, :].shape[0]):
							if 'coil' not in domains.loc[prot, 5].iloc[i].lower() and '-' not in domains.loc[prot, 12].iloc[i].lower():
								temp = domains.loc[prot, 5].iloc[i]
								break
					except:
						temp = domains.loc[prot, 5]
					x = re.findall('(Gp\d{2,}[^,\d -]|Gp\d{1}[^,\d -])', temp)
					if temp != '-' and 'unknown' not in temp and 'UCP' not in temp and len(temp) > 3 and not x:
						dic[prot][0] = temp
		return dic

	def fillDomainsBLAST(self):
		'''
		Using the NCBIWWW package, it searches for domains with BLAST. Domains are saved in the protdomains variable.
		:return: phageDomains, a dictionary that, for each protein in a given species, has domains associated
		'''
		print('Finding functions/domains with BLAST')
		from Bio.Blast import NCBIWWW
		from Bio.Blast import NCBIXML
		import pickle
		from pathlib import Path
		my_file = Path("files/phage_list_blast")
		if my_file.is_file():
			with open('files/phage_list_blast', 'rb') as f:
				list_done = pickle.load(f)
		else:
			list_done = []
		for spec in self.phagesProteins:
			if spec not in list_done:
				for prot in self.phagesProteins[spec]:
					if 'hypothetical' in self.phagesProteins[spec][prot][0].lower() or 'uncharacterized' in self.phagesProteins[spec][prot][0].lower() or 'unknown' in self.phagesProteins[spec][prot][0].lower():
					# if not self.phageDomains[bac][prot]:
						result_handle = NCBIWWW.qblast('blastp', 'nr', self.phagesProteins[spec][prot][1], entrez_query='Acinetobacter baumannii (taxid:470), Escherichia coli (taxid:562), Klebsiella pneumonia (taxid:573)')
						blastout = NCBIXML.read(result_handle)
						for ali in blastout.alignments:
							if 'hypothetical' not in ali.hit_def.lower() and 'uncharacterized' not in ali.hit_def.lower():
								print(ali.hit_def[:ali.hit_def.find(' [')])
								self.phagesProteins[spec][prot][0] = ali.hit_def[:ali.hit_def.find(' [')]
								break
				list_done.append(spec)
				with open('files/phage_list_blast', 'wb') as f:
					pickle.dump(list_done, f)
				self.saveDomains()

	def find_domains_blast(self, dic):
		from Bio.Blast import NCBIWWW
		from Bio.Blast import NCBIXML

		for prot in dic.keys():
			if 'hypothetical' in dic[prot][0].lower() or 'uncharacterized' in dic[prot][0].lower() or 'unknown' in dic[prot][0].lower():
				result_handle = NCBIWWW.qblast('blastp', 'nr', prot, entrez_query='Acinetobacter baumannii (taxid:470), Escherichia coli (taxid:562), Klebsiella pneumonia (taxid:573)')
				blastout = NCBIXML.read(result_handle)
				for ali in blastout.alignments:
					if 'hypothetical' not in ali.hit_def.lower() and 'uncharacterized' not in ali.hit_def.lower():
						print(ali.hit_def[:ali.hit_def.find(' [')])
						self.phagesProteins[spec][prot][0] = ali.hit_def[:ali.hit_def.find(' [')]
						break
		return dic

	def fillDomainsUniProt(self):
		'''
		Using the UniProt website, similar sequences are obtained and the ones with function assigned are saved into the domains. Domains are saved in the protdomains variable.
		:return: phageDomains, a dictionary that, for each protein in a given species, has domains associated
		'''
		print('Finding functions/domains with UniProt')
		import requests
		import pickle
		from pathlib import Path
		my_file = Path("files/phage_list_uniprot")
		if my_file.is_file():
			with open('files/phage_list_uniprot', 'rb') as f:
				list_done = pickle.load(f)
		else:
			list_done = []
		for phage in self.phagesProteins:
			if phage not in list_done:
				for accID in self.phagesProteins[phage]:
					if 'hypothetical' in self.phagesProteins[phage][accID][0].lower() or 'uncharacterized' in self.phagesProteins[phage][accID][0].lower() or 'unknown' in self.phagesProteins[phage][accID][0].lower():
					# if not self.phageDomains[phage][accID]:
						fullURL = ('https://www.uniprot.org/uniprot/?query=' + accID + '&sort=score&format=list')
						result = requests.get(fullURL)
						uniprot_acc = result.text.strip()
						fullURL = ('https://www.uniprot.org/uniprot/?query=cluster:(uniprot:' + uniprot_acc + '* identity:1.0) not id:' + uniprot_acc + '&format=txt')
						result = requests.get(fullURL)
						listResults = result.text.split('\n')
						for entry in listResults:
							if entry[:2] == 'DE':
								start_pos = entry.find('Full=') + 5
								end_pos = entry.find(' {ECO')
								domain = entry[start_pos:end_pos]
								if not any(z in domain.lower() for z in ['uncharacterized', 'flags', 'domain', 'bacteriophage protein', 'family protein', 'phage-like', 'phage protein', 'unassigned', 'orf', 'gene']) and len(domain) > 5:
									print(domain)
									self.phagesProteins[phage][accID][0] = domain
									break
				list_done.append(phage)
				with open('files/phage_list_uniprot', 'wb') as f:
					pickle.dump(list_done, f)
				self.saveDomains()

	def find_domains_uniprot(self, dic):
		import requests
		for accID in dic.keys():
			if 'hypothetical' in dic[accID][0].lower() or 'uncharacterized' in dic[accID][0].lower() or 'unknown' in dic[accID][0].lower():
				fullURL = ('https://www.uniprot.org/uniprot/?query=' + accID + '&sort=score&format=list')
				result = requests.get(fullURL)
				uniprot_acc = result.text.strip()
				fullURL = ('https://www.uniprot.org/uniprot/?query=cluster:(uniprot:' + uniprot_acc + '* identity:1.0) not id:' + uniprot_acc + '&format=txt')
				result = requests.get(fullURL)
				listResults = result.text.split('\n')
				for entry in listResults:
					if entry[:2] == 'DE':
						start_pos = entry.find('Full=') + 5
						end_pos = entry.find(' {ECO')
						domain = entry[start_pos:end_pos]
						if not any(z in domain.lower() for z in ['uncharacterized', 'flags', 'domain', 'bacteriophage protein', 'family protein', 'phage-like', 'phage protein', 'unassigned', 'orf', 'gene']) and len(domain) > 5:
							dic[accID][0] = domain
							break
		return dic

	def cdHit(self):
		import os
		from pathlib import Path
		my_file = Path('files/phagesProteins.fasta')
		if not my_file.is_file():
			self._create_fasta(self.phagesProteins, 'phagesProteins.fasta')
		my_file = Path('files/complete_cdhit.clstr')
		if not my_file.is_file():
			os.system('cd-hit -i files/phagesProteins.fasta -d 50 -o files/complete_cdhit')
		# clusters = {}
		temp_cluster = []
		list_found = []
		found = False
		with open('files/complete_cdhit.clstr', 'r') as f:
			for line in f.readlines():
				if '>Cluster' in line:
					if temp_cluster and found:
						if len(list_found) == 1:
							function = list_found[0]
						else:
							x = int(input(str(list_found) + '\nChoose from 1 to ' + str(len(list_found)) + ': ')) - 1
							function = list_found[x]
						for clust in temp_cluster:
							self.phagesProteins[clust[clust.find('-') + 1:]][clust[:clust.find('-')]][0] = function

					temp_cluster = []
					list_found = []
					found = False
				else:
					pos_i = line.find('>') + 1
					pos_f = line.find('...')
					pos_m = line.find('-')
					prot = line[pos_i:pos_m]
					phage = line[pos_m + 1:pos_f]
					if prot in self.known_function[phage].keys() and not found:
						function = self.known_function[phage][prot][0]
						list_found.append(function)
						found = True
					elif prot in self.known_function[phage].keys() and found:
						if function != self.known_function[phage][prot][0] and self.known_function[phage][prot][0] not in list_found:
							function = self.known_function[phage][prot][0]
							list_found.append(function)
					elif prot in self.unknown_function[phage].keys():
						temp_cluster.append(line[pos_i:pos_f])

	def create_blast_db(self):
		import os
		self._create_fasta(self.known_function, 'database_phages.fasta')
		os.system('makeblastdb -in files/database_phages.fasta -dbtype prot -title PhageProts -parse_seqids -out files/database_phages')
		self._create_fasta(self.unknown_function, 'unknown_phages.fasta')
		os.system('blastp -db files/database_phages -query files/unknown_phages.fasta -out files/test_blast -num_threads 2 -outfmt 6')

	def process_blastdb(self, blastdb):
		import pandas as pd
		blast_domains = pd.read_csv('files/' + blastdb, sep='\t', header=None)
		for phage in self.unknown_function.keys():
			for prot in self.unknown_function[phage]:
				evalue = []
				bitscore = []
				pred = blast_domains[blast_domains[0] == phage + '-' + prot]
				if pred.shape[0] == 0: break
				for i in pred[10]:
					evalue.append(float(i))
				for i in pred[11]:
					bitscore.append(float(i))
				if min(evalue) < 1.0 and max(bitscore) > 30.0:
					ind = evalue.index(min(evalue))
					if ind != bitscore.index(max(bitscore)):
						ind = bitscore.index(max(bitscore))
					temp = pred.iloc[ind,1]
					known_phage = temp[:temp.find('-')]
					known_prot = temp[temp.find('-')+1:]
					if self.known_function[known_phage][known_prot]:
						new_func = self.known_function[known_phage][known_prot][0]
					for j in self.known_function.keys():
						if pred.iloc[ind,1] in self.known_function[j].keys():
							new_func = self.known_function[j][pred.iloc[ind,1]][0]
							break
					self.phagesProteins[phage][prot][0] = new_func
		self.saveDomains()

	def extract_bact_location(self):
		import pandas as pd
		import ast
		import requests
		import re
		from pathlib import Path
		data = pd.read_csv('files/NCBI_Phage_Bacteria_Data.csv', header=0, index_col=0)
		all_bact = []
		for i in data.index:
			for bact in ast.literal_eval(data.loc[i, 'Host_ID']):
				if bact[:-2] not in all_bact:
					all_bact.append(bact[:-2])
		fullURL = ('https://db.psort.org/downloads/precomputed?version=3.00')
		result = requests.get(fullURL)
		psort = result.text.strip()
		urls = re.findall('https?://(?:[-\w.]|(?:%[\da-fA-F]{2}))+/[a-z]+/\S+\"{1}', psort)
		i = 1
		while i < len(urls):
			temp = urls[i]
			bact = temp[temp.rfind('=') + 1:temp.find('"')]
			if bact not in all_bact:
				i += 3
			else:
				my_file = Path('files/psort/' + bact + ".faa.out")
				if not my_file.is_file():
					temp_url = urls[i+1].strip('"')
					r = requests.get(temp_url)
					with open('files/psort/' + bact + ".faa.out", 'wb') as f:
						f.write(r.content)
				i += 3

	def create_fasta_psort(self):
		from pathlib import Path
		import json
		import os
		for bact in os.listdir('files/bacteria'):
			my_file = Path('files/psort/' + bact[:-5] + '.faa.out')
			if not my_file.is_file():
				with open('files/bacteria/' + bact, encoding='utf-8') as F:
					bact_prots = json.loads(F.read())
				self._create_fasta(bact_prots, 'psort/' + bact[:-5] + '.fasta')
			os.system('./psortb -n -i /home/pedro-linux/OneDrive/UMinho/Cenas_de_tese_idk/test_tese_process/files/psort/' + bact[:-5] + '.fasta -r . -o long')
			os.listdir('./psortb')
			os.replace('', '/home/pedro-linux/OneDrive/UMinho/Cenas_de_tese_idk/test_tese_process/files/psort/' + bact[:-5] + '.faa.out')  # move and rename output
			os.remove('files/psort/' + bact[:-5] + '.fasta')

	def saveDomains(self):
		'''
		Saves the protdomain variable in a file.
		:return: SearchedDomains.json
		'''
		import json
		with open('files/phagesProteins.json', 'w') as f:
			json.dump(self.phagesProteins, f)
		# with open('files/phagesProteins.fasta', 'w') as F:
		# 	for phage in self.phagesProteins.keys():
		# 		for prot in self.phagesProteins[phage]:
		# 			F.write('>' + prot + '\n' + self.phagesProteins[phage][prot][1] + '\n')


if __name__ == '__main__':
	test = DomainSearch()

	test.extract_bact_location()
	test.create_fasta_psort()

	test.create_blast_db()
	test.process_blastdb('test_blast')

	test.cdHit()
	test.scanInterPro()
	test.processInterPro()

	test.fillDomainsBLAST()
	test.fillDomainsUniProt()
	test.saveDomains()
