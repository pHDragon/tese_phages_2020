
class PredictInteraction:

	def __init__(self, data = 'FeatureDataset'):
		import pickle
		from sklearn.preprocessing import LabelEncoder
		with open('files/' + data, 'rb') as f:
			self.dataset = pickle.load(f)
		self.dataset = self.dataset.dropna()
		le = LabelEncoder()
		le.fit(['Yes', 'No'])
		self.output = le.transform(self.dataset['Infects'].values)
		self.dataset = self.dataset.drop('Infects', 1)
		self.__standardize()
		self.__split_train_test()

	def __standardize(self):
		from sklearn.preprocessing import StandardScaler
		self.scaler = StandardScaler()
		self.scaler.fit(self.dataset)
		self.data_z = self.scaler.transform(self.dataset)

	def __cross_validation(self, method):
		from sklearn.model_selection import cross_val_score
		from sklearn.model_selection import StratifiedKFold
		# from sklearn.model_selection import ShuffleSplit
		# cv = ShuffleSplit(n_splits=4, test_size=0.3) # cross validation normal
		skf = StratifiedKFold(5, shuffle=True)
		scores = cross_val_score(method, self.data_z, self.output, cv=skf, scoring='f1_weighted')
		print(scores)
		return skf

	def __split_train_test(self):
		from sklearn.model_selection import train_test_split
		self.data_z, self.X_val, self.output, self.y_val = train_test_split(self.data_z, self.output, test_size=0.2)
		self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.data_z, self.output, test_size=0.3)

	def run_knn(self, n=2):
		from sklearn.neighbors import KNeighborsClassifier
		from sklearn.metrics import confusion_matrix
		neigh = KNeighborsClassifier(n_neighbors=2, weights='distance', p=1, algorithm='auto', leaf_size=5)
		neigh.fit(self.X_train, self.y_train)
		# print(neigh.score(self.X_test, self.y_test))
		self.__score_metrics(neigh)
		print(confusion_matrix(self.y_test, neigh.predict(self.X_test)))
		# import time
		# start_time = time.time()
		# cv = self.__cross_validation(neigh)
		# print("--- %s seconds ---" % (time.time() - start_time))
		# self.__plot_roc_cv(neigh, cv)
		# self.__auc_curve(neigh)
		# self.__permutation_importance(self._hyperparameters_knn(neigh))

	def _hyperparameters_knn(self, method):
		from sklearn.model_selection import GridSearchCV
		parameters = {'leaf_size': [5, 15, 30, 50], 'n_neighbors': [2, 3, 5, 7], 'weights': ['uniform', 'distance'], 'algorithm': ['auto', 'ball_tree', 'kd_tree', 'brute'], 'p': [1, 2]}
		clf = GridSearchCV(method, parameters)
		clf.fit(self.X_val, self.y_val)
		self.__score_metrics(clf)
		print(clf.best_params_)
		return clf

	def run_random_forest(self):
		from sklearn.ensemble import RandomForestClassifier
		from sklearn.metrics import confusion_matrix
		clf = RandomForestClassifier(n_estimators=200, bootstrap=False, criterion='gini', min_samples_leaf=2, min_samples_split=4, oob_score=False)
		clf = clf.fit(self.X_train, self.y_train)
		# print(clf.score(self.X_test, self.y_test))
		self.__score_metrics(clf)
		print(confusion_matrix(self.y_test, clf.predict(self.X_test)))
		# self.recursive_feature_elimination(clf)
		# import time
		# start_time = time.time()
		# cv = self.__cross_validation(clf)
		# print("--- %s seconds ---" % (time.time() - start_time))
		# self.__plot_roc_cv(clf, cv)
		# self.__auc_curve(clf)
		# self.__permutation_importance(clf)

	def _hyperparameters_rf(self, method):
		from sklearn.model_selection import GridSearchCV
		parameters = {'n_estimators': [50, 100, 150, 200, 250], 'criterion': ['gini', 'entropy'], 'min_samples_split': [2, 4, 6], 'min_samples_leaf': [2, 4, 6], 'bootstrap': [True, False], 'oob_score': [True, False]}
		clf = GridSearchCV(method, parameters)
		clf.fit(self.X_val, self.y_val)
		self.__score_metrics(clf)
		print(clf.best_params_)
		return clf

	def run_svm(self, c=0.1):
		from sklearn.svm import SVC
		from sklearn.metrics import confusion_matrix
		svm = SVC(C=10, degree=2, gamma='auto', kernel='rbf')
		svm = svm.fit(self.X_train, self.y_train)
		# print(svm.score(self.X_test, self.y_test))
		self.__score_metrics(svm)
		print(confusion_matrix(self.y_test, svm.predict(self.X_test)))
		# import time
		# start_time = time.time()
		# cv = self.__cross_validation(svm)
		# print("--- %s seconds ---" % (time.time() - start_time))
		# self.recursive_feature_elimination(svm)
		# self.__plot_roc_cv(svm, cv)
		# self.__auc_curve(svm)
		# self.__permutation_importance(svm)

	def _hyperparameters_svm(self, method):
		from sklearn.model_selection import GridSearchCV
		parameters = {'C': [0.01, 0.1, 1, 10], 'kernel': ['linear','rbf','poly','sigmoid', 'precomputed'], 'degree': [2, 3, 4], 'gamma': ['scale', 'auto']}
		clf = GridSearchCV(method, parameters)
		clf.fit(self.X_val, self.y_val)
		self.__score_metrics(clf)
		print(clf.best_params_)
		return clf

	def run_neural_networks(self, alpha=1):
		from sklearn.neural_network import MLPClassifier
		from sklearn.metrics import confusion_matrix
		clf = MLPClassifier(alpha=0.0001, activation='tanh', hidden_layer_sizes=200, learning_rate='adaptive', solver='adam')
		clf.fit(self.X_train, self.y_train)
		self.__score_metrics(clf)
		print(confusion_matrix(self.y_test, clf.predict(self.X_test)))
		# import time
		# start_time = time.time()
		# cv = self.__cross_validation(clf)
		# print("--- %s seconds ---" % (time.time() - start_time))
		# self.__plot_roc_cv(clf, cv)
		# self.__auc_curve(clf)
		# self.__permutation_importance(clf)

	def _hyperparameters_ann(self, method):
		from sklearn.model_selection import GridSearchCV
		parameters = {'hidden_layer_sizes': [50, 100, 200], 'activation': ['identity', 'logistic', 'tanh', 'relu'], 'solver': ['lbfgs', 'sgd', 'adam'], 'alpha': [0.0001, 0.05], 'learning_rate': ['constant', 'invscaling', 'adaptive']}
		clf = GridSearchCV(method, parameters)
		clf.fit(self.X_val, self.y_val)
		self.__score_metrics(clf)
		print(clf.best_params_)
		return clf

	def run_logistic_reg(self, c=1):
		from sklearn.linear_model import LogisticRegression
		from sklearn.metrics import confusion_matrix
		clf = LogisticRegression(C=10, penalty='l2', solver='liblinear')
		clf.fit(self.X_train, self.y_train)
		self.__score_metrics(clf)
		print(confusion_matrix(self.y_test, clf.predict(self.X_test)))
		# import time
		# start_time = time.time()
		# cv = self.__cross_validation(clf)
		# print("--- %s seconds ---" % (time.time() - start_time))
		# self.__plot_roc_cv(clf, cv)
		# self.__auc_curve(clf)
		# self.__permutation_importance(clf)

	def _hyperparameters_lr(self, method):
		from sklearn.model_selection import GridSearchCV
		parameters = {'penalty': ['l1', 'l2', 'elasticnet', 'none'], 'C': [0.01, 0.1, 1, 10], 'solver': ['lbfgs', 'liblinear', 'newton-cg', 'sag', 'saga']}
		clf = GridSearchCV(method, parameters)
		clf.fit(self.X_val, self.y_val)
		self.__score_metrics(clf)
		print(clf.best_params_)
		return clf

	def hyperparameter_tuning(self):  # Best Params:  {'bootstrap': False, 'ccp_alpha': 0.0, 'criterion': 'gini', 'max_features': 'sqrt', 'n_estimators': 250}
		from sklearn.model_selection import GridSearchCV
		from sklearn.ensemble import RandomForestClassifier
		clf = RandomForestClassifier()
		n_estimators = [50, 100, 150, 200, 250]
		criterion = ['gini', 'entropy']
		max_features = ['auto', 'sqrt', 'log2']
		bootstrap = [True, False]
		ccp_alpha = [0.0, 0.01, 0.02]
		param_grid = dict(n_estimators=n_estimators, criterion=criterion, max_features=max_features, bootstrap=bootstrap, ccp_alpha=ccp_alpha)
		grid = GridSearchCV(clf, param_grid, n_jobs=-1)
		grid_result = grid.fit(self.X_val, self.y_val)
		print('Best Score: ', grid_result.best_score_)
		print('Best Params: ', grid_result.best_params_)

	def __score_metrics(self, method):
		from sklearn.metrics import matthews_corrcoef, f1_score, precision_score, recall_score
		print(matthews_corrcoef(self.y_test, method.predict(self.X_test)))
		print(f1_score(self.y_test, method.predict(self.X_test), average='binary'))
		print(precision_score(self.y_test, method.predict(self.X_test)))
		print(recall_score(self.y_test, method.predict(self.X_test)))

	def __auc_curve(self, method):
		import matplotlib.pyplot as plt
		from sklearn import metrics
		metrics.plot_roc_curve(method, self.X_test, self.y_test)
		plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
		plt.xlim(0, 0.2)
		plt.show()

	def __plot_roc_cv(self, method, cv):
		import numpy as np
		import matplotlib.pyplot as plt
		from sklearn import datasets
		from sklearn.metrics import auc
		from sklearn.metrics import plot_roc_curve
		tprs = []
		aucs = []
		mean_fpr = np.linspace(0, 1, 100)
		fig, ax = plt.subplots()
		for i, (train, test) in enumerate(cv.split(self.data_z, self.output)):
			method.fit(self.data_z[train], self.output[train])
			viz = plot_roc_curve(method, self.data_z[test], self.output[test], name='ROC fold {}'.format(i), alpha=0.3, lw=1, ax=ax)
			interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
			interp_tpr[0] = 0.0
			tprs.append(interp_tpr)
			aucs.append(viz.roc_auc)
		ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)
		mean_tpr = np.mean(tprs, axis=0)
		mean_tpr[-1] = 1.0
		mean_auc = auc(mean_fpr, mean_tpr)
		std_auc = np.std(aucs)
		ax.plot(mean_fpr, mean_tpr, color='b', label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc), lw=2, alpha=.8)
		std_tpr = np.std(tprs, axis=0)
		tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
		tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
		ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2, label=r'$\pm$ 1 std. dev.')
		ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="Receiver operating characteristic example")
		ax.legend(loc="lower right")
		plt.show()

	def __permutation_importance(self, method):
		from sklearn.inspection import permutation_importance
		r = permutation_importance(method, self.X_test, self.y_test, n_repeats=5)
		for i in r.importances_mean.argsort()[::-1]:
			if r.importances_mean[i] - 2 * r.importances_std[i] > 0.001:
				print(f"{self.dataset.columns[i]:<8}"
					  f" {r.importances_mean[i]:.3f}"
					  f" +/- {r.importances_std[i]:.3f}")

	def recursive_feature_elimination(self, method):
		from sklearn.feature_selection import RFECV
		selector = RFECV(method, cv=5)#, min_features_to_select=200)
		selector.fit(self.data_z, self.output)
		print(selector.ranking_)
		self.data_reduced = selector.transform(self.data_z)

	def predict_interaction(self, phage, bacteria):
		from sklearn.svm import LinearSVC
		from sklearn.ensemble import RandomForestClassifier
		import numpy as np
		from feature_construction import FeatureConstruction

		phageProts = self.__find_phage_proteins(phage)  # dictionary
		bactProts = self.__find_bact_proteins(bacteria)
		if not phageProts or not bactProts:
			print('oops')
			return None
		list_carb = {}
		list_prot = {}
		for prot in phageProts.keys():
			if any(z in phageProts[prot][0].lower() for z in ['lysin', 'collagen', 'glyco', 'galac', 'chitin', 'wall', 'pectin', 'glycan', 'sialidase', 'neuramin', 'amid', 'lysozyme', 'murami', 'pectate', 'sgnh']):
				list_carb[prot] = phageProts[prot]
			else:
				list_prot[prot] = phageProts[prot]
		inter = np.array([])
		fc = FeatureConstruction()
		grouping = fc.get_grouping(phage=list_prot, phage_carb=list_carb, bacteria=bactProts)  # takes a list of protein sequences
		inter = np.append(inter, grouping)
		composition = fc.get_composition(phage=list_prot, phage_carb=list_carb, bacteria=bactProts)
		inter = np.append(inter, composition)
		kmers = fc.get_kmers(phage=list_prot, phage_carb=list_carb, bacteria=bactProts)
		inter = np.append(inter, kmers)
		inter = inter.reshape(1, -1)
		inter = self.scaler.transform(inter)
		# svm = LinearSVC(C=0.01, tol=0.010, dual=False)
		clf = RandomForestClassifier(n_estimators=250, bootstrap=False, ccp_alpha=0.0, max_features='sqrt')
		clf = clf.fit(self.data_z, self.output)
		pred = clf.predict(inter)[0]
		print(pred)
		return pred

	def __find_phage_proteins(self, phage):
		import json
		with open('files/phageTails.json', encoding='utf-8') as F:
			phageTails = json.loads(F.read())
		phageProts = {}
		if phage in phageTails.keys():
			for prot in phageTails[phage]:
				phageProts[prot] = [phageTails[phage][prot][0], phageTails[phage][prot][1]]
		else:
			from domain_search import DomainSearch
			phageProts = self.__find_proteins(phage)
			ds = DomainSearch()
			phageProts = ds.find_domains_interpro(phageProts)
			phageProts = ds.find_domains_blast(phageProts)
			phageProts = ds.find_domains_uniprot(phageProts)
		return phageProts

	def __find_bact_proteins(self, bacteria):
		import os
		import json
		if bacteria + '.json' in os.listdir('files/bacteria'):
			with open('files/bacteria/' + bacteria + '.json', encoding='utf-8') as F:
				bactProts = json.loads(F.read())
		else:
			pass
			# bactProts = self.__find_proteins(bacteria)
			# Implementar previsão de localização celular
		return bactProts

	def __find_proteins(self, id):
		from Bio import Entrez
		from Bio import SeqIO
		Entrez.email = 'pedro_araujo97@hotmail.com'
		prots = {}
		with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=id) as handle:
			genomeBac = SeqIO.read(handle, "gb")
		for feat in genomeBac.features:
			if feat.type == 'CDS':
				try: prots[feat.qualifiers['protein_id'][0]] = [feat.qualifiers['product'][0], feat.qualifiers['translation'][0]]
				except: pass
		if len(genomeBac.features) <= 5:
			with Entrez.efetch(db="nucleotide", rettype="gbwithparts", retmode="text", id=id) as handle:
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
					prots[protKey] = [product, protSeq[:protSeq.find('"')]]
		return prots


if __name__ == '__main__':
	ml = PredictInteraction('dataset_reduced')  # feature_dataset
	# ml.predict_interaction('NC_050143', 'NC_020524.1')  # NC_010468.1 NC_013941.1 NZ_CP029060.1 NZ_CP027394.1 NZ_CP025089.1
	ml.predict_interaction('KM607000', 'NC_020524')  # NC_010468.1 NC_013941.1 NZ_CP029060.1 NZ_CP027394.1 NZ_CP025089.1
	ml.run_knn(2)
	ml.run_random_forest()
	ml.run_svm(0.001)
	ml.run_neural_networks(0.0001)
	ml.run_logistic_reg(0.01)

	import pandas as pd
	import ast
	data = pd.read_csv('files/NCBI_Phage_Bacteria_Data.csv', header=0, index_col=0)
	abaumannii = {}
	for phage in data.index:
		name = data.loc[phage, 'Host']
		if 'acinetobacter' in name.lower():
			for bact in ast.literal_eval(data.loc[phage, 'Host_ID']):
				abaumannii[bact] = 0
	list_yes = {}
	list_yes['KT588074'] = []
	for bact in abaumannii.keys():
		predict = ml.predict_interaction('KT588074', bact)
		if predict == 'Yes':
			list_yes['KT588074'].append(bact)
