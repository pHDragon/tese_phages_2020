
**PhageHostPrediction**

Predict interactions between phages and bacterial strains.

PhageHostPrediction is a python script that predicts phage-host interactions for *E. coli*, *K. pneumoniae* and *A. baumannii* phages, using supervised machine learning models. The models were built from a dataset containing 252 features and 23 987 entries with balanced outputs of 'Yes' and 'No'. The positive cases of interaction predicted are described in the file "NCBI_Phage_Bacteria_Data.csv", contained within this tool, while the negative were randomly assigned by pairing phages with bacteria of different species.

The prediction resorts to complete host proteome and to phage tail proteins, that are inferred within the tool. This inference is made with a locally created database of phage protein functions, available in the file "phagesProteins.json". Unknown proteins are predicted against this database. To help with this prediction, the use of InterProScan is made optional.

**Inputs:**

* phage/bacteria genome format: ID vs fasta; 
* ID: must be a GenBank ID, with the proteome described;
* fasta file: must contain the whole proteome of the organism;
* machine learning model: random forests have better predictive power, while SVM can be slightly faster to run;
* interpro search: should predict tails with higher confidence, but it significantly increases time to run.

**Outputs:**
This tool outputs a tabular file in which phage-host pairs are present in the first column and the prediction result in the second.

**Requirements:**

* Biopython
* Scikit-learn 
* Numpy
* Pandas 
* Scikit-bio
* BLAST_ - must be installed locally and available globally as an environment variable
* InterProScan_ (optional) - must be installed locally and available globally as an environment variable

.. _BLAST: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
.. _InterProScan: http://www.ebi.ac.uk/interpro/download/
