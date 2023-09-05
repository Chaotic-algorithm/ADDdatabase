# ADD: a comprehensive database and benchmarking tool for specificity prediction of adenylation domains in nonribosomal peptide synthetases.
Nonribosomal peptide synthetases (NRPS) are modular enzymes that produce many important secondary metabolites, including antibiotics. Adenylation (A) domains within NRPS determine the final product by recognizing and activating its building blocks - specific amino acids. A-domain specificity prediction is vital for exploring metabolites and engineering NRPS pathways. However, the existing prediction tools have limited accuracy. The software developers lack a comprehensive resource with confirmed A-domain specificities and train their tools on ad hoc datasets.

## About the database
To address this gap, we present ADD, a database encompassing A-domain sequences, specificities, neighboring domains, biosynthetic gene clusters (BGCs), and producers' taxonomies. With 3459 entries, our database is the largest of its kind and includes both bacterial (3063) and fungal (396) A-domains. ADD incorporates and unifies records from previously published A-domain specificity datasets and MIBiG, the largest collection of experimentally validated links between secondary metabolite BGCs and their products.

## About benchmarking
We complement our database with a benchmarking utility for assessing the quality of specificity prediction algorithms. For a more detailed description of used scripts and files please visit the `benchmarking` folder.

### Disclaimer:
The code provided in this repository is a **sample** for benchmarking A-domain specificity prediction algorithms.
Here it contains only benchmarking of the SANDPUMA predictions (that is, SANDPUMA (ensemble), NRPSsp (pHMM), NRPSpredictor2 (ASM and SVM)) using a dataset of all A-domains.
A similar benchmarking procedure was applied for other tools and separately for both bacterial and fungal datasets.
Please note that this code is still in development to enhance its usability for a broader user base. Updates are coming soon!
