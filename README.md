# HSS-PPI
***************************************************************************************************
## HSS-PPI tool for "A Hierarchical Representation for PPIs Sites Prediction" - version 0.1
***************************************************************************************************

***************************************************************************************************
## REQUIREMENTS
***************************************************************************************************
Python 3 
singularity version >= 3.4.2 (see https://sylabs.io/guides/3.4/user-guide/installation.html for information)


***************************************************************************************************
## INSTALLATION
***************************************************************************************************

Download the zip file of last version of HSS-PPI from folder "download" at

https://gitlab.com/sebastiandaberdaku/hss-ppi.git

and put it in any position of your drive. 



***************************************************************************************************
## RUNNING THE SIMULATION 
****************************************************************************************************


To run the simulation for the bound receptors of protein class X (X can be any of the following protein classes: A, AB, EI, ER, ES, OG, OR, OX) use the following command: 
```bash
singularity exec singularity-python.sif python3 train.py --protein_class $X --is_bound --is_receptor
```

To run the simulation for the bound ligands of protein class X, use the following command: 
```bash
singularity exec singularity-python.sif python3 train.py --protein_class $X --is_bound
```


To run the simulation for the unbound receptors of protein class X, use the following command: 
```bash
singularity exec singularity-python.sif python3 train.py --protein_class $X --is_receptor
```

To run the simulation for the unbound ligands of protein class X, use the following command: 
```bash
singularity exec singularity-python.sif python3 train.py --protein_class $X
```


The user can specify a set of parameters: the learning rate, dropout, hidden, weight_decay, epochs, distance threshold and representation type. Their default values are 0.001, 0.5, [35,32], 5e-4, 500,12 and hierarchical, respectively. For any modification, please check out *utils/config.py*. Once the dataset of a given protein class is loaded, it is saved as a pickle in the data folder. Users can use the *--load_pickle* flag to speed up subsequent computation times once the pickle file has been created.


***************************************************************************************************
## EXAMPLES 
***************************************************************************************************
You use the following command:  

```bash
singularity exec singularity-python.sif python3 train.py --protein_class A --is_bound --is_receptor --epochs 1390 --distance_threshold 6 --representation_type contact_map
```

to run the simulation for the bound receptors of protein class A training the model with 1390 epochs and abstracting the proteins as contact map representation and setting the distance threshold (distance between C_alpha - C_alpha carbons) to 6 \AA. The omitted parameters (learning rate, dropout, hidden, weight_decay) take the default values.  


You use the following command:  

```bash
singularity exec singularity-python.sif python3 train.py --protein_class AB --is_bound --epochs 805 --distance_threshold 8 --representation_type contact_map
```

to run the simulation for the bound ligands of protein class OG training the model with 905 epochs and abstracting the proteins as contact map representation and setting the distance threshold (distance between C_alpha - C_alpha carbons) to 10 \AA. The omitted parameters (learning rate, dropout, hidden, weight_decay) take the default values.  

You use the following command:  

```bash
singularity exec singularity-python.sif python3 train.py --protein_class EI --is_receptor --representation_type hierarchical --distance_threshold 8 --epochs 835
```

To run the simulation for the unbound receptors of protein class EI training the model with 835 epochs and abstracting the proteins as hierarchical representation and setting the distance threshold (distance between C_alpha - C_alpha carbons) to 8 \AA. The omitted parameters (learning rate, dropout, hidden, weight_decay) take the default values.  


You use the following command:  

```bash
singularity exec singularity-python.sif python3 train.py --protein_class AB --epochs 580 --distance_threshold 10 --representation_type hierarchical
```
to run the unbound ligands of protein class AB training the model with 580 epochs and abstracting the proteins as hierarchical representation and setting the distance threshold (distance between C_alpha - C_alpha carbons) to 10 \AA. The omitted parameters (learning rate, dropout, hidden, weight_decay) take the default values.  
 

 
***************************************************************************
## DATA FORMAT and PREPROCESSING 
**************************************************************************
To prepare the training, validation and testing datasets for the interface prediction task, it was necessary to determine the distance matrix and  the sequence matrix for each protein, the amino acids of each proteins that interaction interface with the other molecule,  and the ASA, the angle Phi and Psi, the relative ASA  and the number of each contacts. All the preprocessing scripts are available in the *preprocessing* folder.

The following preprocessing steps were used

### determine_interface.py

To compute the interaction interface for each monomer, the *determine_interface.py* script can be used. The script can be executed with the following command from the root of the project:
```bash
singularity exec singularity-python.sif python3 preprocessing/determine_interface.py
```

### aa_distance_matrix.py 
To compute the amino acid distance matrix for each monomer, the *aa_distance_matrix.py* script can be used. The script can be executed with the following command from the root of the project: 
```bash
singularity exec singularity-python.sif python3 preprocessing/aa_distance_matrix.py
```
The computation of all distance matrices is a bit slow because the script evaluate all atom pairs from each residue pair in all PDF files, but it needs to be computed only once. The script will also try to use all the available cores on the host machine to speed up the computation in parallel.

### aa_sequence_matrix.py
To compute the amino acid sequence matrix for each monomer, the *aa_sequence_matrix.py* script can be used. The script can be executed with the following command from the root of the project: 
```bash
singularity exec singularity-python.sif python3 preprocessing/aa_sequence_matrix.py
```

### ASA.py
To compute the ASA for each amino acid, the *ASA.py* script can be used. The script can be executed with the following command from the root of the project: 
```bash
singularity exec singularity-python.sif python3 preprocessing/ASA.py
```

### PHI-PSI.py
To compute the PHI-PSI angles for each amino acid, the *PHI-PSI.py* script can be used. The script can be executed with the following command from the root of the project: 
```bash
singularity exec singularity-python.sif python3 preprocessing/PHI-PSI.py
```
### relative_ASA.py
To compute the relative ASA for each amino acid, the *relative_ASA.py* script can be used. The script can be executed with the following command from the root of the project: 
```bash
singularity exec singularity-python.sif python3 preprocessing/relative_ASA.py
```

### residue_contacts.py
To compute the residue contacts for each monomer, the *residue_contacts.py* script can be used. The script can be executed with the following command from the root of the project: 
```bash
singularity exec singularity-python.sif python3 preprocessing/residue_contacts.py
```

***************************************************************************
CONTACT INFORMATION
***************************************************************************

Please report any issue to michela.quadrini@unicam.it or to Michela Quadrini, Polo
Informatico, via Madonna delle Carceri 9, 62032 Camerino (MC) Italy.
