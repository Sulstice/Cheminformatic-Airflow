# Task 1: To Convert the SMILES to SDF

# Imports
# -------
import os
import pandas as pd

# RDKit & Configurations
# ----------------------
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import DrawingOptions

# SMILES to SDF
# -------------

def smiles_to_sdf(row):

    '''

    Writes the smiles into SDF

    Arguments:
        row (str): row of data to process that contains the smiles and the name.
        database (str): name of the database you would want to process.

    '''
    molecule = Chem.MolFromSmiles(row['smiles'])
    molecule_with_hs = Chem.AddHs(molecule)
    molecule_with_hs.SetProp('smiles', row['smiles'])

    sdwriter = Chem.SDWriter(os.path.join( '/home/sulstice/airflow/dags/task_1/data/sdf_files', str(row.name) + '.sdf'))
    sdwriter.write(molecule_with_hs)

def step_1():
    smiles_dataframe = pd.read_csv('/home/sulstice/airflow/dags/task_1/data/smiles.txt', sep='\n', header=None, names=['smiles'])
    _ = smiles_dataframe.apply(smiles_to_sdf, axis=1)
