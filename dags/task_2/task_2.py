# Pipeline Configurations
# -----------------------
morgan_radius = 1
bit_representation = 512

# Imports 
# -------
import sys

# Scientific Imports
# ------------------
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

# RDkit Imports 
# -------------
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Draw
from rdkit.Chem import DataStructs
rdDepictor.SetPreferCoordGen(True)

# Graphing Imports 
# ----------------

from bokeh.plotting import ColumnDataSource, figure, output_notebook, output_file, show, save
from bokeh.io import output_notebook, export_png
from bokeh.layouts import gridplot
output_notebook()

# Global Configs
# --------------

TOOLTIPS = """<div>\nMolID: @ids<br>\n@img{safe}\n</div>\n"""
colormaps = { 0: '#e6194b', 1: '#3cb44b',  2: '#ffe119', 3: '#4363d8', 4: '#f58231',  5: '#911eb4'}

# Standard Functions
# ------------------

def mol2svg(mol):
    AllChem.Compute2DCoords(mol)
    d2d = rdMolDraw2D.MolDraw2DSVG(200,100)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()

def mol2fparr(mol):

    global morgan_radius
    global bit_representation

    arr = np.zeros((0,))
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, morgan_radius, nBits=bit_representation)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def step_4():
    
    smiles_dataframe = pd.read_csv('/home/sulstice/airflow/dags/task_1/data/smiles.txt', sep='\n', header=None, names=['smiles'])
    smiles_list = smiles_dataframe['smiles'].to_list()

    molecules_list = [Chem.MolFromSmiles(i) for i in smiles_list]
    fingerprints_list = np.array([mol2fparr(m) for m in molecules_list])

    pca = PCA(n_components=0.95)
    chemicalspace = pca.fit_transform(fingerprints_list)
    kmean = KMeans(n_clusters=5, random_state=0)
    kmean.fit(fingerprints_list)
    kmeanc = [colormaps[i] for i in kmean.labels_]

    kmean_data = dict(
        x=chemicalspace[:,0],
        y=chemicalspace[:,2],
        img=[mol2svg(m) for m in molecules_list],
        ids=[str(i) for i in range(0, len(smiles_list))],
        fill_color=kmeanc,
    )

    source = ColumnDataSource(kmean_data)
    plot = figure(plot_width=500, plot_height=500, tooltips=TOOLTIPS, title='error_compounds')
    plot.circle('x', 'y',color='fill_color', size=10, fill_alpha=0.2,source=source)
    
    plot = gridplot([
        [plot]
    ])

        
    output_file(filename="/home/sulstice/airflow/plugins/plotting_plugins/templates/pca_analysis.html", title="Static HTML file")

    save(plot)
