# Cheminformatic SMILES to CGenFF Pipeline
# author: Suliman Sharif
# email: sharifsuliman1@gmail.com
# notes: Please email me if you don't know how to run this :)

# Standard Python Internal Packages
# ---------------------------------
import os, sys
import pandas as pd 
import datetime as dt

sys.path.insert(0,os.path.abspath(os.path.dirname(__file__)))

# Airflow & Configurations
# ------------------------
from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from airflow.operators.python_operator import PythonOperator

default_args = {
    'owner': 'sulstice',
    'start_date': dt.datetime(2018, 9, 24, 10, 00, 00),
    'concurrency': 1,
    'retries': 0
}

# Logging Configurations
# ----------------------

# Airflow Steps
# -------------

from task_1.task_1 import * 
from task_2.task_2 import *

with DAG('smiles_to_sdf',
         default_args=default_args,
         schedule_interval='*/10 * * * *',
         catchup=False,
         ) as dag:

    step1 = PythonOperator(
        task_id='smiles_to_sdf',
        python_callable=step_1
    )

step1 
