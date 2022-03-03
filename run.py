# Imports 

import os
import json
import argparse
import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import deque
import gzip, sys
from tqdm import tqdm
from collections import defaultdict
from typing import Any, Callable, Dict, Optional, Text, Union, Iterable
import threading



def get_vcf_names(vcf_path):
    """
    This function return the name of the collumn in a given vcf file
    """
    with gzip.open(vcf_path, "rt") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line.split('\t')]
                break
                
    ifile.close()
    return vcf_names

def get_dict_range(df):
    """
    This function create a maping between chromosome id and known
    SNP position related to PD desease
    """
    dict_chrom_range = defaultdict(lambda : set())
    for index, row in tqdm(df.iterrows(), \
        desc="looping on filtered file to create dictionary", total = len(df)):
        for number in range(row['start_pos'], row['end_pos']+1):
            dict_chrom_range[row['chrm_id']].add(number)
        
    return dict_chrom_range

def clean_and_make(vcf_chunk, df_clean, df_created, lock):
    
    lock.acquire()
    current_df = vcf_chunk.get_chunk()
    lock.release()
    
    current_df = current_df[~current_df["ALT"].isin(['<NON_REF>'])]
            
    for index, row in tqdm(current_df.iterrows(), \
                desc="looping a chunk of the dataframe", total = len(current_df)):

        ALT = row['ALT'].split("<NON_REF>")[0][:-1] if row['ALT'].split(",")[-1]=="<NON_REF>" else row['ALT'] 

        new_row = {'#CHROM':row['#CHROM'], 'POS':f"{row['POS']}", 'ID':row['ID'], 'REF':row['REF'],
                   'ALT':ALT, 'QUAL':row['QUAL'], 'FILTER':row['FILTER'], 'INFO':row['INFO']}
        df_clean = df_clean.append(new_row, ignore_index=True)

        if int(row["POS"]) in dict_chrom_range[row["#CHROM"]]:        
            df_created = df_created.append(new_row, ignore_index=True)
    
        
#     return df_created, current_df


def read_db_by_chunk(file_path, dict_chrom_range):
        
    names = get_vcf_names(file_path)
    vcf_chunk = pd.read_csv(file_path, 
                compression='gzip', comment='#', 
                chunksize=5000000, 
                delim_whitespace=True, header=None, names=names)
    
    df_clean = pd.DataFrame(columns=['#CHROM', 'POS', 'ID', 'REF', \
                                            'ALT', 'QUAL', 'FILTER', 'INFO'])
    
    df_created = pd.DataFrame(columns=['#CHROM', 'POS', 'ID', 'REF', \
                                            'ALT', 'QUAL', 'FILTER', 'INFO'])
    
    while True:
        try:

            # creating a lock
            lock = threading.Lock()
            
            # creating threads
            tasks = list(range(threads))
            
            # creating threads
            for i in range(len(tasks)):
                tasks[i] = threading.Thread(target=clean_and_make, args=(vcf_chunk, df_clean, df_created, lock,))
            # t1 = threading.Thread(target=clean_and_make, args=(vcf_chunk, df_clean, df_created, lock,))
            # t2 = threading.Thread(target=clean_and_make, args=(vcf_chunk, df_clean, df_created, lock,))
                        
            # start threads
            for i in range(len(tasks)):
                tasks[i].start()
            # t1.start()
            # t2.start()
                                  
            # wait until threads finish their job
            for i in range(len(tasks)):
                tasks[i].join() 
            # t1.join()
            # t2.join()
                    
        except:
            print("Done looping through data")
            extype, value, tb = sys.exc_info()
            print(f"Error details : {value}")
            
            # wait until threads finish their job
            for i in range(len(tasks)):
                tasks[i].join() 
              
            break
                
    return df_clean, df_created              


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Run .vcf cleaning pipeline")
    
    parser.add_argument("-fd", "--files_directory", default="/storage0/psychiatric_disorders", help="Path to the directoty containing our .vcf files")
    parser.add_argument("-gf", "--genome_filter", default="/storage0/psychiatric_disorders/results/PD_associated_genes_position_data_10thJan2022.csv", help="Path to the filtering files")
    parser.add_argument("-o", "--output", default="/storage0/psychiatric_disorders/data_proccessed", help="Output Path to processed data")
    parser.add_argument("-T", "--threads", default="5", help="Number of thread for the code")

    # for terminal arg
    # args = parser.parse_args(args=[])
    args = parser.parse_args()

    files_directory = args.files_directory
    output = args.output
    threads = int(args.threads)
    
    if not os.path.exists(f"{output}"):
        os.makedirs(f"{output}")
        
    genome_filter = args.genome_filter
    
    PD_associated_genes = genome_filter
    gene_filters = pd.read_csv(PD_associated_genes)
    dict_chrom_range = get_dict_range(gene_filters)
    
    
    dir_list = os.listdir(files_directory)
    for file in dir_list:
        # 1.6G <> 9.7G    
    
        # file = "LRRK2AJ6290001.g.vcf.gz"
        start_time = time.time()    

        if len(file.split(".")) != 4: 
            print(file.split("."))
    #         continue 
        else:
            file_name, _, _, _ = file.split(".")

        df_clean, df_created = read_db_by_chunk(f"{files_directory}/{file}", dict_chrom_range)

        df_clean.to_csv(f"{output}/{file_name}_filtered.csv", index=False)
        df_created.to_csv(f"{output}/{file_name}_clean.csv", index=False)

        runtime = round(time.time() - start_time, 3)

        print("#####################################")
        print(f"file {file_name} runtime: {(runtime/60)} minutes")
        print("#####################################")