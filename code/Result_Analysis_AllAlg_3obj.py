"""
Analysis the pareto optimal solutions for optimisation letter paper
"""

import pandas as pd
import numpy as np
import glob
import os
from pymoo.indicators.hv import HV
from pymoo.indicators.gd import GD
from pymoo.indicators.gd_plus import GDPlus
from pymoo.indicators.igd import IGD
from pymoo.indicators.igd_plus import IGDPlus

# Function to find Pareto-optimal points from a given points
def is_pareto_efficient(costs):
    """
    Determine Pareto efficiency. Points are Pareto-efficient if no other point
    is strictly better considering (NPC, GHE, and ERC).
    """
    is_efficient = np.ones(costs.shape[0], dtype=bool)  # All points are initially considered efficient
    for i, c in enumerate(costs):
        if is_efficient[i]:
            # Keep points that are not dominated by this one
            is_efficient[is_efficient] = np.any(costs[is_efficient] < c, axis=1) | np.all(costs[is_efficient] == c, axis=1)
    return is_efficient

###----------------------------------------------------------------------------------
##### Function for reading the simulation results for all algorithms

def read_and_combine_csv_files(names_X, base_path, file_pattern, Pareto_Flag):
    """
    Reads and combines CSV files for each name in names_X based on the provided base path and file pattern.

    Parameters:
    - names_X (list): List of names (X) used in the CSV file names.
    - base_path (str): The base directory where the CSV files are stored.
    - file_pattern (str): The pattern to match the CSV file names (e.g., "AllPop_{X}_noBAT_seedno_*.csv").
    - Pareto_Flag: The flag which indicates the csv files are for all pop or pareto (True or False ).

    Returns:
    - combined_dataframes_dict (dict): A dictionary containing combined DataFrames for each name in names_X.
    """
    # Dictionary to store combined DataFrames, keys are X
    combined_dataframes_dict = {}
    
    # Loop through each name in names_X
    for X in names_X:
        # Create the complete file path pattern using the provided pattern
        complete_pattern = os.path.join(base_path, file_pattern.format(X=X))
        
        # Get the list of files matching the pattern
        file_list = glob.glob(complete_pattern)
        
        # Initialize an empty dataframe to store selected columns
        combined_df = pd.DataFrame()
        
        # Flag to check if "Gen(#)" column has been added
        gen_column_added = False
        
        # Loop through each file and read the required columns into a dataframe
        for i, file in enumerate(file_list):
            
            if not Pareto_Flag:
            
                # Read the specific columns from the CSV file
                if not gen_column_added:
                    # Include the "Gen(#)" column only the first time
                    df = pd.read_csv(file, usecols=["Gen (#)", "NPC ($/kWh)", "GHE (kg CO2/kWh)", "ERC (GWh/yr)", "ENS (%)"])
                    gen_column_added = True  # Mark that the "Gen(#)" column is added
                    df.columns = ["Gen (#)",f"NPC_{i}", f"GHE_{i}", f"ERC_{i}", f"ENS_{i}"]
                else:
                    # Exclude the "Gen(#)" column for subsequent files
                    df = pd.read_csv(file, usecols=["NPC ($/kWh)", "GHE (kg CO2/kWh)", "ERC (GWh/yr)", "ENS (%)"])
                
                # Rename the columns (except "Gen(#)") to distinguish between files
                if i > 0:
                    df.columns = [f"NPC_{i}", f"GHE_{i}", f"ERC_{i}", f"ENS_{i}"]
                
                # Concatenate the dataframe column-wise (axis=1)
                combined_df = pd.concat([combined_df, df], axis=1)
            else:

                df = pd.read_csv(file, usecols=["NPC ($/kWh)", "GHE (kg CO2/kWh)", "ERC (GWh/yr)", "ENS (%)"])
                
                # Concatenate the dataframe row-wise (axis=0)
                combined_df = pd.concat([combined_df, df], axis=0, ignore_index=True)
                
            # Add combined DataFrame to the dictionary
            combined_dataframes_dict[X] = combined_df

    return combined_dataframes_dict

# ###----------------------------------------------------------------------------------
# ##### Reading the simulation results of ALL pop and pareto solutions

# Usage
names_X = ["NSGAII", "NSGA3", "UNSGA3", "SMSEMOA", "AGEMOEA", "AGEMOEA2"]
base_path = "Results//"
file_pattern_Allpop = "AllPop_{X}_noBAT_seedno_*.csv"
file_pattern_Pareto = "Pareto_{X}_noBAT_seedno_*.csv"

# Call the function with the list of names, the base path, and the file pattern
AllPop_dataframes_dict = read_and_combine_csv_files(names_X, base_path, file_pattern_Allpop, Pareto_Flag=False)

# Call the function with the list of names, the base path, and the file pattern
Pareto_dataframes_dict = read_and_combine_csv_files(names_X, base_path, file_pattern_Pareto, Pareto_Flag=True)

df_Pareto = pd.DataFrame()
for X in names_X:
    df = Pareto_dataframes_dict[X]
    df_Pareto = pd.concat([df_Pareto, df], axis=0, ignore_index=True)
    
df_Pareto = df_Pareto.loc[df_Pareto["ENS (%)"] == 0].reset_index(drop=True)
df_Pareto =df_Pareto.drop("ENS (%)", axis=1)

costs = np.column_stack((df_Pareto["NPC ($/kWh)"], df_Pareto["GHE (kg CO2/kWh)"], df_Pareto["ERC (GWh/yr)"]))  # Shape (n, 3)
pareto_mask = is_pareto_efficient(costs)
Final_pareto = df_Pareto.iloc[pareto_mask]


###----------------------------------------------------------------------------------
#### Defining performance indicators for each generation

pf = np.array(Final_pareto) 

## 1- HV for reference value of DG only case
ref_point = np.array([0.1255411, 0.6844676, 0]) # reference point for only DG case 11 MW (NPC,GHE)
ind_HV = HV(ref_point=ref_point)

### 2- GD, GD+, IGD, IGD+ 
ind_GD = GD(pf)
ind_GDplus = GDPlus(pf)
ind_IGD = IGD(pf)
ind_IGDplus = IGDPlus(pf)


###----------------------------------------------------------------------------------
##### Determine the Pareto for each iteration

# Dictionary to store Pareto solutions for each (Alg, Gen, i)
pareto_solutions_dict = {}

HyperVal = {}
GDVal = {}
GDplusVal = {}
IGDVal = {}
IGDplusVal = {}

HV_Avg = {}
GD_Avg = {}
GDplus_Avg = {}
IGD_Avg = {}
IGDplus_Avg = {}

for X in names_X:
    
    df = AllPop_dataframes_dict[X]
    # Extract NPC and GHE column names
    npc_cols = [col for col in df.columns if col.startswith("NPC")]
    ghe_cols = [col for col in df.columns if col.startswith("GHE")]
    erc_cols = [col for col in df.columns if col.startswith("ERC")]
    ens_cols = [col for col in df.columns if col.startswith("ENS")]
    
    # Ensure that we have the same number of NPC and GHE columns
    assert len(npc_cols) == len(ghe_cols), "Mismatch in number of NPC and GHE columns."
    assert len(npc_cols) == len(erc_cols), "Mismatch in number of NPC and ERC columns."

    genMax = np.max(df.loc[:,"Gen (#)"])+1
    
    ## 
    HVal = np.zeros([genMax,len(npc_cols)])
    GVal = np.zeros([genMax,len(npc_cols)])
    GpVal = np.zeros([genMax,len(npc_cols)])
    IGVal = np.zeros([genMax,len(npc_cols)])
    IGpVal = np.zeros([genMax,len(npc_cols)])
    
    # Loop over each Gen value from 0 to 1000
    for gen_value in range(genMax):
        # Filter data corresponding to this specific Gen value
        gen_data = df[df["Gen (#)"] == gen_value]
        
        # Loop over each i (from 0 to the number of NPC columns)
        for i in range(len(npc_cols)):

            # Extract NPC_i, GHE_i, and ENS_i for the current i
            npc_i = gen_data[npc_cols[i]].values
            ghe_i = gen_data[ghe_cols[i]].values
            erc_i = gen_data[erc_cols[i]].values
            ens_i = gen_data[ens_cols[i]].values

            mask = ens_i == 0
            gg = gen_data.iloc[mask][[npc_cols[i], ghe_cols[i], erc_cols[i]]]
            
            # Stack NPC and GHE values as a cost matrix
            costs = np.column_stack((npc_i[mask], ghe_i[mask], erc_i[mask]))  # Shape (n, 3)
            
            # Determine Pareto-efficient points
            pareto_mask = is_pareto_efficient(costs)
            
            # Retrieve the Pareto-optimal solutions for the current (Gen, i)
            pareto_solutions = gg.iloc[pareto_mask][[npc_cols[i], ghe_cols[i], erc_cols[i]]]
            
            # Store the Pareto solutions in the dictionary using the (gen_value, i) as the key
            pareto_solutions_dict[(X, gen_value, i)] = pareto_solutions
            
            # calculating the HV for each Gen(#) and seed_no
            HVal[gen_value,i] = ind_HV(np.array(pareto_solutions))
            GVal[gen_value,i] = ind_GD(np.array(pareto_solutions))
            GpVal[gen_value,i] = ind_GDplus(np.array(pareto_solutions))
            IGVal[gen_value,i] = ind_IGD(np.array(pareto_solutions))
            IGpVal[gen_value,i] = ind_IGDplus(np.array(pareto_solutions))
        
    HyperVal[X]= HVal
    GDVal[X] = GVal
    GDplusVal[X]= GpVal
    IGDVal[X]= IGVal
    IGDplusVal[X]= IGpVal

    HV_Avg[X]=np.mean(HVal[2:-1],axis=1)
    GD_Avg[X]=np.mean(GVal[2:-1],axis=1)
    GDplus_Avg[X]=np.mean(GpVal[2:-1],axis=1)
    IGD_Avg[X]=np.mean(IGVal[2:-1],axis=1)
    IGDplus_Avg[X]=np.mean(IGpVal[2:-1],axis=1)



