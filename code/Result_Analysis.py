"""
Analysis the pareto optimal solutions for optimisation letter paper
"""

import pandas as pd
import numpy as np
import glob
from pymoo.indicators.hv import HV
from pymoo.indicators.gd import GD
from pymoo.indicators.gd_plus import GDPlus
from pymoo.indicators.igd import IGD
from pymoo.indicators.igd_plus import IGDPlus

# Function to find Pareto-optimal points from a given points
def is_pareto_efficient(costs):
    """
    Determine Pareto efficiency. Points are Pareto-efficient if no other point
    is strictly better in both dimensions (NPC and GHE).
    """
    is_efficient = np.ones(costs.shape[0], dtype=bool)  # All points are initially considered efficient
    for i, c in enumerate(costs):
        if is_efficient[i]:
            # Keep points that are not dominated by this one
            is_efficient[is_efficient] = np.any(costs[is_efficient] <= c, axis=1)
            is_efficient[i] = True  # Keep this point
    return is_efficient



###----------------------------------------------------------------------------------
##### Reading the simulation results for all generations

# Specify the file path pattern
file_pattern = "Results//AllPop_NSGAII_noBAT_seedno_*.csv"

# Use glob to find all files matching the pattern
file_list = glob.glob(file_pattern)

# Initialize an empty dataframe to store selected columns
combined_df = pd.DataFrame()

# Flag to check if "Gen(#)" column has been added
gen_column_added = False

# Loop through each file and read the required columns into a dataframe
for i, file in enumerate(file_list):
    
    # Read the specific columns from the CSV file
    if not gen_column_added:
        # Include the "Gen(#)" column only the first time
        df = pd.read_csv(file, usecols=["Gen (#)", "NPC ($/kWh)", "GHE (kg CO2/kWh)"])
        gen_column_added = True  # Mark that the "Gen(#)" column is added
        df.columns = ["Gen (#)",f"NPC_{i}", f"GHE_{i}"]
    else:
        # Exclude the "Gen(#)" column for subsequent files
        df = pd.read_csv(file, usecols=["NPC ($/kWh)", "GHE (kg CO2/kWh)"])
    
    # Rename the columns (except "Gen(#)") to distinguish between files
    if i > 0:
        df.columns = [f"NPC_{i}", f"GHE_{i}"]
    
    # Concatenate the dataframe column-wise (axis=1)
    combined_df = pd.concat([combined_df, df], axis=1)



###----------------------------------------------------------------------------------
##### Reading the simulation results for pareto solutions

# Specify the file path pattern
file_pattern = "Results//Pareto_NSGAII_noBAT_seedno_*.csv"

# Use glob to find all files matching the pattern
file_list = glob.glob(file_pattern)

# Initialize an empty dataframe to store selected columns
df_Pareto = pd.DataFrame()

# Loop through each file and read the required columns into a dataframe
for i, file in enumerate(file_list):
    
    df = pd.read_csv(file, usecols=["NPC ($/kWh)", "GHE (kg CO2/kWh)"])
    
    # Concatenate the dataframe row-wise (axis=0)
    df_Pareto = pd.concat([df_Pareto, df], axis=0, ignore_index=True)

costs = np.column_stack((df_Pareto["NPC ($/kWh)"], df_Pareto["GHE (kg CO2/kWh)"]))  # Shape (n, 2)
pareto_mask = is_pareto_efficient(costs)
Final_pareto = df_Pareto.iloc[pareto_mask]



###----------------------------------------------------------------------------------
#### Defining performance indicators for each generation

pf = np.array(Final_pareto) 

## 1- HV for reference value of DG only case
ref_point = np.array([0.1255411, 0.6844676]) # reference point for only DG case 11 MW (NPC,GHE)
ind_HV = HV(ref_point=ref_point)

### 2- GD, GD+, IGD, IGD+ 
ind_GD = GD(pf)
ind_GDplus = GDPlus(pf)
ind_IGD = IGD(pf)
ind_IGDplus = IGDPlus(pf)



###----------------------------------------------------------------------------------
##### Determine the Pareto for each iteration

# Dictionary to store Pareto solutions for each (Gen, i) pair
pareto_solutions_dict = {}

# Extract NPC and GHE column names
npc_cols = [col for col in combined_df.columns if col.startswith("NPC")]
ghe_cols = [col for col in combined_df.columns if col.startswith("GHE")]

# Ensure that we have the same number of NPC and GHE columns
assert len(npc_cols) == len(ghe_cols), "Mismatch in number of NPC and GHE columns."

## 
HyperVal = np.zeros([1000,31])
GDVal = np.zeros([1000,31])
GDplusVal = np.zeros([1000,31])
IGDVal = np.zeros([1000,31])
IGDplusVal = np.zeros([1000,31])

# Loop over each Gen value from 0 to 1000
for gen_value in range(1000):
    # Filter data corresponding to this specific Gen value
    gen_data = combined_df[combined_df["Gen (#)"] == gen_value]
    
    # Loop over each i (from 0 to the number of NPC columns)
    for i in range(len(npc_cols)):
        # Extract NPC_i and GHE_i for the current i
        npc_i = gen_data[npc_cols[i]].values
        ghe_i = gen_data[ghe_cols[i]].values
        
        # Stack NPC and GHE values as a cost matrix
        costs = np.column_stack((npc_i, ghe_i))  # Shape (n, 2)
        
        # Determine Pareto-efficient points
        pareto_mask = is_pareto_efficient(costs)
        
        # Retrieve the Pareto-optimal solutions for the current (Gen, i)
        pareto_solutions = gen_data.iloc[pareto_mask][[npc_cols[i], ghe_cols[i]]]
        
        # Store the Pareto solutions in the dictionary using the (gen_value, i) as the key
        pareto_solutions_dict[(gen_value, i)] = pareto_solutions
        
        # calculating the HV for each Gen(#) and seed_no
        HyperVal[gen_value,i] = ind_HV(np.array(pareto_solutions))
        GDVal[gen_value,i] = ind_GD(np.array(pareto_solutions))
        GDplusVal[gen_value,i] = ind_GDplus(np.array(pareto_solutions))
        IGDVal[gen_value,i] = ind_IGD(np.array(pareto_solutions))
        IGDplusVal[gen_value,i] = ind_IGDplus(np.array(pareto_solutions))
        

HV_Avg=np.mean(HyperVal[2:-1],axis=1)
GD_Avg=np.mean(GDVal[2:-1],axis=1)
GDplus_Avg=np.mean(GDplusVal[2:-1],axis=1)
IGD_Avg=np.mean(IGDVal[2:-1],axis=1)
IGDplus_Avg=np.mean(IGDplusVal[2:-1],axis=1)



