

'''
The Multi-objective main code for MG optimization using BLAST_lite

'''

import pandas as pd
import numpy as np
# import argparse

from SB_My_Fun import MG_Model_rulebased
from input_data.MG_Input_Data import L_mg

from pymoo.algorithms.moo.sms import SMSEMOA
from pymoo.algorithms.moo.age import AGEMOEA
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.moo.age2 import AGEMOEA2
from pymoo.algorithms.moo.ctaea import CTAEA
from pymoo.algorithms.moo.unsga3 import UNSGA3
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem
from pymoo.termination import get_termination


import warnings
warnings.filterwarnings("ignore")

### Function for running the scheduling problem

def get_Results_Run(df, DT, P_WT, P_PV, P_bat, E_bat, P_FF):
    
    MG_Object = MG_Model_rulebased(df=df, DT=DT, P_WT=P_WT, P_PV=P_PV, P_bat = P_bat, E_bat = E_bat, P_FF=P_FF)  
        
    output = dict(); 
    output = {
        'ENS':0,'ERC':0,'NPC':0,'GHE':0
    }
    
    if P_bat==0:
        
        Res_MG = MG_Object.get_Result_noBat() # Scheduling problem for noBAT case
        output['NPC'] = round(Res_MG[0],4)    ## NPC for DG+Li Million AUD
        output['ERC'] = round(Res_MG[1],4)    ## RES Curtailment in GWh/yr
        output['ENS'] = round(Res_MG[2],0)    ## ENergy not served in GWh/yr
        output['GHE'] = round(Res_MG[3],4)    ## kton/yr CO2-e from DG
 
    else:
        
        Res_MG = MG_Object.get_Result_MG_rulebased_Li_BLAST() # Scheduling problem for Li+DG with BLAST
                    
        output['NPC'] = round(Res_MG[0],4)    ## NPC for DG+Li Million AUD
        output['ERC'] = round(Res_MG[1],4)    ## RES Curtailment in GWh/yr
        output['ENS'] = round(Res_MG[2],0)    ## ENergy not served in GWh/yr
        output['GHE'] = round(Res_MG[3],4)    ## kton/yr CO2-e from DG

    return output


#### Input Data

df = pd.read_csv("input_data//RES_Power_IGO_hourly.csv")

t1 = pd.to_datetime(df.loc[0,'timestamp'])
t2 = pd.to_datetime(df.loc[1,'timestamp'])
DT = pd.Timedelta(t2 - t1).seconds/3600


PL = np.ceil(df["PD"].max()) ##MW
ENL = DT*df["PD"].sum()/1e3 ## Total annual load energy GWh/yr

# RES = get_Results_Run(df, DT, P_WT=0, P_PV=0, P_bat=0, E_bat=0, P_FF=11)


## Using the NSGAII in PYMOO

class MyProblem(ElementwiseProblem):

    def __init__(self):

        super().__init__(n_var=5,
                          n_obj=3,
                          n_ieq_constr=3,
                          xl=np.array([0, 0, 0, 0, 0]),
                          xu=np.array([3*PL, 3*PL, 3*PL, 3*PL, PL]))

    def _evaluate(self, x, out, *args, **kwargs):
        P_WT = x[0]
        P_PV = x[1]
        P_bat = x[2]
        E_bat = x[3]
        P_FF = x[4]
        
        RES = get_Results_Run(df, DT, P_WT, P_PV, P_bat, E_bat, P_FF)
        
        g1 = x[3] - 4*x[2]
        g2 = 0.5*x[2] - x[3]
        g3 = RES['ENS']/ENL - 0.0 ## ENS should be less than 5%
        # out["F"] = [RES['NPC'], RES['ERC'], RES['ENS'], RES['GHE']]
        out["F"] = [RES['NPC'], RES['GHE'], RES['ERC']]
        out["G"] = [g3, g2, g1]
        
problem = MyProblem()


####

class MyProblem_noBAT(ElementwiseProblem):

    def __init__(self):

        super().__init__(n_var=3,
                          n_obj=3,
                          n_ieq_constr=1,
                          xl=np.array([0, 0, 0]),
                          xu=np.array([3*PL, 3*PL, PL]))

    def _evaluate(self, x, out, *args, **kwargs):
        P_WT = x[0]
        P_PV = x[1]
        P_FF = x[2]
        
        RES = get_Results_Run(df, DT, P_WT, P_PV, 0, 0, P_FF)
        
        g3 = RES['ENS']/ENL - 0.0 ## ENS should be less than 0%
        # out["F"] = [RES['NPC'], RES['ERC'], RES['ENS'], RES['GHE']]
        out["F"] = [RES['NPC'], RES['GHE'], RES['ERC']]
        out["G"] = [g3] 

problem_1 = MyProblem_noBAT()



# F, G = problem.evaluate(X=np.array([0,0,0,0,12]))



def res_2_dataframe(res):    
    xval = res.X
    fval = res.F
    ENS_val = res.G
    
    if DesVar==1:
        sol = pd.DataFrame({'PWT (MW)':xval[:,0],
                        'PV (MW)':xval[:,1],
                        'PBAT (MW)':0,
                        'EBAT (MWh)':0,
                        'PDG (MW)': xval[:,2],
                        'NPC ($/kWh)':fval[:,0]/(ENL*L_mg),
                        'GHE (kg CO2/kWh)':fval[:,1]/ENL,
                        'ERC (GWh/yr)':fval[:,2],
                        'ENS (%)':ENS_val[:,0]})
    else:
        sol = pd.DataFrame({'PWT (MW)':xval[:,0],
                        'PV (MW)':xval[:,1],
                        'PBAT (MW)':xval[:,2],
                        'EBAT (MWh)':xval[:,3],
                        'PDG (MW)': xval[:,4],
                        'NPC ($/kWh)':fval[:,0]/(ENL*L_mg),
                        'GHE (kg CO2/kWh)':fval[:,1]/ENL,
                        'ERC (GWh/yr)':fval[:,2],
                        'ENS (%)':ENS_val[:,0]})
    
    hist = res.history
    
    for ii in range(len(hist)):
        a = hist[ii]
        b = a.pop
        x_vg = b.get('X')
        f_vg = b.get('F')
        g_vg = b.get('G')
        if ii==0:
            if DesVar==1:
                sol_g = pd.DataFrame({            
                    'Gen (#)':[ii]*nPop,
                    'PWT (MW)':x_vg[:,0],
                    'PV (MW)':x_vg[:,1],
                    'PBAT (MW)':0,
                    'EBAT (MWh)':0,
                    'PDG (MW)': x_vg[:,2],
                    'NPC ($/kWh)':f_vg[:,0]/(ENL*L_mg),
                    'GHE (kg CO2/kWh)':f_vg[:,1]/ENL,
                    'ERC (GWh/yr)':f_vg[:,2],
                    'ENS (%)':g_vg[:,0]})
            else:    
                sol_g = pd.DataFrame({            
                    'Gen (#)':[ii]*nPop,
                    'PWT (MW)':x_vg[:,0],
                    'PV (MW)':x_vg[:,1],
                    'PBAT (MW)':x_vg[:,2],
                    'EBAT (MWh)':x_vg[:,3],
                    'PDG (MW)': x_vg[:,4],
                    'NPC ($/kWh)':f_vg[:,0]/(ENL*L_mg),
                    'GHE (kg CO2/kWh)':f_vg[:,1]/ENL,
                    'ERC (GWh/yr)':f_vg[:,2],
                    'ENS (%)':g_vg[:,0]})
        else:
            if DesVar==1:
                c = pd.DataFrame({            
                    'Gen (#)':[ii]*nPop,
                    'PWT (MW)':x_vg[:,0],
                    'PV (MW)':x_vg[:,1],
                    'PBAT (MW)':0,
                    'EBAT (MWh)':0,
                    'PDG (MW)': x_vg[:,2],
                    'NPC ($/kWh)':f_vg[:,0]/(ENL*L_mg),
                    'GHE (kg CO2/kWh)':f_vg[:,1]/ENL,
                    'ERC (GWh/yr)':f_vg[:,2],
                    'ENS (%)':g_vg[:,0]})
            else:
                c = pd.DataFrame({            
                    'Gen (#)':[ii]*nPop,
                    'PWT (MW)':x_vg[:,0],
                    'PV (MW)':x_vg[:,1],
                    'PBAT (MW)':x_vg[:,2],
                    'EBAT (MWh)':x_vg[:,3],
                    'PDG (MW)': x_vg[:,4],
                    'NPC ($/kWh)':f_vg[:,0]/(ENL*L_mg),
                    'GHE (kg CO2/kWh)':f_vg[:,1]/ENL,
                    'ERC (GWh/yr)':f_vg[:,2],
                    'ENS (%)':g_vg[:,0]})
                    
        
            sol_g = pd.concat([sol_g,c],ignore_index=True)
    
    return (sol,sol_g)


def NSGAII_run(nPop, nFeval, DesVar, seed_no):
    algorithm = NSGA2(pop_size=nPop)
    termination = get_termination("n_eval",nFeval)

    if DesVar==1:
        res = minimize(problem_1,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)

    else:        
        res = minimize(problem,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)
    
    Sol = res_2_dataframe(res=res)
    
    return Sol

def AGEMOEA_run(nPop, nFeval, DesVar, seed_no):
    algorithm = AGEMOEA(pop_size=nPop)
    termination = get_termination("n_eval",nFeval)

    if DesVar==1:
        res = minimize(problem_1,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)

    else:        
        res = minimize(problem,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)

    Sol = res_2_dataframe(res=res)

    return Sol

def SMSEMOA_run(nPop, nFeval, DesVar, seed_no):
    algorithm = SMSEMOA(pop_size=nPop)
    termination = get_termination("n_eval",nFeval)

    if DesVar==1:
        res = minimize(problem_1,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)

    else:        
        res = minimize(problem,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)
    Sol = res_2_dataframe(res=res)
    
    return Sol

def AGEMOEA2_run(nPop, nFeval, DesVar, seed_no):
    algorithm = AGEMOEA2(pop_size=nPop)
    termination = get_termination("n_eval",nFeval)

    if DesVar==1:
        res = minimize(problem_1,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)

    else:        
        res = minimize(problem,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)

    Sol = res_2_dataframe(res=res)
    
    return Sol

def CTAEA_run(nPop, nFeval, DesVar, seed_no):
    termination = get_termination("n_eval",nFeval)

    if DesVar==1:
        ref_dirs = get_reference_directions("das-dennis", problem_1.n_obj, n_partitions=nPop, n_points=nPop)
        algorithm = CTAEA(ref_dirs=ref_dirs)
        
        res = minimize(problem_1,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)

    else:
        ref_dirs = get_reference_directions("das-dennis", problem.n_obj, n_partitions=nPop, n_points=nPop)
        algorithm = CTAEA(ref_dirs=ref_dirs)
        
        res = minimize(problem,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)

    Sol = res_2_dataframe(res=res)
    
    return Sol

def NSGA3_run(nPop, nFeval, DesVar, seed_no):
    termination = get_termination("n_eval",nFeval)

    if DesVar==1:
        ref_dirs = get_reference_directions("energy", problem_1.n_obj, n_points=nPop)
        algorithm = NSGA3(ref_dirs=ref_dirs, pop_size=nPop)
        
        res = minimize(problem_1,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)

    else:
        ref_dirs = get_reference_directions("energy", problem.n_obj, n_points=nPop)
        algorithm = NSGA3(ref_dirs=ref_dirs, pop_size=nPop)
        
        res = minimize(problem,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)

    Sol = res_2_dataframe(res=res)
    
    return Sol

def UNSGA3_run(nPop, nFeval, DesVar, seed_no):
    termination = get_termination("n_eval",nFeval)

    if DesVar==1:
        ref_dirs = get_reference_directions("energy", problem_1.n_obj, n_points=nPop)
        algorithm = UNSGA3(ref_dirs=ref_dirs, pop_size=nPop)
        
        res = minimize(problem_1,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)

    else:
        ref_dirs = get_reference_directions("energy", problem.n_obj, n_points=nPop)
        algorithm = UNSGA3(ref_dirs=ref_dirs, pop_size=nPop)
        
        res = minimize(problem,
                    algorithm,
                    termination,
                    verbose=False,
                    save_history=True,
                    seed=seed_no)

    Sol = res_2_dataframe(res=res)
    
    return Sol

nPop = 100
nFeval = 100000
DesVar = 1

seed_vect = [4,34,42,76,	87,121,146,170,200,240,270,319,344,542,617,625,706,814,
             	938,	1281,1748,1846,2206,3065,3383,3742,7253,12697,13391,15760,128765]

Alg_Flag =  1    # NSGAII: 1, AGE-MOEA: 2, SMS-EMOA: 3, AGE-MOEA2: 4, NSGA-III: 5, UNSGAIII: 6, C-TAEA: 7


for seed_no in range(len(seed_vect)):
    print('Seed_no = {}'.format(seed_no))
    if Alg_Flag==1:    
        Solution = NSGAII_run(nPop, nFeval, DesVar, seed_vect[seed_no])
        Pareto_Sol = Solution[0]
        Gen_Sol = Solution[1]
        if DesVar==1:
            Pareto_Sol.to_csv("./Results//Pareto_NSGAII_noBAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results//AllPop_NSGAII_noBAT_seedno_{}.csv".format(seed_no))
        else:
            Pareto_Sol.to_csv("./Results_BAT//Pareto_NSGAII_BAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results_BAT//AllPop_NSGAII_BAT_seedno_{}.csv".format(seed_no))

    elif Alg_Flag==2:    
        Solution = AGEMOEA_run(nPop, nFeval, DesVar, seed_vect[seed_no])
        Pareto_Sol = Solution[0]
        Gen_Sol = Solution[1]
        if DesVar==1:
            Pareto_Sol.to_csv("./Results//Pareto_AGEMOEA_noBAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results//AllPop_AGEMOEA_noBAT_seedno_{}.csv".format(seed_no))
        else:
            Pareto_Sol.to_csv("./Results_BAT//Pareto_AGEMOEA_BAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results_BAT//AllPop_AGEMOEA_BAT_seedno_{}.csv".format(seed_no))

    elif Alg_Flag==3:    
        Solution = SMSEMOA_run(nPop, nFeval, DesVar, seed_vect[seed_no])
        Pareto_Sol = Solution[0]
        Gen_Sol = Solution[1]
        if DesVar==1:
            Pareto_Sol.to_csv("./Results//Pareto_SMSEMOA_noBAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results//AllPop_SMSEMOA_noBAT_seedno_{}.csv".format(seed_no))
        else:
            Pareto_Sol.to_csv("./Results_BAT//Pareto_SMSEMOA_BAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results_BAT//AllPop_SMSEMOA_BAT_seedno_{}.csv".format(seed_no))
            
    elif Alg_Flag==4:    
        Solution = AGEMOEA2_run(nPop, nFeval, DesVar, seed_vect[seed_no])
        Pareto_Sol = Solution[0]
        Gen_Sol = Solution[1]
        if DesVar==1:
            Pareto_Sol.to_csv("./Results//Pareto_AGEMOEA2_noBAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results//AllPop_AGEMOEA2_noBAT_seedno_{}.csv".format(seed_no))
        else:
            Pareto_Sol.to_csv("./Results_BAT//Pareto_AGEMOEA2_BAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results_BAT//AllPop_AGEMOEA2_BAT_seedno_{}.csv".format(seed_no))

    elif Alg_Flag==5:    
        Solution = NSGA3_run(nPop, nFeval, DesVar, seed_vect[seed_no])
        Pareto_Sol = Solution[0]
        Gen_Sol = Solution[1]
        if DesVar==1:
            Pareto_Sol.to_csv("./Results//Pareto_NSGA3_noBAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results//AllPop_NSGA3_noBAT_seedno_{}.csv".format(seed_no))
        else:
            Pareto_Sol.to_csv("./Results_BAT//Pareto_NSGA3_BAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results_BAT//AllPop_NSGA3_BAT_seedno_{}.csv".format(seed_no))
    
    elif Alg_Flag==6:    
        Solution = UNSGA3_run(nPop, nFeval, DesVar, seed_vect[seed_no])
        Pareto_Sol = Solution[0]
        Gen_Sol = Solution[1]
        if DesVar==1:
            Pareto_Sol.to_csv("./Results//Pareto_UNSGA3_noBAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results//AllPop_UNSGA3_noBAT_seedno_{}.csv".format(seed_no))
        else:
            Pareto_Sol.to_csv("./Results_BAT//Pareto_UNSGA3_BAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results_BAT//AllPop_UNSGA3_BAT_seedno_{}.csv".format(seed_no))

    elif Alg_Flag==7:    
        Solution = CTAEA_run(nPop, nFeval, DesVar, seed_vect[seed_no])
        Pareto_Sol = Solution[0]
        Gen_Sol = Solution[1]
        if DesVar==1:
            Pareto_Sol.to_csv("./Results//Pareto_CTAEA_noBAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results//AllPop_CTAEA_noBAT_seedno_{}.csv".format(seed_no))
        else:
            Pareto_Sol.to_csv("./Results_BAT//Pareto_CTAEA_BAT_seedno_{}.csv".format(seed_no))
            Gen_Sol.to_csv("./Results_BAT//AllPop_CTAEA_BAT_seedno_{}.csv".format(seed_no))
