
import numpy as np
import pandas as pd

from SB_functions import PV_cost_model, Wind_cost_model, DG_cost_model, LiBattery_cost_model_BLAST

from input_data.MG_Input_Data import L_mg


class MG_Model_rulebased():
       
    def __init__(self, df, DT, P_WT, P_PV, P_bat, E_bat, P_FF):
        self.DT = DT
        self.PWT = P_WT
        self.PSF = P_PV
        self.PWind = df.loc[:,'NPW'].mul(P_WT)
        self.PSolar = df.loc[:,'NPV'].mul(P_PV)
        self.PLoad = df.loc[:,'PD']
        self.PBat = P_bat
        self.EBat = E_bat
        self.PFF = P_FF
        self.temp = df.loc[:,'Temp']
        

### Function for All technologies_rule based_Lithium-ion with Blast
    def get_Result_MG_rulebased_Li_BLAST(self):
        
        P_nom_PV = self.PSF*1000     # Installed capacity of PV in kW-AC
        P_nom_WT = self.PWT*1000     # Installed capacity of WT in kW
        PV_cost = PV_cost_model(P_nom_PV)
        WF_cost = Wind_cost_model(P_nom_WT)
        
        eta_ch =0.95
        eta_dis = 0.95
        SOC_Max = 0.95
        SOC_Min = 0.10
        SOC_ini = 0.5

        Pload=self.PLoad
        PW= [0]*len(Pload)
        PS = [0]*len(Pload)
        PCh = [0]*len(Pload)
        PDis = [0]*len(Pload)
        PNG = [0]*len(Pload)
        Pnet = [0]*len(Pload)
        SOC = [0]*len(Pload)
        P_cur = [0]*len(Pload)
        P_NS = [0]*len(Pload)
        time_array = [0]*len(Pload)

        for t in range(len(Pload)):
            
            time_array[t] = round(t*self.DT*3600)
            PW[t]=self.PWind[t]
            PS[t]=self.PSolar[t]
            Pnet[t] = Pload[t] - PS[t] - PW[t]
            
            if Pnet[t]>0:
                PCh[t]=0
                PDis[t]=min(self.PBat,Pnet[t])
            elif Pnet[t]<0:
                PDis[t]=0
                PCh[t]=-max(-self.PBat, Pnet[t])
            else:
                PCh[t]=0
                PDis[t]=0
                
            if t==0:
                SOC[t]=SOC_ini*self.EBat + self.DT*PCh[t]*eta_ch - self.DT*PDis[t]/eta_dis
                if SOC[t]>=SOC_Max*self.EBat:
                    SOC[t]=SOC_Max*self.EBat
                    PCh[t]=(SOC_Max-SOC_ini)*self.EBat/(self.DT*eta_ch)
                elif SOC[t]<=SOC_Min*self.EBat:
                    SOC[t]=SOC_Min*self.EBat
                    PDis[t]=(SOC_ini-SOC_Min)*self.EBat*eta_dis/self.DT
            else:
                SOC[t]=SOC[t-1] + self.DT*PCh[t]*eta_ch - self.DT*PDis[t]/eta_dis
                if SOC[t]>=SOC_Max*self.EBat:
                    SOC[t]=SOC_Max*self.EBat
                    PCh[t]=(SOC_Max*self.EBat-SOC[t-1])/(self.DT*eta_ch)
                elif SOC[t]<=SOC_Min*self.EBat:
                    SOC[t]=SOC_Min*self.EBat
                    PDis[t]=(SOC[t-1]-SOC_Min*self.EBat)*eta_dis/self.DT
            
            Pnet[t] = Pload[t] - PS[t] - PW[t] - PDis[t] + PCh[t]
            
            if Pnet[t]>=0:
                PNG[t]=Pnet[t]
                P_cur[t]=0
                if PNG[t]>self.PFF:
                    PNG[t] = self.PFF
                    P_NS[t] = Pnet[t] - self.PFF
            else:
                PNG[t] = 0
                P_cur[t] = -Pnet[t]
                
       
        P_DG = self.PFF*1000  ## Capcity of NG unit [kw]
        P_DGt = self.DT*pd.Series(PNG)*1000    ##Intervally production of NG generators unit [kwh]
        
        DE_cost = DG_cost_model(P_DG, P_DGt)
        
        SOC_Val = np.array(SOC)
        time_sec = np.array(time_array)
        Temp_C = np.array(self.temp)
        
        if self.EBat!=0:
        
            P_battLi = self.PBat*1000
            E_battLi = self.EBat*1000
                    
            Storage_cost_Li = LiBattery_cost_model_BLAST(P_battLi, E_battLi, SOC_Val, Temp_C, time_sec)
    
            NPC_DG_Li = (
                DE_cost.get_C_DGcapex()+DE_cost.get_C_DGopex()-DE_cost.get_C_DGsal()+
                PV_cost.get_C_PVcapex()+PV_cost.get_C_PVopex()+PV_cost.get_C_PVrep()-PV_cost.get_C_PVsal()+
                WF_cost.get_C_WTcapex()+WF_cost.get_C_WTopex()-WF_cost.get_C_WTsal()+
                Storage_cost_Li.get_C_battCAPEXLi()+Storage_cost_Li.get_C_battOPEXLi()+Storage_cost_Li.get_C_battREPLiP()+
                Storage_cost_Li.get_C_battREPLiE()-Storage_cost_Li.get_C_battSALLiE()-Storage_cost_Li.get_C_battSALLiP()
                )/1e6
        
        else:
                
            NPC_DG_Li = (
                DE_cost.get_C_DGcapex()+DE_cost.get_C_DGopex()-DE_cost.get_C_DGsal()+
                PV_cost.get_C_PVcapex()+PV_cost.get_C_PVopex()+PV_cost.get_C_PVrep()-PV_cost.get_C_PVsal()+
                WF_cost.get_C_WTcapex()+WF_cost.get_C_WTopex()-WF_cost.get_C_WTsal()
                )/1e6
##------------------------------------------------------        
        
        EC = self.DT*sum(P_cur)/1e3 ## RES Curtailment in GWh
        ENS = self.DT*sum(P_NS)/1e3 ## ENS in GWh/yr
        EM_DG = DE_cost.get_EM_DG()/(1e6*L_mg)                  ## kton/yr CO2-e from DG
        

        return (NPC_DG_Li, EC, ENS, EM_DG)


### Function for no Battery 
    def get_Result_noBat(self):        
        
        P_nom_PV = self.PSF*1000     # Installed capacity of PV in kW-AC
        P_nom_WT = self.PWT*1000     # Installed capacity of WT in kW
        PV_cost = PV_cost_model(P_nom_PV)
        WF_cost = Wind_cost_model(P_nom_WT)
       
            
        PNET = self.PLoad-self.PWind-self.PSolar
        PNET[PNET<0]=0
        
        P_DG = self.PFF*1000  ## Capcity of NG unit [kw]
        P_DGt = self.DT*np.minimum(PNET,self.PFF)*1000    ## Intervally production of DE generators unit [kwh]            
            
        DE_cost = DG_cost_model(P_DG, P_DGt)
    
        Pn1 = PNET - np.minimum(PNET,self.PFF)
        
        ## Calculating NPC in M$
        
        NPC_DG = (
            PV_cost.get_C_PVcapex()+PV_cost.get_C_PVopex()+PV_cost.get_C_PVrep()-PV_cost.get_C_PVsal()+
            WF_cost.get_C_WTcapex()+WF_cost.get_C_WTopex()-WF_cost.get_C_WTsal()+
            DE_cost.get_C_DGcapex()+DE_cost.get_C_DGopex()-DE_cost.get_C_DGsal()
            )/1e6        
        
        
        Pn = self.PLoad-self.PWind-self.PSolar
        EC = self.DT*abs(Pn.mul(~Pn.gt(0)).sum().sum())/1e3         ## RES Curtailment in GWh/yr
        ENS = self.DT*(Pn1.mul(Pn1.gt(0)).sum().sum())/1e3        ## Energy not Served in GWh/yr
        EM_DG = DE_cost.get_EM_DG()/(1e6*L_mg)                  ## kton/yr CO2-e from DG
        
        
        return (NPC_DG, EC, ENS, EM_DG)
