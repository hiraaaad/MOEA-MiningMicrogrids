
import numpy as np
import pandas as pd
from input_data.MG_Input_Data import *
from BLAST_Python.nmc111_gr_Kokam75Ah_2017 import Nmc111_Gr_Kokam75Ah_Battery

cell = Nmc111_Gr_Kokam75Ah_Battery()


class PV_cost_model():

    def __init__(self, P_nom_AC):
        self.P_nom_DC = P_nom_AC*R_DCtAC/eta_inv    #P_nom_DC = Nominal PV system capacity calculated based on AC power inputed by user

    def get_C_mod (self):
        N_mod = np.ceil(self.P_nom_DC/P_mod)        #N_mod = Number of PV modules
        C_mod = c_mod * N_mod                       #c_mod = Cost of one PV module / C_mod = Cost of PV modules
        return C_mod

    def get_C_inv(self):
        P_nom_AC = eta_inv * (self.P_nom_DC/R_DCtAC)     #P_nom_AC = Total nominal AC power output / eta_inv = Inverter efficieny / R_DCtAC = DC to AC ratio of PV system
        P_inv = S_inv * pf                               #P_inv = Nominal inverter AC power output / S_inv = Apparant power of off-grid inverter / pf = power factor
        N_inv = np.ceil(P_nom_AC/P_inv)                  #N_inv = Number of inverters
        C_inv = c_inv * N_inv                            #C_inv = Cost of inverters / #c_inv = Cost of one inverter
        return C_inv

    def get_C_land(self):
        c_land = c_uland + c_uprot                  #c_land = Unit cost of land area and lightning protection / c_uland = Unit cost of land area / c_uprot = Unit cost of lightning protection
        N_mod = np.ceil(self.P_nom_DC / P_mod)      #A_Lmod = Land area covered by the modules
        A_Lmod = N_mod * L_mod * W_mod              #L_mod = Length of PV module #W_mod = Width of PV module
        A_land = alpha_land * A_Lmod                #A_land = Land area of PV plant
        C_land = c_land * A_land                    #C_land = Cost of land area
        return C_land

    def get_C_struct(self):
        C_struct = c_struct * self.P_nom_DC              #C_struct = Total installation cost of support structure / c_struct = Cost of support structure per nominal DC power
        return C_struct

    def get_C_tsys(self):
        C_tsys = c_tsys * self.P_nom_DC                  #C_tsys = Cost of tracking system / c_tsys = Cost of tracking system per nominal DC power
        return C_tsys

    def get_C_DCcable(self):
        l_DCcable = l_uDCcable * (self.P_nom_DC/1000)    #l_uDCcable = Length of DC cables per MW
        C_DCcable = l_DCcable * c_uDCcable               #l_DCcable = Total length of DC cable / c_uDCcable = Cost of DC cable per meter
        return C_DCcable

    def get_C_ACcable(self):
        l_ACcable = l_uACcable * (self.P_nom_DC / 1000)  #l_uACcable = Length of AC cables per MW
        C_ACcable = l_ACcable * c_uACcable               #l_ACcable = Total length of AC cable / c_uACcable = Cost of AC cable per meter
        return C_ACcable

    def get_C_tr(self):
        P_nom_AC = eta_inv * (self.P_nom_DC/R_DCtAC)
        C_tr = c_tr * P_nom_AC                           #C_tr = Cost of AC components / c_tr = Unit cost of AC component per nominal AC power
        return C_tr

    def get_C_PVcapex (self):
        C_PVcapex = self.get_C_mod() + self.get_C_inv() + self.get_C_land() + self.get_C_struct() \
                    + self.get_C_tsys() + self.get_C_DCcable() + self.get_C_ACcable() + self.get_C_tr()         #C_PVcapex = Total installation capital cost of PV plant
        return C_PVcapex

    def get_C_PVopex(self):                                                                                     #L_pv = Life of PV system / alpha_opex = Operation and maintanance cost factor / intr = Interest rate
        y = np.arange(1, L_mg)                                                                               #L_pv = Life of PV system / alpha_opex = Operation and maintanance cost factor / intr = Interest rate
        C_PVopex1 = np.sum((alpha_opex * self.get_C_PVcapex())/(1+intr)**(y))                                   #C_PVopex = Operation and maintenance cost / y = Year count
        C_PVopex2 = np.sum((gamma_ins * self.get_C_PVcapex())/(1+intr)**(y))                                    #gamma_ins = Insurance cost factor
        C_PVopex = C_PVopex1 + C_PVopex2
        return C_PVopex

    def get_C_PVmodrep(self):
        y = np.arange(Wr_mod + 1, L_mg)                                                                      # Wr_mod = Module warranty period
        C_PVmodrep = np.sum((f_pv * self.get_C_mod()) / (1 + intr) ** (y))                                      # f_pv = Failure rate of PV module
        return C_PVmodrep                                                                                       # C_PVinvrep = Inverter replacement cost

    def get_C_PVinvrep(self):
        y = np.arange(1, int(np.floor(L_mg / L_inv) + 1))                                                    # L_inv = Life of inverter
        C_PVinvrep = np.sum(self.get_C_inv() / ((1 + intr) ** (y * L_inv)))                                     #C_PVmodrep = PV module replacement cost
        return C_PVinvrep

    def get_C_PVrep(self):
        C_PVrep = self.get_C_PVmodrep() + self.get_C_PVinvrep()                                                 #C_PVRep = Replacement cost
        return C_PVrep

    def get_C_modSAL (self):
        d = 1/(1+intr)**L_mg                                                                                    #d = Discount factor
        C_modSAL = self.get_C_mod() * (1-(2/L_pv))**L_mg * d                                                    #C_modSAL = Salvage value of PV modules
        return C_modSAL

    def get_C_invSAL (self):                                                                                    #L_uinv = Number of years last inverter has been used
        d = 1 / (1 + intr) ** L_mg
        C_invSAL = self.get_C_inv() * (1-(2/L_inv))**(L_mg%L_inv) * d                                        #C_invSAL = Salvage value of inverters
        return C_invSAL

    def get_C_landSAL(self):                                                                                    #This function is only needed when the land is bought or owned
        C_landSAL1 = self.get_C_land()
        y = np.arange(1, L_mg+1)
        C_landSAL2 = np.sum((gamma_land * self.get_C_land())/(1+intr)**(y-1))                                   #gamma_land = Yearly increment factor of land price
        C_landSAL = C_landSAL1 + C_landSAL2
        return C_landSAL

    def get_C_structSAL(self):
        d = 1 / (1 + intr) ** L_mg
        C_structSAL = self.get_C_struct() * ((1-(2/L_struct))**L_mg) * d                                        #L_struct = Estimated life of support structure
        return C_structSAL

    def get_C_DCcabSAL(self):
        d = 1 / (1 + intr) ** L_mg
        C_DCcabSAL = self.get_C_DCcable() * ((1-(2/L_DCcab))**L_mg) * d                                         #L_DCcab = Estimated life of DC cables
        return C_DCcabSAL

    def get_C_ACcabSAL(self):
        d = 1 / (1 + intr) ** L_mg
        C_ACcabSAL = self.get_C_ACcable() * ((1-(2/L_ACcab))**L_mg) * d                                         #L_ACcab = Estimated life of AC cables
        return C_ACcabSAL

    def get_C_trSAL(self):
        d = 1 / (1 + intr) ** L_mg
        C_trSAL = self.get_C_tr() * ((1-(2/L_tr))**L_mg) * d                                                    #L_tr = Estimated average life of AC components
        return C_trSAL

    def get_C_tsysSAL(self):
        d = 1 / (1 + intr) ** L_mg
        C_tsysSAL = self.get_C_tsys() * ((1 - (2 / L_tsys)) ** L_mg) * d
        return C_tsysSAL

    def get_C_PVsal(self):
        C_PVsal = self.get_C_modSAL() + self.get_C_invSAL() + self.get_C_structSAL() + self.get_C_tsysSAL() \
                  + self.get_C_DCcabSAL() + self.get_C_ACcabSAL() + self.get_C_trSAL()                          #C_PVsal = Salvage value of PV system
        return C_PVsal                                                                                          #Salvage value of land is not considered here, add when necessary


################ Wind Farm Cost model##############################
class Wind_cost_model():
    def __init__(self, P_WT):
        self.P_WT = P_WT                                                                                        #P_WT = Nominal wind turbine capacity input by user
        
    def get_C_rt(self):
        C_rt = self.P_WT * c_rt                                                                                 #C_rt = Cost of rotor / c_rt = Per unit cost of rotor
        return C_rt

    def get_C_nac(self):
        C_nac = self.P_WT * c_nac                                                                               #C_nac = Cost of nacelle / c_nac = Per unit cost of nacelle
        return  C_nac

    def get_C_twr(self):
        C_twr = self.P_WT * c_twr                                                                               #C_twr = Cost of tower / c_twr = Per unit cost of tower
        return C_twr

    def get_C_WTturb(self):
        C_WTturb = self.get_C_rt() + self.get_C_nac() + self.get_C_twr()                                        #C_WTturb = Cost of turbine
        return C_WTturb

    def get_C_eng(self):
        C_eng = self.P_WT * c_eng                                                                               #C_eng = Cost of engineering / c_eng = Per unit cost of engineering
        return C_eng

    def get_C_pm(self):
        C_pm = self.P_WT * c_pm                                                                                 #C_pm = Cost of project management / c_pm = Per unit cost of project management
        return C_pm

    def get_C_fdn(self):
        C_fdn = self.P_WT * c_fdn                                                                               #C_fdn = Cost of foundation / c_fdn = Per unit cost of foundation
        return C_fdn

    def get_C_ssf(self):
        C_ssf = self.P_WT * c_ssf                                                                               #C_ssf = Cost of site access / c_ssf = Per unit cost of site access
        return C_ssf

    def get_C_inst(self):
        C_inst = self.P_WT * c_inst                                                                             #C_inst = Cost of installation / c_inst = Per unit cost of installation
        return C_inst

    def get_C_ei(self):
        C_ei = self.P_WT * c_ei                                                                                 #C_ei = Cost of electrical infrastructure / c_ei = Per unit cost of electrical infrastructure
        return C_ei

    def get_C_WTbos(self):
        C_WTbos = self.get_C_eng() + self.get_C_pm() + self.get_C_fdn() + self.get_C_ssf() + self.get_C_inst() + self.get_C_ei()  #C_WTbos = Cost of balance of system
        return C_WTbos

    def get_C_ctfn(self):
        C_ctfn = self.P_WT * c_ctfn                                                                             #C_ctfn = Cost of construction and financing / c_ctfn = Per unit cost of construction and financing
        return C_ctfn

    def get_C_cont(self):
        C_cont = self.P_WT * c_cont                                                                             #C_cont = Cost of contingency / c_cont = Per unit cost of contingency
        return C_cont

    def get_C_WTconstrt(self):
        C_WTconstrt = self.get_C_ctfn() + self.get_C_cont()                                                     #C_WTconstrt = Cost of construction
        return C_WTconstrt

    def get_C_WTcapex(self):
        C_WTcapex = self.get_C_WTturb() + self.get_C_WTbos() + self.get_C_WTconstrt()                           #C_WTcapex = Total installation capital cost
        return C_WTcapex

    def get_C_WTopex(self):
        y = np.arange(1, L_mg)
        C_WTopex = np.sum((self.P_WT * c_uWTOaM) / (1 + intr) ** (y))                                           #C_WTopex = Total operation and maintenance cost / c_uWTOaM = Per unit cost of operation and maintenance
        return C_WTopex

    def get_C_WTsal(self):
        d = (1 / ((1 + intr) ** L_mg))
        C_WTsal = self.get_C_WTturb() * ((1 - (2 / L_wt)) ** L_mg) * d                                          #C_WTsal = Total salvage value of wind turbine
        return C_WTsal
    
################ Litium-Ion Battery Cost Model ####################

class LiBattery_cost_model():
    def __init__(self, P_battLi, E_battLi, DoD, total_cycle):
        self.P_battLi = P_battLi                                                                                                                #P_battLi = Battery inverter power capacity input by user
        self.E_battLi = E_battLi    
        self.DoD = DoD
        self.cycle = total_cycle
        self.Li_cost=Li_data.loc[(Li_data['htype']-E_battLi/P_battLi).abs().argsort()[:1]]
        self.Li_cost=self.Li_cost.reset_index(drop=True)


    def get_C_battECLi(self):
        C_battECLi = self.E_battLi * self.Li_cost.loc[0,'c_uECLi']
        return C_battECLi

    def get_C_battBOSLi(self):
        C_battBOSLi = self.E_battLi * self.Li_cost.loc[0,'c_uBOSLi']
        return C_battBOSLi

    def get_C_battPCSLi(self):
        C_battPCSLi = self.P_battLi * c_uPCSLi                                                                                                    #C_battPCSLi = Cost of power conversion system / c_uPCSLi = Per unit cost of power conversion system
        return C_battPCSLi

    def get_C_battCnCLi(self):
        C_battCnCLi = self.P_battLi * c_uCnCLi                                                                                                     #C_battCnCLi = Cost of control and communication / c_uCnCLi = Per unit cost of constrol and communication
        return C_battCnCLi

    def get_C_battIntLi(self):
        C_battIntLi = self.E_battLi * self.Li_cost.loc[0,'c_uIntLi']                                                                                                       #C_battIntLi = Cost of system integration
        return C_battIntLi

    def get_C_battEPCLi(self):
        C_battEPCLi = self.E_battLi * self.Li_cost.loc[0,'c_uEPCLi']                                                                                                     #C_battEPCLi = Cost of engineering, procurement and construction / c_uEPCLi = Per unit cost of EPC
        return C_battEPCLi

    def get_C_battPDLi(self):
        hr = self.E_battLi/self.P_battLi
        C_battPDLi = self.E_battLi * self.Li_cost.loc[0,'c_uPDLi']                                                                                                         #C_battPDLi = Cost of project development / c_uPDLi = Per unit cost of project development
        return C_battPDLi

    def get_C_battCAPEXLi(self):
        C_battCAPEXLi = self.get_C_battECLi() + self.get_C_battPCSLi() + self.get_C_battBOSLi() + self.get_C_battCnCLi() + self.get_C_battIntLi() + self.get_C_battEPCLi() + self.get_C_battPDLi()      #C_battCAPEXLi = Total capital cost of batteries
        return C_battCAPEXLi

    def get_C_battOPEXLi(self):
        y = np.arange(1, L_mg)
        C_battOPEXLi = np.sum((self.P_battLi * self.Li_cost.loc[0,'c_uOaMLi']) / (1 + intr) ** (y))                                                                       #c_uOaMLi = Per unit cost of operation and maintenance
        return C_battOPEXLi

    def get_C_battREPLiE(self):
        e_lossyr = beta1 * (np.exp(beta2 * (self.P_battLi/self.E_battLi))) * (self.cycle * self.DoD * E_cell)
        L_battLi = max(0.5,min(L_battcLi, (1 - EoL)/(e_lossyr/100)))                                                                                                           #EoL = End of life fraction of battery / e_lossyr = Energy loss percentage of battery
        y = np.arange(1, np.floor(L_mg/L_battLi)+1)
        C_battREPLiE = np.sum((self.get_C_battECLi()) / ((1 + intr)**(y*L_battLi)))                                         #C_battREPLiE = Replacement cost of battery energy capacity
        
        return C_battREPLiE

    def get_C_battREPLiP(self):
        y = np.arange(1, np.floor(L_mg/L_battinv)+1)
        C_battREPLiP = np.sum((self.get_C_battPCSLi()) / ((1 + intr)**(y*L_battinv)))                                                           #C_battREPLiP = Replacement cost of battery power capacity
        
        return C_battREPLiP


    def get_C_battREPLiE_high(self):
        e_lossyr = beta1 * (np.exp(beta2 * (self.P_battLi/self.E_battLi))) * (self.cycle * self.DoD * E_cell)
        L_battLi = max(0.5,min(L_battcLi, (1 - EoL)/(e_lossyr/100)))                                                                                                           #EoL = End of life fraction of battery / e_lossyr = Energy loss percentage of battery
        y = np.arange(1, np.floor(L_mg/L_battLi)+1)
        C_battREPLiE = np.sum((self.get_C_battECLi())*(1+0.01*Li_cost_proj_E.iloc[y*L_battLi-1,1]) / ((1 + intr)**(y*L_battLi)))     #C_battREPLiE_h = High forecast scnarios Replacement cost of battery energy capacity
        
        return C_battREPLiE

    def get_C_battREPLiP_high(self):
        y = np.arange(1, np.floor(L_mg/L_battinv)+1)
        C_battREPLiP = np.sum((self.get_C_battPCSLi())*(1+0.01*Li_cost_proj_P.iloc[y*L_battinv-1,1]) / ((1 + intr)**(y*L_battinv)))     #C_battREPLiE_h = High forecast scnarios Replacement cost of battery energy capacity
        
        return C_battREPLiP
    
    def get_C_battREPLiE_mid(self):
        e_lossyr = beta1 * (np.exp(beta2 * (self.P_battLi/self.E_battLi))) * (self.cycle * self.DoD * E_cell)
        L_battLi = max(0.5,min(L_battcLi, (1 - EoL)/(e_lossyr/100)))                                                                                                           #EoL = End of life fraction of battery / e_lossyr = Energy loss percentage of battery
        y = np.arange(1, np.floor(L_mg/L_battLi)+1)
        C_battREPLiE = np.sum((self.get_C_battECLi())*(1+0.01*Li_cost_proj_E.iloc[y*L_battLi-1,2]) / ((1 + intr)**(y*L_battLi)))     #C_battREPLiE_h = High forecast scnarios Replacement cost of battery energy capacity
        
        return C_battREPLiE

    def get_C_battREPLiP_mid(self):
        y = np.arange(1, np.floor(L_mg/L_battinv)+1)
        C_battREPLiP = np.sum((self.get_C_battPCSLi())*(1+0.01*Li_cost_proj_P.iloc[y*L_battinv-1,2]) / ((1 + intr)**(y*L_battinv)))     #C_battREPLiE_h = High forecast scnarios Replacement cost of battery energy capacity
        
        return C_battREPLiP

    def get_C_battREPLiE_low(self):
        e_lossyr = beta1 * (np.exp(beta2 * (self.P_battLi/self.E_battLi))) * (self.cycle * self.DoD * E_cell)
        L_battLi = max(0.5,min(L_battcLi, (1 - EoL)/(e_lossyr/100)))                                                                                                           #EoL = End of life fraction of battery / e_lossyr = Energy loss percentage of battery
        y = np.arange(1, np.floor(L_mg/L_battLi)+1)
        C_battREPLiE = np.sum((self.get_C_battECLi())*(1+0.01*Li_cost_proj_E.iloc[y*L_battLi-1,3]) / ((1 + intr)**(y*L_battLi)))     #C_battREPLiE_h = High forecast scnarios Replacement cost of battery energy capacity
        
        return C_battREPLiE

    def get_C_battREPLiP_low(self):
        y = np.arange(1, np.floor(L_mg/L_battinv)+1)
        C_battREPLiP = np.sum((self.get_C_battPCSLi())*(1+0.01*Li_cost_proj_P.iloc[y*L_battinv-1,3]) / ((1 + intr)**(y*L_battinv)))     #C_battREPLiE_h = High forecast scnarios Replacement cost of battery energy capacity
        
        return C_battREPLiP  


    def get_C_battSALLiE(self):
        e_lossyr = beta1 * (np.exp(beta2 * (self.P_battLi/self.E_battLi))) * (self.cycle * self.DoD * E_cell)                                             #L_battLi = Life of lithium-ion battery
        L_battLi = max(0.5,min(L_battcLi, (1 - EoL)/(e_lossyr/100)))  
        d = 1 / ((1 + intr) ** L_mg)                                                                                                              #d = Discount factor
        L_ubattLi = L_battLi - ((L_battLi * (np.ceil(L_mg / L_battLi))) - L_mg)
        C_battSALLiE = (self.get_C_battECLi()) * ((1 - (2/L_battLi))**L_ubattLi) * d
        C_battSALLiE = np.nan_to_num(C_battSALLiE, nan=0)
        return C_battSALLiE

    def get_C_battSALLiP(self):
        d = 1 / ((1 + intr) ** L_mg)
        L_ubattinv = L_battinv - ((L_battinv * (np.ceil(L_mg / L_battinv))) - L_mg)
        C_battSALLiP1 = (self.get_C_battPCSLi()) * ((1 - (2/L_ubattinv))**L_ubattinv) * d                                                         #C_battSALLiP = Salvage value of battery power capacity
        C_battSALLiP2 = (self.get_C_battBOSLi()) * ((1 - (2/L_battBOS))**L_mg) * d
        C_battSALLiP = C_battSALLiP1 + C_battSALLiP2
        return C_battSALLiP
    
    def get_EM_battLi(self):
        e_lossyr = beta1 * (np.exp(beta2 * (self.P_battLi/self.E_battLi))) * (self.cycle * self.DoD * E_cell)                                             #L_battLi = Life of lithium-ion battery
        L_battLi = max(0.5,min(L_battcLi, (1 - EoL)/(e_lossyr/100)))
        EM_battLiI = self.E_battLi * em_battLi                                                                                                       #EM_battLi = Life cycle emission of lithium-ion battery
        EM_battLi = (EM_battLiI + np.floor(L_mg/L_battLi) * EM_battLiI) - ((np.ceil(L_mg/L_battLi) - (L_mg/L_battLi)) * EM_battLiI)
        return EM_battLi
    
################ Redox Flow Battery Cost Model ####################    

class RfBattery_cost_model():
    def __init__(self, P_battRf, E_battRf):
        self.P_battRf = P_battRf                                                                                                                #P_battRf = Battery inverter power capacity input by user
        self.E_battRf = E_battRf                                                                                                                #E_battRf = Battery energy capacity input by user
        self.Rf_cost=Rf_data.loc[(Rf_data['htype']-E_battRf/P_battRf).abs().argsort()[:1]]
        self.Rf_cost=self.Rf_cost.reset_index(drop=True)

    def get_C_battECRf(self):
        C_battECRf = self.E_battRf * self.Rf_cost.loc[0,'c_uECRf']                                                                                                   #C_battECRf = Cost of energy capacity / c_uECRf = Per unit cost of energy capacity
        return C_battECRf

    def get_C_battBOSRf(self):
        C_battBOSRf = self.E_battRf * self.Rf_cost.loc[0,'c_uBOSRf']                                                                                                    #C_battBOSRf = Cost of balance of system / c_uBOSRf = Per unit cost of balance of system
        return C_battBOSRf

    def get_C_battPCSRf(self):
        C_battPCSRf = self.P_battRf * c_uPCSRf                                                                                                    #C_battPCSRf = Cost of power conversion system / c_uPCSRf = Per unit cost of power conversion system
        return C_battPCSRf

    def get_C_battCnCRf(self):
        C_battCnCRf = self.P_battRf * c_uCnCRf                                                                                                     #C_battCnCRf = Cost of control and communication / c_uCnCRf = Per unit cost of constrol and communication
        return C_battCnCRf

    def get_C_battIntRf(self):
        C_battIntRf = self.E_battRf * self.Rf_cost.loc[0,'c_uIntRf']                                                                                                       #C_battIntRf = Cost of system integration
        return C_battIntRf

    def get_C_battEPCRf(self):
        C_battEPCRf = self.E_battRf * self.Rf_cost.loc[0,'c_uEPCRf']                                                                                                      #C_battEPCRf = Cost of engineering, procurement and construction / c_uEPCRf = Per unit cost of EPC
        return C_battEPCRf

    def get_C_battPDRf(self):
        C_battPDRf = self.E_battRf * self.Rf_cost.loc[0,'c_uPDRf']                                                                                                        #C_battPDRf = Cost of project development / c_uPDRf = Per unit cost of project development
        return C_battPDRf

    def get_C_battCAPEXRf(self):
        C_battCAPEXRf = self.get_C_battECRf() + self.get_C_battPCSRf() + self.get_C_battBOSRf() + self.get_C_battCnCRf() + self.get_C_battIntRf() + self.get_C_battEPCRf() + self.get_C_battPDRf()      #C_battCAPEXRf = Total capital cost of batteries
        return C_battCAPEXRf

    def get_C_battOPEXRf(self):
        y = np.arange(1, L_mg)
        C_battOPEXRf = np.sum((self.P_battRf * self.Rf_cost.loc[0,'c_uOaMRf']) / (1 + intr) ** (y))                                                           #c_uOaMRf = Per unit cost of operation and maintenance
        return C_battOPEXRf

    def get_C_battREPRfE(self):
        y = np.arange(1, np.floor(L_mg/L_battRf)+1)
        C_battREPRfE = np.sum((self.get_C_battECRf()) / ((1 + intr)**(y*L_battRf)))                                                             #C_battREPRfE = Replacement cost of battery energy capacity
        return C_battREPRfE

    def get_C_battREPRfP(self):
        y = np.arange(1, np.floor(L_mg/L_battinv)+1)
        C_battREPRfP = np.sum((self.get_C_battPCSRf()) / ((1 + intr)**(y*L_battinv)))                                                           #C_battREPRfP = Replacement cost of battery power capacity
        return C_battREPRfP

    def get_C_battSALRfE(self):
        d = 1 / ((1 + intr) ** L_mg)                                                                                                              #d = Discount factor / L_battRf = Life of lithium-ion battery
        L_ubattRf = L_battRf - ((L_battRf * (np.ceil(L_mg / L_battRf))) - L_mg)
        C_battSALRfE = (self.get_C_battECRf()) * ((1 - (2/L_battRf))**L_ubattRf) * d                                                              #C_battSALRfE = Salvage value of battery energy capacity
        return C_battSALRfE

    def get_C_battSALRfP(self):
        d = 1 / ((1 + intr) ** L_mg)
        L_ubattinv = L_battinv - ((L_battinv * (np.ceil(L_mg / L_battinv))) - L_mg)
        C_battSALRfP1 = (self.get_C_battPCSRf()) * ((1 - (2/L_ubattinv))**L_ubattinv) * d                                                         #C_battSALRfP = Salvage value of battery power capacity
        C_battSALRfP2 = (self.get_C_battBOSRf()) * ((1 - (2/L_battBOS))**L_mg) * d
        C_battSALRfP = C_battSALRfP1 + C_battSALRfP2
        return C_battSALRfP
    
    def get_EM_battRf(self):
        EM_battRf = self.E_battRf * em_battRf                                                                                                       #EM_battLi = Life cycle emission of redox-flow battery
        return EM_battRf

################ Gas Turbine Cost Model ####################

class GasTurbine_cost_model():
    def __init__(self, P_GT, P_GThr):
        self.P_GT = P_GT                                                                                                                    #P_GT = Gas turbine power capacity input by user
        self.P_GThr = P_GThr                                                                                                                  #LF_GT = Load factor of gas turbine

    def get_C_GTcapex(self):
        C_GTcapex = self.P_GT * c_GT                                                                                                        #c_GT = Per unit capital cost of gas turbine
        return C_GTcapex

    def get_C_GTopex(self):
        y = np.arange(1, L_mg)
        C_GTopex1 = np.sum((self.P_GT * c_MGT) / (1 + intr) ** (y))                                                                         #c_MGT = Per unit maintenance cost of gas turbine per year
        g_GThr = (A * self.P_GThr) + B
        g_GThr = g_GThr.replace(to_replace = B, value=0)                                                                                                           #A,B = Gas consumption factor
        g_GTyr = np.sum(g_GThr)* u_m3toGJ                                                                                                             #g_GThr/g_GTyr = Amount of natural gas required per hour/per year
        y = np.arange(1, L_mg)                                                                                                              #g_GT = Amount of natural gas required per year
        C_GTopex2 = np.sum((g_GTyr * c_g) / (1 + intr) ** (y))                                                                              #c_g = Per unit cost of natural gas
        C_GTopex = C_GTopex1 + C_GTopex2
        return C_GTopex

    def get_C_GTsal(self):
        d = 1 / ((1 + intr) ** L_mg)
        C_GTsal = self.get_C_GTcapex() * ((1 - (2/L_GT))**L_mg) * d                                                                           #L_GT = Life of gas turbine / L_mg = Life of microgrid
        return C_GTsal

    def get_EM_GT(self):
        g_GThr = (A * self.P_GThr) + B
        g_GThr = g_GThr.replace(to_replace = B, value=0)
        g_GTyr = np.sum(g_GThr)
        EM_GT = (g_GTyr * L_mg) * u_m3toGJ * em_GT                                                                                         #EM_GT = Emission of natural gas turbine / em_GT = Carbon emission of natural gas per GJ unit
        return EM_GT
    
    def get_EMlc_GT(self):
        Emlc_GT = self.get_EM_GT() + (np.sum(self.P_GThr)) * L_mg * emlc_GT                                                                 ## Total life cycle emission of GT 
        return Emlc_GT    

################ Diesel Generator Cost Model ####################
class DG_cost_model():
    def __init__(self, P_DG, P_Dhr):
        self.P_DG = P_DG                                                                                                                    #P_DG = Diesel generator power capacity input by user
        self.P_Dhr = P_Dhr                                                                                                                  #LF_DG = Load factor of diesel generator

    def get_C_DGcapex(self):
        C_DGcapex = self.P_DG * c_DG                                                                                                        #c_DG = Per unit capital cost of diesel generator
        return C_DGcapex

    def get_C_DGopex(self):
        y = np.arange(1, L_mg)
        C_DGopex1 = np.sum((self.P_DG * c_MDG) / (1 + intr) ** (y))                                                                           #c_MDG = Per unit maintenance cost of diesel generator per year
        l_DGhr = (C * self.P_Dhr) + D
        l_DGhr = l_DGhr.replace(to_replace = D, value=0)                                                                                                              #C,D = Diesel consumption factor
        l_DGyr = np.sum(l_DGhr)                                                                                                                #l_DGhr/l_DGyr = Amount of diesel required per hour/per year
        y = np.arange(1, L_mg)                                                                                                                #l_DG = Amount of diesel required per year
        C_DGopex2 = np.sum((l_DGyr * c_l) / (1 + intr) ** (y))                                                                                #c_l = Per unit cost of diesel
        C_DGopex = C_DGopex1 + C_DGopex2
        return C_DGopex

    def get_C_DGsal(self):
        d = 1 / ((1 + intr) ** L_mg)
        C_DGsal = self.get_C_DGcapex() * ((1 - (2/L_DG))**L_mg) * d                                                                         #L_DG = Life of diesel generator
        return C_DGsal

    def get_EM_DG(self):
        l_DGhr = (C * self.P_Dhr) + D
        l_DGhr = l_DGhr.replace(to_replace = D, value=0)
        l_DGyr = np.sum(l_DGhr)
        EM_DG = (l_DGyr * L_mg) * u_LtoGJ * em_DG                                                                                          #EM_DG = Emission of diesel generator / em_DG = Carbon emission of diesel generator per GJ unit
        return EM_DG
    
    def get_EMlc_DG(self):
        Emlc_DG = self.get_EM_DG() + (np.sum(self.P_Dhr)) * L_mg * emlc_DG                                                                 ## Total life cycle emission of GT 
        return Emlc_DG


################ Litium-Ion Battery Cost Model -- BLAST Degradation ####################

class LiBattery_cost_model_BLAST():
    def __init__(self, P_battLi, E_battLi, SOC, temp, time_sec):
        self.P_battLi = P_battLi                                                                                                                #P_battLi = Battery inverter power capacity input by user
        self.E_battLi = E_battLi    
        self.SOC = SOC
        self.Temp = temp
        self.Time = time_sec
        self.Li_cost=Li_data.loc[(Li_data['htype']-E_battLi/P_battLi).abs().argsort()[:1]]
        self.Li_cost=self.Li_cost.reset_index(drop=True)


    def get_C_battECLi(self):
        C_battECLi = self.E_battLi * self.Li_cost.loc[0,'c_uECLi']
        return C_battECLi

    def get_C_battBOSLi(self):
        C_battBOSLi = self.E_battLi * self.Li_cost.loc[0,'c_uBOSLi']
        return C_battBOSLi

    def get_C_battPCSLi(self):
        C_battPCSLi = self.P_battLi * c_uPCSLi                                                                                                    #C_battPCSLi = Cost of power conversion system / c_uPCSLi = Per unit cost of power conversion system
        return C_battPCSLi

    def get_C_battCnCLi(self):
        C_battCnCLi = self.P_battLi * c_uCnCLi                                                                                                     #C_battCnCLi = Cost of control and communication / c_uCnCLi = Per unit cost of constrol and communication
        return C_battCnCLi

    def get_C_battIntLi(self):
        C_battIntLi = self.E_battLi * self.Li_cost.loc[0,'c_uIntLi']                                                                                                       #C_battIntLi = Cost of system integration
        return C_battIntLi

    def get_C_battEPCLi(self):
        C_battEPCLi = self.E_battLi * self.Li_cost.loc[0,'c_uEPCLi']                                                                                                     #C_battEPCLi = Cost of engineering, procurement and construction / c_uEPCLi = Per unit cost of EPC
        return C_battEPCLi

    def get_C_battPDLi(self):
        hr = self.E_battLi/self.P_battLi
        C_battPDLi = self.E_battLi * self.Li_cost.loc[0,'c_uPDLi']                                                                                                         #C_battPDLi = Cost of project development / c_uPDLi = Per unit cost of project development
        return C_battPDLi

    def get_C_battCAPEXLi(self):
        C_battCAPEXLi = self.get_C_battECLi() + self.get_C_battPCSLi() + self.get_C_battBOSLi() + self.get_C_battCnCLi() + self.get_C_battIntLi() + self.get_C_battEPCLi() + self.get_C_battPDLi()      #C_battCAPEXLi = Total capital cost of batteries
        return C_battCAPEXLi

    def get_C_battOPEXLi(self):
        y = np.arange(1, L_mg)
        C_battOPEXLi = np.sum((self.P_battLi * self.Li_cost.loc[0,'c_uOaMLi']) / (1 + intr) ** (y))                                                                       #c_uOaMLi = Per unit cost of operation and maintenance
        return C_battOPEXLi

    def get_C_battREPLiE(self):
        
        cell_used = cell
        cell_used.update_battery_state(t_secs=self.Time, soc=self.SOC, T_celsius=self.Temp)
        e_lossyr = 1 - cell.outputs['q']
        L_battLi = max(0.5,min(L_battcLi, (1 - EoL)/(e_lossyr[1])))                                                                                                           #EoL = End of life fraction of battery / e_lossyr = Energy loss percentage of battery
        y = np.arange(1, np.floor(L_mg/L_battLi)+1)
        C_battREPLiE = np.sum((self.get_C_battECLi()) / ((1 + intr)**(y*L_battLi)))                                         #C_battREPLiE = Replacement cost of battery energy capacity
        
        return C_battREPLiE

    def get_C_battREPLiP(self):
        y = np.arange(1, np.floor(L_mg/L_battinv)+1)
        C_battREPLiP = np.sum((self.get_C_battPCSLi()) / ((1 + intr)**(y*L_battinv)))                                                           #C_battREPLiP = Replacement cost of battery power capacity
        
        return C_battREPLiP


    def get_C_battREPLiE_high(self):
        
        cell_used = cell
        cell_used.update_battery_state(t_secs=self.Time, soc=self.SOC, T_celsius=self.Temp)
        e_lossyr = 1 - cell.outputs['q']
        L_battLi = max(0.5,min(L_battcLi, (1 - EoL)/(e_lossyr[1])))                        
        y = np.arange(1, np.floor(L_mg/L_battLi)+1)
        C_battREPLiE = np.sum((self.get_C_battECLi())*(1+0.01*Li_cost_proj_E.iloc[y*L_battLi-1,1]) / ((1 + intr)**(y*L_battLi)))     #C_battREPLiE_h = High forecast scnarios Replacement cost of battery energy capacity
        
        return C_battREPLiE

    def get_C_battREPLiP_high(self):
        y = np.arange(1, np.floor(L_mg/L_battinv)+1)
        C_battREPLiP = np.sum((self.get_C_battPCSLi())*(1+0.01*Li_cost_proj_P.iloc[y*L_battinv-1,1]) / ((1 + intr)**(y*L_battinv)))     #C_battREPLiE_h = High forecast scnarios Replacement cost of battery energy capacity
        
        return C_battREPLiP
    
    def get_C_battREPLiE_mid(self):
        
        cell_used = cell
        cell_used.update_battery_state(t_secs=self.Time, soc=self.SOC, T_celsius=self.Temp)
        e_lossyr = 1 - cell.outputs['q']
        L_battLi = max(0.5,min(L_battcLi, (1 - EoL)/(e_lossyr[1])))                        
        y = np.arange(1, np.floor(L_mg/L_battLi)+1)
        C_battREPLiE = np.sum((self.get_C_battECLi())*(1+0.01*Li_cost_proj_E.iloc[y*L_battLi-1,2]) / ((1 + intr)**(y*L_battLi)))     #C_battREPLiE_h = High forecast scnarios Replacement cost of battery energy capacity
        
        return C_battREPLiE

    def get_C_battREPLiP_mid(self):
        y = np.arange(1, np.floor(L_mg/L_battinv)+1)
        C_battREPLiP = np.sum((self.get_C_battPCSLi())*(1+0.01*Li_cost_proj_P.iloc[y*L_battinv-1,2]) / ((1 + intr)**(y*L_battinv)))     #C_battREPLiE_h = High forecast scnarios Replacement cost of battery energy capacity
        
        return C_battREPLiP

    def get_C_battREPLiE_low(self):
        
        cell_used = cell
        cell_used.update_battery_state(t_secs=self.Time, soc=self.SOC, T_celsius=self.Temp)
        e_lossyr = 1 - cell.outputs['q']
        L_battLi = max(0.5,min(L_battcLi, (1 - EoL)/(e_lossyr[1])))                        
        y = np.arange(1, np.floor(L_mg/L_battLi)+1)
        C_battREPLiE = np.sum((self.get_C_battECLi())*(1+0.01*Li_cost_proj_E.iloc[y*L_battLi-1,3]) / ((1 + intr)**(y*L_battLi)))     #C_battREPLiE_h = High forecast scnarios Replacement cost of battery energy capacity
        
        return C_battREPLiE

    def get_C_battREPLiP_low(self):
        y = np.arange(1, np.floor(L_mg/L_battinv)+1)
        C_battREPLiP = np.sum((self.get_C_battPCSLi())*(1+0.01*Li_cost_proj_P.iloc[y*L_battinv-1,3]) / ((1 + intr)**(y*L_battinv)))     #C_battREPLiE_h = High forecast scnarios Replacement cost of battery energy capacity
        
        return C_battREPLiP  


    def get_C_battSALLiE(self):
        
        cell_used = cell
        cell_used.update_battery_state(t_secs=self.Time, soc=self.SOC, T_celsius=self.Temp)
        e_lossyr = 1 - cell.outputs['q']
        L_battLi = max(0.5,min(L_battcLi, (1 - EoL)/(e_lossyr[1])))                        
        d = 1 / ((1 + intr) ** L_mg)                                                                                                              #d = Discount factor
        L_ubattLi = L_battLi - ((L_battLi * (np.ceil(L_mg / L_battLi))) - L_mg)
        C_battSALLiE = (self.get_C_battECLi()) * ((1 - (2/L_battLi))**L_ubattLi) * d
        C_battSALLiE = np.nan_to_num(C_battSALLiE, nan=0)
        return C_battSALLiE

    def get_C_battSALLiP(self):
        d = 1 / ((1 + intr) ** L_mg)
        L_ubattinv = L_battinv - ((L_battinv * (np.ceil(L_mg / L_battinv))) - L_mg)
        C_battSALLiP1 = (self.get_C_battPCSLi()) * ((1 - (2/L_ubattinv))**L_ubattinv) * d                                                         #C_battSALLiP = Salvage value of battery power capacity
        C_battSALLiP2 = (self.get_C_battBOSLi()) * ((1 - (2/L_battBOS))**L_mg) * d
        C_battSALLiP = C_battSALLiP1 + C_battSALLiP2
        return C_battSALLiP
    
    def get_EM_battLi(self):
        
        cell_used = cell
        cell_used.update_battery_state(t_secs=self.Time, soc=self.SOC, T_celsius=self.Temp)
        e_lossyr = 1 - cell.outputs['q']
        L_battLi = max(0.5,min(L_battcLi, (1 - EoL)/(e_lossyr[1])))                        
        EM_battLiI = self.E_battLi * em_battLi                                                                                                       #EM_battLi = Life cycle emission of lithium-ion battery
        EM_battLi = (EM_battLiI + np.floor(L_mg/L_battLi) * EM_battLiI) - ((np.ceil(L_mg/L_battLi) - (L_mg/L_battLi)) * EM_battLiI)
        return EM_battLi
