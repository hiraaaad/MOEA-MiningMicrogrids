"""
Microgrid Input Data
"""
import pandas as pd

L_mg = 25           #Life of microgrid
intr = 0.031        #interest rate

### --------------------PV System Data-------------------------------------------

##PV System Components life
L_pv = 25       #life of PV system [years]
L_inv = 10      #life of inverter [years]
L_uinv = 5      #number of years the last inverter is used [years]
L_struct = 60   #estimated life of support structures [years]
L_tsys = 25     #estimated life of tracking system [years]
L_DCcab = 25    #estimated life of DC cable [years]
L_ACcab = 25    #estimated life of AC cable [years]
L_tr = 40       #estimated average life of AC components [years]

##PV Module parameters
P_mod = 0.415   #power of each module [kW]
c_mod = 117     #cost of each module [AUD]
L_mod = 1.77    #length of each module [meter]
W_mod = 1.09    #width of each module [meter]
Wr_mod = 12     #warranty period of the module [years]
T_coef = -0.3   #Temperature coeficient [%]

## PV Inverter parameters
S_inv = 15      #apparant power of inverter [kVA]
c_inv = 6648    #cost of each inverter [AUD]

##Fixed unit cost parameters of PV System
c_uland = 0.0638    #unit cost of land area for lease [AUD/m2]
c_uprot = 9.6       #unit cost of lightning protection [AUD/m2]
c_struct = 126.8    #cost of support structure per nominal DC power [AUD/kWp]
c_tsys = 400        #cost of tracking system [AUD/kWp]
c_uDCcable = 3.5    #cost of DC cables per meter [AUD/meter]
c_uACcable = 3.7    #cost of AC cables per meter [AUD/meter]
c_tr = 66.3         #cost of AC component per nominal AC power [AUD/kW]

##Constant input parameters of PV Cost Model
eta_inv = 0.95      #inverter effciency
pf = 0.9            #power factor
R_DCtAC = 1.3       #DC to AC ratio of PV system
alpha_land = 1.8    #land area extension factor
alpha_opex = 0.02   #operation and maintenance cost factor
gamma_ins = 0.005   #insurance cost factor
f_pv = 0.0005       #PV module failure rate
gamma_land = 0.08   #land price yearly increment factor
T_mod = 75          #degree celsius
l_uDCcable = 7000   #meters/MW
l_uACcable = 2200   #meters/MW

##Life cycle emission of PV system
em_PV = 0.0456      #Life cycle emission of solar PV system [kg C02-e/kWh]


### --------------------Wind System Data-------------------------------------------

##Wind Turbine Components life
L_wt = 25   #Life of turbine

##Wind turbine parameters
h = 110     #height of the wind turbine [m]
vci = 2.5   #cut-in speed [m/s] 
vco = 20    #cut-out speed [m/s]
vr = 10.2   #rated speed [m/s]
pr = 3.57   #rated power [MW]
z0_BHP = 0.03   #the surface roughness length [m]
z0_IGO = 0.80   #the surface roughness length [m]
                #see here for more detail. https://www.homerenergy.com/products/pro/docs/3.11/wind_resource_variation_with_height.html
rho_0 = 1.225 # Air density in standard condition[kg/m3]                
R = 287.05    # Ideal gas constant [J/m3.K]  
T0 = 273.15   # Kelvin Degree

##Fixed unit cost parameters of wind turbine
c_rt = 470      #per unit cost of rotor [AUD/kW]
c_nac = 768     #per unit cost of nacelle [AUD/kW]
c_twr = 306     #per unit cost of tower [AUD/kW]
c_eng = 34.5    #per unit cost of engineering [AUD/kW]
c_pm = 15       #per unit cost of project management [AUD/kW]
c_fdn = 112.5   #per unit cost of foundation [AUD/kW]
c_ssf = 60      #per unit cost of site access, staging and facilities [AUD/kW]
c_inst = 61.5   #per unit cost of assembly and installation [AUD/kW]
c_ei = 198      #per unit cost of electrical infrastructure [AUD/kW]
c_ctfn = 34.5   #per unit cost of construction financing [AUD/kW]
c_cont = 135    #per unit cost of contingency [AUD/kW]
c_uWTOaM = 60   #per unit cost of operation and maintenance per year [AUD/kW]

#Constant input parameters of Wind Turbine Cost Model
beta_OaM = 0.04 #operation and maintenance cost factor
beta_ins = 0.01 #insurance cost factor


##Life cycle emission of Wind Turbine
em_WT = 0.015      #Life cycle emission of solar Wind Turbine [kg C02-e/kWh]


### --------------------Battery Energy Storage Data-------------------------------------------

Li_data = pd.DataFrame({'htype' : [10,4,2,1,0.5],
          'c_uECLi' :[268.2, 273.4 , 277.41, 313.47, 341.68],
          'c_uBOSLi':[60.00, 63.57 , 68.30 , 77.17 , 84.11],
          'c_uIntLi':[68.80, 75.24 , 84.84 , 95.86 , 104.5],
          'c_uEPCLi':[83.16, 91.8  , 104.8 , 118.44, 129.09],
          'c_uPDLi' :[99.80, 110.14, 125.77, 142.12, 154.91],
          'c_uOaMLi':[15.88, 7.57  , 4.75  , 5.36  , 5.85]})

Rf_data = pd.DataFrame({'htype' : [10,4,2,1,0.5],
          'c_uECRf' :[331.4, 414.88, 554.00, 626.00, 682.36],
          'c_uBOSRf':[66.28, 83.00 , 110.80, 125.20, 136.40],
          'c_uIntRf':[63.60, 84.50 , 119.38, 179.07, 195.18],
          'c_uEPCRf':[73.75, 98.31 , 139.50, 157.68, 171.80],
          'c_uPDRf' :[84.60, 113.07, 160.40, 181.25, 197.56],
          'c_uOaMRf':[19.75, 11.25 , 8.40  , 9.50  , 10.35]})

## Other cost terms of Li and Rf
c_uPCSLi = 127      #per unit cost of power conversion system [AUD/kW]
c_uCnCLi = 60       #per unit cost of control and communication [AUD/kW]

c_uPCSRf = 232.3    #per unit cost of power conversion system [AUD/kW]
c_uCnCRf = 60       #per unit cost of control and communication [AUD/kW]

  
##Lithium-ion battery parameters
EoL = 0.8           #end of life fraction of battery
# beta1 = 0.0045      #Pre-exponential factor of temperature function in 46 degC
# beta2 = 0.1826      #Exponential factor of temperature function in 46 degC
beta1 = 0.0008      #Pre-exponential factor of temperature function in 20 degC
beta2 = 0.3903      #Exponential factor of temperature function in 20 degC
E_cell = 1.5        #Energy capacity of each cell of the batteries [kWh]
eta_ch_Li =0.95
eta_dis_Li = 0.95

eta_ch_Rf =0.85
eta_dis_Rf = 0.85


## Lithium-ion Cost projections

Li_cost_proj_E = pd.read_excel('input_data//Cost Projection Values.xlsx', sheet_name='Energy Components',usecols= 'M:Q',nrows=29)
Li_cost_proj_P = pd.read_excel('input_data//Cost Projection Values.xlsx', sheet_name='Power Components',usecols= 'M:Q',nrows=29)


##Battery life
L_battcLi = 12      #Calendar life of lithium ion batteries
L_battRf = 20       #Life of Redox flow battery
L_battinv = 10      #Life of battery inverter capacity
L_battBOS = 40      #Life of balance of system


##Life cycle emission
em_battLi = 120     #Life cycle carbon emission of lithium-ion battery [kg CO2-e/kWh]
em_battRf = 72      #Life cycle carbon emission of redox-flow battery [kg C02-e/kWh]


### --------------------Back-up Generators Data-------------------------------------------

#Fixed unit cost parameters of natural gas turbine
c_GT = 1500 #per unit capital cost of natural gas turbine [AUD/kW]
c_MGT = 52.5 #per unit cost of maintenance for whole life [AUD/kW]
c_g = 12 #per unit cost of natural gas [AUD/Gj]

#Fixed unit cost parameters of diesel generator
c_DG = 1200 #per unit capital cost of diesel generator [AUD/kW]
c_MDG = 52.5 #per unit cost of maintenance for whole life [AUD/kW]
c_l = 1.7/2 #per unit cost of diesel [AUD/L]

#Constant input parameter of natural gas turbine consumption curve
A = 0.2504
B = 40.8

#Constant input parameter of diesel generator consumption curve
C = 0.2009
D = 17.55

#Gas Turbine life
L_GT = 35 #years

#Diesel Generator life
L_DG = 25 #years

#Input parameters for emission model
em_GT = 51.53 #Carbon emission of natural gas per GJ unit [kg CO2-e/GJ]
u_LtoGJ = 0.0386 #Unit conversion factor from L to GJ [GJ/L]
u_m3toGJ = 0.0393 #Unit conversion factor from m3 to GJ [GJ/m3]
em_DG = 87.5 #Carbon emission of diesel generator per GJ unit [kg CO2-e/GJ]
u_kWhtoGJ = 0.0081 #Unit conversion factor from kWh to GJ [GJ/kWh]
u_kWhtoL = 0.3 #Unit conversion factor from kWh to L [L/kWh]


#Life cycle emissions
emlc_GT = 0.1269 #Life cycle carbon emission of natural gas turbine [kg CO2-e/kWh]
emlc_DG = 0.0747 #Life cycle carbon emission of diesel generator [kg CO2-e/kWh]


