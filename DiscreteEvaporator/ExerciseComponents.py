from __future__ import division, print_function, absolute_import
from math import pi,log,exp

#from scipy.optimize import brentq #solver to find roots (zero points) of functions
#import numpy as np
import pandas as pd

#import CoolProp as CP
from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI

from ACHP.convert_units import *
from extra_functions import TXPtoHP, HPtoTXP, PropertyTXPth

from EVAP import Evaporator
    
def ExerciseEvaporator():
    
    Ref = 'R410A'
    #design conditions
    SH = 23.004 #superheat (F)
    #SC = 10.9512 #subcooling F (***no effect)
    ET_SL = 48.3264-23 # F, evaporation temperature at the suction line
    ET_DT = 45+4 # F, add 4F due to evaporator coild pressure drop
    #CTOA_DL = 29.511+4 # F, add 4F due to pressure drop in cond coil, 18F is at the liquid line (** no changes)
    #CTOA_LL = 29.511 # F, at the liquid line  (** no changes)
    AMB = 95 # F
    RA = 83.048#80 # indoor room dry bulb temperature (F)
    RARH = 0.51 # indoor room humidity (-)
    RAWB = 67.856 #  indoor room wetbulb temoperature (F)
    MR = 955.572 # lbm/hr, refrigerant mass flow rate
    #GOA = 3.82 # kg/s/m^2, outdoor air mass flux (this is for condenser component)
    GIA = 5.6 # kg/s/m^2, indoor air mass flux (0.9237/0.3623)
    #DSH = 67.13 # F, discharge superheat (***nothing changes)
    
    
    #initialize the values
    hpi = {'H':0.0,'P':0.0}
    m = {'m':0.0,'V':0.0}
    TPo = {'T':0.0,'P':0.0}
    Aflow = 0.0
    Evap_prms = [1,1,1,1,1] #correction factor (see Evaporator model for more details)
    NSeg = -1 #keep it to -1 to use the default number of segment (10 segments)
    Evap_struc={'Di':0.0,'L':0.0,'xp':0.0,'Df':0.0,'z':0.0,'th':0.0,'vsp':0.0,'Brad':0.0,
                'NSeg':int(0),
                'Dh':0.0,'Do':0.0,'y':0.0,'Ls':0.0,'N':0.0,'Ax':0.0,'Api':0.0,'Apo':0.0,'Af':0.0,'Aflow':0.0,'Vs':0.0,'Blen':0.0,'BVs':0.0,
                'Ro:':0.0,
                'Gr':0.0,'Ga':0.0,
                'HPo':{'H':0.0,'P':0.0},
                'TXPo':{'T':0.0,'X':0.0,'P':0.0},
                'TPi':{'T':0.0,'P':0.0},
                'hAirAdj':0.0,'hRefAdj':0.0,'PRefAdj':0.0,'WAirAdj':0.0,
                'type':int(0),
                'Nrows':int(0),'Ndeep':int(0),
                'NBranchs':int(0), 'NBraTube':int(0),
                'P_l':0.0, #spacing between tubs in the longitudual direction (m)
                'microfin':int(0), #new parameters for micro-fin tubes, and specially configured fins #decide whether is micro-fin tubes microfin=1
                'w_b':0.0, #base length for single micro-fin
                'w_e':0.0, #top length for single micro-fin
                'w_z':0.0,    #width between the neighboring micro-fins
                'finH':0.0,    #micro-fin height
                'beta':0.0,#micro-fin helical angle
                'gama':0.0,    #micro-fin apex angle
                'finN':0.0,    #total micro-fin number
                'Acs':0.0,    # micro-fin cross-sectional area
                'P_H':0.0,    #micro-fin hydraulic circumference
                'Dh_i':0.0,    #micro-fin tube hydraulic diameter
                'K_T':0.0,    #conductance factor of tube wall
                'K_F':0.0,    #conductance factor of fin material
                'L_F':0.0,    #reference fin length for Schimidt fin efficiency calculation
                'D_b':0.0,'D_m':0.0,    #fin base diameter, and mean diameter
                'airFin':int(0),
                'sub1':0.0,'sub2':0.0,'sub3':0.0,'sub4':0.0,'sub5':0.0,#for inputing sub-structures of fin surface
                'ho':0.0, 'wetadj':0.0,#airside heat transfer coefficient
                'Frontal_A':0.0,#frontal area
                'GetP':0.0,#calculate the airside pressure drop
                #variables for generating the simple evaporator function
                'V_TOT':0.0, 'V_TP':0.0, 'V_Liq':0.0, 'V_Vap':0.0,#inner volume of different phase
                'L_TOT':0.0,'LiqL':0.0,'VapL':0.0,'TPL':0.0,#tube length of different phase
                'L_dry':0.0, 'L_wet':0.0,#tube length of dry and wet heat transfer
                'A_TOT':0.0,'A_Liq':0.0,'A_Vap':0.0,'A_TP':0.0,#heat transfer surface area of different phase
                'm_TOT':0.0, 'm_TP':0.0, 'm_Liq':0.0, 'm_Vap':0.0,#airflow rate across different phase
                'rho_TOT':0.0, 'rho_TP':0.0, 'rho_Liq':0.0, 'rho_Vap':0.0,#average density of different phase
                'U_TP':0.0, 'U_Liq':0.0, 'U_Vap':0.0,#averge dry heat conductance of different phase 
                'Uw_TP':0.0, 'Uw_Liq':0.0, 'Uw_Vap':0.0,#average wet heat conductance of different phase
                'DP_TOT':0.0, 'DP_TP':0.0, 'DP_Liq':0.0, 'DP_Vap':0.0,#average pressure gradient of different phase
                'UA_TOT':0.0, 'UA_Liq':0.0,'UA_Vap':0.0,'UA_TP':0.0,#overall dry heat conductance of different phase
                'UAw_TOT':0.0, 'UAw_Liq':0.0,'UAw_Vap':0.0,'UAw_TP':0.0,#overall wet heat conductance of different phase
                'mr':0.0, 'ma_TOT':0.0,#overall refrigerant and air mass flow rate
                'Ga_meanL':0.0,#average air flow rate per tube length
                'HP_out':{'H':0.0,'P':0.0}, 'HP_dry':{'H':0.0,'P':0.0}, 'HP_wet':{'H':0.0,'P':0.0}, 'HP_in':{'H':0.0,'P':0.0},'HP_TP1':{'H':0.0,'P':0.0},'HP_TP2':{'H':0.0,'P':0.0},#state parameters at important locations
                'count1':0.0,'count2':0.0,#count the state points of two-phase flow begin point and end point 
                'Qtp_dry':0.0, 'Qtp_wet':0.0,#two-phase dry heat transfer and wet heat transfer amount
                'r_v':0.0, 'r_l':0.0,'r_tp':0.0,#parameters for adjusting the theoretical heat transfer effectiveness of each phase
                'H2_residual':0.0,#for the consistency of the moving boundary model analysis
                #------------------------------B.S.
                'q_flux':0.0, 'T_w':0.0,#segment heat flux and inside tube wall temperature
                'wet':int(0),#wet=0 to calculate dry heat transfer, wet=1 to calculate wet heat transfer
                'REV':int(0),
                'cfma':0.0,
                'Hout8':[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],'DPr':[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],#superheat and pressure drop at each evaporator branch (array of 10 elements)
                }

    #inlet air state
    TPi = {'T':F2K(RA), 'P':RARH} #T: air inlet (return air) dry bulb temperature [K], P: air inlet (return air) relative humidity [-]
    #air mass flux (kg/s/m^2)
    Ga = GIA
    #refrigerant mass flow rate [kg/s]
    mr = lbh2kgs(MR)
    #outlet refrigerant state
    ETBRA = DeltaF2K(RA-ET_SL) #evaporating temperature below return air [K]
    Tsh = DeltaF2K(SH) #suction line superheat [K]
    
    Tsat_sl = TPi['T']-ETBRA #Temperature of saturated liquid [K]
    P_sl = PropsSI("P", "T", Tsat_sl, "Q", 1, Ref)
    txpo = {'T':Tsat_sl+Tsh,'X':1,'P':P_sl} #initially assume outlet temperature equals to saturated vapor pressure (txpo_P)
    hpo = TXPtoHP(txpo, Ref)  #outlet refrigerant state
    Twbi = HAPropsSI("Twb", "T", TPi['T'], "P", 101325, "R", TPi['P']) #inlet air wet bulb temperature 

    #*Run the evaporator model*#
    hpo, hpi, TPo, m, Aflow, Evap_struc = Evaporator(Ref,'../InputDoc/Evap_Data.xlsx',mr,hpo,Ga,TPi,hpi,TPo,m,Aflow,Evap_struc,Evap_prms,NSeg,evapInit=1)
    
    #inlet refrigerant state and   
    #saturation temperature computed from the inlet refrigerant state
    txpi = HPtoTXP(hpi)
    Tsat_dt = PropertyTXPth('T',txpi,Ref) #[K]
     
    Twbo = HAPropsSI("Twb", "T", TPo['T'], "P", 101325, "R", TPo['P']) #outlet air wet bulb temperature 
     
    #Calculate the total heat transfer and then the UA = Qdot/DT
    Qdot = W2BTUh(mr*(hpo['H']-hpi['H'])) #[BTU/hr] - for all 4 sections of the heat exchanger
    DT = DeltaK2F(TPi['T'] - Tsat_sl) #[F]
    UA = Qdot/DT #[BTU/hr/F] - for all 4 sections of the heat exchanger 
         
    #print inputs
    print("Tai (F), RHai (%)", C2F(TPi['T']),TPi['P']*100.0)
    print("CFM (cfm)",Evap_struc['cfma'])
    print("mr (lbm/hr)", MR)
    print("ETBRA_sl (F)",DeltaK2F(TPi['T']-Tsat_sl))
    print("Tsh_sl (F)",DeltaK2F(txpo['T']-Tsat_sl))
    print("X_sl (0-1)",txpo['X'])
    if (NSeg>0):
        print("NSeg",NSeg)
    else:
        print("NSeg",'Default')
     
    #print outputs
    print("ETBRA_dt (F)",DeltaK2F(TPi['T']-Tsat_dt)) #ETBRA_dt
    print("X_dt (0-1)",txpi['X'])
    print("dP (psi)",kPa2psi(txpi['P']-txpo['P'])) ###BEWARE that txpo_P was set to saturatewd liquid pressure, check if the code modifies it
    print("DTa (F)",DeltaK2F(TPi['T']-TPo['T']))
    print("DTwb (F)",DeltaK2F(Twbi-Twbo))
    print("mass (kg)",m['m']) #charge
    print("Qdot (BTU/hr)",Qdot)
    print("UA (BTU/hr/F)",UA)
    print("Hin (kJ/kg)",hpi['H']/1000.0)
    print("Hout (kJ/kg)",hpo['H']/1000.0)
    print("Tsat_sl (C)",K2C(Tsat_sl))
    print("Tsat_dt (C)",K2C(Tsat_dt))

    return 0
    
def ExerciseComponents():
    ExerciseEvaporator()

if __name__=='__main__':    
#     ExerciseComponents()
    df = pd.read_excel("../InputDoc/EvapStruc.xlsx",sheetname='tube', header = 0)
    print (df)
