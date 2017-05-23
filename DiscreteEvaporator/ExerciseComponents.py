from __future__ import division, print_function, absolute_import
from math import pi,log,exp

#from scipy.optimize import brentq #solver to find roots (zero points) of functions
#import numpy as np
import pandas as pd

from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI

from ACHP.convert_units import *
from extra_functions import TXPtoHP, HPtoTXP, PropertyTXPth, ETdim

from EVAP import Evaporator
    
def ExerciseEvaporator():
    
    Ref = 'R407C'
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
    Evap_struc= ETdim() #initialized evaporator staructure

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
