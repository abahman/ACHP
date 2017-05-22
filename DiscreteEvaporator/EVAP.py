from __future__ import division, print_function, absolute_import
from math import pi,log,exp

#from scipy.optimize import brentq #solver to find roots (zero points) of functions
#import numpy as np
import pandas as pd

#import CoolProp as CP
#from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI

from extra_functions import HPtoTXP
from EvapCirc import EvapCircuit
from PRESSURE import dPelbow, dPmom
from VOLUME import VolumeALL
from CORR import Circuit_DP_EVAP

def Evaporator(Ref, #refrigerant string
               filename, #file nane string
               mr,#refrigerant mass flow rate
               HPo,#refrigerant outlet state
               GaI,#air mass flux
               TPi,#inlet air state
               HPi,#inlet refrigerant state of the segment
               TPo,#outlet air state
               Sm,#mass and volume in the evaporator
               Aflow,#air flow area
               Evap_struc,
               Prms, #adjustment factors
               NSeg,
               evapInit):
   
    '''    
    Model of a direct expansion evaporator  coil. The model is comprised 
    of small pieces of single finned tubes put together to form rows.  
    Rows of finned tubes are put together using elbows and manifolds to 
    make an entire coil. This function defines the coil's circuit pattern.
    Model parameters are stored in the file "evap.dat".
       
    INPUTS
    mr - refrigerant mass flow rate (kg/s).
    HPo - outlet refrigerant thermodynamic state (J/kg,kPa). 
    Ga - air mass flux (kg/s*m^2).
    TPi - inlet air thermodynamic state (K,-).
    Prms - array of three multiplicative correlation correction factors.
        Prms[0] - corrects air side convection coefficient (dry coil).
        Prms[1] - corrects refrigerant side convection coefficient.
        Prms[2] - corrects refrigerant side pressure drop correlation.
        Prms[3] - corrects air side convection coefficient (wet coil).
        Prms[4] - reduction in air flow rate (0-1) to simulate fouling.
        Set these parameters to 1 (Prms={1,1,1,1}) for no corrections.
       
    OUTPUTS
    HPi - inlet refrigerant thermodynamic state (J/kg,kPa).
    TPo - outlet air thermodynamic state (K,-).
    Sm - mass of charge in evaporator (kg) and volume of charge (m^3).
    Aflow = air flow area (ma=Ga*Aflow) (m^2).
    status - error flag.  Set to 1 to indicate error has occurred.
    '''
    D = Evap_struc.copy()
     
 
    if (evapInit > 0 or NSeg>0):
         
        df = pd.read_excel(filename,header = 0) #open the file and all the data are under column 'Evaporator'
         
        D['type'] = df.Evaporator[0]
        D['Di'] = df.Evaporator[1]             # tube inside diameter
        D['L'] = df.Evaporator[2]              # tube length 
        #printf("evap tube length = %lf\n",D.L);
        D['xp'] = df.Evaporator[3]             # tube thickness
        D['Df'] = df.Evaporator[4]              #outside fin diameter (m)
        D['z'] = df.Evaporator[5]            #it is used as the fin pitch
        D['th'] = df.Evaporator[6]          #spacing between tubes in bank (m)
        D['vsp'] = df.Evaporator[7]       #actually it is P_t
        D['P_l'] = df.Evaporator[8]         #spacing between tubs in the longitudual direction (m)
        D['NSeg'] = int(df.Evaporator[9])   #number of segments per tube in finite difference model
        if(NSeg>0): 
            D['NSeg'] = NSeg
        D['Brad']=df.Evaporator[10]     #radius of return bend (m)
         
        GaNom = df.Evaporator[11]      #nominal air mass flux (kg/s/m^2), correponding to the maximum air flux
         
        D['NBranchs'] = int(df.Evaporator[12]) # number of equivalent branches
        D['NBraTube'] = int(df.Evaporator[13]) # number of tube in each branch
        D['Nrows'] = int(df.Evaporator[14])        # number of rows high
        D['Ndeep'] = int(df.Evaporator[15])        # number of rows deep
        D['Frontal_A'] = df.Evaporator[16]    # frontal area
         
        #correction factors
        hAirAdj = df.Evaporator[17]
        hRefAdj = df.Evaporator[18]
        PRefAdj = df.Evaporator[19]
        WAirAdj = df.Evaporator[20]
        FoulFac = df.Evaporator[21]
   
        #B.S. ------------------------
        #new input parameter for evaporator
        D['microfin'] = int(df.Evaporator[22])     #microfin type, 0=smooth tube, 1=helical, 2=cross-grooved, 3=herringbone
        D['finN'] = df.Evaporator[23]              #fin number in a micro-fin tube
        D['gama'] = df.Evaporator[24]              #fin apex angle in a micro-fin tube
        D['beta'] = df.Evaporator[25]              #fin helix angle in a micro-fin tube
        D['finH'] = df.Evaporator[26]              #fin height in a micro-fin tube
        D['w_b'] = df.Evaporator[27]                #base width of a single fin
        D['w_e'] = df.Evaporator[28]                #top width of a single fin
        D['w_z'] = df.Evaporator[29]                #base distance between two neighboring fins
        # airside fin type
        D['airFin'] = int(df.Evaporator[30])       #1-plain, 2-corrugated, 3-slit, 4-louvered, 5-convex louvered, 6-smooth wavy, 7-spine
        D['sub1'] = df.Evaporator[31]              #fin substructures, for the possible use later, such as the waffle height of the wavy fin, louver number of the louverd fin
        D['sub2'] = df.Evaporator[32]               #fin substructures
        D['sub3'] = df.Evaporator[33]                #fin substructures
        D['sub4'] = df.Evaporator[34]                #fin substructures
        D['sub5'] = df.Evaporator[35]                #fin substructures
        D['K_T'] = df.Evaporator[36]                #400-- tube wall conductance which is the conductance factor of copper, for calculating tube resistance
        D['K_F'] = df.Evaporator[37]                    #237-- fin conductance
        #------------------------------B.S.
       
        endSign=int(df.Evaporator[38])      #-256
        if(endSign != -256):
            print ("Evaporator parameter input wrong")
            raise
   
        evapInit = 0
         
        D['Ls']=D['L']/float(D['NSeg'])
           
        #B.S.-----------------------------------
        # tube inside geometry
        if (D['microfin'] > 0):
            D['D_b'] = D['Di'] + 2*D['finH']        #tube diameter at the base of the fin
            D['Do'] = D['D_b'] + 2*D['xp']           #Pipe outside diameter.
            D['D_m'] = (D['Di'] + D['D_b'])/2       #mean diameter of the micro-fin tube
            D['P_H'] = D['finN']*(D['w_z']+2*pow((pow(D['finH'],2.0)+pow((0.5*(D['w_b']-D['w_e'])),2.0)),0.5)+D['w_e'])# the hydraulical circumference
            D['Acs'] = pi*pow((D['D_b']/2.0),2.0)-D['finN']*0.5*D['finH']*(D['w_b']+D['w_e'])#cross area of the micro-fin tube, this is the actual cross-section area
            D['Dh_i'] = 4*D['Acs']/D['P_H']               #Inside hydraulical diameter 
            D['Ax'] = pi/4.0*D['Di']*D['Di']              #Inside pipe cross sectional area, based on fin tips
            D['Api'] = pi*D['Di']*D['Ls']                 #Inside pipe surface area (one tube segment), based on fin tips
        else:
            D['finH']=0                     #tube without fins
            D['D_b'] = D['Di']              #tube diameter at the base of the fin
            D['Do'] = D['D_b']+2*D['xp']          #Pipe outside diameter.
            D['D_m'] = D['Di']                 #averge diameter of the micro-fin tube
            D['Ax'] = pi/4.0*D['Di']*D['Di']      # Inside pipe cross sectional area
            D['Api'] = pi*D['Di']*D['Ls']         # Inside pipe surface area (one tube segment)
            D['P_H'] = pi*D['D_b']             # the hydraulical circumference
            D['Acs'] = D['Ax']                 #cross area of the micro-fin tube, this is the actual cross-section area
            D['Dh_i'] = D['Di']                #inside hydraulical diameter
        #-----------------------------------------B.S.
        D['Vs'] = D['Acs']*D['Ls']
        D['Blen'] = pi*D['Brad']
        D['BVs'] = D['Blen']*D['Ax']
   
        ##*External Tube Parameters*##
        #B.S.-----------------------------
        #get the fin reference length, according to Schimidt equation
        M = D['vsp']/2e0
        R_O = D['Do']/2
        L = 0.5e0*pow((pow(M,2e0)+pow(D['P_l'],2e0)),0.5)
        BETA = L/M
        PHI = M/R_O
        R_EFF = R_O*1.27e0*PHI*pow((BETA-0.3e0),0.5e0)
        D['Df'] = 2.0*R_EFF                             #confirm the outside reference fin diameter
        D['L_F'] = (R_EFF-R_O)*(1e0+0.35e0*log(R_EFF/R_O))   #get the fin reference length for calculating the fin efficiency, according to Schimidt equation
        #--------------------------------------B.S.
   
        D['y'] = (D['Df']-D['Do'])/2        #Distance from outside of pipe to edge of fin. 
        D['N'] = int(D['Ls']/D['z'])        #B.S., change the original, (D.Ls/(D.z+D.th))    
        D['Apo'] = pi*D['Do']*(D['Ls']-D['N']*D['th'])       #B.S., correct pi*D.Do*(D.Ls)
   
        x = D['Df'] + D['th']/2.0
        D['Af'] = D['N']*pi/2.0*(x*x-D['Do']*D['Do'])                   #fin has faces
        D['Aflow'] = (D['vsp']-D['Do'])*(D['Ls']-D['N']*D['th'])        #minimum airflow area per segment
        # printf("Minimum free flow area = %lf (m^2)\n",D.Aflow)
         
        #* hydraulic diameter *#
        D['Dh']=4*D['Aflow']/(D['Af']+D['Apo'])*D['vsp']     
     
   
    D['hAirAdj'] = Prms[0]*hAirAdj
    D['hRefAdj'] = Prms[1]*hRefAdj
    D['PRefAdj'] = Prms[2]*PRefAdj
    D['WAirAdj'] = Prms[3]*WAirAdj
   
    if (GaI <= 0):
        Ga = FoulFac*GaNom*Prms[4]
    else:
        Ga = GaI
       
    # Calculate air cfm after Ga is established.
    A_total = D['Aflow']*D['NSeg']*float(D['Nrows']) #total area air flows throguh for entire coil (m^2)
    v_air = HAPropsSI("V", "T", TPi['T'], "P", 101325, "R", TPi['P']) # specific volume of air [m^3/kg dry air]
    vfra = Ga*A_total*v_air # nominal volumetric flow rate air
    D['cfma'] = vfra*(60.0/((12.0*0.0254)*(12.0*0.0254)*(12.0*0.0254))) # nominal cfm of air flow
   
    #B.S., prepare to set up the simple model (A.B., this might not be necessary because they are already initialized with 0 values in Evap_Sturc dictionary
    D['V_TOT']=0; D['V_TP']=0; D['V_Liq']=0; D['V_Vap']=0;
    D['m_TOT']=0; D['m_TP']=0; D['m_Liq']=0; D['m_Vap']=0;
    D['mr']=0; D['ma_TOT']=0;
    D['UA_TOT']=0; D['UA_Liq']=0; D['UA_Vap']=0; D['UA_TP']=0;
    D['UAw_TOT']=0; D['UAw_Liq']=0; D['UAw_Vap]']=0; D['UAw_TP']=0;
    D['LiqL']=0; D['TPL']=0; D['VapL']=0;
    D['L_dry']=0; D['L_wet']=0;
    D['count1']=0; D['count2']=0;
    D['HP_TP1']['H']=0; D['HP_TP1']['P']=0;
    D['HP_TP2']['H']=0; D['HP_TP2']['P']=0;
    D['Qtp_dry']=0; D['Qtp_wet']=0;
    #-----------------------------------------B.S.
       
    D['GetP'] = 0     #don't calculate the airside pressure drop in the following procedure
       
    D['REV'] = Evap_struc['REV']    #iteration loop direction
     
    HPo, HPi, TPo, Sm, Aflow, D = EvapCircuit(Ref,D['type'],mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,D)
   
    DP_circuit, TPo, D = Circuit_DP_EVAP(Ga,TPi,TPo,D)#B.S., calculate the airside pressure drop across the circuit
    #gOutData.Eva_AirDP = DP_circuit#for output
   
    #B.S.---------------------------------------
    D['HP_out'] = HPo      #suciton site of the evaporator
    D['HP_in'] = HPi       #entrance of the evaporator
    D['mr'] = mr
    D['Ga'] = Ga
    D['TPi'] = TPi
    D['ma_TOT'] = Ga*D['Aflow']*D['NSeg']*D['Nrows']
    #Build_Eva(&D)        #correlate the moving boundary and lumped evaporator model
    Evap_struc = D          #output the evaporator construct
    #---------------------------------------B.S.
  
    return (HPo, HPi, TPo, Sm, Aflow, Evap_struc)



def EvapTubeBend(Ref,Gr,HPo,m,P):
    '''
    /********************************************************************
    Models an evaporator tube return bend. (Reverse)
    ********************************************************************/
    '''
    
    TXPo = HPtoTXP(HPo, Ref)

    vo = VolumeALL(TXPo,Gr,P['Di'],0, Ref)

    dP = dPelbow(TXPo,Gr,P['Di'],P['Brad'], Ref)

    HPm = {'H':0.0,'P':0.0};
    HPm['H']=HPo['H'];
    HPm['P']=HPo['P']-dP;
    
    TXPm=HPtoTXP(HPm, Ref)
    
    vm = VolumeALL(TXPm,Gr,P['Di'],0, Ref)

    dP = dPmom(vo,vm,Gr)

    HPo['H']=HPm['H'];
    HPo['P']=HPm['P']-dP;

    m['V']=P['BVs']; 
    m['m']=P['BVs']/vm;
    
    #B.S.----------------------------------------------------
    if (TXPm['X']==1):      #B.S., the following caculate the inner volume and mass in the tube bend    
        P['m_Vap'] = P['m_Vap'] + m['V']/vm;      #B.S. record the vapor mass
        P['V_Vap'] = P['V_Vap'] + m['V'];         #B.S. record the vapor volume
    elif (TXPm['X']<1 and TXPm['X']>0.00):
        P['m_TP'] = P['m_TP'] + m['V']/vm;        #B.S. record the two-phase mass
        P['V_TP'] = P['V_TP'] + m['V'];           #B.S. record the two-phase volume
    else:
        P['m_Liq'] = P['m_Liq'] + m['V']/vm;    #B.S. record the liquid mass
        P['V_Liq'] = P['V_Liq'] + m['V'];       #B.S. record the liquid volume
    #---------------------------------------------B.S.
    return (HPo, m, P)

def EvapTubeBend_Fwd(Ref,Gr,HPo,m,P):
    '''
    /********************************************************************
    Models an evaporator tube return bend.
    ********************************************************************/
    '''

    TXPo = HPtoTXP(HPo, Ref)

    vo = VolumeALL(TXPo,Gr,P['Di'],0, Ref)
    
    dP = dPelbow(TXPo,Gr,P['Di'],P['Brad'], Ref)

    HPm = {'H':0.0,'P':0.0};
    HPm['H']=HPo['H'];
    HPm['P']=HPo['P']+dP;
    
    TXPm=HPtoTXP(HPm, Ref)

    vm = VolumeALL(TXPm,Gr,P['Di'],0, Ref)

    dP = dPmom(vo,vm,Gr)

    HPo['H']=HPm['H'];
    HPo['P']=HPm['P']-dP;

    m['V']=P['BVs']; 
    m['m']=P['BVs']/vm;

    #B.S.----------------------------------------------------
    if (TXPm['X']==1):      #B.S., the following caculate the inner volume and mass in the tube bend    
        P['m_Vap'] = P['m_Vap'] + m['V']/vm;    #B.S. record the vapor mass
        P['V_Vap'] = P['V_Vap'] + m['V'];       #B.S. record the vapor volume
    elif(TXPm['X']<1 and TXPm['X']>0.00):
        P['m_TP'] = P['m_TP'] + m['V']/vm;      #B.S. record the two-phase mass
        P['V_TP'] = P['V_TP'] + m['V'];         #B.S. record the two-phase volume
    else:
        P['m_Liq'] = P['m_Liq'] + m['V']/vm;    #B.S. record the liquid mass
        P['V_Liq'] = P['V_Liq'] + m['V'];       #B.S. record the liquid volume
    #---------------------------------------------B.S.

    return (HPo, m, P)
    
