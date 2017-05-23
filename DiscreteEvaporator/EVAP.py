from __future__ import division, print_function, absolute_import
from math import pi,log,exp

from scipy.optimize import brentq #solver to find roots (zero points) of functions
#import numpy as np
import pandas as pd

#import CoolProp as CP
#from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI

from extra_functions import HPtoTXP, EVA_Get_Q, toTXP, PreAcc
from EvapCirc import EvapCircuit
from PRESSURE import dPelbow, dPmom, GET_PreAcc
from VOLUME import VolumeALL
from CORR import Circuit_DP_EVAP, ConvCoeffAir_EVA, FinEffect_Schmidt, ConvCoeffSP
from CMIN import CmineCrossFlow_dry

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
    

# def EvapTubeL_Rev(Gr,#refrigerant mass flux
#                   HPo,#refrigerant outlet and inlet state
#                   Ga,#air mass flux
#                   TPi,#air inlet state
#                   WHo,#air outlet state
#                   m,#charge and inner volume in the evaporator
#                   P):#evaporator struct
#     '''
#     /*********************************************************************
#     Evaporator tube segment model which neglects pressure drops
#     Inputs:
#         Gr = refrigerant mass flux (kg/s/m^2)
#         HPi = refrigerant inlet state (h,P)
#         Ga = air mass flux (kg/s/m^2)
#         Tai = air inlet temperature (C)
#     Outputs:
#         HPi = refrigerant outlet state (h,P)
#         hao = air outlet enthalpy (J/kg)
#         m = mass of charge in return bend (kg)
#     *********************************************************************/
#     '''
# 
#     double v,hai,ma,mr;
#     TXP TXPo;
#     EVA_Get_Q EVA_Get_Q;#B.S., this struct is for storing some parameters for iteration
#     double y=0;
#     double hi=0;
#     double Ri=0,R=0;
#     double q=0;
# 
#     /* Calculate air side resistance */
#     P->wet=0;
#     const double ho = P->hAirAdj*ConvCoeffAir_EVA(TPi,Ga,P);#B.S., the airside dry heat transfer coefficient can be got for different fin types
#     if(errorLog.IsError()) {
#         errorLog.Add("EvapTubeL", "get ho");
#         return 0;
#     }
# 
#     const double phi = FinEffect_Schmidt(ho,233,P->th,P->y,P->Do);#B.S. calculate the fin efficiency with the Schmidt equation
#     if(errorLog.IsError()) {
#         errorLog.Add("EvapTubeL", "FinEffect_Schmidt");
#         return 0;
#     }
# 
#     P->Ro = 1/(ho*(P->Apo+phi*P->Af));
#     P->ho=ho;
#     P->Ga=Ga;
# 
#     
#     # convert inlet state into TXP format
#     TXPo = HPtoTXP(*HPo);
#     if(errorLog.IsError()) {
#         errorLog.Add("EvapTubeL_new", "(TXPo,1)");
#         return 0;
#     }
# 
#     #B.S.--------------------------------------
#     TXP TXP_bak;
#     TXP_bak = TXPo;#backup the outlet state of this segment
#     #---------------------------------------B.S.
# 
#     # Mass flow rates of air and refrigerant
#     ma=Ga*P->Aflow;#B.S. Ga is the maximum airflow flux, P->Aflow is the minimun air flow cross-sectional area
#     mr=Gr*P->Ax;
# 
#     # Calculate heat transfered per unit mass of refrigerant.
# 
#     const double hi_max=10e4;#largest possible refrigerant side heat tranfer coeffcient
#     const double hi_min=100;#minimum possible refrigerant side heat transfer coeffcient
# 
#     #B.S.-----------------------------
#     #store the parameters for iteration
#     EVA_Get_Q.ma=ma;
#     EVA_Get_Q.mr=mr;
#     EVA_Get_Q.TPi=TPi;
#     EVA_Get_Q.Gr=Gr;
#     EVA_Get_Q.P=*P;
#     EVA_Get_Q.q=0;
#     EVA_Get_Q.W=wair.HumidityRatio(TPi.T,TPi.P);
#     if(errorLog.IsError()) {
#         errorLog.Add("EvapTubeL_new", "(W)");
#         return 0;
#     }
#     #---------------------------------------B.S.
# 
#     const double X1=0.05,X2=0.95;#    const double X1=0.01,X2=0.99;
#     #for calculating the heat transfer and pressure drop in single-phase or approximated single-phase region
# 
#     if(TXPo.X>X2||TXPo.X<X1) {
# 
#         if(TXPo.X>=1 || TXPo.X<=0) {#Pure single phase region
# 
#             y, P = ConvCoeffSP(TXPo,Gr,P, Ref);#get the single-phase heat transfer coefficient
#             if(errorLog.IsError()) {
#                 errorLog.Add("EvapTubeLNew","(hi_SP)");
#                 return 0;
#             }
#         } else if(TXPo.X<X1) {
#     
#             TXP TXP1 = toTXP(TXPo.T,0.0,TXPo.P);
#             TXP TXP2 = toTXP(TXPo.T,X1,TXPo.P);
# 
#             y1, P = ConvCoeffSP(TXP1,Gr,P, Ref);
#             if(errorLog.IsError()) {
#                 errorLog.Add("EvapTubeLNew","(hi_SP)");
#                 return 0;
#             }
#         
#             EVA_Get_Q.TXPo=TXP2;#local refrigerant state
#             hi = Zbrent(hi_max,hi_min,Get_Q_EVA,1e-5,&EVA_Get_Q);#B.S., get the heat transfer coefficient at the two-phase region
#             if(errorLog.IsError()) {
#                 errorLog.Add("EvapTubeLNew","Zbrent1");
#                 return 0;
#             }
# 
#             double y2 = hi/P->hRefAdj;#B.S., calculate the two-phase heat transfer coefficient at this region
#             
#             if(y1>y2) y1=y2;#B.S., make the single phase flow less than two phase
#             y = y1+TXPo.X*(y2-y1)/X1;#B.S., get the heat transfer coefficient in this region with intropolation
#         } else {
# 
#             TXP TXP1 = toTXP(TXPo.T,X2,TXPo.P);
#             TXP TXP2 = toTXP(TXPo.T,1,TXPo.P);
#     
#             EVA_Get_Q.TXPo = TXP1;
#             hi = Zbrent(hi_max,hi_min,Get_Q_EVA,1e-5,&EVA_Get_Q);#two-phase refrigerant heat transfer coefficient
#             #B.S., get the heat transfer coefficient at the two-phase region
#             if(errorLog.IsError()) {
#                 char msg[256];
#                 sprintf(msg,"Zbrent2: TXPo.T=%.2lf, T_max=%.2lf, T_min=%.2lf",TXPo.T,hi_max,hi_min);
#                 errorLog.Add("EvapTubeLNew",msg);
#                 return 0;
#             }
# 
#             double y1 = hi/P->hRefAdj;#ConvCoeffEvapTP_microfin(TXP1,Gr,P);#B.S. get two-phase heat transfer coefficient at this region    
#         
#             y2, P = ConvCoeffSP(TXP2,Gr,P, Ref);#B.S., single-phase heat transfer coefficient
# 
#             if(y2>y1) y2=y1; #B.S., for making the single phase flow less than two phase
#             y = y2-(1-TXPo.X)*(y2-y1)/(1-X2);#B.S., get the heat transfer coefficient in this region with intropolation
#         } 
# 
#         hi = P->hRefAdj*y;#B.S., heat transfer coefficient at this region
#         Ri = 1/(hi*P->Api);#B.S., inside thermal resistance
#         const double R_W=log(P->Do/(P->Do-2.0*P->xp))/(2.0*3.1415*P->K_T*P->Ls);
#         R = P->Ro+R_W+Ri;#overall thermal resistance
#     
#         #B.S., prepare the parameters for iteration
#         EVA_Get_Q.TXPo = TXPo;
#         EVA_Get_Q.Cmin = CmineCrossFlow_dry(R,mr,ma,TXPo,TPi, Ref);
#         EVA_Get_Q.HPo = *HPo;
#         
#         if(errorLog.IsError()) {
#                 errorLog.Add("EvapTubeL_new", "T_sat");
#                 return 0;
#             }
# 
#         if(0/*TXPo.T<T_sat+1*/) {
#             # TXPo.T is actually the outlet refrigerant temperature
#             # however TXPo.T = TXPi.T in the two phase region and
#             # TXPo.T ~ to TXPi.T in the superheated region if steps are small
# 
#             q = (TPi.T-TXPo.T)/mr*CmineCrossFlow_dry(R,mr,ma,TXPo,TPi, Ref);
# 
#             if(errorLog.IsError()) {
#                 errorLog.Add("EvapTubeL_new", "(q_singlephase)");
#                 return 0;
#             }    
# 
#         } else {#when it is close to the saturated state, the function zbrent can not converge
# 
#             const double H_max=HPo->H;
#             const double H_min=HPo->H-4e4;
#             EVA_Get_Q.hi=hi;#input the refrigerant side heat transfer coefficient
#             Zbrent(H_max,H_min,Get_Q_Single,1e-7,&EVA_Get_Q);#B.S., to calculate the single-phase heat transfer and pressure drop at this region
#             q=EVA_Get_Q.q;
#         }
# 
#         if(errorLog.IsError()) {
#             errorLog.ClearError("EvapTubeL_new", "(q_singlephase)");
#             q = (TPi.T-TXPo.T)/mr*CmineCrossFlow_dry(R,mr,ma,TXPo,TPi, Ref);
#         }/**#/don't delete this
# 
#         #when the superheat is setted too high, there might exist a problem
#     /*    if(q<=0) {#B.S., to remove the wrong result
#             q=0;
#         }*#/don't delete this
# 
#         if(errorLog.IsError()) {
#             errorLog.Add("EvapTubeL_new", "(q_singlephase)");
#             return 0;
#         }    
# 
#     
#         #B.S.---------------------------------------------
#     
#         if(TXPo.X>=0.9999999) {#B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.    
#             R = (P->Ro+Ri);
#             P->UA_Vap=P->UA_Vap+1/R;#B.S., vapor heat transfer conductance
#             P->VapL = P->VapL+ P->Ls;#B.S., vapor length
#         } else if(TXPo.X<=0.00000001) {
#             R = (P->Ro + Ri);
#             P->LiqL = P->LiqL+ P->Ls;#B.S., keep the liquid length in the evaporator
#             P->UA_Liq=P->UA_Liq+1/R;#B.S. liquid heat transfer conductance    
#         } else {
#             R = P->Ro+Ri;
#             P->UA_TP=P->UA_TP+1/R;#B.S., dry thermal resistance of the two-phase heat transfer
#             P->TPL = P->TPL+ P->Ls;#B.S., length of two-phase flow
#         }
#     #-----------------------------------------B.S.
# 
#     } else {#calculating the heat transfer in the two-phase region
# 
# 
#         EVA_Get_Q.TXPo=TXPo;
#         hi = Zbrent(hi_max,hi_min,Get_Q_EVA,1e-5,&EVA_Get_Q);#iterate the heat transfer coefficient
# 
#         if(errorLog.IsError()) {
#             # plot error vs. T to see whgat function looks like that caused error
#         #    ZbrentPlot(hi_max,hi_min,Get_Q_EVA,1e-4,&EVA_Get_Q);#Haorong change from -7 to-2
#             char msg[2048];
#             sprintf(msg,"(Zbrent3) T=%.4lf, T_min=%.4lf, T_max=%.4lf",hi,hi_min,hi_max);
#             errorLog.Add("EvapTubeLnew", msg);
#             return 0;
#         }; 
# 
#         q=EVA_Get_Q.q;
# 
#         #B.S.-----------------------------------
#         P->Qtp_wet = P->Qtp_wet+EVA_Get_Q.P.Qtp_wet;#wet heat transfer amount
#         P->Qtp_dry =P->Qtp_dry+EVA_Get_Q.P.Qtp_dry;#dry heat transfer amount
#         P->UAw_TP = EVA_Get_Q.P.UAw_TP+P->UAw_TP;#heat conductance of two-phase wet heat transfer
#         P->UA_TP = EVA_Get_Q.P.UA_TP+P->UA_TP;#heat conductance of two-phase dry heat transfer
#         P->L_wet = EVA_Get_Q.P.L_wet + P->L_wet;#wet heat transfer length in the two-phase region
#         P->L_dry = EVA_Get_Q.P.L_dry + P->L_dry;#wet heat transfer length in the two-phase region
#         P->TPL = P->TPL+ P->Ls;#length of the two-phase heat transfer
#         #------------------------------------------B.S.
#     }
# 
#     HPo->H-=q;
# 
#     # calculate pressure drop
#     double P_out;
# 
#     Gr = mr/P->Acs;#this mass flux is for calculating the pressure drop
#     const double DP_FR = FricDP(TXPo, Gr, q, P);
#     
#     if(errorLog.IsError()) {
#         errorLog.Add("EvapTubeL_new", "DP_FR");
#         return 0;
#     }
#     
# 
#     if(TXPo.X>0.05&&TXPo.X<0.99)
#     {
#     #B.S., prepared to calculate the two-phase acceleration pressure drop
#     PreAcc Preacc;
#     Preacc.DP_FR=-1*DP_FR;
#     Preacc.G=Gr;
#     Preacc.H_OUT=HPo->H;
#     Preacc.P_IN=TXPo.P;
#     Preacc.X_IN=TXPo.X;
#     
#     const double DP_ACC=Zbrent(10,-10,GET_PreAcc,1e-7,&Preacc);#B.S., calculate the acceleration pressure drop 
#     if(errorLog.IsError()) {
#         errorLog.Add("EvapTubeL_new", "(Preacc)");
#         return 0;
#     }
#     P_out=TXPo.P+DP_FR-DP_ACC;
#     }
#     else {
#     P_out=TXPo.P+DP_FR;
#     }
# 
#     HPo->P=P_out;
# 
#     # Determine mass of charge.  It is based on the inlet state.
#     # The specific volume is recalculated so that a different model
#     # can be used from the one used to calculate the pressure drop.
# 
#     TXPo = HPtoTXP(*HPo);
#     if(errorLog.IsError()) {
#         errorLog.Add("EvapTubeLnew", "(TXPo,2)");
#         return 0;
#     }
# 
#     v = VolumeALL(TXPo,Gr,P->Di,mr*q/P->Api);        # seperate flow model
#     # v = PropertyTXP(VOL,TXPo);        # homogeneous flow model
#     if(errorLog.IsError()) {
#         errorLog.Add("EvapTubeLnew", "(v)");
#         return 0;
#     }
# 
#     m->V=P->Ls*P->Acs;#B.S., use the real cross-sectional area to calculate the inner volume
#     m->m=m->V/v;
# 
#     # Calculate output air state
#     hai = wair.h(TPi.T,TPi.P);
#     if(errorLog.IsError()) {
#         errorLog.Add("EvapTubeLnew", "(hai)");
#         return 0;
#     }
# 
#     const double W_I=wair.HumidityRatio(TPi.T,TPi.P);
#     if(errorLog.IsError()) {
#         errorLog.Add("EvapTubeLnew", "(W_I)");
#         return 0;
#     }
# 
#     const double HF_water=(EVA_Get_Q.T_S_O*4.1877+0.0594)*1e3;
# 
#     WHo->H=hai-q*mr/ma-HF_water*(W_I-EVA_Get_Q.W);
# 
#     WHo->W = EVA_Get_Q.W;#B.S. the outlet air humidity is calculated with the fuction Get_Q_EVA below.
#     if(errorLog.IsError()) {
#         errorLog.Add("EvapTubenew");
#         return 0;
#     }
# 
#     
#     #B.S.----------------------------------------
#     if(TXPo.X>=0.999999999999)#B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.    
#     {
#     P->m_Vap = P->m_Vap+m->V/v;#B.S.vapor mass
#     P->V_Vap = P->V_Vap + m->V;#B.S.vapor volume
#     }
#     else if(TXPo.X<=0.0000001)
#     {
#     P->m_Liq = P->m_Liq+m->V/v;#B.S., liquid mass
#     P->V_Liq = P->V_Liq + m->V;#B.S., liquid volume
#     
#     if(TXP_bak.X>0.0000001)
#     {
#     P->HP_TP1.H=P->HP_TP1.H+mr*HPo->H;#B.S., corresponding to the first segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
#     P->HP_TP1.P=P->HP_TP1.P+mr*HPo->P;#B.S., corresponding to the first segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
#     P->count1 = P->count1+mr;#B.S., this addition is for counting how many evaporator circuits that involved with heat transfer calculation
#     }
# 
#     }
#     else 
#     {
#     P->m_TP = P->m_TP+m->V/v;#B.S., mass of two-phase flow
#     P->V_TP = P->V_TP + m->V;#B.S., volume of two-phase flow
# 
#     if(TXP_bak.X>=0.999999999999)
#     {
#     P->HP_TP2.H=P->HP_TP2.H+mr*HPo->H;#B.S. corresponding to the last segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
#     P->HP_TP2.P=P->HP_TP2.P+mr*HPo->P;#B.S. corresponding to the last segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
#     P->count2 = P->count2+mr;#B.S., this addition is for counting how many evaporator circuits that involved with heat transfer calculation
#     }
#     }
#     #---------------------------------------------B.S.
# 
#     return (HPo, WHo, m, P)


def EvapTubeL_Fwd(Ref, #refrigerant 
                  Gr,#refrigerant mass flux
                  HPo,#refrigerant outlet and inlet state
                  Ga,#air mass flux
                  TPi,#air inlet state
                  WHo,#air outlet state
                  m,#charge and inner volume in the evaporator
                  P):#evaporator struct
    '''
    /*********************************************************************
    Evaporator tube segment model which neglects pressure drops
    Inputs:
        Gr = refrigerant mass flux (kg/s/m^2)
        HPi = refrigerant inlet state (h,P)
        Ga = air mass flux (kg/s/m^2)
        Tai = air inlet temperature (C)
    Outputs:
        HPi = refrigerant outlet state (h,P)
        hao = air outlet enthalpy (J/kg)
        m = mass of charge in return bend (kg)
    *********************************************************************/
    '''
    
    TXPo = {'T':0.0,'X':0.0,'P':0.0};
    
    #B.S., this dictionary struct is for storing some parameters for iteration
    EVA_Get_Q = EVA_Get_Q()
    
    y=0;
    hi=0;
    Ri=0;R=0;
    q=0;

    ### Calculate air side resistance ###
    P['wet']=0;#without considering wet heat tranfer adjustment

    ho, P = P['hAirAdj']*ConvCoeffAir_EVA(TPi,Ga,P);#B.S., the airside dry heat transfer coefficient can be got for different fin types

    phi = FinEffect_Schmidt(ho,233,P['th'],P['y'],P['Do']);#B.S. calculate the fin efficiency with the Schmidt equation

    P['Ro'] = 1/(ho*(P['Apo']+phi*P['Af']));
    P['ho']=ho;
    P['Ga']=Ga;

    # convert inlet state into TXP format
    TXPo = HPtoTXP(HPo);

    TXP_bak = TXPo.copy();#backup the inlet state of this segment

    # Mass flow rates of air and refrigerant
    ma=Ga*P['Aflow'];#B.S. Ga is the maximum airflow flux, P['Aflow'] is the minimun air flow cross-sectional area
    mr=Gr*P['Ax'];

    # Calculate heat transfered per unit mass of refrigerant.

    hi_max=10e4;#largest possible refrigerant side heat tranfer coeffcient
    hi_min=100;#minimum possible refrigerant side heat transfer coeffcient

    #store the parameters for iteration
    EVA_Get_Q['ma']=ma;
    EVA_Get_Q['mr']=mr;
    EVA_Get_Q['TPi']=TPi;
    EVA_Get_Q['Gr']=Gr;
    EVA_Get_Q['P']=P;
    EVA_Get_Q['q']=0;
    EVA_Get_Q['W']= HAPropsSI('W','T',TPi['T'],'P',101325,'R',TPi['P']) #wair.HumidityRatio(TPi.T,TPi.P); #[kg water/kg dry air]


    X1=0.05; X2=0.95;   # X1=0.01,X2=0.99;
    #for calculating the heat transfer and pressure drop in single-phase or approximated single-phase region
    if(TXPo['X']>X2 or TXPo['X']<X1):
        if (TXPo['X']>=1 or TXPo['X']<=0):#Pure single phase region (subcooled or superheated)
            y, P = ConvCoeffSP(TXPo,Gr,P, Ref) #get the single-phase heat transfer coefficient
        elif (TXPo['X']<X1): #transition from saturated liquid and subcooled
            TXP1 = toTXP(TXPo['T'],0.0,TXPo['P']);
            TXP2 = toTXP(TXPo['T'],X1,TXPo['P']);
            y1, P = ConvCoeffSP(TXP1,Gr,P, Ref);
            EVA_Get_Q['TXPo']=TXP2;#local refrigerant state
            #two-phase refrigerant heat transfer coefficient
            hi = brentq(Get_Q_EVA,hi_max,hi_min,args=(EVA_Get_Q),xtol=1e-5,rtol=6e-8,maxiter=40) #xtol, rtol and maxiter are changed to match "Zbrent" solver in ACMODEL
            y2 = hi/P['hRefAdj'];#ConvCoeffEvapTP_microfin(TXP2,Gr,P)
            if (y1>y2):
                y1=y2;#B.S., make the single phase flow less than two phase
            y = y1+TXPo['X']*(y2-y1)/X1;#B.S., get the heat transfer coefficient in this region with intropolation
        else: #transition from saturated vapor and superheat
            TXP1 = toTXP(TXPo['T'],X2,TXPo['P']);
            TXP2 = toTXP(TXPo['T'],1,TXPo['P']);
            EVA_Get_Q['TXPo'] = TXP1;
            #two-phase refrigerant heat transfer coefficient
            hi = brentq(Get_Q_EVA,hi_max,hi_min,args=(EVA_Get_Q),xtol=1e-5,rtol=6e-8,maxiter=40) #xtol, rtol and maxiter are changed to match "Zbrent" solver in ACMODEL
            #B.S., get the heat transfer coefficient at the two-phase region
            y1 = hi/P['hRefAdj'];#B.S. get two-phase heat transfer coefficient at this region    
            y2, P = ConvCoeffSP(TXP2,Gr,P, Ref);#B.S., single-phase heat transfer coefficient
            if (y2>y1):
                y2=y1; #B.S., for making the single phase flow less than two phase
            y = y2-(1-TXPo['X'])*(y2-y1)/(1-X2);#B.S., get the heat transfer coefficient in this region with intropolation
    
        hi = P['hRefAdj']*y;#B.S., heat transfer coefficient at this region
        Ri = 1/(hi*P['Api']);#B.S., inside thermal resistance
        R_W=log(P['Do']/(P['Do']-2.0*P['xp']))/(2.0*3.1415*P['K_T']*P['Ls']);
        R = P['Ro']+R_W+Ri;#overall thermal resistance
    
        #B.S., prepare the parameters for iteration
        EVA_Get_Q['TXPo'] = TXPo;
        EVA_Get_Q['Cmin'] = CmineCrossFlow_dry(R,mr,ma,TXPo,TPi, Ref);
        EVA_Get_Q['HPo'] = HPo;
    
        if (0):#TXPo['T']<T_sat+1
            # TXPo['T'] is actually the outlet refrigerant temperature
            # however TXPo['T'] = TXPi['T'] in the two phase region and
            # TXPo['T'] ~ to TXPi['T'] in the superheated region if steps are small
            q = (TPi['T']-TXPo['T'])/mr*CmineCrossFlow_dry(R,mr,ma,TXPo,TPi, Ref);
        else: #when it is close to the saturated state, the function zbrent can not converge
            EVA_Get_Q['hi']=hi;#input the refrigerant side heat transfer coefficient
            HPo, EVA_Get_Q = Get_Q_Single_For(HPo, EVA_Get_Q);   
            q=EVA_Get_Q['q'];
    
    
        if (TXPo['X']>=0.9999999): #B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.    
            R = (P['Ro']+R_W+Ri);
            P['UA_Vap']=P['UA_Vap']+1/R;#B.S., vapor heat transfer conductance
            P['VapL'] = P['VapL']+ P['Ls'];#B.S., vapor length
        elif (TXPo['X']<=0.00000001):
            R = (P['Ro'] + R_W+ Ri);
            P['LiqL'] = P['LiqL']+ P['Ls'];#B.S., keep the liquid length in the evaporator
            P['UA_Liq']=P['UA_Liq']+1/R;#B.S. liquid heat transfer conductance    
        else:
            R = P['Ro']+ R_W+ Ri;
            P['UA_TP']=P['UA_TP']+1/R;#B.S., dry thermal resistance of the two-phase heat transfer
            P['TPL'] = P['TPL']+ P['Ls'];#B.S., length of two-phase flow


    
    else: #calculating the heat transfer in the two-phase region
        EVA_Get_Q['TXPo']=TXPo;
        #iterate the heat transfer coefficient
        try:
            hi = brentq(Get_Q_EVA,hi_max,hi_min,args=(EVA_Get_Q),xtol=1e-5,rtol=6e-8,maxiter=40) #xtol, rtol and maxiter are changed to match "Zbrent" solver in ACMODEL
        except: #if failed, try again
            EVA_Get_Q['TXPo']['X'] = EVA_Get_Q['TXPo']['X']+0.05
            hi = brentq(Get_Q_EVA,hi_max,hi_min,args=(EVA_Get_Q),xtol=1e-5,rtol=6e-8,maxiter=40)
            EVA_Get_Q['TXPo']['X'] = EVA_Get_Q['TXPo']['X']-0.05;
            dhi = Get_Q_EVA(hi,EVA_Get_Q)
            ## plot error vs. T to see whgat function looks like that caused error
            #ZbrentPlot(hi_max,hi_min,Get_Q_EVA,1e-4,&EVA_Get_Q);#Haorong change from -7 to-2
            print("EVAP::EvapTubeL_Fwd (brentq) dhi= ",str(dhi))
            print("EVAP::EvapTubeL_Fwd (brentq) h={:f}, h_min={:f}, h_max={:f}".format(hi,hi_min,hi_max));


        q=EVA_Get_Q['q'];

        P['Qtp_wet'] = P['Qtp_wet']+EVA_Get_Q['P']['Qtp_wet'];#wet heat transfer amount
        P['Qtp_dry'] =P['Qtp_dry']+EVA_Get_Q['P']['Qtp_dry'];#dry heat transfer amount
        P['UAw_TP'] = EVA_Get_Q['P']['UAw_TP']+P['UAw_TP'];#heat conductance of two-phase wet heat transfer
        P['UA_TP'] = EVA_Get_Q['P']['UA_TP']+P['UA_TP'];#heat conductance of two-phase dry heat transfer
        P['L_wet'] = EVA_Get_Q['P']['L_wet'] + P['L_wet'];#wet heat transfer length in the two-phase region
        P['L_dry'] = EVA_Get_Q['P']['L_dry'] + P['L_dry'];#wet heat transfer length in the two-phase region
        P['TPL'] = P['TPL']+ P['Ls'];#length of the two-phase heat transfer


    HPo['H']+=q;

    #===========================================================================
    # calculate pressure drop
    #===========================================================================
    Gr = mr/P['Acs'];#this mass flux is for calculating the pressure drop
    DP_FR, P = FricDP(TXPo, Gr, q, P);

    if (TXPo['X']>0.05 and TXPo['X']<0.99):
        #B.S., prepared to calculate the two-phase acceleration pressure drop 
        Preacc = PreAcc()
        Preacc['DP_FR']=DP_FR;
        Preacc['G']=Gr;
        Preacc['H_OUT']=HPo['H'];
        Preacc['P_IN']=TXPo['P'];
        Preacc['X_IN']=TXPo['X'];
        
        #B.S., calculate the acceleration pressure drop
        DP_ACC = brentq(GET_PreAcc,10,-10,args=(Preacc, Ref),xtol=1e-7,rtol=6e-8,maxiter=40) 
        P_out=TXPo['P']-(DP_FR+DP_ACC)*P['PRefAdj'];
    else:
        P_out=TXPo['P']-DP_FR*P['PRefAdj'];

    HPo['P']=P_out;

    # Determine mass of charge.  It is based on the inlet state.
    # The specific volume is recalculated so that a different model
    # can be used from the one used to calculate the pressure drop.

    v = VolumeALL(TXPo,Gr,P['Di'],mr*q/P['Api'], Ref);        # seperate flow model
    #v = 1/PropertyTXP('D',TXPo, Ref);        # homogeneous flow model

    m['V']=P['Ls']*P['Acs'];#B.S., use the real cross-sectional area to calculate the inner volume
    m['m']=m['V']/v;

    TXPo = HPtoTXP(HPo);

    # Calculate output air state
    hai = HAPropsSI('Hha','T',TPi['T'],'P',101325,'R',TPi['P']) #wair.h(TPi.T,TPi.P); #[J/kg humid air]
    W_I = HAPropsSI('W','T',TPi['T'],'P',101325,'R',TPi['P']) #wair.HumidityRatio(TPi.T,TPi.P); #[kg water/kg dry air]
    HF_water=(EVA_Get_Q['T_S_O']*4.1877+0.0594)*1e3;

    WHo['H']=hai-q*mr/ma-HF_water*(W_I-EVA_Get_Q['W']);

    WHo['W'] = EVA_Get_Q['W'];#B.S. the outlet air humidity is calculated with the fuction Get_Q_EVA below.

    
    #B.S.----------------------------------------
    if (TXP_bak['X']>=0.999999999999):#B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.    
        P['m_Vap'] = P['m_Vap']+m['V']/v;#B.S.vapor mass
        P['V_Vap'] = P['V_Vap'] + m['V'];#B.S.vapor volume
    elif (TXP_bak['X']<=0.0000001):
        P['m_Liq'] = P['m_Liq']+m['V']/v;#B.S., liquid mass
        P['V_Liq'] = P['V_Liq'] + m['V'];#B.S., liquid volume
        if (TXPo['X']>0.0000001):
            P['HP_TP1']['H']=P['HP_TP1']['H']+mr*HPo['H'];#B.S., corresponding to the first segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
            P['HP_TP1']['P']=P['HP_TP1']['P']+mr*HPo['P'];#B.S., corresponding to the first segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
            P['count1'] = P['count1']+mr;#B.S., this addition is for counting how many evaporator circuits that involved with heat transfer calculation
    else: 
        P['m_TP'] = P['m_TP']+m['V']/v;#B.S., mass of two-phase flow
        P['V_TP'] = P['V_TP'] + m['V'];#B.S., volume of two-phase flow
        if (TXPo['X']>=0.999999999999):
            P['HP_TP2']['H']=P['HP_TP2']['H']+mr*HPo['H'];#B.S. corresponding to the last segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
            P['HP_TP2']['P']=P['HP_TP2']['P']+mr*HPo['P'];#B.S. corresponding to the last segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
            P['count2'] = P['count2']+mr;#B.S., this addition is for counting how many evaporator circuits that involved with heat transfer calculation
    #---------------------------------------------B.S.

    return (HPo, WHo, m, P)