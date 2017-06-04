from __future__ import division, print_function, absolute_import
from math import pi,log,exp,tanh

from scipy.optimize import brentq #solver to find roots (zero points) of functions
#import numpy as np
import pandas as pd

from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI

from ACHP.convert_units import *

from extra_functions import CGP, TXPtoHP, HPtoTXP, PropertyTXPth, PropertyTXPtr
# from PRESSURE import 
# from VOLUME import 
from CORR import Circuit_DP_COND
# from CMINE import 
            
        
        
def EvapTubeBend_Rev(Gr,HPo,m,P, Ref):
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
    return 0

def EvapTubeBend_Fwd(Gr,HPo,m,P, Ref):
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

    return 0
    
def EvapTubeL_Rev(Gr,#refrigerant mass flux
                  HPo,#refrigerant outlet or inlet state
                  Ga,#air mass flux
                  TPi,#air inlet state
                  WHo,#air outlet state
                  m,#charge and inner volume in the evaporator
                  P,#evaporator struct
                  Ref):#refrigerant
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
 
    EVA_Get_Q = EVA_Get_Q_dic();#B.S., this struct is for storing some parameters for iteration
    y=0;
    hi=0;
    Ri=0;
    R=0;
    q=0;
    
    ### Calculate air side resistance ###
    P['wet']=0;
    
    ho = ConvCoeffAir_EVA(TPi,Ga,P);#B.S., the airside dry heat transfer coefficient can be got for different fin types
    ho = P['hAirAdj']*ho
 
    phi = FinEffect_Schmidt(ho,233,P['th'],P['y'],P['Do']);#B.S. calculate the fin efficiency with the Schmidt equation
    
    P['Ro'] = 1/(ho*(P['Apo']+phi*P['Af']));
    P['ho']=ho;
    P['Ga']=Ga;
     
    # convert inlet state into TXP format
    TXPo = HPtoTXP(HPo,Ref);
    
    #B.S.--------------------------------------
    TXP_bak = TXPo.copy();#backup the outlet state or inlet of this segment
    #---------------------------------------B.S.
 
    # Mass flow rates of air and refrigerant
    ma=Ga*P['Aflow'];#B.S. Ga is the maximum airflow flux, P->Aflow is the minimun air flow cross-sectional area
    mr=Gr*P['Ax'];
 
    #===========================================================================
    # Calculate heat transfered per unit mass of refrigerant
    #===========================================================================
    hi_max=100000;#largest possible refrigerant side heat transfer coefficient
    hi_min=100;#minimum possible refrigerant side heat transfer coefficient
 
    #B.S.-----------------------------
    #store the parameters for iteration
    EVA_Get_Q['ma']=ma;
    EVA_Get_Q['mr']=mr;
    EVA_Get_Q['TPi']=TPi.copy();
    EVA_Get_Q['Gr']=Gr;
    EVA_Get_Q['P']=P.copy();
    EVA_Get_Q['q']=0;
    EVA_Get_Q['W']=HAPropsSI('W','T',TPi['T'],'P',101325,'R',TPi['P']) #[kg water/kg dry air]
    #---------------------------------------B.S.
 
    X1=0.05; X2=0.95;       
    #X1=0.01; X2=0.99;
    #for calculating the heat transfer and pressure drop in single-phase or approximated single-phase region
    if(TXPo['X']>X2 or TXPo['X']<X1):
        if(TXPo['X']>=1 or TXPo['X']<=0): #Pure single phase region
            y = ConvCoeffSP(TXPo,Gr,P, Ref);#get the single-phase heat transfer coefficient
        elif (TXPo['X']<X1):
            TXP1 = toTXP(TXPo['T'],0.0,TXPo['P']);
            TXP2 = toTXP(TXPo['T'],X1,TXPo['P']);
            y1 = ConvCoeffSP(TXP1,Gr,P, Ref);
            EVA_Get_Q['TXPo']=TXP2.copy();#local refrigerant state
            #B.S., get the heat transfer coefficient at the two-phase region
            hi = brentq(Get_Q_EVA,hi_max,hi_min,args=(Ref,EVA_Get_Q),xtol=1e-5,rtol=6e-8,maxiter=40) #xtol, rtol and maxiter are changed to match "Zbrent" solver in ACMODEL
            y2 = hi/P['hRefAdj'];#B.S., calculate the two-phase heat transfer coefficient at this region
             
            if(y1>y2):
                y1=y2;#B.S., make the single phase flow less than two phase
            y = y1+TXPo['X']*(y2-y1)/X1;#B.S., get the heat transfer coefficient in this region with intropolation
        else:
            TXP1 = toTXP(TXPo['T'],X2,TXPo['P']);
            TXP2 = toTXP(TXPo['T'],1,TXPo['P']);
            EVA_Get_Q['TXPo'] = TXP1.copy();
            try:
                #two-phase refrigerant heat transfer coefficient
                hi = brentq(Get_Q_EVA,hi_max,hi_min,args=(Ref,EVA_Get_Q),xtol=1e-5,rtol=6e-8,maxiter=40) #xtol, rtol and maxiter are changed to match "Zbrent" solver in ACMODEL
            except:
                print("EVAP::EvapTubeL_Rev (brentq) TXPo['T']={:f}, h_max={:f}, h_min={:f}".format(TXPo['T'],hi_max,hi_min));

            #B.S., get the heat transfer coefficient at the two-phase region
            y1 = hi/P['hRefAdj'];#ConvCoeffEvapTP_microfin(TXP1,Gr,P);#B.S. get two-phase heat transfer coefficient at this region   
            y2 = ConvCoeffSP(TXP2,Gr,P, Ref);#B.S., single-phase heat transfer coefficient
            if (y2>y1):
                y2=y1; #B.S., for making the single phase flow less than two phase
            y = y2-(1-TXPo['X'])*(y2-y1)/(1-X2);#B.S., get the heat transfer coefficient in this region with intropolation
 
        
        hi = P['hRefAdj']*y;#B.S., heat transfer coefficient at this region
        Ri = 1/(hi*P['Api']);#B.S., inside thermal resistance
        R_W=log(P['Do']/(P['Do']-2.0*P['xp']))/(2.0*pi*P['K_T']*P['Ls']);
        R = P['Ro']+R_W+Ri;#overall thermal resistance
        
        #B.S., prepare the parameters for iteration
        EVA_Get_Q['TXPo'] = TXPo.copy();
        EVA_Get_Q['Cmin'] = CmineCrossFlow_dry(R,mr,ma,TXPo,TPi, Ref);
        EVA_Get_Q['HPo'] = HPo.copy();
        
        try:
            if(0): #(TXPo.T<T_sat+1)
                # TXPo.T is actually the outlet refrigerant temperature
                # however TXPo.T = TXPi.T in the two phase region and
                # TXPo.T ~ to TXPi.T in the superheated region if steps are small
                q = (TPi['T']-TXPo['T'])/mr*CmineCrossFlow_dry(R,mr,ma,TXPo,TPi, Ref);
            else:#when it is close to the saturated state, the function zbrent can not converge
                H_max=HPo['H'];
                H_min=HPo['H']-4e4;
                EVA_Get_Q['hi']=hi;#input the refrigerant side heat transfer coefficient
                NN = brentq(Get_Q_Single,H_max,H_min,args=(Ref,EVA_Get_Q),xtol=1e-7,rtol=6e-8,maxiter=40)#B.S., to calculate the single-phase heat transfer and pressure drop at this region
                q=EVA_Get_Q['q'];
        except:
            print('EvapTubeL_Rev:: Exception :: Get_Q_Single')
            q = (TPi['T']-TXPo['T'])/mr*CmineCrossFlow_dry(R,mr,ma,TXPo,TPi, Ref); ####don't delete this
            
        #when the superheat is setted too high, there might exist a problem
        #if(q<=0): #B.S., to remove the wrong result
        #    q=0;
        ####don't delete this
   
        #B.S.---------------------------------------------
        if (TXPo['X']>=0.9999999): #B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.    
            R = (P['Ro']+Ri);
            P['UA_Vap']=P['UA_Vap']+1/R;#B.S., vapor heat transfer conductance
            P['VapL'] = P['VapL']+ P['Ls'];#B.S., vapor length
        elif (TXPo['X']<=0.00000001):
            R = (P['Ro'] + Ri);
            P['LiqL'] = P['LiqL']+ P['Ls'];#B.S., keep the liquid length in the evaporator
            P['UA_Liq']=P['UA_Liq']+1/R;#B.S. liquid heat transfer conductance    
        else:
            R = P['Ro']+Ri;
            P['UA_TP']=P['UA_TP']+1/R;#B.S., dry thermal resistance of the two-phase heat transfer
            P['TPL'] = P['TPL']+ P['Ls'];#B.S., length of two-phase flow
        #-----------------------------------------B.S.
 
    else: #calculating the heat transfer in the two-phase region
        EVA_Get_Q['TXPo']=TXPo.copy();
        try:
            hi = brentq(Get_Q_EVA,hi_max,hi_min,args=(Ref,EVA_Get_Q),xtol=1e-5,rtol=6e-8,maxiter=40) #iterate the heat transfer coefficient
        except:
            #plot error vs. T to see whgat function looks like that caused error
            #ZbrentPlot(hi_max,hi_min,Get_Q_EVA,1e-4,&EVA_Get_Q);#Haorong change from -7 to-2
            print("EVAP::EvapTubeL_Rev (brentq) h={:f}, h_min={:f}, h_max={:f}".format(hi,hi_min,hi_max))

        
        q=EVA_Get_Q['q'];
 
        #B.S.-----------------------------------
        P['Qtp_wet'] = P['Qtp_wet']+EVA_Get_Q['P']['Qtp_wet'];#wet heat transfer amount
        P['Qtp_dry'] =P['Qtp_dry']+EVA_Get_Q['P']['Qtp_dry'];#dry heat transfer amount
        P['UAw_TP'] = EVA_Get_Q['P']['UAw_TP']+P['UAw_TP'];#heat conductance of two-phase wet heat transfer
        P['UA_TP'] = EVA_Get_Q['P']['UA_TP']+P['UA_TP'];#heat conductance of two-phase dry heat transfer
        P['L_wet'] = EVA_Get_Q['P']['L_wet'] + P['L_wet'];#wet heat transfer length in the two-phase region
        P['L_dry'] = EVA_Get_Q['P']['L_dry'] + P['L_dry'];#wet heat transfer length in the two-phase region
        P['TPL'] = P['TPL']+ P['Ls'];#length of the two-phase heat transfer
        #------------------------------------------B.S.
 
 
    HPo['H']=HPo['H']-q;
 
    #===========================================================================
    # calculate pressure drop
    #===========================================================================
    Gr = mr/P['Acs'];#this mass flux is for calculating the pressure drop
    DP_FR = FricDP(TXPo, Gr, q, P, Ref);
  
    if (TXPo['X']>0.05 and TXPo['X']<0.99):
    #B.S., prepared to calculate the two-phase acceleration pressure drop
        Preacc=PreAcc()
        Preacc['DP_FR']=-1*DP_FR;
        Preacc['G']=Gr;
        Preacc['H_OUT']=HPo['H'];
        Preacc['P_IN']=TXPo['P'];
        Preacc['X_IN']=TXPo['X'];
        DP_ACC=brentq(GET_PreAcc,10000,-10000,args=(Ref, Preacc),xtol=1e-7,rtol=6e-8,maxiter=40)#B.S., calculate the acceleration pressure drop 
        P_out=TXPo['P']+DP_FR-DP_ACC;   
    else:
        P_out=TXPo['P']+DP_FR;
    
    HPo['P']=P_out;
    
    #===========================================================================
    # Determine mass of charge.  
    #===========================================================================
    # It is based on the inlet state.
    # The specific volume is recalculated so that a different model
    # can be used from the one used to calculate the pressure drop.
    TXPo = HPtoTXP(HPo,Ref);

    v = VolumeALL(TXPo,Gr,P['Di'],mr*q/P['Api'], Ref); #seperate flow model
    #v = PropertyTXP(VOL,TXPo); #homogeneous flow model
    
    m['V']=P['Ls']*P['Acs'];#B.S., use the real cross-sectional area to calculate the inner volume
    m['m']=m['V']/v;
 
    #===========================================================================
    # Calculate output air state
    #===========================================================================
    hai = HAPropsSI('H','T',TPi['T'],'P',101325,'R',TPi['P']) #[J/kg humid air]
    W_I = HAPropsSI('W','T',TPi['T'],'P',101325,'R',TPi['P']) #[kg water/kg dry air]
    #HF_water=(K2C(EVA_Get_Q['T_S_O'])*4.1877+0.0594)*1e3; #correction for saturated liquid water enthalpy using T_S_O 
    HF_water = HAPropsSI('H','T',EVA_Get_Q['T_S_O'],'P',101325,'R',1) #A.B much faster convergence using HAProps to find HF_water
    
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
 
    return 0

def EvapTubeL_Fwd(Gr,#refrigerant mass flux
                  HPo,#refrigerant outlet and inlet state
                  Ga,#air mass flux
                  TPi,#air inlet state
                  WHo,#air outlet state
                  m,#charge and inner volume in the evaporator
                  P,#evaporator struct
                  Ref): #refrigerant
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
    
    #B.S., this dictionary struct is for storing some parameters for iteration
    EVA_Get_Q = EVA_Get_Q_dic()

    ### Calculate air side resistance ###
    P['wet']=0;#without considering wet heat tranfer adjustment

    ho = ConvCoeffAir_EVA(TPi,Ga,P);#B.S., the airside dry heat transfer coefficient can be got for different fin types
    ho = P['hAirAdj']*ho
    
    phi = FinEffect_Schmidt(ho,233,P['th'],P['y'],P['Do']);#B.S. calculate the fin efficiency with the Schmidt equation

    P['Ro']=1/(ho*(P['Apo']+phi*P['Af']));
    P['ho']=ho;
    P['Ga']=Ga;

    # convert inlet state into TXP format
    TXPo = HPtoTXP(HPo, Ref);

    TXP_bak = TXPo.copy();#backup the inlet state of this segment

    # Mass flow rates of air and refrigerant
    ma=Ga*P['Aflow'];#B.S. Ga is the maximum airflow flux, P['Aflow'] is the minimun air flow cross-sectional area
    mr=Gr*P['Ax'];

    #===========================================================================
    # Calculate heat transfered per unit mass of refrigerant
    #===========================================================================
    hi_max=10e4;#largest possible refrigerant side heat tranfer coefficient
    hi_min=100;#minimum possible refrigerant side heat transfer coefficient

    #store the parameters for iteration
    EVA_Get_Q['ma']=ma;
    EVA_Get_Q['mr']=mr;
    EVA_Get_Q['TPi']=TPi.copy();
    EVA_Get_Q['Gr']=Gr;
    EVA_Get_Q['P']=P.copy();
    EVA_Get_Q['q']=0;
    EVA_Get_Q['W']= HAPropsSI('W','T',TPi['T'],'P',101325,'R',TPi['P']) #wair.HumidityRatio(TPi.T,TPi.P); #[kg water/kg dry air]

    #X1=0.01; X2=0.99;
    X1=0.05; X2=0.95;
    #for calculating the heat transfer and pressure drop in single-phase or approximated single-phase region
    if(TXPo['X']>X2 or TXPo['X']<X1):
        if (TXPo['X']>=1 or TXPo['X']<=0):#Pure single phase region (subcooled or superheated)
            y = ConvCoeffSP(TXPo,Gr,P, Ref) #get the single-phase heat transfer coefficient
        elif (TXPo['X']<X1): #transition from saturated liquid and subcooled
            TXP1 = toTXP(TXPo['T'],0.0,TXPo['P']);
            TXP2 = toTXP(TXPo['T'],X1,TXPo['P']);
            y1 = ConvCoeffSP(TXP1,Gr,P, Ref);
            EVA_Get_Q['TXPo']=TXP2.copy();#local refrigerant state
            #two-phase refrigerant heat transfer coefficient
            hi = brentq(Get_Q_EVA,hi_max,hi_min,args=(Ref,EVA_Get_Q),xtol=1e-5,rtol=6e-8,maxiter=40) #xtol, rtol and maxiter are changed to match "Zbrent" solver in ACMODEL
            y2 = hi/P['hRefAdj'];#ConvCoeffEvapTP_microfin(TXP2,Gr,P)
            if (y1>y2):
                y1=y2;#B.S., make the single phase flow less than two phase
            y = y1+TXPo['X']*(y2-y1)/X1;#B.S., get the heat transfer coefficient in this region with intropolation
        else: #transition from saturated vapor and superheat
            TXP1 = toTXP(TXPo['T'],X2,TXPo['P']);
            TXP2 = toTXP(TXPo['T'],1,TXPo['P']);
            EVA_Get_Q['TXPo'] = TXP1.copy();
            #two-phase refrigerant heat transfer coefficient
            hi = brentq(Get_Q_EVA,hi_max,hi_min,args=(Ref,EVA_Get_Q),xtol=1e-5,rtol=6e-8,maxiter=40) #xtol, rtol and maxiter are changed to match "Zbrent" solver in ACMODEL
            #B.S., get the heat transfer coefficient at the two-phase region
            y1 = hi/P['hRefAdj'];#B.S. get two-phase heat transfer coefficient at this region    
            y2 = ConvCoeffSP(TXP2,Gr,P, Ref);#B.S., single-phase heat transfer coefficient
            if (y2>y1):
                y2=y1; #B.S., for making the single phase flow less than two phase
            y = y2-(1-TXPo['X'])*(y2-y1)/(1-X2);#B.S., get the heat transfer coefficient in this region with intropolation
    
        hi = P['hRefAdj']*y;#B.S., heat transfer coefficient at this region
        Ri = 1/(hi*P['Api']);#B.S., inside thermal resistance
        R_W=log(P['Do']/(P['Do']-2.0*P['xp']))/(2.0*3.1415*P['K_T']*P['Ls']);
        R = P['Ro']+R_W+Ri;#overall thermal resistance
    
        #B.S., prepare the parameters for iteration
        EVA_Get_Q['TXPo'] = TXPo.copy();
        EVA_Get_Q['Cmin'] = CmineCrossFlow_dry(R,mr,ma,TXPo,TPi, Ref);
        EVA_Get_Q['HPo'] = HPo.copy();
    
        if (0):#TXPo['T']<T_sat+1
            # TXPo['T'] is actually the outlet refrigerant temperature
            # however TXPo['T'] = TXPi['T'] in the two phase region and
            # TXPo['T'] ~ to TXPi['T'] in the superheated region if steps are small
            q = (TPi['T']-TXPo['T'])/mr*CmineCrossFlow_dry(R,mr,ma,TXPo,TPi, Ref);
        else: #when it is close to the saturated state, the function zbrent can not converge
            EVA_Get_Q['hi']=hi;#input the refrigerant side heat transfer coefficient
            Get_Q_Single_For(HPo, Ref, EVA_Get_Q);   
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
        EVA_Get_Q['TXPo']=TXPo.copy();
        #iterate the heat transfer coefficient
        try:
            hi = brentq(Get_Q_EVA,hi_max,hi_min,args=(Ref,EVA_Get_Q),xtol=1e-5,rtol=6e-8,maxiter=40) #xtol, rtol and maxiter are changed to match "Zbrent" solver in ACMODEL
        except: #if failed, try again
            EVA_Get_Q['TXPo']['X'] = EVA_Get_Q['TXPo']['X']+0.05
            hi = brentq(Get_Q_EVA,hi_max,hi_min,args=(Ref,EVA_Get_Q),xtol=1e-5,rtol=6e-8,maxiter=40)
            EVA_Get_Q['TXPo']['X'] = EVA_Get_Q['TXPo']['X']-0.05;
            dhi = Get_Q_EVA(hi,Ref,EVA_Get_Q)
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


    HPo['H'] = HPo['H']+q;

    #===========================================================================
    # calculate pressure drop
    #===========================================================================
    Gr = mr/P['Acs'];#this mass flux is for calculating the pressure drop
    DP_FR = FricDP(TXPo, Gr, q, P, Ref);

    if (TXPo['X']>0.05 and TXPo['X']<0.99):
        #B.S., prepared to calculate the two-phase acceleration pressure drop 
        Preacc = PreAcc()
        Preacc['DP_FR']=DP_FR;
        Preacc['G']=Gr;
        Preacc['H_OUT']=HPo['H'];
        Preacc['P_IN']=TXPo['P'];
        Preacc['X_IN']=TXPo['X'];
        
        #B.S., calculate the acceleration pressure drop [Pa]
        DP_ACC = brentq(GET_PreAcc,10000,-10000,args=(Ref, Preacc),xtol=1e-7,rtol=6e-8,maxiter=40) 
        P_out=TXPo['P']-(DP_FR+DP_ACC)*P['PRefAdj'];
    else:
        P_out=TXPo['P']-DP_FR*P['PRefAdj'];

    HPo['P']=P_out;

    #===========================================================================
    # Determine mass of charge.  
    #===========================================================================
    #It is based on the inlet state.
    # The specific volume is recalculated so that a different model
    # can be used from the one used to calculate the pressure drop.

    v = VolumeALL(TXPo,Gr,P['Di'],mr*q/P['Api'], Ref);        # seperate flow model
    #v = PropertyTXP(VOL,TXPo);        # homogeneous flow model

    m['V']=P['Ls']*P['Acs'];#B.S., use the real cross-sectional area to calculate the inner volume
    m['m']=m['V']/v;

    TXPo = HPtoTXP(HPo, Ref);

    #===========================================================================
    # Calculate output air state
    #===========================================================================
    hai = HAPropsSI('H','T',TPi['T'],'P',101325,'R',TPi['P']) #[J/kg humid air]
    W_I = HAPropsSI('W','T',TPi['T'],'P',101325,'R',TPi['P']) #[kg water/kg dry air]
    #HF_water=(K2C(EVA_Get_Q['T_S_O'])*4.1877+0.0594)*1e3; #correction for saturated liquid water enthalpy using T_S_O 
    HF_water = HAPropsSI('H','T',EVA_Get_Q['T_S_O'],'P',101325,'R',1) #A.B much faster convergence using HAProps to find HF_water
    
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

    return 0

def Get_Q_EVA(hi_0,#this function is for getting the refrigerant side heat transfer coefficent 
              Ref,
              Params=None):#this struct contains the necessary parameters for iteration
    '''
    /********************************************************************
    B.S., add for iteration get evaporative heat transfer at two-phase region
    wet coil calculation from 
    Braun, J. E., Klein, S.A., and Michell, J.W., 1989,"effectiveness models for cooling towers
    and cooling coils", ASHRAE Transactions, Vol. 95-2, pp. 164-174
    ***********************************************************************/
    '''
    
    if (Params == None):
        EVA_Q = EVA_Get_Q_dic()
    else:    
        EVA_Q = Params #dictionary that keep updating
    
    #B.S.--------------------------
    EVA_Q['P']['Qtp_wet'] = 0;#new
    EVA_Q['P']['Qtp_dry'] = 0;#new
    EVA_Q['P']['UAw_TP'] = 0;
    EVA_Q['P']['UA_TP'] = 0;
    EVA_Q['P']['L_wet'] = 0;
    EVA_Q['P']['L_dry'] = 0;
    #---------------------------B.S.
    
    hi=hi_0;

    W_I= HAPropsSI('W','P',101325,'T',EVA_Q['TPi']['T'],'R',EVA_Q['TPi']['P']) # inlet air humidity #[kg water/kg dry air]
        
    #Eva_dim = ETdim();
    Eva_dim = EVA_Q['P'].copy();#B.S. get the evaporator struct parameters

    #===========================================================================
    # dry condition
    #===========================================================================
    Ri = 1/(hi*EVA_Q['P']['Api']);#B.S., inside thermal resistance
    R_W=log(EVA_Q['P']['Do']/(EVA_Q['P']['Do']-2.0*EVA_Q['P']['xp']))/(2.0*pi*EVA_Q['P']['K_T']*EVA_Q['P']['Ls']);
    R = EVA_Q['P']['Ro']+R_W+Ri;#B.S., external thermal resistance under dry condition

    CP_M=HAPropsSI('cp','P',101325,'T',EVA_Q['TPi']['T'],'R',EVA_Q['TPi']['P']) #[J/kg humid air/K]
    
    #===========================================================================
    # begin with wet condition
    #===========================================================================
    H_A_I=HAPropsSI('H','P',101325,'T',EVA_Q['TPi']['T'],'R',EVA_Q['TPi']['P']) #inlet air enthalpy #[J/kg humid air]
    H_SAT=HAPropsSI('H','P',101325,'T',EVA_Q['TXPo']['T'],'R',0.999) #[J/kg humid air]
    enthal1=HAPropsSI('H','P',101325,'T',EVA_Q['TXPo']['T'],'R',0.999)  #[J/kg humid air]
    enthal2=HAPropsSI('H','P',101325,'T',EVA_Q['TXPo']['T']-0.1,'R',0.999) #[J/kg humid air]

    C_S=(enthal1-enthal2)/0.1e0;

    M_DOT_A=EVA_Q['ma'];
    EVA_Q['P']['wet']=1;
    
    H_W = ConvCoeffAir_EVA(EVA_Q['TPi'],EVA_Q['P']['Ga'],EVA_Q['P']);
    H_W=EVA_Q['P']['hAirAdj']*H_W
    EVA_Q['P']['wet']=0;
    
    #===========================================================================
    # calculate wet fin efficiency
    #===========================================================================
    M_F_W=pow((2e0*H_W*C_S/(EVA_Q['P']['K_F']*EVA_Q['P']['th']*CP_M)),0.5e0);#for calculating the wet fin efficiency, including both the heat and mass transfer
    ETA_F_W=tanh(M_F_W*EVA_Q['P']['L_F'])/(M_F_W*EVA_Q['P']['L_F']);#wet fin efficiency
    
    A_F=EVA_Q['P']['Af'];#fin surface area
    A_T=EVA_Q['P']['Af']+EVA_Q['P']['Apo'];#the whole external area
    ETA_O_W=1e0-A_F/A_T*(1e0-ETA_F_W);#overall external fin efficiency
    R_SURF_W=1/(ETA_O_W*H_W*A_T);#airside thermal resistance under wet condition
    NTU_O_W=1/(R_SURF_W*M_DOT_A*CP_M);#for calculating the effective saturated air enthalpy at T_S_O

    NTU_O=1/(EVA_Q['P']['Ro']*M_DOT_A*CP_M);#for calculating the outlet air temperature, which is involved with sensible heat transfer

    UA_W=1e0/(C_S*Ri+C_S*R_W+CP_M*R_SURF_W);#UA under wet condition, considering both the heat and mass transfer
    NTU_W=1*UA_W/(M_DOT_A);#NTU under wet condition, for calculating heat transfer
    EPSILON_W=1e0-exp(-NTU_W);#epsilon under wet condition

    Q=EPSILON_W*M_DOT_A*(H_A_I-H_SAT);#heat transfer under wet condition
    EVA_Q['q'] =Q/EVA_Q['mr'];

    #get outlet humidity and effective surface temperature
    H_A_O=H_A_I-Q/(M_DOT_A);#outlet air enthalpy, actually the energy balance should include the water flowing away, like -(83.84e3+(T_S-20)*4.183e3)*(W_I-W_O), since this part is small and need iteration, so it is ignored
    H_S_S_O=H_A_I-(H_A_I-H_A_O)/(1-exp(-NTU_O_W));#saturated air enthalpy at T_S_S_O, effective surface temperature
    

#     H_zero=wair.h(2,0.99);#fixed the minimum possible H_S_S_O
#     if (H_S_S_O<=H_zero): 
#         H_S_S_O=H_zero;
    #TP_S_O = {'T':0.0,'P':0.0}
    TP_S_O = HPtoTP(H_S_S_O,0.999);
    T_S_O=TP_S_O['T'];#effective surface temperature
    T_O=T_S_O+(EVA_Q['TPi']['T']-T_S_O)*exp(-NTU_O);#outlet air temperature
    EVA_Q['T_S_O']= T_S_O;

    #TP_dew = {'T':0.0,'P':0.0}
    TP_dew = WPtoTP(W_I,0.999);#get the dew point corresponding to the inlet humidity ratio
    T_dew=TP_dew['T'];#dew temperature corresponding to the inlet air enthalpy
    
    #===========================================================================
    # dry condition
    #===========================================================================
    if (T_S_O>=T_dew):
        # TXPo['T'] is actually the outlet refrigerant temperature
        # however TXPo['T'] = TXPi['T'] in the two phase region and
        # TXPo['T'] ~ to TXPi['T'] in the superheated region if steps are small
        q_dry=(EVA_Q['TPi']['T']-EVA_Q['TXPo']['T'])/EVA_Q['mr']*CmineCrossFlow_wet(R,EVA_Q['mr'],EVA_Q['ma'],EVA_Q['TXPo'],EVA_Q['TPi']['T'],CP_M, Ref);
        if (T_S_O>=T_dew+0.2):
            EVA_Q['q'] = (EVA_Q['TPi']['T']-EVA_Q['TXPo']['T'])/EVA_Q['mr']*CmineCrossFlow_wet(R,EVA_Q['mr'],EVA_Q['ma'],EVA_Q['TXPo'],EVA_Q['TPi']['T'],CP_M, Ref);
        else:
            EVA_Q['q']=(q_dry-Q/EVA_Q['mr'])/(0.2)*(T_S_O-T_dew) + Q/EVA_Q['mr'];

        Q=EVA_Q['q']*EVA_Q['mr'];#heat transfer amount under dry condition
        EVA_Q['W']=W_I;#without dehumidifying
        T_w=Q*Ri+EVA_Q['TXPo']['T'];
        EVA_Q['T_w']=T_w;
        #B.S.---------------------------
        EVA_Q['P']['Qtp_dry'] = Q;#dry heat transfer
        EVA_Q['P']['UA_TP'] = 1/R;#dry heat transfer conductance of this segment
        EVA_Q['P']['L_dry'] = EVA_Q['P']['Ls'];#dry heat transfer tube length of this segment
        #--------------------------------B.S.
        Eva_dim['T_w']=T_w;
        Eva_dim['q_flux']=Q/EVA_Q['P']['Api'];
        
        hi = ConvCoeffEvapTP_microfin(EVA_Q['TXPo'],EVA_Q['Gr'],Eva_dim, Ref);#B.S., based on the tube wall temperature to get the refrigerant side heat transfer coefficient
        hi = EVA_Q['P']['hRefAdj']*hi
        dhi_dry = (hi-hi_0)/(hi+hi_0);
        
        return dhi_dry

    TP_O = {'T':0.0,'P':0.0}
    TP_O['T']=T_O;
    HH_min = HAPropsSI('H','P',101325,'T',T_O,'R',0.05) #[J/kg humid air]
    HH_max = HAPropsSI('H','P',101325,'T',T_O,'R',0.999) #[J/kg humid air]
    
    if (H_A_O>=HH_max):
        TP_O['T']=T_O;
        TP_O['P']=0.999;
    elif (H_A_O<=HH_min):
        TP_O['T']=T_O;
        TP_O['P']=0.05;
    else:
        try:
            TP_O=THtoTP(T_O,H_A_O);#outlet air state
        except:
            print('Get_Q_EVA:: state above saturation')
            TP_O['T']=T_O;
            TP_O['P']=0.999;
    
    try:
        W_O=HAPropsSI('W','P',101325,'T',T_O,'R',TP_O['P']) #outlet air humidity ratio #[-]
    except:
        print('Get_Q_EVA:: W_O = W_I')
        W_O=W_I;

    EVA_Q['W']=W_O;
    T_w=Q*Ri+EVA_Q['TXPo']['T'];
    EVA_Q['T_w']=T_w;

    #B.S.---------------------------
    EVA_Q['P']['Qtp_wet'] = EVA_Q['mr']*EVA_Q['q'];#wet heat transfer
    EVA_Q['P']['UAw_TP'] = UA_W;#wet heat transfer conductance of this segment
    EVA_Q['P']['L_wet'] = EVA_Q['P']['Ls'];#wet heat transfer tube length of this segment
    #--------------------------------B.S.

    Eva_dim['T_w']=T_w;
    Eva_dim['q_flux']=Q/EVA_Q['P']['Api'];
    hi =ConvCoeffEvapTP_microfin(EVA_Q['TXPo'],EVA_Q['Gr'],Eva_dim, Ref);#B.S., based on the tube wall temperature to get the refrigerant side heat transfer coefficient
    hi = EVA_Q['P']['hRefAdj']*hi
    dhi=(hi-hi_0)/(hi+hi_0);
    
    return dhi

def Get_Q_Single(H_in, #inlet refrigerant temperature in the segment
                 Ref,
                 Params=None):#this struct stores the parameters for iteration
    '''
    #for iteration to get evaporative heat transfer at single-phase region, use the inlet refrigerant state as the reference state
    '''
    if (Params==None):
        EVA_Q=EVA_Get_Q_dic()
    else:
        EVA_Q=Params #dictionary that keep updating
    
    
    #TPi={'T':0.0,'P':0.0};
    #HPo={'H':0.0,'P':0.0};
    #TXPo={'T':0.0,'X':0.0,'P':0.0};
    
    #mr=EVA_Q['mr'];
    #q=0;
    #Q=0;
    #Ri=0;
    Gr=EVA_Q['Gr'];
    
    #Eva_dim=ETdim();
    Eva_dim = EVA_Q['P'].copy(); #B.S. get the evaporator struct parameters

    HPo=EVA_Q['HPo'].copy(); #refrigerant outlet state
    H_out=HPo['H']; #backup the outlet enthalpy
    TPi = EVA_Q['TPi'].copy(); #air inlet state
    TXPo = EVA_Q['TXPo'].copy();#refrigerant outlet state
    
    
    q=HPo['H']-H_in;
    
    DP_SP = FricDP(TXPo,Gr,q,Eva_dim, Ref);#single-phase pressure drop
    P_in=TXPo['P']+DP_SP;
    HPo['P'] = P_in;
    HPo['H'] = H_in;
    
    TXPo = HPtoTXP(HPo,Ref);
    
    W_I=HAPropsSI('W','P',101325,'T',EVA_Q['TPi']['T'],'R',EVA_Q['TPi']['P']) #B.S., inlet air humidity #[kg water/kg dry air]
       
    #===========================================================================
    # dry condition
    #===========================================================================
    Ri = 1/(EVA_Q['hi']*EVA_Q['P']['Api']);#B.S., inside thermal resistance
    R_W=log(EVA_Q['P']['Do']/(EVA_Q['P']['Do']-2.0*EVA_Q['P']['xp']))/(2.0*pi*EVA_Q['P']['K_T']*EVA_Q['P']['Ls']);
    R = EVA_Q['P']['Ro']+ R_W +Ri;#B.S., external thermal resistance under dry condition

    CP_M=HAPropsSI('cp','P',101325,'T',EVA_Q['TPi']['T'],'R',EVA_Q['TPi']['P']) # [J/kg humid air/K]
    
    #===========================================================================
    # begin with wet condition
    #===========================================================================
    H_A_I=HAPropsSI('H','P',101325,'T',EVA_Q['TPi']['T'],'R',EVA_Q['TPi']['P']) #inlet air enthalpy #[J/kg humid air]
    H_SAT=HAPropsSI('H','P',101325,'T',TXPo['T'],'R',0.999) #[J/kg humid air]
    enthal1=HAPropsSI('H','P',101325,'T',TXPo['T'],'R',0.999) #[J/kg humid air]
    enthal2=HAPropsSI('H','P',101325,'T',TXPo['T']-0.1,'R',0.999) #[J/kg humid air]

    C_S=(enthal1-enthal2)/0.1e0;
   
    M_DOT_A=EVA_Q['ma'];
    EVA_Q['P']['wet']=1;
    H_W = ConvCoeffAir_EVA(EVA_Q['TPi'],EVA_Q['P']['Ga'],EVA_Q['P']);
    H_W = EVA_Q['P']['hAirAdj']*H_W
    EVA_Q['P']['wet']=0;
    
    #===========================================================================
    # calculate wet fin efficiency
    #===========================================================================
    M_F_W=pow((2e0*H_W*C_S/(EVA_Q['P']['K_F']*EVA_Q['P']['th']*CP_M)),0.5e0);#for calculating the wet fin efficiency, including both the heat and mass transfer
    ETA_F_W=tanh(M_F_W*EVA_Q['P']['L_F'])/(M_F_W*EVA_Q['P']['L_F']);#wet fin efficiency
    
    A_F=EVA_Q['P']['Af'];#fin surface area
    A_T=EVA_Q['P']['Af']+EVA_Q['P']['Apo'];#the whole external area
    ETA_O_W=1e0-A_F/A_T*(1e0-ETA_F_W);#overall external fin efficiency
    R_SURF_W=1/(ETA_O_W*H_W*A_T);#airside thermal resistance under wet condition
    NTU_O_W=1/(R_SURF_W*M_DOT_A*CP_M);#for calculating the effective saturated air enthalpy at T_S_O

    NTU_O=1/(EVA_Q['P']['Ro']*M_DOT_A*CP_M);#for calculating the outlet air temperature, which is involved with sensible heat transfer

    UA_W=1e0/(C_S*Ri+C_S*R_W+CP_M*R_SURF_W);#UA under wet condition, considering both the heat and mass transfer
    NTU_W=1*UA_W/(M_DOT_A);#NTU under wet condition, for calculating heat transfer
    EPSILON_W=1e0-exp(-NTU_W);#epsilon under wet condition

    Q=EPSILON_W*M_DOT_A*(H_A_I-H_SAT);#heat transfer under wet condition
    EVA_Q['q'] =Q/EVA_Q['mr'];
    q=Q/EVA_Q['mr'];

    #get outlet humidity and effective surface temperature
    H_A_O=H_A_I-Q/(M_DOT_A);#outlet air enthalpy, actually the energy balance shoudl include the water flowing away, like -(83.84e3+(T_S-20)*4.183e3)*(W_I-W_O), since this part is small and need iteration, so it is ignored
    H_S_S_O=H_A_I-(H_A_I-H_A_O)/(1-exp(-NTU_O_W));#saturated air enthalpy at T_S_S_O, effective surface temperature
    
    
#        H_zero=wair.h(2,0.99);#fixed the minimum possible H_S_S_O
#        if (H_S_S_O<=H_zero):
#              H_S_S_O=H_zero;
    #TP_S_O={'T':0.0,'P':0.0};
    TP_S_O=HPtoTP(H_S_S_O,0.999);
    T_S_O=TP_S_O['T'];#effective surface temperature
    T_O=T_S_O+(EVA_Q['TPi']['T']-T_S_O)*exp(-NTU_O);#outlet air temperature
    EVA_Q['T_S_O']= T_S_O;

    #TP_dew={'T':0.0,'P':0.0};
    TP_dew=WPtoTP(W_I,0.999);#get the dew point corresponding to the inlet humidity ratio
    T_dew=TP_dew['T'];#dew temperature corresponding to the inlet air enthalpy
            
    if(T_S_O>T_dew): #dry condition
        EVA_Q['q'] = (EVA_Q['TPi']['T']-TXPo['T'])/EVA_Q['mr']*CmineCrossFlow_wet(R,EVA_Q['mr'],EVA_Q['ma'],EVA_Q['TXPo'],EVA_Q['TPi']['T'],CP_M, Ref);
        q=EVA_Q['q'];
        Q=EVA_Q['q']*EVA_Q['mr'];#heat transfer amount under dry condition
        EVA_Q['W']=W_I;#without dehumidifying

    TP_O={'T':0.0,'P':0.0};
    TP_O['T']=T_O;
    HH_min = HAPropsSI('H','P',101325,'T',T_O,'R',0.05)   #[J/kg humid air]
    HH_max = HAPropsSI('H','P',101325,'T',T_O,'R',0.999) #[J/kg humid air]
    
    if(H_A_O>=HH_max):
        TP_O['T']=T_O;
        TP_O['P']=0.999;
    elif(H_A_O<=HH_min):
        TP_O['T']=T_O;
        TP_O['P']=0.05;
    else:
        try:
            TP_O=THtoTP(T_O,H_A_O);#outlet air state
        except:
            print('Get_Q_Single:: state above saturation')
            TP_O['T']=T_O;
            TP_O['P']=0.999;

    try:
        W_O=HAPropsSI('W','P',101325,'T',T_O,'R',TP_O['P']) #outlet air humidity ratio #[-]
    except:
        print('Get_Q_Single:: W_O = W_I')
        W_O=W_I;

    EVA_Q['W']=W_O;

    HPo['H']=H_out - q;
    dH=(HPo['H']-H_in)/(HPo['H']+H_in);
    
    return dH

def Get_Q_Single_For(HPi, #inlet refrigerant temperature in the segment
                     Ref,
                     Params=None):#this struct stores the parameters for iteration
    '''
    #for iteration to get evaporative heat transfer at single-phase region, use the inlet refrigerant state as the reference state
    '''
    
    if (Params==None):
        EVA_Q = EVA_Get_Q_dic()
    else:
        EVA_Q = Params #dictionary that keep updating

    #TPi= {'T':0.0,'P':0.0};
    #TXPi={'T':0.0,'X':0.0,'P':0.0};
    #mr=EVA_Q['mr'];
    #Q=0;
    #Ri=0;
    #Gr=EVA_Q['Gr'];

    #Eva_dim = ETdim()
    #Eva_dim = EVA_Q['P'];#B.S. get the evaporator struct parameters

    #TPi = EVA_Q['TPi'].copy();#air inlet state
    
    TXPi = HPtoTXP(HPi, Ref);
    
    W_I=HAPropsSI('W','P',101325,'T',EVA_Q['TPi']['T'],'R',EVA_Q['TPi']['P']) #B.S., inlet air humidity #[kg water/kg dry air]
        
    #===========================================================================
    # dry condition
    #===========================================================================
    Ri = 1/(EVA_Q['hi']*EVA_Q['P']['Api']);#B.S., inside thermal resistance
    R_W=log(EVA_Q['P']['Do']/(EVA_Q['P']['Do']-2.0*EVA_Q['P']['xp']))/(2.0*pi*EVA_Q['P']['K_T']*EVA_Q['P']['Ls']);
    R = EVA_Q['P']['Ro']+ R_W +Ri;#B.S., external thermal resistance under dry condition
    
    CP_M=HAPropsSI('cp','P',101325,'T',EVA_Q['TPi']['T'],'R',EVA_Q['TPi']['P']) #[J/kg humid air/K]

    #===========================================================================
    # begin with wet condition
    #===========================================================================
    H_A_I=HAPropsSI('H','P',101325,'T',EVA_Q['TPi']['T'],'R',EVA_Q['TPi']['P']) #inlet air enthalpy #[J/kg humid air]
    H_SAT=HAPropsSI('H','P',101325,'T',TXPi['T'],'R',0.999) #[J/kg humid air]
    enthal1=HAPropsSI('H','P',101325,'T',TXPi['T'],'R',0.999) #[J/kg humid air]
    enthal2=HAPropsSI('H','P',101325,'T',TXPi['T']-0.1,'R',0.999) #[J/kg humid air]

    C_S=(enthal1-enthal2)/0.1e0;
    
    M_DOT_A=EVA_Q['ma'];
    EVA_Q['P']['wet']=1;
    H_W = ConvCoeffAir_EVA(EVA_Q['TPi'],EVA_Q['P']['Ga'],EVA_Q['P'])
    H_W = EVA_Q['P']['hAirAdj']*H_W
    EVA_Q['P']['wet']=0;
    
    #===========================================================================
    # calculate wet fin efficiency
    #===========================================================================
    M_F_W=pow((2e0*H_W*C_S/(EVA_Q['P']['K_F']*EVA_Q['P']['th']*CP_M)),0.5e0);#for calculating the wet fin efficiency, including both the heat and mass transfer
    ETA_F_W=tanh(M_F_W*EVA_Q['P']['L_F'])/(M_F_W*EVA_Q['P']['L_F']);#wet fin efficiency
    
    A_F=EVA_Q['P']['Af'];#fin surface area
    A_T=EVA_Q['P']['Af']+EVA_Q['P']['Apo'];#the whole external area
    ETA_O_W=1e0-A_F/A_T*(1e0-ETA_F_W);#overall external fin efficiency
    R_SURF_W=1/(ETA_O_W*H_W*A_T);#airside thermal resistance under wet condition
    NTU_O_W=1/(R_SURF_W*M_DOT_A*CP_M);#for calculating the effective saturated air enthalpy at T_S_O

    NTU_O=1/(EVA_Q['P']['Ro']*M_DOT_A*CP_M);#for calculating the outlet air temperature, which is involved with sensible heat transfer

    UA_W=1e0/(C_S*Ri+C_S*R_W+CP_M*R_SURF_W);#UA under wet condition, considering both the heat and mass transfer
    NTU_W=1*UA_W/(M_DOT_A);#NTU under wet condition, for calculating heat transfer
    EPSILON_W=1e0-exp(-NTU_W);#epsilon under wet condition

    Q=EPSILON_W*M_DOT_A*(H_A_I-H_SAT);#heat transfer under wet condition
    EVA_Q['q']=Q/EVA_Q['mr'];

    #get outlet humidity and effective surface temperature
    H_A_O=H_A_I-Q/(M_DOT_A);#outlet air enthalpy, actually the energy balance shoudl include the water flowing away, like -(83.84e3+(T_S-20)*4.183e3)*(W_I-W_O), since this part is small and need iteration, so it is ignored
    H_S_S_O=H_A_I-(H_A_I-H_A_O)/(1-exp(-NTU_O_W));#saturated air enthalpy at T_S_S_O, effective surface temperature
    
    
#        H_zero=wair.h(2,0.99);#fixed the minimum possible H_S_S_O
#        if(H_S_S_O<=H_zero):
#             H_S_S_O=H_zero;
    #TP_S_O={'T':0.0,'P':0.0};
    TP_S_O=HPtoTP(H_S_S_O,0.999);
    T_S_O=TP_S_O['T'];#effective surface temperature
    T_O=T_S_O+(EVA_Q['TPi']['T']-T_S_O)*exp(-NTU_O);#outlet air temperature
    EVA_Q['T_S_O']= T_S_O;

    #TP_dew={'T':0.0,'P':0.0};
    TP_dew=WPtoTP(W_I,0.999);#get the dew point corresponding to the inlet humidity ratio
    T_dew=TP_dew['T'];#dew temperature corresponding to the inlet air enthalpy
    
    #===========================================================================
    # dry condition
    #===========================================================================
    if(T_S_O>T_dew):
        EVA_Q['q'] = (EVA_Q['TPi']['T']-TXPi['T'])/EVA_Q['mr']*CmineCrossFlow_wet(R,EVA_Q['mr'],EVA_Q['ma'],EVA_Q['TXPo'],EVA_Q['TPi']['T'],CP_M, Ref);
        Q=EVA_Q['q']*EVA_Q['mr'];#heat transfer amount under dry condition
        EVA_Q['W']=W_I;#without dehumidifying
        
        return 0


    TP_O={'T':0.0,'P':0.0};
    TP_O['T']=T_O;
    HH_min = HAPropsSI('H','P',101325,'T',T_O,'R',0.05) #[J/kg humid air]
    HH_max = HAPropsSI('H','P',101325,'T',T_O,'R',0.999) #[J/kg humid air]
    
    if(H_A_O>=HH_max):
        TP_O['T']=T_O;
        TP_O['P']=0.999;
    elif(H_A_O<=HH_min):
        TP_O['T']=T_O;
        TP_O['P']=0.05;
    else:
        try:
            TP_O=THtoTP(T_O,H_A_O);#outlet air state
        except:
            print('Get_Q_Single_For:: state above saturation')
            TP_O['T']=T_O;
            TP_O['P']=0.999;

    try:
        W_O=HAPropsSI('W','P',101325,'T',T_O,'R',TP_O['P']) #wair.HumidityRatio(T_O,TP_O.P);#outlet air humidity ratio #[-]
    except:
        print('Get_Q_Single_For:: W_O = W_I')
        W_O=W_I;

    EVA_Q['W']=W_O;

    return 0




def Build_Cond(P,Ref):
    '''
    This function for correlating the moving boundary and lumped heat exchanger models
    '''

    TXP_prop={'T':0.0,'X':0.0,'P':0.0};
    #===========================================================================
    # average density, pressure drop, cross-sectional area of each phase
    #===========================================================================
    if(P['count1']>0):
        P['HP_TP1']['H']=P['HP_TP1']['H']/P['count1'];#average enthalpy at the first two-phase segment of different circuits
        P['HP_TP1']['P']=P['HP_TP1']['P']/P['count1'];#average pressure at the first two-phase segment of different circuits
    else: #otherwise use the inlet state of the condenser
        P['HP_TP1'] = P['HP_in'];

    if(P['count2']>0):
        P['HP_TP2']['H']=P['HP_TP2']['H']/P['count2'];#average enthalpy at the last two-phase segment of different circuits
        P['HP_TP2']['P']=P['HP_TP2']['P']/P['count2'];#average pressure at the last two-phase segment of different circuits
    else:#otherwise use the outlet state of the condenser
        P['HP_TP2'] = P['HP_out'];
        
    
    P['L_TOT'] = P['VapL']+P['TPL']+P['LiqL'];#total tube length of the condenser
    P['V_TOT'] = P['V_Liq']+P['V_Vap']+P['V_TP'];#total inner volume of the condenser
    P['A_TOT']=P['V_TOT']/P['L_TOT'];#average heat transfer area per tube length

    #===========================================================================
    # prepare heat transfer calculation
    #===========================================================================
    P['Ga_meanL'] = P['ma_TOT']/P['L_TOT'];#average air flow rate across per tube length

    ma_v = P['Ga_meanL'] * P['VapL'];#air flow rate across the gas-phase region
    ma_tp = P['Ga_meanL'] * P['TPL'];#air flow rate across the two-phase region
    ma_l = P['Ga_meanL'] * P['LiqL'];#air flow rate across the liquid-phase region
    
    cp_a= HAPropsSI('cp','T',P['Tai']['T'],'P',101325,'R',0) #air.Cp(P->Tai); [J/kg/K]
    
    #TXP1={'T':0.0,'X':0.0,'P':0.0};
    #inlet condenser state
    TXP1= HPtoTXP(P['HP_in'],Ref);
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1;
    Tsat1 = PropertyTXPth('T',TXP_prop,Ref);#saturation refrigerant temperature at the inlet condenser state [K]

    TXP_prop['P']=TXP1['P'];
    TXP_prop['T']=Tsat1;
    TXP_prop['X']=1;
    cp_rv = PropertyTXPtr('C',TXP_prop,Ref); #[J/kg/K]

    #TXP2={'T':0.0,'X':0.0,'P':0.0};
    #state at the first two-phase segment of the condenser
    TXP2= HPtoTXP(P['HP_TP1'],Ref);
    TXP_prop['P']=TXP2['P'];
    TXP_prop['X']=1;
    Tsat2 = PropertyTXPth('T',TXP_prop,Ref); #[K]

    TXP_prop['P']=P['HP_TP1']['P'];
    TXP_prop['T']=Tsat2;
    TXP_prop['X']=1;
    Hsat1= PropertyTXPth('H',TXP_prop,Ref);#[J/kg]
    
    P['H1_residual'] = Hsat1 - P['HP_TP1']['H'];#for maintaining the consistency in the moving boundary model analysis

    #TXP3={'T':0.0,'X':0.0,'P':0.0};
    TXP3= HPtoTXP(P['HP_TP2'],Ref);#state at the last two-phase segment of the condenser
    TXP_prop['P']=TXP3['P'];
    TXP_prop['X']=1;
    Tsat3 = PropertyTXPth('T',TXP_prop,Ref);#[K]

    TXP_prop['P']=P['HP_TP2']['P'];
    TXP_prop['T']=Tsat3;
    TXP_prop['X']=0;
    Hsat2= PropertyTXPth('H',TXP_prop,Ref);#[J/kg]

    if(P['LiqL']>0):
        P['H2_residual'] = Hsat2 - P['HP_TP2']['H'];#for maintaining the consistency in the moving boundary model analysis
    #Note: without "if(P->LiqL>0)" it may cause divengence at some cases
    
    TXP_prop['P']=TXP3['P'];
    TXP_prop['T']=Tsat3;
    TXP_prop['X']=0;
    cp_rl = PropertyTXPtr('C',TXP_prop,Ref);#[J/kg/K]
    
    #TXP4={'T':0.0,'X':0.0,'P':0.0};
    #condenser outlet state
    TXP4= HPtoTXP(P['HP_out'],Ref);
    

    #===========================================================================
    # heat transfer of the vapor phase
    #===========================================================================
    cmin_v=0.0;
    Cr_v = 0.0;

    if(ma_v>0):
        Q_v = P['mr']*(P['HP_in']['H']-P['HP_TP1']['H']);#vapor-phase heat transfer
        Cr_v = (P['mr']*cp_rv)/(ma_v*cp_a);
    
        if(ma_v*cp_a>P['mr']*cp_rv):
            cmin_v = P['mr']*cp_rv;
            Cr_v = Cr_v;
        else:
            cmin_v = ma_v*cp_a;
            Cr_v = 1/Cr_v;
        Qv_max = cmin_v*(TXP1['T'] - P['Tai']['T']);
        e_v = Q_v/Qv_max;#actual heat transfer effectiveness of the vapor-phase
        NTU_v = P['UA_Vap']/cmin_v;
        
        e_sup =0.0;#theoritical heat transfer effectivenss of the vapor-phase
        if(ma_v*cp_a>P['mr']*cp_rv):#cmax unmixed
            l = 1-exp(-Cr_v*NTU_v);
            e_sup = 1 - exp(-l/Cr_v);
        else:#cmin unmixed
            l = -Cr_v*(1-exp(-NTU_v));
            e_sup = (1-exp(l))/Cr_v;
        
        P['r_v'] = e_v/e_sup;#for ajusting the theortical heat transfer effectivenss of vapor-phase
        P['U_Vap'] = P['UA_Vap']/(P['VapL']);#average vapor-phase heat transfer conduntance per tube length
    
        P['rho_Vap'] = P['m_Vap']/P['V_Vap'];#average vapor-phase density
        P['DP_Vap'] = (P['HP_in']['P']-P['HP_TP1']['P'])/P['V_Vap'];#average pressure drop gradient
        P['A_Vap']=P['V_Vap']/P['VapL'];#average heat transfer surface area per tube length

    
    #===========================================================================
    # heat transfer of liquid phase
    #===========================================================================
    cmin_l=0.0;
    Cr_l = 0.0;
        
    if(ma_l>0):
        Q_l = P['mr']*(P['HP_TP2']['H']-P['HP_out']['H']);#liquid-phase heat transfer amount
        Cr_l = (P['mr']*cp_rl)/(ma_l*cp_a);
    
        if(ma_l*cp_a>P['mr']*cp_rl):
            cmin_l = P['mr']*cp_rl;
            Cr_l=Cr_l;
        else:
            cmin_l = ma_l*cp_a;
            Cr_l =1/Cr_l;
        Ql_max = cmin_l*(TXP3['T'] - P['Tai']['T']);
        e_l = Q_l/Ql_max;#actual liquid-phase heat tranfer effectivenss
        NTU_l=P['UA_Liq']/cmin_l;
    
        e_sub =0.0;#theoritical liquid-phase heat transfer effectivenss
    
        if(ma_l*cp_a>P['mr']*cp_rl):#cmax unmixed
            l = 1-exp(-Cr_l*NTU_l);
            e_sub = 1 - exp(-l/Cr_l);
        else:
            l = -Cr_l*(1-exp(-NTU_l));
            e_sub = (1-exp(l))/Cr_l;
    
        P['r_l'] = e_l/e_sub;#parameter for adjusting the theoritical liquid-phase heat transfer effectivess
        P['U_Liq'] = P['UA_Liq']/(P['LiqL']);#averge liquid-phase heat transfer conductance
    
        P['DP_Liq'] = (P['HP_TP2']['P']-P['HP_out']['P'])/P['V_Liq'];#avergage liquid-phase pressure drop gradient
        P['A_Liq']=P['V_Liq']/P['LiqL'];#avergage liquid-phase heat transfer surface area
        P['rho_Liq'] = P['m_Liq']/P['V_Liq'];#average liquid-phase density

    else:#otherwise use the default parameters
        FirCond=0;
        if(FirCond==0):
            P['r_l'] = 1.0;
            P['U_Liq'] = 29;
            P['DP_Liq'] = 102247;
            P['A_Liq']=P['A_TOT'];
            P['rho_Liq'] = 1100;
            FirCond=1;

    #===========================================================================
    # two-phase heat transfer
    #===========================================================================
    Q_tp = P['mr']*(P['HP_TP1']['H']-P['HP_TP2']['H']);#two-phase heat transfer amount
    cmin_tp = ma_tp*cp_a;    
    Qtp_max = cmin_tp*(TXP2['T'] - P['Tai']['T']);
    e_TP = Q_tp/Qtp_max;#actual two-phase heat transfer effectivenss
    NTU_tp = P['UA_TP']/cmin_tp;
    e_sat = 1-exp(-NTU_tp);#theoritical two-phase heat transfer effectivenss

    P['r_tp'] = e_TP/e_sat;#parameter for adjusting the theoritical two-phase heat transfer effectivenss    
    P['U_TP'] = P['UA_TP']/(P['TPL']);#average two-phase heat transfer conductance per tube length

    P['rho_TP'] = P['m_TP']/P['V_TP'];#average two-phase density
    P['A_TP']=P['V_TP']/P['TPL'];#average two-phase heat transfer surface area per tube length
    P['DP_TP'] = (P['HP_TP1']['P']-P['HP_TP2']['P'])/P['V_TP'];#average pressure drop gradient

    #===========================================================================
    # lumped model
    #===========================================================================
    Q_TOT = P['mr']*(P['HP_in']['H']-P['HP_out']['H']);#overall heat transfer amount
    epsilon_TOT = Q_TOT/(P['ma_TOT']*cp_a)/(Tsat1-P['Tai']);#overall effectivenss
    P['DP_TOT'] = (P['HP_in']['P']-P['HP_out']['P']);#overall pressure drop
    NTU_TOT = -log(1-epsilon_TOT);
    P['UA_TOT'] = NTU_TOT*P['ma_TOT']*cp_a;#overall heat transfer conductance

    return 0



class StructEvapClass():
    
    def __init__(self,filename,Ref):#mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,P,Ref):
        '''Initialize'''
        #reading the data file in
        #the followings from   Class StructEvap(filename) in the ACMODLE
        
        self.Rev = 0 # for calculation evaporator solver direction, 0-- for forward evaporator solver, 1--for reversed evaporator solver
        self.Ref = Ref
        #self.mr = mr
        #self.HPo = HPo
        #self.Ga = Ga
        #self.TPi = TPi
        #self.HPi = HPi
        #self.TPo = TPo
        #self.Sm = Sm
        #self.Aflow = Aflow
        #self.P = P      
           
        #=======================================================================
        # filename = "../InputDoc/EvapStruc.xlsx"
        #=======================================================================
        
        df = pd.read_excel(filename,sheetname='description',header = 0)
        df_node = pd.read_excel(filename,sheetname='node',header = 0)
        df_branch = pd.read_excel(filename,sheetname='branch',header = 0)
        df_tube = pd.read_excel(filename,sheetname='tube',header = 0)
           
        self.NodNum = int(df.EvapStruc[0]);          #Node Number
        self.BranNum = int(df.EvapStruc[1]);         #Branch Number
        self.RowNum = int(df.EvapStruc[2]);          #Row Number
        self.TubeNum = int(df.EvapStruc[3]);         #overall tube Number
        self.SegNum = int(df.EvapStruc[4]);          #Segment Number of per tube
        self.AreaFront = df.EvapStruc[5];       #frontal area
        self.Volum = df.EvapStruc[6];           #total nominal air mass flow rate [m^3/s]
       
        self.Nod = [EvapNode() for k in range(self.NodNum)]
        self.Bra = [EvapBranch()for k in range(self.BranNum)]
        self.GaRow = [0.0 for k in range(self.RowNum)]
        self.Tub = [TubeEvap() for k in range(self.TubeNum)]
       
        for i in range(self.NodNum):   
            self.Nod[i]['NodNo'] = df_node.NodNo[i]
            self.Nod[i]['OutNum'] = df_node.OutNum[i]    #get outlet branches at first
            self.Nod[i]['InNum'] = df_node.InNum[i]      #get inlet branches second
            self.Nod[i]['BranOUT'] = list(df_node.ix[i][3 : 3+self.Nod[i]['OutNum']]) #get outlet branches first
            self.Nod[i]['BranIN'] = list(df_node.ix[i][3+self.Nod[i]['OutNum'] : 3+self.Nod[i]['OutNum']+self.Nod[i]['InNum']]) #get inlet branches second
       
        for i in range(self.BranNum):
            self.Bra[i]['BranNo'] = df_branch.BranNo[i]
            self.Bra[i]['EqulNo'] = df_branch.EqulNo[i]
            self.Bra[i]['GrFac'] = df_branch.GrFac[i]
            self.Bra[i]['TubNum'] = df_branch.TubNum[i]
            self.Bra[i]['TubNo'] = list(df_branch.ix[i][4 : 4+self.Bra[i]['TubNum']])#important, the tube number in branch, inputted first from the compressor suction
       
        for i in range(self.TubeNum):
            self.Tub[i]['Seg'] = [TubEvpSeg() for k in range(self.SegNum)]
            self.Tub[i]['TubNo'] = df_tube.TubNo[i]
            self.Tub[i]['RowNo'] = df_tube.RowNo[i]
            self.Tub[i]['Refdownstream'] = df_tube.Refdownstream[i]
            self.Tub[i]['AirUpstreamUpper'] = df_tube.AirUpstreamUpper[i]
            self.Tub[i]['AirUpstreamLower'] = df_tube.AirUpstreamLower[i]
            self.Tub[i]['GaFac'] = df_tube.GaFac[i]
            self.Tub[i]['even'] = df_tube.even[i]
        
    def Update(self):
        """Update the parameters passed in using the dictionary"""
        pass
   
    def DeStructEvap(self):
        '''
        ######This function is NOT in use######
        delete the memory
        '''
        for i in range(self.TubeNum):
            del self.Tub[i]['Seg']
        del self.Tub
         
        for i in range(self.BranNum):
            del self.Bra[i]['TubNo']
        del self.Bra
     
        for i in range(self.NodNum):
            del self.Nod[i]['BranIN']
            del self.Nod[i]['BranOUT']
        del self.Nod
          
    def _EvapCircuit_Fwd(self,mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,P):
        '''
        Forward Evaporator solver
        '''
        
        Gr=0.0;#mass flux
        H_in=0;
        
        if (self.Rev): #reversed calculation
            return self._EvapCircuit_Rev(mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,P)
            
        
        TPo=TPi.copy();#inlet air state
    
        H_in=HPi['H'];#inlet refrigerant enthalpy
        
        #air side initializing
        rho_air=1/HAPropsSI("V", "T", TPi['T'], "P", 101325, "R", TPi['P']) #[kg dry air/ m^3]
    
        #const double Ma =Ga*rho_air*4.719e-4; 
        Ga=Ga/self.AreaFront*P['vsp']*P['Ls']/((P['vsp']-P['Do'])*(P['Ls']-P['N']*P['th']));
    
    
        for i in range(self.RowNum):
            self.GaRow[i]=0 #this step can be deleted because GaRow[i] already initialized with zero
            
        for j in range(self.RowNum):
            N=0;
            for i in range(self.TubeNum):
                if (self.Tub[i]['RowNo']==j):
                    for k in range(self.SegNum):
                        self.Tub[i]['Seg'][k]['TPi']['T']=TPi['T'];
                        self.Tub[i]['Seg'][k]['TPi']['P']=TPi['P'];
                        self.Tub[i]['Seg'][k]['WHo']['W'] = HAPropsSI("W", "T", TPi['T'], "P", 101325, "R", TPi['P']) #[kg water/kg dry air] ####wair.HumidityRatio(TPi.T,TPi.P)
                        self.Tub[i]['Seg'][k]['WHo']['H'] = HAPropsSI("H", "T", TPi['T'], "P", 101325, "R", TPi['P']) #[J/kg humid air] ####wair.h(TPi.T,TPi.P);
        
                    if (self.Tub[i]['RowNo']>0):#not the first row
                        Upper = self.Tub[i]['AirUpstreamUpper']
                        Lower = self.Tub[i]['AirUpstreamLower']
                        if (Upper>=0 and Lower>=0):
                            self.Tub[i]['GaFac']=(self.Tub[Upper]['GaFac']+self.Tub[Lower]['GaFac'])/2
                        elif (Upper>=0):
                            self.Tub[i]['GaFac']=self.Tub[Upper]['GaFac']
                        else:
                            self.Tub[i]['GaFac']=self.Tub[Lower]['GaFac']
                    
                    self.GaRow[j]=self.GaRow[j]+self.Tub[i]['GaFac']
                    N=N+1
                #ifend    
            #i loop end
            self.GaRow[j]=self.GaRow[j]/N
        #j loop end
        
        for i in range(self.TubeNum):
            RowN=self.Tub[i]['RowNo'];
            self.Tub[i]['Ga']=self.Tub[i]['GaFac']/self.GaRow[RowN]*Ga
    
        H_out=0.0;  #intermidiate evaporator exit enthalpy
        Res=0;
        DP = 20000.0; #[Pa]
        P_in=0.0;
        P_out=0.0;
        HPi['P']=HPo['P']+DP;   #temperary inlet refrigerant pressure
        HPi['H']=H_in;          #evaporator inlet refrigerant enthalpy
        #HPi is the intermediate variable
        P_in = HPo['P']+DP;       #temperary inlet refrigerant pressure
        P_out = HPo['P'];         #evaporator outlet refrigerant pressure
        IterN=0;
        
    
        #refrigerant state intialize first time
        for i in range(self.NodNum):#inlet nodes
            if (self.Nod[i]['BranIN'][0]<0):#no inlet branch, only from the distributor
                Gr = mr/(P['Ax']*self.Nod[i]['OutNum']);    

                for j in range(self.Nod[i]['OutNum']):#states flowing out from the node
                    jj = self.Nod[i]['BranOUT'][j];#index of the outlet branches
                    self.Bra[jj]['HPi']['H']=HPi['H'];
                    self.Bra[jj]['HPi']['P']=HPi['P'];
                    self.Bra[jj]['HPo']['H']=HPo['H'];
                    self.Bra[jj]['HPo']['P']=HPo['P'];
                    self.Bra[jj]['Gr']=Gr*self.Bra[jj]['GrFac'];
                    self.Bra[jj]['Ini']=1;
                #end j loop
            #endif
        #end i loop
    
        while True:
            
            P['VapL'] =0;
            P['TPL'] = 0;
            P['LiqL'] = 0;
            P['V_Vap'] =0;
            P['V_TP'] =0;
            P['V_Liq'] = 0;
            P['m_Vap'] = 0;
            P['m_TP'] = 0;
            P['m_Liq'] = 0;    
            P['UA_Vap'] = 0;
            P['UA_TP'] = 0;
            P['UA_Liq'] = 0;
    
            IterN=IterN+1;
            
            HPi['H']=H_in; #evaporator inlet refrigerant enthalpy
            HPi['P']=P_in;
            
            for i in range(self.BranNum):
                self.Bra[i]['Ini']=0
        
            for i in range(self.NodNum): #inlet nodes
                if (self.Nod[i]['BranIN'][0]<0):#from the distributor
                    Gr = mr/(P['Ax']*self.Nod[i]['OutNum']);    
                    for j in range(self.Nod[i]['OutNum']): #states flowing out from the node
                        jj = self.Nod[i]['BranOUT'][j];     #index of the outlet branches
                        self.Bra[jj]['HPi']['H']=HPi['H'];
            
                        if (self.Nod[i]['OutNum']==self.BranNum):#no middle nodes
                            DDP = (self.Bra[jj]['HPi']['P']-self.Bra[jj]['HPo']['P']);
                            self.Bra[jj]['HPi']['P']=HPo['P']+DDP;
                            P['DPr'][jj]=DDP;#output the pressure drop of each branch
                        else:
                            self.Bra[jj]['HPi']['P']=HPi['P'];
                        self.Bra[jj]['Gr']=Gr*self.Bra[jj]['GrFac'];
                        self.Bra[jj]['Ini']=1;
                    #end j circle
                #endif
            #end i loop
            
            self.Cal_HP(P,0,HPi);    #Heat transfer and pressure drop calculation #A.B. return updated values for P and HPi
            
            for i in range(self.BranNum):
                self.Bra[i]['Ini']=0;
        
            Gr=0;
            HPi['H']=0;
            HPi['P']=0;
            
            #nodes in the middle
            for i in range(self.NodNum):
                if (self.Nod[i]['BranIN'][0]>=0 and self.Nod[i]['BranOUT'][0]>=0): #nodes in the middle
                    for j in range(self.Nod[i]['InNum']):   #node inlet state
                        jj= self.Nod[i]['BranIN'][j];
                        Gr=Gr+self.Bra[jj]['Gr'];
                        HPi['H']=self.Bra[jj]['HPo']['H']*self.Bra[jj]['Gr']+HPi['H'];
                        HPi['P']=self.Bra[jj]['HPo']['P']*self.Bra[jj]['Gr']+HPi['P'];
        
                    HPi['H']=HPi['H']/Gr;
                    HPi['P']=HPi['P']/Gr;
                    Gr=Gr/self.Nod[i]['InNum'];
                    Gr=Gr*self.Nod[i]['InNum']/self.Nod[i]['OutNum'];
        
                    for j in range(self.Nod[i]['OutNum']):
                        jj = self.Nod[i]['BranOUT'][j];   #index of outlet branches
                        self.Bra[jj]['HPi']['H']=HPi['H'];
                        self.Bra[jj]['HPi']['P']=HPi['P'];
                        self.Bra[jj]['Gr']=Gr*self.Bra[jj]['GrFac'];
                        self.Bra[jj]['Ini']=1;
                    #end j circle
                #endif
            #end i circle
        
        
            self.Cal_HP(P,1,HPi);   #Heat transfer and pressure drop calculation #A.B. return updated values for P and HPi
            
            Gr=0;
            HPi['H']=0;
            HPi['P']=0;
            
            #exit nodes
            for i in range(self.NodNum):
                if (self.Nod[i]['BranOUT'][0]<0):#no outlet branch except the suction line
                    for j in range(self.Nod[i]['InNum']):
                        jj=self.Nod[i]['BranIN'][j];
                        Gr=Gr+self.Bra[jj]['Gr'];
                        HPi['H']=self.Bra[jj]['HPo']['H']*self.Bra[jj]['Gr']+HPi['H'];
                        HPi['P']=self.Bra[jj]['HPo']['P']*self.Bra[jj]['Gr']+HPi['P'];
                        P['Hout8'][jj] = self.Bra[jj]['HPo']['H'];  #output the exit enthalpy of each branch
                    
                    HPi['H']=HPi['H']/Gr;
                    HPi['P']=HPi['P']/Gr;
                    Gr=Gr/self.Nod[i]['InNum'];
                    Gr=Gr*self.Nod[i]['InNum']/self.Nod[i]['OutNum'];
                #endif
            #end i circle
            
            if (self.RowNum==1):
                Res=0;
            else:
                Res=2*(HPi['H']-H_out)/(HPi['H']+H_out);
            H_out=HPi['H'];
            P_out=HPi['P'];
            P_in=HPo['P']+(P_in-P_out);
            
            #print the res and iteration no.
            print('Res {}, IterN {}'.format(Res, IterN))
            
            if (abs(Res)<1e-6 or IterN>20): #condition to break the while loop
                break

        #end while loop
            
        if (abs(Res)>1e-5):
            print()
            print('######################################################')
            print('_EvapCircuit_Fwd, Can NOT reach the required tolerance, Res= '+str(Res))
            print('######################################################')
            print()
    
        Sm['m']=0;
        Sm['V']=0;
    
        for i in range(self.BranNum):
            Sm['m']=Sm['m']+self.Bra[i]['m']['m'];
            Sm['V']=Sm['V']+self.Bra[i]['m']['V'];
    
        WH_out={'W':0,'H':0};
        Ma_out=0.0;
    
        #air outputs
        for i in range(self.BranNum):
            if (self.Bra[i]['EqulNo']<0): #no equivalent circuits
                for j in range(self.Bra[i]['TubNum']):
                    Tubj=0;
                    if (self.Rev):   #counted from the point connected to compressor suction
                        Tubj=j;
                    else:
                        Tubj=self.Bra[i]['TubNum']-1-j;#counted from the point of evaporator entrance
                    TubeN=self.Bra[i]['TubNo'][Tubj];
    
                    if (self.Tub[TubeN]['RowNo']==self.RowNum-1):#outlet row tube at air side
                        for k in range(self.SegNum):
                            WH_out['W']=self.Tub[TubeN]['Ga']*self.Tub[TubeN]['Seg'][k]['WHo']['W']+WH_out['W'];
                            WH_out['H']=self.Tub[TubeN]['Ga']*self.Tub[TubeN]['Seg'][k]['WHo']['H']+WH_out['H'];
                            Ma_out=Ma_out+self.Tub[TubeN]['Ga'];
    
        WH_out['W'] = WH_out['W']/Ma_out;
        WH_out['H'] = WH_out['H']/Ma_out;
        
        try:
            #TPo_temp = {'T':0.0,'P':0.0}
            #TPo_temp['T']=TPi['T']-5;
            #TPo_temp['P']=0.8;
            #TPo = WHtoTP(WH_out,TPo_temp)
            TPo['T'] = HAPropsSI('T','P',101325,'H',WH_out['H'],'W',WH_out['W'])
            TPo['P'] = HAPropsSI('R','P',101325,'H',WH_out['H'],'W',WH_out['W'])
        except:
            print('_EvapCircuit_Fwd :: check TPo that need to be numerically solved!')
            print('Other source of exception might be due to relative humidity ~100%, try to solve TPo with relative humidity of 99.5%')
            TPo['T'] = HAPropsSI('T','P',101325,'H',WH_out['H'],'R',0.995)
            TPo['P'] = 0.995

    
        #refrigerant outputs
        HPo['H']=H_out;
        HPo['P']=HPo['P'];
        HPi['H']=H_in;
        HPi['P']=P_in;
    
        #Sm->m=Sm->m + (2.87-2.766)+(3.10-3.05602)/(2.0792-7.955)*(P->VapL-7.955);#for three points charge tuning
        
        return 0
    
    def _EvapCircuit_Rev(self,mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,P):
        '''
        reversed evaporator solver
        '''
        
        Gr=0.0;#mass flux
        
        self.Rev=1;
        #air side initializing
        rho_air=1/HAPropsSI("Vda", "T", TPi['T'], "P", 101325, "R", 0) #[kg dry air/ m^3]#
        Ma =self.Volum*rho_air; 
        GaNom=Ma/self.AreaFront*P['vsp']*P['Ls']/((P['vsp']-P['Do'])*(P['Ls']-P['N']*P['th']));
        if(Ga<0):
            Ga = GaNom;
    
        for i in range(self.RowNum):
            self.GaRow[i]=0
        
        
        for j in range(self.RowNum):
            N=0;
            for i in range(self.TubeNum):
                if (self.Tub[i]['RowNo']==j):
                    for k in range(self.SegNum):
                        self.Tub[i]['Seg'][k]['TPi']['T']=TPi['T'];
                        self.Tub[i]['Seg'][k]['TPi']['P']=TPi['P'];
                        self.Tub[i]['Seg'][k]['WHo']['W'] = HAPropsSI("W", "T", TPi['T'], "P", 101325, "R", TPi['P']) #[kg water/kg dry air]
                        self.Tub[i]['Seg'][k]['WHo']['H'] = HAPropsSI("H", "T", TPi['T'], "P", 101325, "R", TPi['P']) #[J/kg humid air]
                        
                    if(self.Tub[i]['RowNo']>0): #not the first row
                        Upper = self.Tub[i]['AirUpstreamUpper'];
                        Lower = self.Tub[i]['AirUpstreamLower'];
                        if (Upper>=0 and Lower>=0):
                            self.Tub[i]['GaFac']=(self.Tub[Upper]['GaFac']+self.Tub[Lower]['GaFac'])/2
                        elif(Upper>=0):
                            self.Tub[i]['GaFac']=self.Tub[Upper]['GaFac']
                        else:
                            self.Tub[i]['GaFac']=self.Tub[Lower]['GaFac']
                    
                    self.GaRow[j]=self.GaRow[j]+self.Tub[i]['GaFac'];
                    N=N+1;
                #ifend    
            #i circle end
            self.GaRow[j]=self.GaRow[j]/N;
        #j circle end
        
        for i in range(self.TubeNum):
            RowN=self.Tub[i]['RowNo'];
            self.Tub[i]['Ga']=self.Tub[i]['GaFac']/self.GaRow[RowN]*Ga;
        
    
        # Refrigerant side initialize
        H_in=0;
        Res=0;
        P_in=0;
        IterN=0;
    
        while True:
            P['VapL'] =0;
            P['TPL'] = 0;
            P['LiqL'] = 0;
            P['V_Vap'] =0;
            P['V_TP'] =0;
            P['V_Liq'] = 0;
            P['m_Vap'] = 0;
            P['m_TP'] = 0;
            P['m_Liq'] = 0;    
            P['UA_Vap'] = 0;
            P['UA_TP'] = 0;
            P['UA_Liq'] = 0;
    
            IterN=IterN+1;
            
            HPi['H']=HPo['H'];
            HPi['P']=HPo['P'];
            
            for i in range(self.BranNum):
                self.Bra[i]['Ini']=0;
    
            for i in range(self.NodNum):#Outlet nodes
                if(self.Nod[i]['BranOUT'][0]<0):#to compressor suction
                    Gr = mr/(P['Ax']*self.Nod[i]['InNum']);
                    for j in range(self.Nod[i]['InNum']):#states flowing into the node
                        jj = self.Nod[i]['BranIN'][j];#index of the inlet branches
                        self.Bra[jj]['HPo']['H']=HPi['H'];
                        self.Bra[jj]['HPo']['P']=HPi['P'];
                        self.Bra[jj]['Gr']=Gr*self.Bra[jj]['GrFac'];
                        self.Bra[jj]['Ini']=1;
                    #end j circle
                #endif
            #end i circle
            
            self.Cal_HP(P,0,HPi);   #Heat transfer and pressure drop calculation #A.B. return updated values for P and HPi
            
            for i in range(self.BranNum):
                self.Bra[i]['Ini']=0;
            
            Gr=0;
            HPi['H']=0;
            HPi['P']=0;
        
            #nodes in the middle
            for i in range(self.NodNum):
                if(self.Nod[i]['BranOUT'][0]>=0 and self.Nod[i]['BranIN'][0]>=0):
                    for j in range(self.Nod[i]['OutNum']):#node outlet state
                        jj= self.Nod[i]['BranOUT'][j];
                        Gr=Gr+self.Bra[jj]['Gr'];
                        HPi['H']=self.Bra[jj]['HPi']['H']*self.Bra[jj]['Gr']+HPi['H'];
                        HPi['P']=self.Bra[jj]['HPi']['P']*self.Bra[jj]['Gr']+HPi['P'];
                        
                    HPi['H']=HPi['H']/Gr;
                    HPi['P']=HPi['P']/Gr;
                    Gr=Gr/self.Nod[i]['OutNum'];
                    Gr=Gr*self.Nod[i]['OutNum']/self.Nod[i]['InNum'];
                
                    for j in range(self.Nod[i]['InNum']):
                        jj = self.Nod[i]['BranIN'][j];#index of outlet branches
                        self.Bra[jj]['HPo']['H']=HPi['H'];
                        self.Bra[jj]['HPo']['P']=HPi['P'];
                        self.Bra[jj]['Gr']=Gr*self.Bra[jj]['GrFac'];
                        self.Bra[jj]['Ini']=1;
                    #end j circle
                #endif
            #end i circle
            
            self.Cal_HP(P,1,HPi); #Heat transfer and pressure drop calculation #A.B. return updated values for P and HPi

            Gr=0;
            HPi['H']=0;
            HPi['P']=0;

            #inlet nodes
            for i in range(self.NodNum):
                if(self.Nod[i]['BranIN'][0]<0):#no inlet branch except the distribution line
                    for j in range(self.Nod[i]['OutNum']):
                        jj= self.Nod[i]['BranOUT'][j];
                        Gr=Gr+self.Bra[jj]['Gr'];
                        HPi['H']=self.Bra[jj]['HPi']['H']*self.Bra[jj]['Gr']+HPi['H'];
                        HPi['P']=self.Bra[jj]['HPi']['P']*self.Bra[jj]['Gr']+HPi['P'];
                        
                    HPi['H']=HPi['H']/Gr;
                    HPi['P']=HPi['P']/Gr;
                    Gr=Gr/self.Nod[i]['OutNum'];
                    Gr=Gr*self.Nod[i]['OutNum']/self.Nod[i]['InNum'];
                #endif
            #end i circle
                
            if (self.RowNum==1):
                Res=0;
            else:
                Res=2*(HPi['H']-H_in)/(HPi['H']+H_in)
            
            H_in=HPi['H'];
            
            #print Res and iteration no.
            print('Res {}, IterN {}'.format(Res, IterN))
            
            if (abs(Res)<1e-6 or IterN>20): #condition to break the while loop
                break

        #end while loop
    
        if (abs(Res)>1e-5):
            print()
            print('######################################################')
            print('_EvapCircuit_Rev, Can NOT reach the required tolerance, Res= '+str(Res))
            print('######################################################')
            print()
            
            
        Sm['m']=0; 
        Sm['V']=0;
        
        for i in range(self.BranNum):
            Sm['m']=Sm['m']+self.Bra[i]['m']['m'];
            Sm['V']=Sm['V']+self.Bra[i]['m']['V'];
    
        WH_out={'W':0,'H':0};
        Ma_out=0;
    
        #air outputs
        for i in range(self.TubeNum):
            if(self.Tub[i]['RowNo']==self.RowNum-1):#airside outlet row
                for j in range(self.SegNum):
                    WH_out['W']=self.Tub[i]['Ga']*self.Tub[i]['Seg'][j]['WHo']['W']+WH_out['W'];
                    WH_out['H']=self.Tub[i]['Ga']*self.Tub[i]['Seg'][j]['WHo']['H']+WH_out['H'];
                    Ma_out=Ma_out+self.Tub[i]['Ga'];
    
        WH_out['W'] = WH_out['W']/Ma_out;
        WH_out['H'] = WH_out['H']/Ma_out;
        
        try:
            #TPo_temp = {'T':0.0,'P':0.0}
            #TPo_temp['T']=TPi['T']-5;
            #TPo_temp['P']=0.8;
            #TPo = WHtoTP(WH_out,TPo_temp)
            TPo['T'] = HAPropsSI('T','P',101325,'H',WH_out['H'],'W',WH_out['W'])
            TPo['P'] = HAPropsSI('R','P',101325,'H',WH_out['H'],'W',WH_out['W'])        
        except:
            print('_EvapCircuit_Rev :: check TPo that need to be numerically solved!')
            print('Other source of exception might be due to relative humidity ~100%, try to solve TPo with relative humidity of 99.5%')
            TPo['T'] = HAPropsSI('T','P',101325,'H',WH_out['H'],'R',0.995)
            TPo['P'] = 0.995
        
        
        return 0
    
    def Cal_HP(self,P,Pos,HPi):
        ''''function to calculate the heat transfer and pressure drop'''
        
        #heat transfer and pressure drop of the two-phase region
        Logic=0;
        WHo = {'W':0,'H':0};
        Sm = {'m':0,'V':0};
        mi = {'m':0,'V':0};
        #Bak = ETdim()
        
        Bak = P.copy() #keeping the information of the evaporator
        
        if(Pos<0 or Pos>1):
            print("Cal_HP, Wrong position")
        
    
        for i in range(self.BranNum):
            #branch parameter clear zero
            P['VapL'] =0;
            P['TPL'] = 0;
            P['LiqL'] = 0;
            P['V_Vap'] =0;
            P['V_TP'] =0;
            P['V_Liq'] = 0;
            P['m_Vap'] = 0;
            P['m_TP'] = 0;
            P['m_Liq'] = 0;    
            P['UA_Vap'] = 0;
            P['UA_TP'] = 0;
            P['UA_Liq'] = 0;
            
            if Pos == 0: 
                if (self.Bra[i]['Ini']==1):#this branch has been initialized
                    Logic=1;
                else:
                    Logic=0;  
            elif Pos == 1:
                if(self.Bra[i]['Ini']==1):#this branch has been initialized
                    Logic=1;
                else:
                    Logic=0;
            
            if(Logic):
                
                if(self.Rev):
                    HPi['H']=self.Bra[i]['HPo']['H']; #opposite to flow direction
                    HPi['P']=self.Bra[i]['HPo']['P'];
                else:
                    HPi['H']=self.Bra[i]['HPi']['H']; #paralell to flow direction
                    HPi['P']=self.Bra[i]['HPi']['P'];
                
                if(self.Bra[i]['EqulNo']<0):    #no equivalent branch

                    for j in range(self.Bra[i]['TubNum']):
                        Tubj=0;
                        if (self.Rev):   #counted from the point connected to compressor suction
                            Tubj=j;
                        else:
                            Tubj=self.Bra[i]['TubNum']-1-j;#counted from the point of evaporator entrance
            
                        TubeN=self.Bra[i]['TubNo'][Tubj];
                        
                        self.Tub[TubeN]['HPo']['H']=HPi['H'];
                        self.Tub[TubeN]['HPo']['P']=HPi['P'];
                        self.Tub[TubeN]['m']['m']=0;
                        self.Tub[TubeN]['m']['V']=0;
            
                        if (self.Tub[TubeN]['RowNo']>0): #not the first row, to get the air state  
                            Upper = self.Tub[TubeN]['AirUpstreamUpper'];
                            Lower = self.Tub[TubeN]['AirUpstreamLower'];
                            for k in range(self.SegNum):
                                if(Upper>=0 and Lower>=0):
                                    WHo['W']=(self.Tub[Upper]['Seg'][k]['WHo']['W']*self.Tub[Upper]['Ga']+self.Tub[Lower]['Seg'][k]['WHo']['W']*self.Tub[Lower]['Ga'])/(self.Tub[Upper]['Ga']+self.Tub[Lower]['Ga']);
                                    WHo['H']=(self.Tub[Upper]['Seg'][k]['WHo']['H']*self.Tub[Upper]['Ga']+self.Tub[Lower]['Seg'][k]['WHo']['H']*self.Tub[Lower]['Ga'])/(self.Tub[Upper]['Ga']+self.Tub[Lower]['Ga']);
                                elif (Upper>=0):
                                    WHo['H']=self.Tub[Upper]['Seg'][k]['WHo']['H'];
                                    WHo['W']=self.Tub[Upper]['Seg'][k]['WHo']['W'];
                                else:
                                    WHo['H']=self.Tub[Lower]['Seg'][k]['WHo']['H'];
                                    WHo['W']=self.Tub[Lower]['Seg'][k]['WHo']['W']
            
                                try:
                                    #self.Tub[TubeN]['Seg'][k]['TPi'] = WHtoTP(WHo,self.Tub[TubeN]['Seg'][k]['TPi']);
                                    self.Tub[TubeN]['Seg'][k]['TPi']['T'] = HAPropsSI('T','P',101325,'H',WHo['H'],'W',WHo['W'])
                                    self.Tub[TubeN]['Seg'][k]['TPi']['P'] = HAPropsSI('R','P',101325,'H',WHo['H'],'W',WHo['W'])
                                except:
                                    print('Cal_HP :: check self.Tub that need to be numerically solved!')
                                    print('Other source of exception might be due to relative humidity ~100%, try to solve TPo with relative humidity of 99.5%')
                                    raise
                                    self.Tub[TubeN]['Seg'][k]['TPi']['T'] = HAPropsSI('T','P',101325,'H',WHo['H'],'R',0.995)
                                    self.Tub[TubeN]['Seg'][k]['TPi']['P'] = 0.995
                                    
                            #end k circle
                        #endif
                        
                        for k in range(self.SegNum):
                            #for debugging
                            #if(k==9 and TubeN==13):    
                            #    shenb=0;
                
                            realk=0;
                            if (self.Tub[TubeN]['even']):
                                realk=(self.SegNum-1-k);
                            else:
                                realk=k;
                            
                            if (self.Rev):
                                EvapTubeL_Rev(self.Bra[i]['Gr'],HPi,self.Tub[TubeN]['Ga'],self.Tub[TubeN]['Seg'][realk]['TPi'],WHo,mi,P, self.Ref) #return updated (HPi, WHo, mi, P) for segment tubeL model
                            else:
                                EvapTubeL_Fwd(self.Bra[i]['Gr'],HPi,self.Tub[TubeN]['Ga'],self.Tub[TubeN]['Seg'][realk]['TPi'],WHo,mi,P, self.Ref) #return updated (HPi, WHo, mi, P) for segment tubeL model
                
                            self.Tub[TubeN]['Seg'][realk]['WHo']['H']=WHo['H'];
                            self.Tub[TubeN]['Seg'][realk]['WHo']['W']=WHo['W'];
                            self.Tub[TubeN]['m']['m']=self.Tub[TubeN]['m']['m']+mi['m'];
                            self.Tub[TubeN]['m']['V']=self.Tub[TubeN]['m']['V']+mi['V'];
                            Sm['m']=Sm['m']+mi['m'];
                            Sm['V']=Sm['V']+mi['V'];
                        #end k circle
        
                        
                        if (self.Rev):
                            EvapTubeBend_Rev(self.Bra[i]['Gr'],HPi,mi,P, self.Ref) #return updated (HPi, mi, P) for bend model
                        else:
                            EvapTubeBend_Fwd(self.Bra[i]['Gr'],HPi,mi,P, self.Ref) #return updated (HPi, mi, P) for bend model
                        
                        self.Tub[TubeN]['m']['m']=self.Tub[TubeN]['m']['m']+mi['m'];
                        self.Tub[TubeN]['m']['V']=self.Tub[TubeN]['m']['V']+mi['V'];
                        Sm['m']=Sm['m']+mi['m'];
                        Sm['V']=Sm['V']+mi['V'];
                        if (self.Rev):
                            self.Tub[TubeN]['HPi']['H']=HPi['H'];
                            self.Tub[TubeN]['HPi']['P']=HPi['P'];
                        else:
                            self.Tub[TubeN]['HPo']['H']=HPi['H'];
                            self.Tub[TubeN]['HPo']['P']=HPi['P'];
                    #end j circle
            
                    #output of this branch
                    if(self.Rev):
                        self.Bra[i]['HPi']['H']=HPi['H'];
                        self.Bra[i]['HPi']['P']=HPi['P'];
                    else:
                        self.Bra[i]['HPo']['H']=HPi['H'];
                        self.Bra[i]['HPo']['P']=HPi['P'];
                    self.Bra[i]['m']['m']=Sm['m'];
                    self.Bra[i]['m']['V']=Sm['V'];
                    Sm['m']=0;
                    Sm['V']=0;
                    self.Bra[i]['Para_Struc'][0]=P['VapL'];
                    self.Bra[i]['Para_Struc'][1]=P['TPL'];
                    self.Bra[i]['Para_Struc'][2]=P['LiqL'];
                    self.Bra[i]['Para_Struc'][3]=P['V_Vap'];
                    self.Bra[i]['Para_Struc'][4]=P['V_TP'];
                    self.Bra[i]['Para_Struc'][5]=P['V_Liq'];
                    self.Bra[i]['Para_Struc'][6]=P['m_Vap'];
                    self.Bra[i]['Para_Struc'][7]=P['m_TP'];
                    self.Bra[i]['Para_Struc'][8]=P['m_Liq'];    
                    self.Bra[i]['Para_Struc'][9]=P['UA_Vap'];
                    self.Bra[i]['Para_Struc'][10]=P['UA_TP'];
                    self.Bra[i]['Para_Struc'][11]=P['UA_Liq'];
                    
                else: #there is an equivalent branch for it
                    NoBra=self.Bra[i]['EqulNo'];
                    if (self.Rev):
                        self.Bra[i]['HPi']['H']=self.Bra[NoBra]['HPi']['H'];
                        self.Bra[i]['HPi']['P']=self.Bra[NoBra]['HPi']['P'];
                    else:
                        self.Bra[i]['HPo']['H']=self.Bra[NoBra]['HPo']['H'];
                        self.Bra[i]['HPo']['P']=self.Bra[NoBra]['HPo']['P'];
                    self.Bra[i]['m']['m']=self.Bra[NoBra]['m']['m'];
                    self.Bra[i]['m']['V']=self.Bra[NoBra]['m']['V'];
                    for NN in range(12):
                        self.Bra[i]['Para_Struc'][0]=self.Bra[NoBra]['Para_Struc'][0];
                        self.Bra[i]['Para_Struc'][1]=self.Bra[NoBra]['Para_Struc'][1];
                        self.Bra[i]['Para_Struc'][2]=self.Bra[NoBra]['Para_Struc'][2];
                        self.Bra[i]['Para_Struc'][3]=self.Bra[NoBra]['Para_Struc'][3];
                        self.Bra[i]['Para_Struc'][4]=self.Bra[NoBra]['Para_Struc'][4];
                        self.Bra[i]['Para_Struc'][5]=self.Bra[NoBra]['Para_Struc'][5];
                        self.Bra[i]['Para_Struc'][6]=self.Bra[NoBra]['Para_Struc'][6];
                        self.Bra[i]['Para_Struc'][7]=self.Bra[NoBra]['Para_Struc'][7];
                        self.Bra[i]['Para_Struc'][8]=self.Bra[NoBra]['Para_Struc'][8];
                        self.Bra[i]['Para_Struc'][9]=self.Bra[NoBra]['Para_Struc'][9];
                        self.Bra[i]['Para_Struc'][10]=self.Bra[NoBra]['Para_Struc'][10];
                        self.Bra[i]['Para_Struc'][11]=self.Bra[NoBra]['Para_Struc'][11];
                    #end NN loop
                #end else
                Bak['VapL'] = Bak['VapL']+self.Bra[i]['Para_Struc'][0];
                Bak['TPL'] = Bak['TPL']+self.Bra[i]['Para_Struc'][1];
                Bak['LiqL'] = Bak['LiqL']+self.Bra[i]['Para_Struc'][2];
                Bak['V_Vap'] = Bak['V_Vap']+self.Bra[i]['Para_Struc'][3];
                Bak['V_TP'] = Bak['V_TP']+self.Bra[i]['Para_Struc'][4];
                Bak['V_Liq'] = Bak['V_Liq']+self.Bra[i]['Para_Struc'][5];
                Bak['m_Vap'] = Bak['m_Vap']+self.Bra[i]['Para_Struc'][6];
                Bak['m_TP'] = Bak['m_TP']+self.Bra[i]['Para_Struc'][7];
                Bak['m_Liq'] = Bak['m_Liq']+self.Bra[i]['Para_Struc'][8];    
                Bak['UA_Vap'] = Bak['UA_Vap']+self.Bra[i]['Para_Struc'][9];
                Bak['UA_TP'] = Bak['UA_TP']+self.Bra[i]['Para_Struc'][10];
                Bak['UA_Liq'] = Bak['UA_Liq']+self.Bra[i]['Para_Struc'][11];
            
            #endif
            
        #end i loop
    
        #Return all backup values to the original P dictionary
        for key, value in Bak.iteritems():
            P[key] = value
        
        return 0
        





def CondCircuit(Type,mr,HPi,Tai,Ga,HPo,Tao,m,P, Ref):
    '''
    /***********************************************************************
    Condenser circuitry pattern.
    ***********************************************************************/
    '''

    CondInfo = StructCondClass("/InputDoc/CondStruc.xlsx", Ref);
    
    if (Type == 201): #user defined
        CondInfo._CondCircuit(mr,HPi,Tai,Ga,HPo,Tao,m,P) #returns updated (HPo, Tao, m, P)
#     elif type == 202: #Carrier RTU at Herrick
#         CondCircuitMultipleTube(mr,HPi,Tai,Ga,HPo,Tao,m,P);
#     elif type == 203: #Single finned tube condenser coil
#         CondCircuitSingleTube(mr,HPi,Tai,Ga,HPo,Tao,m,P);
#     elif type == 204: #user defined
#         CondCircuitMultipleTube_204(mr,HPi,Tai,Ga,HPo,Tao,m,P); 
    else:
        print("Condenser model "+str(type)+" is not found")
    
    return 0


def Condenser(Ref,
              mr,
              HPi,
              Tai,
              GaI,
              HPo,
              Tao,
              m,
              Cond_struc,#new
              Prms,
              NSeg,
              condInit):
    '''
    /**********************************************************************
    Model of an air cooled condenser coil.  The model is comprised of
    small pieces of single finned tubes put together to form rows.
    Rows of finned tubes are put together using elbows and manifolds
    to make an entire coil. This function defines the coil's circuit
    pattern.  Model parameters are stored in the file "cond.dat".
    
    Inputs:
        mr - refrigerant mass flow rate (kg/s).
        HPi - inlet refrigerant thermodynamic state (J/kg,Pa).
        Tai - inlet air temperature (K), relative humidity (-).
        Ga - air mass flux (kg/s*m^2).
            If Ga<0, then nominal mass flux is used from data file.
        Prms - array of three multiplicative correlation correction factors.
            Prms[0] - corrects air side convection coefficient.
            Prms[1] - corrects refrigerant side convection coefficient.
            Prms[2] - corrects refrigerant side pressure drop correlation.
            Prms[3] - fouling factor (1=no fouling,0=completely fouled).
            Set these parameters to 1 (Prms={1,1,1,1}) for no corrections.
            NSeg - number of segments (-1 default - use value from file)
    
    Outputs:
        HPo - outlet refrigerant thermodynamic state (J/kg,Pa).
        Tao - outlet air temperature (K), relative humidity (-).
        m - mass of charge in condenser (kg) and volume of charge (m^3).
    ***********************************************************************/
    '''

    P = CGP();

    if(condInit or NSeg>0):
        
        df = pd.read_excel('InputDoc/Cond_Data.xlsx',sheetname='definition',header = 0) #open the file and all the data are under column 'Condenser'
        
        P['type'] = df.Condenser[0]
        P['Di']=df.Condenser[1];           # tube inside diameter (m)
        P['L']=df.Condenser[2];            # tube length (m)
        P['xp']=df.Condenser[3];           # pipe wall thickness (m)
        P['Df']=df.Condenser[4];           # fin diameter (m)
        P['z']=df.Condenser[5];            # B.S., originally it is space between fins (m), here it is used as fin pitch, which includes a fin thickness
        P['th']=df.Condenser[6];           # fin thickness (m)
        P['vsp']=df.Condenser[7];          # vertical tube spacing (m)
        P['P_l']=df.Condenser[8];          # B.S. tube spacing along the airflow direction (m)
        P['NSeg']=int(df.Condenser[9]);    # number of lengthwise segments (-)
        if(NSeg>0):
            P['NSeg'] = NSeg;
        P['Brad']=df.Condenser[10];         # radius of return bend (m)
        
        GaNom=df.Condenser[11];          # nominal air mass flux (kg/s/m^2)

        P['Nbranchs']=int(df.Condenser[12]);        # number of branchs in main section
        P['Nmaintubes']=int(df.Condenser[13]);      # number of tubes in main section
        P['Nsubtubes']=df.Condenser[14];            # number of tubes in subsooled section
        P['Ndeep']=df.Condenser[15];                # B.S. number of high rows
        P['Frontal_A']=df.Condenser[16];            # B.S. frontal area

        hAirAdj=df.Condenser[17];        # hAirAdj
        hRefAdj=df.Condenser[18];
        hRefAdj_Sub=df.Condenser[19];   #B.S., it is a newly added adjustment factor for correcting the inaccurate ratio of two-phase heat transfer and single-phase heat transfer.
        PRefAdj=df.Condenser[20];
        FoulFac=df.Condenser[21];

        #B.S. ------------------------
        #new input parameter for condenser
        P['Microfin']=int(df.Condenser[22]); #microfin type, 0=smooth tube, 1=helical, 2=cross-grooved, 3=herringbone
        P['finN'] = df.Condenser[23]; #fin number in a micro-fin tube
        P['gama'] = df.Condenser[24];#fin apex angle in a micro-fin tube
        P['beta'] = df.Condenser[25];    #fin helix angle in a micro-fin tube
        P['finH'] = df.Condenser[26]; #fin height in a micro-fin tube
        P['w_b'] = df.Condenser[27]; #base width of a single fin, in the axial direction
        P['w_e'] = df.Condenser[28]; #top width of a single fin, in the axial direction
        P['w_z'] = df.Condenser[29]; #base distance between two neighboring fins, in the axial direction
        # airside fin type
        P['airFin'] = int(df.Condenser[30]);#1-plain, 2-corrugated, 3-slit, 4-louvered, 5-convex louvered, 6-smooth wavy, 7-spine
        P['sub1'] = df.Condenser[31];#fin substructures, for the possible use later, such as the waffle height of the wavy fin, louver number of the louverd fin
        P['sub2'] = df.Condenser[32];#fin substructures
        P['sub3'] = df.Condenser[33];#fin substructures
        P['sub4'] = df.Condenser[34];#fin substructures
        P['sub5'] = df.Condenser[35];#fin substructures
        P['K_T'] = df.Condenser[36];#400, this is the conductance factor of copper
        P['K_F'] = df.Condenser[37];#237, conductance factor of fin
        #------------------------------B.S.

        if(P['type']==202):
            P['cfmA']=df.Condenser[38];
            P['cfmB']=df.Condenser[39];
            P['cfmC']=df.Condenser[40];

        condInit=0
        
        # Calculate physical sizes needed in calculations
        # from input dimensions.
        # Length of a tube segment.
        P['Ls']=P['L']/P['NSeg'];

        #B.S.-----------------------------------
        # tube inside geometry
        if(P['Microfin']):
            P['D_b'] = P['Di'] + 2*P['finH'];#tube diameter at the base of the fin
            P['Do'] = P['D_b']+2*P['xp']; #Pipe outside diameter.
            P['D_m'] = (P['Di'] + P['D_b'])/2; #mean diameter of the micro-fin tube
            P['P_H'] = P['finN']*(P['w_z']+2*pow((pow(P['finH'],2.0)+pow((0.5*(P['w_b']-P['w_e'])),2.0)),0.5)+P['w_e']);# the hydraulical circumference
            P['Acs'] = pi*pow((P['D_b']/2.0),2.0)-P['finN']*0.5*P['finH']*(P['w_b']+P['w_e']);#cross area of the micro-fin tube, this is the actual cross-section area
            P['Dh_i'] = 4*P['Acs']/P['P_H'];#inside hydraulical diameter 
            P['Ax'] = pi/4.0*P['Di']*P['Di'];# Inside pipe cross sectional area, based on fin tips
            P['Api'] = pi*P['Di']*P['Ls'] ;# Inside pipe surface area (one tube segment), based on fin tips. All the two-phase heat transfer are calculated with Cavallini correlations, which are based on the heat transfer surface area on the tip of the fins. 
        else:
            P['finH'] = 0;#without fin height
            P['D_b'] = P['Di'];#tube diameter at the base of the fin
            P['Do'] = P['D_b']+2*P['xp']; #Pipe outside diameter.
            P['D_m'] = P['Di'] ; #averge diameter of the micro-fin tube
            P['Ax'] = pi/4.0*P['Di']*P['Di'];# Inside pipe cross sectional area
            P['Api'] = pi*P['Di']*P['Ls'] ;# Inside pipe surface area (one tube segment)
            P['P_H'] = pi*P['D_b'];# the hydraulical circumference
            P['Acs'] = P['Ax'];#cross area of the micro-fin tube, this is the actual cross-section area
            P['Dh_i'] = P['Di'];#inside hydraulical diameter
        #-----------------------------------------B.S.
        
        # Volume inside tube segment.
        P['Vs']=P['Acs']*P['Ls'];
        # Length of travel through bend.
        P['Blen']=pi*P['Brad'];
    
        
        #tube outside geometry
        #B.S. --------------------------------------------------
        #get the fin reference length, according to Schimidt equation
        M=P['vsp']/2e0;
        R_O=P['Do']/2;
        L = 0;

        if (P['Ndeep']>1):
            L=0.5e0*pow((pow(M,2e0)+pow(P['P_l'],2e0)),0.5);
        else:
            L = P['P_l']/2.0;

        BETA=L/M;
        PHI=M/R_O;
        R_EFF=R_O*1.27e0*PHI*pow((BETA-0.3e0),0.5e0);
        P['Df'] = R_EFF*2;# confirm the reference fin outside diameter
        #-------------------------------------------------------B.S.

        # Distance from outside of pipe to edge of fin. 
    
        P['y']=(P['Df']-P['Do'])/2 ;#B.S. correct the original term: P.y=(P.Df-P.Do)
        # Number of fins along one tube segment.
        P['N']= P['Ls']/(P['z']);#B.S. correct P.Ls/(P.z+P.th); here P.z is inputted as the fin pith while not the fin spacing
        # Outside pipe surface area (one tube segment).
        P['Apo']=pi*P['Do']*(P['Ls']-P['N']*P['th']);#B.S.

        # Fin wetted area (one tube segment).
        x = P['Df'] + P['th']/2.0;
        P['Af']=P['N']*pi/2.0*(x*x-P['Do']*P['Do']);# the fin has double faces
        # Air flow area (one tube segment).
        P['Aflow']=(P['vsp']-P['Do'])*(P['Ls']-P['N']*P['th']) ;#minimum airflow cross-sectional area
        # Hydrolic diameter
        P['Dh']=4*P['Aflow']/(P['Af']+P['Apo'])*P['vsp'];#hydraulic diameter of the fin channel


    #/* Process scaling factors */
    P['hAirAdj']=hAirAdj*Prms[0];
    P['hRefAdj']=hRefAdj*Prms[1];
    P['PRefAdj']=PRefAdj*Prms[2];
    #B.S. the following is for this is for adjusting the single-phase heat transfer
    P['hRefAdj_Sub']=hRefAdj_Sub*Prms[4];#B.S. Prms[9] is for avoiding being confusing with other adjust parameters #A.B. changes from Prms[9] to Prms[4]

    if(P['type']==202):
        # Date: Tue, 07 Mar 2000 09:36:46 -0500
        # To: Todd Rossi <rossi@mail.fielddiagnostics.com>
        # From: "Natascha S. Castro" <natascha@nist.gov>
        #
        # The polynomial that I extracted from experimental
        # data is : y = 16.677x2 - 7567.8x + 861160
        # I chose a polynomial fit.  It would be good to have the model read
        # coefficients from a file.  

        # Date: Thu, 23 Mar 2000 15:10:49 -0500
        # To: Todd Rossi <rossi@mail.fielddiagnostics.com>
        # From: "Natascha S. Castro" <natascha@nist.gov>
        #
        # y=cfm
        # x=psia

        # Date: Fri, 24 Mar 2000 09:41:39 -0500
        # To: Todd Rossi <rossi@mail.fielddiagnostics.com>
        # From: "Natascha S. Castro" <natascha@nist.gov>
        #    
        # The range for normal data is 220 psia to 245 psia, and the equation holds
        # for that range.  I think selecting extreme constant values is a good idea,
        # particularly because the data looks like a step, with a lower values level
        # and upper values level.  Based on the normal data, set discharge pressures
        # less than 225 psia and below to 2400 cfm and discharge pressures greater
        # than 245 to 3500 cfm.
        Pdis_kPa = HPi['P']/1000;
        Pdis_psia = Pdis_kPa/101.3*14.7;
        af_cfm=0.0;
        if(Pdis_psia<225.0):
            af_cfm = 2400.0;
        elif(Pdis_psia>245.0):
            af_cfm = 3500.0;
        else:
            af_cfm = P['cfmA']+Pdis_psia*(P['cfmB']+Pdis_psia*P['cfmC']);
            
        af_m3s = af_cfm/60.0*(12.0*0.0254)*(12.0*0.0254)*(12.0*0.0254);
        af_kgs = af_m3s/HAPropsSI('Vda','T',Tai['T'],'P',101325,'R',0)
        Ntubes = P['Nsubtubes']+P['Nbranchs']*P['Nmaintubes'];
        Aflow_seg = P['Aflow']*P['NSeg'];
        Aflow_tot = Aflow_seg*Ntubes;
        GaNom = af_kgs/Aflow_tot;

        Ga=GaNom*FoulFac*Prms[3];

    else:
        if(GaI<0):
            Ga=GaNom*FoulFac*Prms[3];
        else:
            Ga=GaI*Prms[3];#shenbo correct the original Ga=GaI

    # TR Calculate the air cfm after GA is established
    Ntubes = P['Nsubtubes']+P['Nbranchs']*P['Nmaintubes'];
    A_total = P['Aflow']*P['NSeg']*Ntubes; #total area air flows throguh for entire coil (m^2)
    v_air = HAPropsSI('Vda','T',Tai['T'],'P',101325,'R',0) #[m^3/kg dry air] # specific volume of air
    vfra = Ga*A_total*v_air; # nominal volumetric flow rate air
    P['cfma'] = vfra*(60.0/((12.0*0.0254)*(12.0*0.0254)*(12.0*0.0254))); # nominal cfm of air flow
    
    #B.S., prepare to set up the simple model
    P['V_TOT'] =0; P['V_TP']=0; P['V_Liq'] =0; P['V_Vap'] = 0;
    P['m_TOT'] = 0; P['m_TP'] =0; P['m_Liq'] =0; P['m_Vap'] =0;
    P['mr'] = 0; P['ma_TOT'] =0;
    P['UA_TOT']=0; P['UA_Liq'] =0;P['UA_Vap'] = 0; P['UA_TP'] =0;
    P['LiqL'] = 0; P['TPL'] = 0; P['VapL'] = 0;
    P['count1']=0; P['count2']=0;
    P['HP_TP1']['H']=0; P['HP_TP1']['P']=0;
    P['HP_TP2']['H']=0; P['HP_TP2']['P']=0;
    #-----------------------------------------B.S.

    P['GetP']=0;#don't calculate the airside pressure drop in the following procedure
    
    CondCircuit(P['type'],mr,HPi,Tai,Ga,HPo,Tao,m,P, Ref) #Return updated (HPo,Tao,m,P)
    
    DP_circuit = Circuit_DP_COND(Ga,Tai,Tao,P);#B.S. calculate the airside pressure drop cross the circuit
    #gOutData.Cond_AirDP = DP_circuit;#for output

    #B.S.---------------------------------------
    P['HP_in']['H'] = HPi['H'];
    P['HP_in']['P'] = HPi['P']; 
    P['HP_out']['H'] = HPo['H'];
    P['HP_out']['P'] = HPo['P'];
    P['mr'] = mr; 
    P['Ga'] = Ga; 
    P['Tai']=Tai; 
    P['ma_TOT'] = Ga*P['Aflow']*P['NSeg']*(P['Nbranchs']*P['Nmaintubes']+P['Nsubtubes']);
    P['Cond_AirDP']=DP_circuit       #A.B.
    Build_Cond(P,Ref);
    for key, value in P.iteritems(): #A.B. output the condenser construct (copy back dictionary "P" to dictionary "Cond_struc")
        Cond_struc[key] = value
    #---------------------------------------B.S.

    return 0

def ExerciseCondenser():
    
    Ref = 'R410A'
    #design conditions
    SC = 10.9512 #subcooling F
    CTOA_DL = 29.511+4 # F, add 4F due to pressure drop in cond coil, 18F is at the liquid line
    CTOA_LL = 29.511 # F, at the liquid line
    AMB = 95 #[F]
    AMB_RH = 0.2 #[-] #A.B.
    MR = 955.572 # lbm/hr, refrigerant mass flow rate
    GOA = 3.82 # kg/s/m^2, outdoor air mass flux
    DSH = 67.13 # F, discharge superheat 
    
    #initialize the values
    hpo = {'H':0.0,'P':0.0}
    m = {'m':0.0,'V':0.0}
    Tao = {'T':0.0, 'P':0.0} #A.B.
    CondPrms = [1,1,1,1,1] #correction factor (see Condenser model for more details)
    NSeg = -1 #keep it to -1 to use the default number of segment (10 segments)
    Cond_struc= CGP() #initialized condenser staructure

    print('Discretize Condenser model ')
    print()
    
    #inlet air state [K]
    Tai = {'T':F2K(AMB), 'P':AMB_RH}
    #air mass flux (kg/s/m^2)
    Goa = GOA;
    #refrigerant mass flow rate [kg/s]
    mr = lbh2kgs(MR);
    #inlet refrigerant state
    ctoa = DeltaF2K(CTOA_DL); # condensing temperature over ambient
    Tsh = DeltaF2K(DSH); #discharge superheat
      
    Tsat_dl = Tai['T']+ctoa #Temperature of saturated vapor [K]
    P_dl = PropsSI("P", "T", Tsat_dl, "Q", 1, Ref) #CHECK THIS it might be either liquid or vapor 
    txpi = {'T':Tsat_dl+Tsh,'X':1,'P':P_dl} #initially assume outlet temperature equals to saturated vapor pressure (txpo_P)
    hpi = TXPtoHP(txpi, Ref)  #inlet refrigerant state
    
    #*Run the evaporator model*#
    Condenser(Ref,mr,hpi,Tai,Goa,hpo,Tao,m,Cond_struc,CondPrms,NSeg,condInit=1); #A.B. return updated values for  (hpo, Tao, m, Cond_struc)
              
    txpo = HPtoTXP(hpo, Ref)
    
    Tsat_ll = PropsSI("T", "P", txpo['P'], "Q",1, Ref) #reftplthP.Tsat(txpo.P) #[K] #CHECK this it might be vapor of liquid
    
    # Calculate the total heat transfer and then the UA = Qdot/DT
    Qdot = W2BTUh(mr*(hpi['H']-hpo['H'])); #[BTU/hr ]- for all 4 sections of the heat exchanger
    DT = DeltaK2F(Tsat_ll - Tai['T']); #[F]
    UA = Qdot/DT; #[BTU/hr/F] - for all 4 sections of the heat exchanger
         
    #print inputs
    print()
    print("Tai (F)", K2F(Tai['T']))
    print("CFM (cfm): ",Cond_struc['cfma'])
    print("mr (lbm/hr): ", MR)
    print("CTOA_dl (F): ",DeltaK2F(Tsat_dl-Tai['T']))
    print("Tsh_dl (F): ",DeltaK2F(txpi['T']-Tsat_dl))
    print("X_dl (0-1): ",txpi['X'])
    #print outputs
    print()
    print("CTOA_ll (F): ",DeltaK2F(Tsat_ll-Tai['T'])) #CTOA_ll
    print("Tsc_ll (F): ",DeltaK2F(Tsat_ll-txpo['T'])) # Tsc_ll
    print("X_ll (0-1): ",txpo['X'])
    print("dP (psi): ",kPa2psi((txpi['P']-txpo['P'])/1000))
    print("DTa (F): ",DeltaK2F(Tao['T']-Tai['T']))
    print("mass (kg): ",m['m']) #charge
    print("Qdot (BTU/hr): ",Qdot)
    print("UA (BTU/hr/F): ",UA)
    print("Hin (kJ/kg): ",hpi['H']/1000)
    print("Hout (kJ/kg): ",hpo['H']/1000)
    print()
    return 0


    
if __name__=='__main__': 
    print('Hello world')
    print()
    from time import time
    t1=time()
    ExerciseCondenser()
    print ('Took '+str(time()-t1)+' seconds to run discretize Condenser model')