from __future__ import division, print_function, absolute_import
from math import pi,log,exp,tanh

from scipy.optimize import brentq #solver to find roots (zero points) of functions
#import numpy as np
import pandas as pd

from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI

from ACHP.convert_units import *

from extra_functions import HPtoTXP, TXPtoHP, PropertyTXPth, EVA_Get_Q_dic, toTXP, PreAcc, ETdim, HPtoTP, WPtoTP, THtoTP, WHtoTP, EvapNode, EvapBranch, TubeEvap, TubEvpSeg
from PRESSURE import dPelbow, dPmom, GET_PreAcc
from VOLUME import VolumeALL
from CORR import Circuit_DP_EVAP, ConvCoeffAir_EVA, FinEffect_Schmidt, ConvCoeffSP, FricDP, ConvCoeffEvapTP_microfin
from CMINE import CmineCrossFlow_dry, CmineCrossFlow_wet
            
        
        
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
    return HPo, m, P

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

    return HPo, m, P
    
def EvapTubeL_Rev(Gr,#refrigerant mass flux
                  HPo,#refrigerant outlet and inlet state
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
    TXP_bak = TXPo.copy();#backup the outlet state of this segment
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
    Eva_dim = EVA_Q['P'];#B.S. get the evaporator struct parameters

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
    Eva_dim = EVA_Q['P']; #B.S. get the evaporator struct parameters

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

    CP_M=HAPropsSI('cp','P',101325,'T',EVA_Q['TPi']['T'],'R',EVA_Q['TPi']['P']) #inlet air enthalpy [J/kg humid air/K]
    
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
    
    CP_M=HAPropsSI('cp','P',101325,'T',EVA_Q['TPi']['T'],'R',EVA_Q['TPi']['P']) #inlet air enthalpy [J/kg humid air/K]

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






class StructEvapClass():
    
    def __init__(self,filename,Ref):#mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,P,Ref):
        '''Initialize'''
        #reading the data file in
        #the followings from   Class StructEvap(filename) in the ACMODLE
        
        self.Rev = 0 # for reverse calculation, 1--for reversed evaporator calculation, 0-- for forward evaporator solver
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
        
        df = pd.read_excel(filename,sheetname=0,header = 0)
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
            
            P, HPi = self.Cal_HP(P,0,HPi);    #heat transfer and pressure drop calculation
            
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
        
        
            P, HPi = self.Cal_HP(P,1,HPi);#heat transfer and pressure drop calculation
            
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
            TPo_temp = {'T':0.0,'P':0.0}
            TPo_temp['T']=TPi['T']-5;
            TPo_temp['P']=0.8;
            TPo = WHtoTP(WH_out,TPo_temp)
            #TPo = {'T': HAPropsSI('T','P',101325,'H',WH_out['H'],'W',WH_out['W']), 'P': HAPropsSI('R','P',101325,'H',WH_out['H'],'W',WH_out['W'])} 
        except:
            print('_EvapCircuit_Fwd :: check TPo that need to be numerically solved!')
            print('Other source of exception might be due to relative humidity ~100%, try to solve TPo with relative humidity of 95%')
            TPo = {'T': HAPropsSI('T','P',101325,'H',WH_out['H'],'R',0.999), 'P': 0.999}

    
        #refrigerant outputs
        HPo['H']=H_out;
        HPo['P']=HPo['P'];
        HPi['H']=H_in;
        HPi['P']=P_in;
    
        #Sm->m=Sm->m + (2.87-2.766)+(3.10-3.05602)/(2.0792-7.955)*(P->VapL-7.955);#for three points charge tuning
        
        return HPo, HPi, TPo, Sm, Aflow, P
    
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
            
            P, HPi = self.Cal_HP(P,0,HPi);#heat transfer and pressure drop calculation
            
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
            
            P, HPi = self.Cal_HP(P,1,HPi);#heat transfer and pressure drop calculation

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
            TPo_temp = {'T':0.0,'P':0.0}
            TPo_temp['T']=TPi['T']-5;
            TPo_temp['P']=0.8;
            TPo = WHtoTP(WH_out,TPo_temp)
            #TPo = {'T': HAPropsSI('T','P',101325,'H',WH_out['H'],'W',WH_out['W']), 'P': HAPropsSI('R','P',101325,'H',WH_out['H'],'W',WH_out['W'])} 
        except:
            print('_EvapCircuit_Rev :: check TPo that need to be numerically solved!')
            print('Other source of exception might be due to relative humidity ~100%, try to solve TPo with relative humidity of 95%')
            TPo = {'T': HAPropsSI('T','P',101325,'H',WH_out['H'],'R',0.999), 'P': 0.999}
        
        
        return HPo, HPi, TPo, Sm, Aflow, P
    
    def Cal_HP(self,P,Pos,HPi):
        ''''function to calculate the heat transfer and pressure drop'''
        
        #heat transfer and pressure drop of the two-phase region
        Logic=0;
        WHo = {'W':0,'H':0};
        Sm = {'m':0,'V':0};
        mi = {'m':0,'V':0};
        #Bak = ETdim()
        
        Bak=P.copy()  ;#keeping the information of the evaporator
        
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
                                    self.Tub[TubeN]['Seg'][k]['TPi'] = WHtoTP(WHo,self.Tub[TubeN]['Seg'][k]['TPi']);
                                    #self.Tub[TubeN]['Seg'][k]['TPi']['T'] = HAPropsSI('T','P',101325,'H',WHo['H'],'W',WHo['W'])
                                    #self.Tub[TubeN]['Seg'][k]['TPi']['P'] = HAPropsSI('R','P',101325,'H',WHo['H'],'W',WHo['W'])
                                except:
                                    print('Cal_HP :: check self.Tub that need to be numerically solved!')
                                    print('Other source of exception might be due to relative humidity ~100%, try to solve TPo with relative humidity of 95%')
                                    raise
                                    #self.Tub[TubeN]['Seg'][k]['TPi'] = {'T': HAPropsSI('T','P',101325,'H',WHo['H'],'R',0.999), 'P': 0.999}
                                    
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
                                HPi, WHo, mi, P = EvapTubeL_Rev(self.Bra[i]['Gr'],HPi,self.Tub[TubeN]['Ga'],self.Tub[TubeN]['Seg'][realk]['TPi'],WHo,mi,P, self.Ref) #return updated (HPi, WHo, mi, P)
                            else:
                                HPi, WHo, mi, P = EvapTubeL_Fwd(self.Bra[i]['Gr'],HPi,self.Tub[TubeN]['Ga'],self.Tub[TubeN]['Seg'][realk]['TPi'],WHo,mi,P, self.Ref) #return updated (HPi, WHo, mi, P)
                
                            self.Tub[TubeN]['Seg'][realk]['WHo']['H']=WHo['H'];
                            self.Tub[TubeN]['Seg'][realk]['WHo']['W']=WHo['W'];
                            self.Tub[TubeN]['m']['m']=self.Tub[TubeN]['m']['m']+mi['m'];
                            self.Tub[TubeN]['m']['V']=self.Tub[TubeN]['m']['V']+mi['V'];
                            Sm['m']=Sm['m']+mi['m'];
                            Sm['V']=Sm['V']+mi['V'];
                        #end k circle
        
                        
                        if (self.Rev):
                            HPi, mi, P = EvapTubeBend_Rev(self.Bra[i]['Gr'],HPi,mi,P, self.Ref) #return updated HPi, mi, P
                        else:
                            HPi, mi, P = EvapTubeBend_Fwd(self.Bra[i]['Gr'],HPi,mi,P, self.Ref) #return updated HPi, mi, P
                        
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
    
        P=Bak
        
        return P, HPi
        





def EvapCircuit(Type,mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D, Ref):
    
    EvapInfo = StructEvapClass("../InputDoc/EvapStruc.xlsx", Ref)
    EvapInfo.Rev = D['REV'];
    
    if (Type == 301): #user defined
        HPo, HPi, TPo, Sm, Aflow, D = EvapInfo._EvapCircuit_Fwd(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D)
#     elif type == 302: #Carrier RTU at Purdue
#         EvapCircuit301(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D); 
#     elif type == 303: #Single finned tube
#         EvapCircuit302(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D);
#     elif type == 304: #this function is for CS30TR
#         EvapCircuit303(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D); 
    else:
        print("Evaporator model "+str(type)+" is not found")
    
    return HPo, HPi, TPo, Sm, Aflow, D

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
    D = ETdim()
     
 
    if (evapInit > 0 or NSeg>0):
         
        df = pd.read_excel(filename,header = 0) #open the file and all the data are under column 'Evaporator'
         
        D['type'] = df.Evaporator[0]
        D['Di'] = df.Evaporator[1]             # tube inside diameter
        D['L'] = df.Evaporator[2]              # tube length 
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
        D['Nrows'] = df.Evaporator[14]        # number of rows high
        D['Ndeep'] = df.Evaporator[15]        # number of rows deep
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
       
        #endSign=int(df.Evaporator[38])      #-256
        #if(endSign != -256):
        #    print ("Evaporator parameter input wrong")
        #    raise
   
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
    A_total = D['Aflow']*D['NSeg']*D['Nrows'] #total area air flows throguh for entire coil (m^2)
    v_air = HAPropsSI("Vda", "T", TPi['T'], "P", 101325, "R", 0) # specific volume of air [m^3/kg dry air]
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
    
    HPo, HPi, TPo, Sm, Aflow, D = EvapCircuit(D['type'],mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,D, Ref)
    
    DP_circuit = Circuit_DP_EVAP(Ga,TPi,TPo,D)#B.S., calculate the airside pressure drop across the circuit
    #gOutData.Eva_AirDP = DP_circuit#for output
    
    #B.S.---------------------------------------
    D['HP_out']['H'] = HPo['H']      #suciton site of the evaporator
    D['HP_out']['P'] = HPo['P']
    D['HP_in']['H'] = HPi['H']      #entrance of the evaporator
    D['HP_in']['P'] = HPi['P'] 
    D['mr'] = mr
    D['Ga'] = Ga
    D['TPi']['T'] = TPi['T']
    D['TPi']['P'] = TPi['P']
    D['ma_TOT'] = Ga*D['Aflow']*D['NSeg']*D['Nrows']
    #Build_Eva(&D)        #correlate the moving boundary and lumped evaporator model
    Evap_struc = D.copy()          #output the evaporator construct
    #---------------------------------------B.S.
    print(DP_circuit, Evap_struc)
    return (HPo, HPi, TPo, Sm, Aflow, Evap_struc)

def ExerciseEvaporator():
    
    Ref = 'R410A'
    #design conditions
    SH = 23.004 #superheat (F)
    #SC = 10.9512 #subcooling F (***no effect)
    ET_SL = 48.3264-23 # F, evaporation temperature at the suction line
    #ET_DT = 45+4 # F, add 4F due to evaporator coild pressure drop
    #CTOA_DL = 29.511+4 # F, add 4F due to pressure drop in cond coil, 18F is at the liquid line (** no changes)
    #CTOA_LL = 29.511 # F, at the liquid line  (** no changes)
    #AMB = 95 # F
    RA = 83.048#80 # indoor room dry bulb temperature (F)
    RARH = 0.51 # indoor room humidity (-)
    #RAWB = 67.856 #  indoor room wetbulb temoperature (F)
    MR = 955.572 # lbm/hr, refrigerant mass flow rate
    #GOA = 3.82 # kg/s/m^2, outdoor air mass flux (this is for condenser component)
    GIA = 2.55 # kg/s/m^2, indoor air mass flux (0.9237/0.3623)
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
    P_sl = PropsSI("P", "T", Tsat_sl, "Q", 0, Ref)
    txpo = {'T':Tsat_sl+Tsh,'X':1,'P':P_sl} #initially assume outlet temperature equals to saturated vapor pressure (txpo_P)
    hpo = TXPtoHP(txpo, Ref)  #outlet refrigerant state
    Twbi = HAPropsSI("Twb", "T", TPi['T'], "P", 101325, "R", TPi['P']) #inlet air wet bulb temperature [K]
    
    #*Run the evaporator model*#
    Evap_struc['REV'] = 1 #0-- forward evaporator solver, 1-- reversed evaporator (backward solver)
    hpo, hpi, TPo, m, Aflow, Evap_struc = Evaporator(Ref,'../InputDoc/Evap_Data.xlsx',mr,hpo,Ga,TPi,hpi,TPo,m,Aflow,Evap_struc,Evap_prms,NSeg,evapInit=1)
    
    #inlet refrigerant state and   
    #saturation temperature computed from the inlet refrigerant state
    txpi = HPtoTXP(hpi, Ref)
    
    Tsat_dt = PropsSI("T", "P", txpi['P'], "Q",0, Ref) #PropertyTXPth('T',txpi,Ref) #[K]
     
    Twbo = HAPropsSI("Twb", "T", TPo['T'], "P", 101325, "R", TPo['P']) #outlet air wet bulb temperature [K]
     
    #Calculate the total heat transfer and then the UA = Qdot/DT
    Qdot = W2BTUh(mr*(hpo['H']-hpi['H'])) #[BTU/hr] - for all 4 sections of the heat exchanger
    DT = DeltaK2F(TPi['T'] - Tsat_sl) #[F]
    UA = Qdot/DT #[BTU/hr/F] - for all 4 sections of the heat exchanger 
         
    #print inputs
    print()
    print("Tai (F), RHai (%): ", K2F(TPi['T']),TPi['P']*100.0)
    print("CFM (cfm): ",Evap_struc['cfma'])
    print("mr (lbm/hr): ", MR)
    print("ETBRA_sl (F): ",DeltaK2F(TPi['T']-Tsat_sl))
    print("Tsh_sl (F): ",DeltaK2F(txpo['T']-Tsat_sl))
    print("X_sl (0-1): ",txpo['X'])
    if (NSeg>0):
        print("NSeg: ",NSeg)
    else:
        print("NSeg: ",'Default') 
    #print outputs
    print()
    print("ETBRA_dt (F): ",DeltaK2F(TPi['T']-Tsat_dt)) #ETBRA_dt
    print("X_dt (0-1): ",txpi['X'])
    print("dP (psi): ",kPa2psi((txpi['P']-txpo['P'])/1000)) ###BEWARE that txpo_P was set to saturated liquid pressure, check if the code modifies it
    print("DTa (F): ",DeltaK2F(TPi['T']-TPo['T']))
    print("DTwb (F): ",DeltaK2F(Twbi-Twbo))
    print("mass (kg): ",m['m']) #charge
    print("Qdot (BTU/hr): ",Qdot)
    print("UA (BTU/hr/F): ",UA)
    print("Hin (kJ/kg): ",hpi['H']/1000.0)
    print("Hout (kJ/kg): ",hpo['H']/1000.0)
    print("Tsat_sl (C): ",K2C(Tsat_sl))
    print("Tsat_dt (C): ",K2C(Tsat_dt))
    print()
    return 0



    
if __name__=='__main__': 
    print('Hello world')
    print()
    from time import time
    t1=time()
    ExerciseEvaporator()
    print ('Took '+str(time()-t1)+' seconds to run discretize Evaporator model')