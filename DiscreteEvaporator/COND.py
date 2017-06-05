from __future__ import division, print_function, absolute_import
from math import pi,log,exp,tanh

from scipy.optimize import brentq, minimize
import numpy as np
import pandas as pd

from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI

from ACHP.convert_units import *

from extra_functions import CGP, TXPtoHP, HPtoTXP, PropertyTXPth, PropertyTXPtr, CondNode, CondBranch, TubeCond, TubCndSeg, HatoTa, ACRP, CLRP, PreAcc
from PRESSURE import LossCoeffReturnBend, GET_PreAcc
from VOLUME import VolumeALL
from CORR import Circuit_DP_COND, ConvCoeffAir_CON, FinEffect_Schmidt, ConvCoeffInside, FricDP
from CMINE import CmineCrossFlow_dry




def CondTubeL_new(Gr,   #refigerant mass flux
                  HPi,  #refrigerant inlet state
                  Ga,   #air mass flux
                  Tai,  #air inlet temperature
                  hao,  #air outlet enthalpy
                  m,    #inner volume and the refrigerant in the segment
                  P,    #condenser struct
                  Ref): #Refrigerant
    '''
    /*********************************************************************
    Condenser tube segment model.
    Inputs:
        Gr = refrigerant mass flux (kg/s/m^2)
        HPi = refrigerant inlet state (H:J/kg,P:Pa)
        Ga = air mass flux (kg/s/m^2)
        Tai = air inlet temperature (K)
    Outputs:
        HPi = refrigerant outlet state (H:J/kg,P:Pa)
        hao = air outlet enthalpy (J/kg)
        m = mass of charge in return bend (kg)
    *********************************************************************/
    '''

    # convert inlet state into TXP format
    TXPi = HPtoTXP(HPi,Ref)
    
    #air heat transfer 
    ho = ConvCoeffAir_CON(Tai,Ga,P);#B.S. ConvCoeffAir_CON(Tai,Ga,P) it is a new function, which can call heat transfer calculation corresponding to different fin type
    ho = P['hAirAdj']*ho
    
    phi = FinEffect_Schmidt(ho,237,P['th'],P['y'],P['Do']);#fin efficiency, using Schmidt equation 
    P['Ro'] = 1/(ho*(P['Apo']+phi*P['Af']));#airside thermal resistance
    
    P['airT'] = Tai['T'];#B.S., this may be used as the local air inlet temperature, while Tai is used as the global inlet air temperature
    P['Ga'] = Ga;#B.S., keep the air state parameters for late use

    # Mass flow rates of air and refrigerant.
    ma = Ga*P['Aflow'];
    mr = Gr*P['Ax'];

    # Calculate heat transfered per unit mass of refrigerant.
    hi = ConvCoeffInside(TXPi,Gr,P['Di'], P, Ref);#B.S., this function can call different heat transfer calculation, correponding to different tube type
    hi = P['hRefAdj']*hi
    
    Ri = 1/(hi*P['Api']);#refrigerant side thermal resistance
    R_W=log(P['Do']/(P['Do']-2.0*P['xp']))/(2.0*pi*P['K_T']*P['Ls']);
    #R = 0.0;

    R = (P['Ro']+R_W+Ri);
    
    #B.S.---------------------------------------------
    if(TXPi['X']>=1):#B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.    
        P['UA_Vap']=P['UA_Vap']+1/R;#B.S, vapor-phase heat transfer conductance
        P['VapL'] = P['VapL']+ P['Ls'];#B.S., vapor-phase tube length
        R = (P['Ro']+R_W+Ri);
    elif(TXPi['X']<=0):
        P['LiqL'] = P['LiqL']+ P['Ls'];#B.S., liquid-phase tube length
        P['UA_Liq']=P['UA_Liq']+1/R;#B.S. liquid-phase heat transfer conductance
        R = (P['Ro']*P['hAirAdj']+R_W+Ri);
        #R = (P['Ro']+R_W+Ri);
    else:
        P['UA_TP']=P['UA_TP']+1/R;#B.S. two-phase heat transfer conductance
        P['TPL'] = P['TPL']+ P['Ls'];#B.S. two-phase heat transfer tube length
        R = (P['Ro']+R_W+Ri);
    #-----------------------------------------B.S.

    q = (TXPi['T']-Tai['T'])/mr*CmineCrossFlow_dry(R,mr,ma,TXPi,Tai,Ref);

    # Calculate outlet state.
    # q = dh (energy balance)

    HPi['H'] = HPi['H'] - q;

    # calculating pressure drop
    D_H=P['Dh_i'];#B.S., get the hydraulic diameter for calculating the pressure drop with Kedzierski correlation

    Gr = mr/P['Acs'];#B.S., this mass flux is for calculating the pressure drop
    DP_FR = FricDP(TXPi, Gr, q, P, Ref);#calculate the fictional pressure drop, no matter single-phase or two-phase flow
    
    if(TXPi['X']>0.05 and TXPi['X']<0.99):
        #B.S. packing up the parameters for calculating the acceleration pressure drop
        Preacc=PreAcc();
        Preacc['DP_FR']=DP_FR;
        Preacc['G']=Gr;
        Preacc['H_OUT']=HPi['H'];
        Preacc['P_IN']=TXPi['P'];
        Preacc['X_IN']=TXPi['X'];
        DP_ACC=brentq(GET_PreAcc,10000,-10000,args=(Ref, Preacc),xtol=1e-7,rtol=6e-8,maxiter=40)#B.S. two-phase acceleration pressure drop
        P_out=TXPi['P']-(DP_FR+DP_ACC)*P['PRefAdj'];
    else:
        P_out=TXPi['P'] - DP_FR*P['PRefAdj'];#assume single-phase flow has no accelerattion pressure drop


    HPi['P']=P_out;
    
    # Determine mass of charge.  It is based on the inlet state.
    #    The specific volume is recalculated so that a different model
    #    can be used from the one used to calculate the pressure drop.

    v = VolumeALL(TXPi,Gr,P['Di'],-mr*q/P['Api'],Ref);        # seperate flow model

    m['V'] = P['Ls']*P['Acs'];#B.S. use the real cross-sectional area to calculate the inner volume of micro-fin tube
    m['m'] = m['V']/v;

    # Calculate output air state.
    hai = HAPropsSI('H','T',Tai['T'],'P',101325,'R',0) #[J/kg dry air]
    
    hao['H'] = hai+q*mr/ma;

    #B.S.----------------------------------------
    TXPo = HPtoTXP(HPi,Ref);
    if(TXPi['X']>=0.999999999):#B.S., the following calculate the whole thermal resistances separately, prepared to adjust them separately.    
        P['m_Vap'] = P['m_Vap'] + m['V']/v;#B.S. vapor-phase mass
        P['V_Vap'] = P['V_Vap'] + m['V'];#B.S. vapor-phase volume
        if(TXPo['X']<0.999999999):
            P['HP_TP1']['H']=HPi['H']*mr+P['HP_TP1']['H'];#corresponding to inlet state of the first segment that involve with two-phase heat transfer calculation
            P['HP_TP1']['P']=HPi['P']*mr+P['HP_TP1']['P'];#corresponding to inlet state of the first segment that involve with two-phase heat transfer calculation
            P['count1'] = P['count1'] + mr;#count the number of circuit
    elif(TXPi['X']<=0.00000000001):
        P['m_Liq'] = P['m_Liq'] + m['V']/v;#B.S. liquid-phase mass
        P['V_Liq'] = P['V_Liq'] + m['V'];#B.S. liquid-phase volume
    else:
        P['m_TP'] = P['m_TP'] + m['V']/v;#B.S. two-phase mass
        P['V_TP'] = P['V_TP'] + m['V'];#B.S. two-phase volume
        if(TXPo['X']<=0.00000000001):
            P['HP_TP2']['H'] = P['HP_TP2']['H']+mr*HPi['H'];#corresponding to the inlet state of the last segment that involve with two-phase heat transfer calculation
            P['HP_TP2']['P'] = P['HP_TP2']['P']+mr*HPi['P'];#corresponding to the inlet state of the last segment that involve with two-phase heat transfer calculation
            P['count2'] = P['count2'] + mr;#count the number of circuit
    #---------------------------------------------B.S.
    return 0
      
def CondReturnBend(G,HPi,m,P,Ref):
    '''
    /*********************************************************************
    Condenser tube return bend model.
    Inputs:
        G = refrigerant mass flux (kg/s/m^2)
        HPi = refrigerant inlet state (h,P)
    Outputs:
        HPi = refrigerant outlet state (h,P)
        m = mass of charge in return bend (kg)
    *********************************************************************/
    '''

    Q=CLRP()

    TXPi=HPtoTXP(HPi,Ref);
    
    #Q['vi']=VolumeALL(TXPi,G,P['Di'],Ref);    # seperate flow model
    Q['vi'] = 1/PropertyTXPth('D',TXPi,Ref);    # homogeneous flow model

    Q['K']=LossCoeffReturnBend(TXPi,G,P['Brad'],P['Di'],Ref);

    Q['D']=P['Di'];
    Q['HPi']['H']=HPi['H'];
    Q['HPi']['P']=HPi['P'];
    Q['G']=G;
    Q['q']=0;          # adiabatic

    # Solve for outlet state.
    try:
        X = (HPi['H'],HPi['P']) #initial guess
        cons = {'type': 'ineq', 'fun': HPLimitsConst, 'args': (Ref,)}
        res = minimize(CompLossResid, X, args=(Q,Ref), method='SLSQP', bounds=(), constraints=cons, options={'eps':1e-4, 'maxiter': 100, 'ftol':1e-6})
        HPi['H'] = res.x[0]
        HPi['P'] = res.x[1]
    except:
        print("CondReturnBend G=[] HPi['H']={} HPi['P']={}".format(G,HPi['H'],HPi['P']));
        raise

    #/* Determine mass of charge.  It is based on the inlet state.
    #    The specific volume is recalculated so that a different model
    #    can be used from the one used to calculate the pressure drop.  */

    v = VolumeALL(TXPi,G,P['Di'],0,Ref);    # seperate flow model
    #v = 1/PropertyTXP('D',TXPi,Ref);    # homogeneous flow model

    m['V']=P['Blen']*P['Ax'];
    m['m']=m['V']/v;

    #B.S.----------------------------------------------------
    if(TXPi['X']==1):#B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.    
        P['m_Vap'] = P['m_Vap'] + m['V']/v;#B.S. vapor-phase mass
        P['V_Vap'] = P['V_Vap'] + m['V'];#B.S. vapor-phase volume
    elif(TXPi['X']<1 and TXPi['X']>0.0):
        P['m_TP'] = P['m_TP'] + m['V']/v;#B.S. two-phase mass
        P['V_TP'] = P['V_TP'] + m['V'];#B.S. two-phase volume
    else:
        P['m_Liq'] = P['m_Liq'] + m['V']/v;#B.S. liquid-phase mass
        P['V_Liq'] = P['V_Liq'] + m['V'];#B.S. liquid-phase volume
    #---------------------------------------------B.S.
    
    return 0

def CondMan(ni,no,D,G,HPi,Ref):
    '''
    /*********************************************************************
    Condnenser manifold model.  Combines tubes together.  Holds no
    refrigerant charge.
    Inputs:
        ni = number of inlet tubes
        no = number of outlet tubes
        D = tube diameters, assumed same (m)
        G = inlet refrigerant mass flux (kg/s/m^2)
        HPi = refrigerant inlet state (h,P)
    Outputs:
        HPi = refrigerant outlet state (h,P)
        G = outlet refrigerant mass flux (kg/s/m^2)
    *********************************************************************/
    '''
    P=ACRP()
    Q=CLRP()
    #TXPi={'T':0.0,'X':0.0,'P':0.0}

    ##/* isentropic area change */
    P['Di']=D;
    P['HPi']['H']=HPi['H'];
    P['HPi']['P']=HPi['P'];
    P['Gi']=G;
    G *= ni/no;
    P['Go']=G;

    TXPi = HPtoTXP(HPi,Ref);
   
    P['si'] = PropertyTXPth('S',TXPi,Ref);

    # specific volume
    # P['vi']=VolumeALL(TXPi,P['Gi'],P['Di'],Ref);    # seperate flow model
    P['vi'] = 1/PropertyTXPth('D',TXPi,Ref);            # homogeneous flow model
    
    try:
        X = (HPi['H'],HPi['P']) #initial guess
        cons = {'type': 'ineq', 'fun': HPLimitsConst, 'args': (Ref,)}
        res = minimize(AreaChangeResid, X, args=(P,Ref,), method='SLSQP', bounds=(), constraints=cons, options={'eps':1e-4, 'maxiter': 100, 'ftol':1e-6})
        HPi['H'] = res.x[0]
        HPi['P'] = res.x[1]
        G = P['Go']
    except:
        print("CondMan: (AreaChangeResid) Gr={} HPi['H']={} HPi['P']={}".format(P['Go'],HPi['H'],HPi['P']));
        raise

    # irreversible enterance loss
    Q['K']=0.5;
    Q['D']=D;
    Q['HPi']['H']=HPi['H'];
    Q['HPi']['P']=HPi['P'];
    Q['G']=G;
    Q['q']=0;

    TXPi=HPtoTXP(HPi,Ref);

    #specific volume
    #Q['vi'] = VolumeALL(TXPi,Q['G'],Q['D']);          # seperate flow model
    Q['vi'] = 1/PropertyTXPth('D',TXPi,Ref);           # homogeneous flow model
    try:
        X = (HPi['H'],HPi['P']) #initial guess
        cons = {'type': 'ineq', 'fun': HPLimitsConst, 'args': (Ref,)}
        res = minimize(CompLossResid, X, args=(Q,Ref,), method='SLSQP', bounds=(), constraints=cons, options={'eps':1e-4, 'maxiter': 100, 'ftol':1e-6})
        HPi['H'] = res.x[0]
        HPi['P'] = res.x[1]
        G = Q['G']
    except:
        print("CondMan: (CompLossResid) Gr={} HPi['H']={} HPi['P']={}".format(Q['G'],HPi['H'],HPi['P']));
        raise

    return G

def AreaChangeResid(X,Params,Ref):
    '''
    /********************************************************************
    Generates two residuals that must be zero at the corrent outlet
    state of the manifold.  Energy conservation requires that the
    stagnation enthalpy across the manifold be zero and the assumption
    that the process is reversible and adiabatic results in and
    isentropic process that forms the second residual (s1=s2).
    ********************************************************************/
    '''
    HPo = {'H':0.0,'P':0.0}
    HPo['H'] = X[0];
    HPo['P'] = X[1];
    
    P = Params;

    # outlet specific volume and entropy
    TXPa = HPtoTXP(HPo,Ref);
    
    so = PropertyTXPth('S',TXPa,Ref);

    # specific volume
    # double vo = VolumeALL(TXPa,P['Go'],P['Di']);    # seperate flow model
    vo = 1/PropertyTXPth('D',TXPa,Ref);                # homogeneous flow model
    
    ##/* inlet and outlet stagnation enthalpy */
    xi = P['vi']*P['Gi'];
    h0i = P['HPi']['H'] + 0.5*xi*xi;
    xo = vo*P['Go'];
    h0o = HPo['H'] + 0.5*xo*xo;

    F = np.zeros(2)
    # stagnation enthalpy residual
    F[0] = (h0o-h0i)/h0i;
    # entropy residual
    F[1] = (so-P['si'])/P['si'];

    return np.dot(F,F)

def CompLossResid(X,Params,Ref):
    '''
    /********************************************************************
    Used to solve for the change of state considering the energy and
    momentum conservation equations for compressible flow.
    ********************************************************************/
    '''

    HPo = {'H':0.0,'P':0.0}
    HPo['H'] = X[0];
    HPo['P'] = X[1];
    
    P = Params;

    # outlet specific volume and entropy
    TXPo = HPtoTXP(HPo,Ref);
   
    # specific volume
    # vo = VolumeALL(TXPo,P->G,P->D);           # seperate flow model
    vo = 1/PropertyTXPth('D',TXPo,Ref);         # homogeneous flow model

    # inlet and outlet stagnation enthalpy
    xi = P['vi']*P['G'];
    h0i = P['HPi']['H']+0.5*xi*xi;        # with kinetic energy terms
    xo = vo*P['G'];
    h0o = HPo['H']+0.5*xo*xo;

    # momentum conservation
    t1 = P['G']*P['G']*(P['vi'] + vo)*P['K']/4.0;
    t2 = (HPo['P']-P['HPi']['P']);
    t3 = P['G']*P['G']*(vo - P['vi']);

    # stagnation enthalpy residual
    F = np.zeros(2)
    F[0] = (h0i-h0o-P['q'])/h0i;
    F[1] = (t1+t2+t3)/t1;

    return np.dot(F,F)

def HPLimitsConst(X,Ref):
    '''
    /********************************************************************
    Provides constraints used in solving for the outlet state of the
    manifold.  The numerical solver updates guesses of enthalpy and
    pressure.  This functions states a new guess to ensure that is
    within the property tables.
    ********************************************************************/
    '''
    PMINth = 20*1000 #[Pa]
    PMAXth = 4200*1000 #[Pa]
    TMIN = -50+273.15 #[K]
    TMAX = 160+273.15 #[K]
    
    HPo = {'H':0.0,'P':0.0}
    HPo['H'] = X[0];
    HPo['P'] = X[1];
    
    TXP_prop={'T':0.0,'X':0.0,'P':0.0};
    
    F = np.zeros(4)
    
    F[0] = PMINth - HPo['P']
    F[1] = HPo['P'] - PMAXth

    TXP_prop['P']=HPo['P'];
    TXP_prop['T']=TMIN;
    TXP_prop['X']=0;
    hmin = PropertyTXPth('H',TXP_prop,Ref);
    
    F[2] = hmin - HPo['H'] 
    
    TXP_prop['P']=HPo['P'];
    TXP_prop['T']=TMAX;
    TXP_prop['X']=1;
    hmax = PropertyTXPth('H',TXP_prop,Ref);
    
    F[3] = HPo['H'] - hmax

    return F




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
    
    cp_a= HAPropsSI('cp','T',P['Tai']['T'],'P',101325,'R',0) #[J/kg/K]
    
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
    epsilon_TOT = Q_TOT/(P['ma_TOT']*cp_a)/(Tsat1-P['Tai']['T']);#overall effectivenss
    P['DP_TOT'] = (P['HP_in']['P']-P['HP_out']['P']);#overall pressure drop
    NTU_TOT = -log(1-epsilon_TOT);
    P['UA_TOT'] = NTU_TOT*P['ma_TOT']*cp_a;#overall heat transfer conductance

    return 0



class StructCondClass():
    
    def __init__(self,filename,Ref):
        '''Initialize'''
        #reading the data file in
        #the followings from  Class StructCond(filename) in the ACMODLE
        
        self.Ref = Ref 
           
        #=======================================================================
        # filename = "../InputDoc/CondStruc.xlsx"
        #=======================================================================
        
        df = pd.read_excel(filename,sheetname='description',header = 0)
        df_node = pd.read_excel(filename,sheetname='node',header = 0)
        df_branch = pd.read_excel(filename,sheetname='branch',header = 0)
        df_tube = pd.read_excel(filename,sheetname='tube',header = 0)
           
        self.NodNum = int(df.CondStruc[0]);          #Node Number
        self.BranNum = int(df.CondStruc[1]);         #Branch Number
        self.RowNum = int(df.CondStruc[2]);          #Row Number
        self.TubeNum = int(df.CondStruc[3]);         #overall tube Number
        self.SegNum = int(df.CondStruc[4]);          #Segment Number of per tube
        self.AreaFront = df.CondStruc[5];       #frontal area
        self.Volum = df.CondStruc[6];           #total nominal air mass flow rate [m^3/s]
       
        self.Nod = [CondNode() for k in range(self.NodNum)]
        self.Bra = [CondBranch()for k in range(self.BranNum)]
        self.GaRow = [0.0 for k in range(self.RowNum)]
        self.Tub = [TubeCond() for k in range(self.TubeNum)]
       
        for i in range(self.NodNum):   
            self.Nod[i]['NodNo'] = df_node.NodNo[i]
            self.Nod[i]['InNum'] = df_node.InNum[i]      #get inlet branches
            self.Nod[i]['OutNum'] = df_node.OutNum[i]    #get outlet branches
            self.Nod[i]['BranIN'] = list(df_node.ix[i][3 : 3+self.Nod[i]['InNum']]) #get outlet branches first
            self.Nod[i]['BranOUT'] = list(df_node.ix[i][3+self.Nod[i]['InNum'] : 3+self.Nod[i]['InNum']+self.Nod[i]['OutNum']]) #get inlet branches second
            
        for i in range(self.BranNum):
            self.Bra[i]['BranNo'] = df_branch.BranNo[i]
            self.Bra[i]['EqulNo'] = df_branch.EqulNo[i]
            self.Bra[i]['GrFac'] = df_branch.GrFac[i]
            self.Bra[i]['TubNum'] = df_branch.TubNum[i]
            self.Bra[i]['TubNo'] = list(df_branch.ix[i][4 : 4+self.Bra[i]['TubNum']])#important, the tube number in branch, inputted first from the compressor suction
        
        for i in range(self.TubeNum):
            self.Tub[i]['Seg'] = [TubCndSeg() for k in range(self.SegNum)]
            self.Tub[i]['TubNo'] = df_tube.TubNo[i]
            self.Tub[i]['RowNo'] = df_tube.RowNo[i]
            self.Tub[i]['RefUpstream'] = df_tube.RefUpstream[i]
            self.Tub[i]['AirUpstreamUpper'] = df_tube.AirUpstreamUpper[i]
            self.Tub[i]['AirUpstreamLower'] = df_tube.AirUpstreamLower[i]
            self.Tub[i]['GaFac'] = df_tube.GaFac[i]
            self.Tub[i]['even'] = df_tube.even[i]

    def Update(self):
        """Update the parameters passed in using the dictionary"""
        pass
   
    def DeStructCond(self):
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
    
    def _CondCircuit(self,mr,HPi,Tai,Ga,HPo,Tao,m,P):
        ''' Condenser solver'''
    
        hai={'H':0.0,'W':0.0};#air inlet state #A.B.
        hao={'H':0.0,'W':0.0};#air outlet state #A.B.
        Gr=0;#mass flux
    
        #air side initializing
        hai['H'] = HAPropsSI('H','T',Tai['T'],'P',101325,'R',0) #[J/kg]
        rho_air=1/HAPropsSI('Vda','T',Tai['T'],'P',101325,'R',0) #[kg/m^3 dry air]
        Ma =self.Volum*rho_air; 
        Ga=Ma/self.AreaFront*P['vsp']*P['Ls']/((P['vsp']-P['Do'])*(P['Ls']-P['N']*P['th']));
    
        for i in range(self.RowNum):
            self.GaRow[i]=0;
        
        for j in range(self.RowNum):
            N=0;
            for i in range(self.TubeNum):
                if(self.Tub[i]['RowNo']==j):
                    for k in range(self.SegNum):
                        self.Tub[i]['Seg'][k]['Tai']['T']=Tai['T'];
                        self.Tub[i]['Seg'][k]['hao']['H']=HAPropsSI('H','T',Tai['T'],'P',101325,'R',0) #[J/kg]
        
                    if(self.Tub[i]['RowNo']>0):
                        Upper = self.Tub[i]['AirUpstreamUpper'];
                        Lower = self.Tub[i]['AirUpstreamLower'];
                        if(Upper>=0 and Lower>=0):
                            self.Tub[i]['GaFac']=(self.Tub[Upper]['GaFac']+self.Tub[Lower]['GaFac'])/2;
                        elif(Upper>=0):
                            self.Tub[i]['GaFac']=self.Tub[Upper]['GaFac'];
                        else:
                            self.Tub[i]['GaFac']=self.Tub[Lower]['GaFac'];
                    self.GaRow[j]=self.GaRow[j]+self.Tub[i]['GaFac'];
                    N=N+1;
                #ifend    
            #i circle end
            self.GaRow[j]=self.GaRow[j]/N;
        #j circle end
        
        for i in range(self.TubeNum):
            RowN=self.Tub[i]['RowNo'];
            self.Tub[i]['Ga']=self.Tub[i]['GaFac']/self.GaRow[RowN]*Ga;
        
        #initialize
        H_out=0;
        Res=0;
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
            Gr=mr/P['Ax'];
            
            HPo['H']=HPi['H'];#HPo is the intermediate variable in the calculation
            HPo['P']=HPi['P'];
            
            for i in range(self.BranNum):
                self.Bra[i]['Ini']=0;
        
            for i in range(self.NodNum):#inlet nodes
                if(self.Nod[i]['BranIN'][0]<0):#no inlet branch, only from the hot gas line
                    
                    Gr = CondMan(self.Nod[i]['InNum'],self.Nod[i]['OutNum'],P['Di'],Gr,HPo,self.Ref) #return updated (HPo and Gr)
                    
                    for j in range(self.Nod[i]['OutNum']):#states flowing out from the node
                        jj = int(self.Nod[i]['BranOUT'][j]);#index of the outlet branches
                        self.Bra[jj]['HPi']['H']=HPo['H'];
                        self.Bra[jj]['HPi']['P']=HPo['P'];
                        self.Bra[jj]['Gr']=Gr*self.Bra[jj]['GrFac'];
                        self.Bra[jj]['Ini']=1;
                    #end j circle
                #endif
            #end i circle
            
            self.Cal_HP(P,0,HPo); #Heat transfer and pressure drop calculation #A.B. return updated values for P and HPo
        
            
            for i in range(self.BranNum):
                self.Bra[i]['Ini']=0;
        
            Gr=0;
            HPo['H']=0;
            HPo['P']=0;
        
            #nodes in the middle
            for i in range(self.NodNum):
                if(self.Nod[i]['BranIN'][0]>=0 and self.Nod[i]['BranOUT'][0]>=0):#nodes in the middle
                    for j in range(self.Nod[i]['InNum']):#node inlet state
                        jj=int(self.Nod[i]['BranIN'][j]);
                        Gr=Gr+self.Bra[jj]['Gr'];
                        HPo['H']=self.Bra[jj]['HPo']['H']*self.Bra[jj]['Gr']+HPo['H'];
                        HPo['P']=self.Bra[jj]['HPo']['P']*self.Bra[jj]['Gr']+HPo['P'];
                
                    HPo['H']=HPo['H']/Gr;
                    HPo['P']=HPo['P']/Gr;
                    Gr=Gr/self.Nod[i]['InNum'];
                
                    Gr = CondMan(self.Nod[i]['InNum'],self.Nod[i]['OutNum'],P['Di'],Gr,HPo,self.Ref) #return updated (HPo and Gr)
                    
                    for j in range(self.Nod[i]['OutNum']):
                        jj = int(self.Nod[i]['BranOUT'][j]);#index of outlet branches
                        self.Bra[jj]['HPi']['H']=HPo['H'];
                        self.Bra[jj]['HPi']['P']=HPo['P'];
                        self.Bra[jj]['Gr']=Gr*self.Bra[jj]['GrFac'];
                        self.Bra[jj]['Ini']=1;
                    #end j circle
                #endif
            #end i circle
        
            self.Cal_HP(P,1,HPo) #Heat transfer and pressure drop calculation #A.B. return updated values for P and HPo
        
            Gr=0;
            HPo['H']=0;
            HPo['P']=0;
        
            for i in range(self.BranNum):
                self.Bra[i]['Ini']=0;
            
            #end nodes
            for i in range(self.NodNum):
                if(self.Nod[i]['BranOUT'][0]<0):#no outlet branch except the liquid line
                    for j in range(self.Nod[i]['InNum']):
                        jj= int(self.Nod[i]['BranIN'][j]);
                        Gr=Gr+self.Bra[jj]['Gr'];
                        HPo['H']=self.Bra[jj]['HPo']['H']*self.Bra[jj]['Gr']+HPo['H'];
                        HPo['P']=self.Bra[jj]['HPo']['P']*self.Bra[jj]['Gr']+HPo['P'];
                    
                    HPo['H']=HPo['H']/Gr;
                    HPo['P']=HPo['P']/Gr;
                    Gr=Gr/self.Nod[i]['InNum'];
                
                    Gr = CondMan(self.Nod[i]['InNum'],self.Nod[i]['OutNum'],P['Di'],Gr,HPo,self.Ref) #return updated (HPo and Gr)
                #endif
            #end i circle
            
            if(self.RowNum==1):
                Res=0;
            else:
                Res=2*(HPo['H']-H_out)/(HPo['H']+H_out);
            
            H_out=HPo['H'];
            
            #print the res and iteration no.
            print('Res {}, IterN {}'.format(Res, IterN))
            
            if (abs(Res)<1e-6 or IterN>20): #condition to break the while loop
                break
        
        #end while loop    
        
        if (abs(Res)>1e-5):
            print()
            print('######################################################')
            print('_CondCircuit, Can NOT reach the required tolerance, Res= '+str(Res))
            print('######################################################')
            print()
    
        m['m']=0;
        m['V']=0;
    
        for i in range(self.BranNum):
            m['m']=m['m']+self.Bra[i]['m']['m'];
            m['V']=m['V']+self.Bra[i]['m']['V'];
        
        hao['H']=mr*(HPi['H']-HPo['H'])/Ma+hai['H'];
        #Tao = HatoTa(hao);
        Tao['T'] = HAPropsSI('T','P',101325,'H',hao['H'],'R',0)
        
        #m['m'] = m['m']+(2.86-2.444)+(2.87-2.785);#one point charge tuning for varied charge
        #m->m=m->m+(2.858-2.533)+(2.953-2.8419)+(3.28-3.211)/(13.981-7.667)*(P->LiqL-7.667);#corresponding to three points charge tuning
        #m->m = m->m+(2.86-2.444)+(2.87-2.785)+(3.10-3.045)/(5.863-3.608)*(P->LiqL-3.608);#For two-point charge tuning
        #m->m = m->m+(2.86-2.444); one point charge tuning for fixed charge
        
        return 0
    
    def Cal_HP(self,P,Pos,HPo):
    
        #heat transfer and pressure drop of the two-phase region
        Logic=0;
        hao={'H':0.0,'W':0.0};
        Sm={'m':0.0,'V':0.0};
        mi={'m':0.0,'V':0.0};
    
        Bak = P.copy() #keeping the information of the condenser
        
        if(Pos<0 or Pos>1):
            print("StructCondClass::Cal_HP, Wrong position")

    
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
            
            if (Pos == 0): 
                if (self.Bra[i]['Ini']==1):#this branch has been initialized
                    Logic=1;
                else:
                    Logic=0;  
            elif (Pos == 1):
                if(self.Bra[i]['Ini']==1):#this branch has been initialized
                    Logic=1;
                else:
                    Logic=0;
            
            if(Logic):
                HPo['H']=self.Bra[i]['HPi']['H'];
                HPo['P']=self.Bra[i]['HPi']['P'];
            
                if(self.Bra[i]['EqulNo']<0):#no equivalent branch
                    for j in range(self.Bra[i]['TubNum']):
                        TubeN=self.Bra[i]['TubNo'][j];
                        self.Tub[TubeN]['HPi']['H']=HPo['H'];
                        self.Tub[TubeN]['HPi']['P']=HPo['P'];
                        self.Tub[TubeN]['m']['m']=0;
                        self.Tub[TubeN]['m']['V']=0;
                
                        if(self.Tub[TubeN]['RowNo']>0):#not the first row, to get the air state
                            Upper = self.Tub[TubeN]['AirUpstreamUpper'];
                            Lower = self.Tub[TubeN]['AirUpstreamLower'];
                            for k in range(self.SegNum):
                                if(Upper>=0 and Lower>=0):
                                    hao['H']=(self.Tub[Upper]['Seg'][k]['hao']['H']*self.Tub[Upper]['Ga']+self.Tub[Lower]['Seg'][k]['hao']['H']*self.Tub[Lower]['Ga'])/(self.Tub[Upper]['Ga']+self.Tub[Lower]['Ga']);
                                elif(Upper>=0):
                                    hao['H']=self.Tub[Upper]['Seg'][k]['hao']['H'];
                                else:
                                    hao['H']=self.Tub[Lower]['Seg'][k]['hao']['H'];
        
                                #self.Tub[TubeN]['Seg'][k]['Tai'] = HatoTa(hao);
                                self.Tub[TubeN]['Seg'][k]['Tai']['T'] = HAPropsSI('T','P',101325,'H',hao['H'],'R',0)
                            #end k circle
                        #endif
                
                        for k in range(self.SegNum):
                            #if(k==9 and TubeN==13):
                            #    shenb=0;#for debugging
                            realk=0;
                            if(self.Tub[TubeN]['even']):
                                realk=(self.SegNum-1-k);
                            else:
                                realk=k;
                            
                            CondTubeL_new(self.Bra[i]['Gr'],HPo,self.Tub[TubeN]['Ga'],self.Tub[TubeN]['Seg'][realk]['Tai'],self.Tub[TubeN]['Seg'][realk]['hao'],mi,P,self.Ref) #return updated(self.Tub[TubeN]['Seg'][realk]['hao']['H'], mi, P, HPo)
                            
                            self.Tub[TubeN]['m']['m']=self.Tub[TubeN]['m']['m']+mi['m'];
                            self.Tub[TubeN]['m']['V']=self.Tub[TubeN]['m']['V']+mi['V'];
                            Sm['m']=Sm['m']+mi['m'];
                            Sm['V']=Sm['V']+mi['V'];
                        #end k circle
                
                        CondReturnBend(self.Bra[i]['Gr'],HPo,mi,P,self.Ref) #return updated (HPo, mi, P)
                        
                        self.Tub[TubeN]['m']['m']=self.Tub[TubeN]['m']['m']+mi['m'];
                        self.Tub[TubeN]['m']['V']=self.Tub[TubeN]['m']['V']+mi['V'];
                        Sm['m']=Sm['m']+mi['m'];
                        Sm['V']=Sm['V']+mi['V'];
                        self.Tub[TubeN]['HPo']['H']=HPo['H'];
                        self.Tub[TubeN]['HPo']['P']=HPo['P'];
                    #end j circle
            
                    #output of this branch
                    self.Bra[i]['HPo']['H']=HPo['H'];
                    self.Bra[i]['HPo']['P']=HPo['P'];
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
        #end i circle
    
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

    CondInfo = StructCondClass("InputDoc/CondStruc.xlsx", Ref);
    
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
    MR = 530 # lbm/hr, refrigerant mass flow rate
    GOA = 3.82 # kg/s/m^2, outdoor air mass flux
    DSH = 72 # F, discharge superheat 
    
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
      
    Tsat_dl = Tai['T']+ctoa #Temperature of saturated liquid [K]
    P_dl = PropsSI("P", "T", Tsat_dl, "Q", 0, Ref) #CHECK THIS it might be either liquid or vapor 
    txpi = {'T':Tsat_dl+Tsh,'X':1,'P':P_dl} #initially assume outlet temperature equals to saturated vapor pressure (txpo_P)
    hpi = TXPtoHP(txpi, Ref)  #inlet refrigerant state
    
    #*Run the evaporator model*#
    Condenser(Ref,mr,hpi,Tai,Goa,hpo,Tao,m,Cond_struc,CondPrms,NSeg,condInit=1); #A.B. return updated values for  (hpo, Tao, m, Cond_struc)
              
    txpo = HPtoTXP(hpo, Ref)
    
    Tsat_ll = PropsSI("T", "P", txpo['P'], "Q",0, Ref) #Temperature of saturated liquid [K]
    
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