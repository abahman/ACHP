from __future__ import division, print_function, absolute_import
from math import log,pi,sqrt,exp,cos,sin,tan,log10
#from scipy.integrate import quad,quadrature,trapz,simps,fixed_quad
from scipy.optimize import brentq

from CoolProp.CoolProp import PropsSI

from extra_functions import PropertyTXPth, PropertyTXPtr, toTXP
from CORR import Eva_FlowPattern


def Xtt(TXPm, Ref):
    '''
    /********************************************************************
    Martinelli parameter. Sqrt of the ratio of the frictional
    pressure drop of liquid to vapor if they were flowing alone
    in a smooth pipe.
    ********************************************************************/
    '''

    TXP_prop={'T':0,'X':0,'P':0};
    
    TXP_prop['P']=TXPm['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    vl = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    mul = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]
    
    TXP_prop['P']=TXPm['P'];
    TXP_prop['X']=1;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    vv  = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    muv = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]

    #TXP_prop['P']=TXPm['P'];
    #TXP_prop['X']=0;
    #TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref)
    #mul = PropertyTXPtr('V',TXP_prop, Ref) #THIS has been moved up to save computational TIME
    
    #TXP_prop['P']=TXPm['P'];
    #TXP_prop['X']=1;
    #TXP_prop['T']=PropertyTXPth('T',TXP_prop,Ref);
    #muv = PropertyTXPtr('V',TXP_prop, Ref); #THIS has been moved up to save computational TIME

    xtt = pow((1-TXPm['X'])/TXPm['X'],0.9)*pow(vl/vv,0.5)*pow(mul/muv,0.1)

    return xtt

def VolumeTPFunc(Alpha,Params=None):
    '''
    /********************************************************************
    Solves for the void fraction (Alpha) in the Hughmark model.
    It takes a guess of Alpha as an argument and returns a
    residual on Alpha that must be zero.
    ********************************************************************/
    '''
    if (Params == None):
        VolParams = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
        P = VolParams
    else:
        P = Params 

    # Z is the "correlating parameter"
    Z1 = pow((P['D']*P['G'])/(P['ul']+Alpha*(P['uv']-P['ul'])),0.16666667);
    Z2 = pow((P['G']*P['x'])/(P['Beta']*(1-P['Beta'])/P['vv']),2);
    Z  = Z1 * pow(Z2/(9.81*P['D']),0.125);

    # K is the flow parameter
    if(Z<=10):
        K = -0.1636+Z*(0.3103+Z*(-0.0352+Z*1.365e-3))
    else:
        K = 0.7257+Z*(6.122e-3+Z*(-6.26e-5+Z*2.396e-7));

    Alpha2 = P['Beta']*K;

    return 2*(Alpha-Alpha2)/(Alpha+Alpha2)


def VolumeALL(TXPi,G,D,q,Ref):
    '''
    /********************************************************************
    Computes the specific volume of refrigerant for single
    and two phase.  It seperates the phase regimes and uses the
    appropriate correlation to get the solution.
        
    TXP is the state of the refrigerant
    returns the specific volume.
    ********************************************************************/
    '''

    X1=0.1;
    X2=0.9;#0.1--original

    if (TXPi['X']>0.999 or TXPi['X']<0.001): #superheated or subcooled
        v = 1/PropertyTXPth('D',TXPi, Ref) #[m^3/kg]
         
    elif (TXPi['X']<X1): #transition region from 0 to 0.1
        TXP1 = toTXP(TXPi['T'],0,TXPi['P']);
        v1 = 1/PropertyTXPth('D',TXP1, Ref) #[m^3/kg]
        TXP2 = toTXP(TXPi['T'],X1,TXPi['P'])
        v2 = VolumeTP(TXP2,G,D,q,Ref)
        v = v1 + TXPi['X']*(v2-v1)/X1;
    
    elif (TXPi['X']>X2): #transition region from 0.9 to 1
        TXP1 = toTXP(TXPi['T'],X2,TXPi['P']);
        v1 = VolumeTP(TXP1,G,D,q, Ref)
        TXP2 = toTXP(TXPi['T'],1,TXPi['P']);
        v2 = 1/PropertyTXPth('D',TXP2, Ref) #[m^3/kg]
        v = v2-(1-TXPi['X'])*(v2-v1)/(1-X2);
    
    else: #two-phase
        v = VolumeTP(TXPi,G,D,q, Ref)
        
    return v


def VolumeTP(TXP,G,D,q,Ref):
    '''
    /********************************************************************
    Specific volume of two phase refrigerant.  Void fraction for two
    phase flow.  This corrlation implements Hugemark model.  ASHRAE
    Transactions ???, pp.309-316.
    
    TXP {temperature(K), quality(-), pressure(Pa)} = thermodynamic state.
    G refrigerant mass flux (kg/m^2/s)
    D internal tube diameter (m).
    returns the specific volume (m^3/kg).
    ********************************************************************/
    '''

    V=0.0;
    try:
        V=VolumeTP_Baroc(TXP,G,D,Ref)
    except:
        print('Baroczy void faild, repeat')
        V=VolumeTP_Baroc(TXP,G,D,Ref)
#     switch(1) {
#             case 1: 
#                 V=VolumeTP_Baroc(TXP,G,D,Ref);
#                 break;
#             case 2:
#                 V=VolumeTP_Zivi(TXP,G,D,Ref);
#                 break;
#             case 3:
#                 V=VolumeTP_Hugh(TXP,G,D,Ref);
#                 break;
#             case 4:
#                 V=VolumeTP_ACRC(TXP,G,D,Ref);
#                 break;
#             case 5:
#                 V=VolumeTP_LM(TXP,G,D,Ref);
#                 break;
#             case 6:
#                 V=VolumeTP_Rigot(TXP,G,D,Ref);
#                 break;
#             case 7:
#                 V=VolumeTP_Smith(TXP,G,D,Ref);
#                 break;
#             case 8:
#                 V=VolumeTP_Tandon(TXP,G,D,Ref);
#                 break;
#             case 9:
#                 V=VolumeTP_Thom(TXP,G,D,Ref);
#                 break;
#             case 10:
#                 V=VolumeTP_Premoli(TXP,G,D,Ref);
#                 break;
#             case 11:
#                 V=VolumeTP_Homo(TXP,G,D,Ref);
#                 break;
#             case 12:
#                 V=VolumeTP_Taitel(TXP,G,D,Ref);
#                 break;
#             case 13:
#                 V=VolumeTP_Rouhani(TXP,G,D,Ref);
#                 break;
#             case 14:
#                 if(fabs(q)<0.00000001): 
#                     {V=VolumeTP_Baroc(TXP,G,D,Ref);
#                      break;}
#                 V=VolumeTP_Thome(TXP,G,D,q,Ref);
#                 break;
#             case 15:
#                 if(fabs(q)<0.00000001): 
#                     {V=VolumeTP_Baroc(TXP,G,D,Ref);
#                     break;}
#                 V=VolumeTP_FlowPattern(TXP,G,D,q,Ref);
#                 break;
#             default: 
#                 V=VolumeTP_Baroc(TXP,G,D,Ref);
#                 break;
#         };

    return V


def VolumeTP_Baroc(TXP1,G,D,Ref):
    '''#B.S.------------------------------------------------
    #Baroczy
    '''
    
    VolParams = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
    P = VolParams #P is a structure of VolParams
    TXP_prop={'T':0,'X':0,'P':0};

    # liquid refrigerant properties
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vl'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    P['ul'] = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]

    #vapor refrigerant properties
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vv'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    P['uv'] = PropertyTXPtr('V',TXP_prop , Ref) #[Pa-s]

    # other parameters that need to be passed into the structure
    P['x']=TXP1['X'];
    P['G']=G; 
    P['D']=D;

    #Baroczy model
    #a is a matrix of 7x7
    a=[[0.49938,-3.9952e-2,1.0212e-2,2.8821e-3,3.6251e-4,2.41e-5,6.6334e-7],
        [-2.5137e-1,-1.5633e-2,7.125e-3,1.7547e-3,1.4822e-4,4.3686e-6,0.0],
        [6.0946e-3,6.5706e-3,1.5014e-4,-6.1361e-5,-1.2632e-6,1.6624e-7,0.0],
        [1.4514e-2,2.4809e-3,-6.9250e-4,-2.2043e-4,-2.1028e-5,-7.3585e-7,0.0],
        [-6.3527e-4,-1.3860e-4,5.3116e-5,2.3740e-5,3.2202e-6,1.5895e-7,0.0],
        [-3.6424e-4,-9.3832e-5,6.0560e-6,2.4383e-6,4.7549e-8,-8.1633e-9,0.0],
        [1.6539e-5,0.0,0.0,0.0,0.0,0.0,0.0]];

    PI2=P['vl']/P['vv']*pow(P['ul']/P['uv'],0.2);
    X_tt=Xtt(TXP1,Ref)
    
    Alpha=0.0;
    for j in range(7):
        for i in range(7):
            Alpha=Alpha+a[i][j]*pow(log(X_tt),i)*pow(log(PI2),j)


    if (Alpha>=1): 
        Alpha=0.9999999999
    if (Alpha<=0): 
        Alpha=0.0000000001
    
    return 1/(Alpha/P['vv']+(1-Alpha)/P['vl'])


def VolumeTP_Zivi(TXP1,G,D,Ref):
    '''
    #B.S.-----------------------------------------------------------
    #Zivi
    '''
    
    VolParams = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
    P = VolParams;

    TXP_prop={'T':0,'X':0,'P':0};

    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vl'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vv'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]

    PI1=P['vl']/P['vv'];
    S=pow(PI1,-0.3333333);
    Alpha=1/(1+(1-TXP1['X'])/TXP1['X']*P['vl']/P['vv']*S);
    
    if (Alpha<=0):
        Alpha=0.000000000001;
    if (Alpha>=1):
        Alpha=0.999999999999;

    return 1/(Alpha/P['vv']+(1-Alpha)/P['vl'])


def VolumeTP_Hugh(TXP1,G,D,Ref):
    '''
    #B.S.----------------------------------------------------
    #Hughmark
    '''

    VolParams = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
    P = VolParams

    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vl'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    P['ul'] = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vv'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    P['uv'] = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]

    #Other parameters that need to be passed into the structure
    P['x']=TXP1['X'];
    P['G']=G;
    P['D']=D;

    #Beta is the prediction of the homogeneous flow model.
    P['Beta'] = (P['x']*P['vv'])/(P['x']*P['vv']+(1-P['x'])*P['vl']);

    # Alpha is the void fraction.
    Alpha = brentq(VolumeTPFunc,0,1,args=(P),xtol=1e-7,rtol=6e-8,maxiter=40) #xtol, rtol and maxiter are changed to match "Zbrent" solver in ACMODEL
    
    if (Alpha<=0):
        Alpha=0.001;
    if (Alpha>=1):
        Alpha=0.999;
    
    return 1/(Alpha/P['vv']+(1-Alpha)/P['vl'])    


def VolumeTP_ACRC(TXP1,G,D,Ref):
    '''
    #B.S.--------------------------------------------------
    #ACRC
    '''

    VolParams = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
    P = VolParams

    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vl'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vv'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]


    j_g=TXP1['X']*G*P['vv'];

    Fr_ACRC=pow((TXP1['X']/(1-TXP1['X'])),0.5)*pow(j_g*j_g/(9.8*D),0.5);

    X_tt=Xtt(TXP1, Ref)

    n=-0.321;

    Alpha=pow((1+1/Fr_ACRC+X_tt),n);

    if (Alpha>=1):
        Alpha=0.9999999999;
    if (Alpha<=0):
        Alpha=0.0000000001;

    return 1/(Alpha/P['vv']+(1-Alpha)/P['vl'])


def VolumeTP_LM(TXP1,G,D,Ref):
    '''
    #B.S.----------------------------------------------------
    #Lockhart-Martinelli
    '''

    VolParams = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
    P = VolParams

    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vl'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vv'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]

    PI1=P['vl']/P['vv']

    X_tt=Xtt(TXP1, Ref)

    if (X_tt<=10):
        Alpha=pow((1+pow(X_tt,0.8)),-0.378)
    else:
        Alpha=0.823-0.157*log(X_tt)

    if (Alpha<=0):
        Alpha=0.000000000001;
    if (Alpha>=1):
        Alpha=0.999999999999;
    
    return 1/(Alpha/P['vv']+(1-Alpha)/P['vl'])


def VolumeTP_Rigot(TXP1,G,D,Ref):
    '''
    #B.S.--------------------------------------------------
    #Rigot
    '''

    VolParams = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
    P = VolParams

    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vl'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vv'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    
    S=2;
    Alpha=1/(1+(1-TXP1['X'])/TXP1['X']*P['vl']/P['vv']*S);

    if (Alpha<=0):
        Alpha=0.000000000001;
    if (Alpha>=1):
        Alpha=0.999999999999;

    return 1/(Alpha/P['vv']+(1-Alpha)/P['vl'])


def VolumeTP_Smith(TXP1,G,D,Ref):
    '''
    #B.S.-------------------------------------------------
    #Smith
    '''
    VolParams = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
    P = VolParams

    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vl'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vv'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]

    
    K=0.4;
    PI1=P['vl']/P['vv'];

    S=K+(1-K)*pow(((1/PI1+K*((1-TXP1['X'])/TXP1['X']))/(1+K*((1-TXP1['X'])/TXP1['X']))),0.5);

    Alpha=1/(1+(1-TXP1['X'])/TXP1['X']*P['vl']/P['vv']*S);
    
    if (Alpha<=0):
        Alpha=0.000000000001;
    if (Alpha>=1):
        Alpha=0.999999999999;

    return 1/(Alpha/P['vv']+(1-Alpha)/P['vl'])


def VolumeTP_Tandon(TXP1,G,D,Ref):
    '''
    #B.S.-----------------------------------------------------
    #Tandon
    '''

    VolParams = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
    P = VolParams

    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vl'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    P['ul'] = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vv'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    

    Re_L=G*D/P['ul'];

    X_tt=Xtt(TXP1, Ref)
    
    F_Xtt=0.15*(1/X_tt+2.85/pow(X_tt,0.476));

    if(Re_L<1125):
        Alpha=1-1.928*pow(Re_L,-0.315)/F_Xtt+0.9293*pow(Re_L,-0.63)/pow(F_Xtt,2);
    else:
        Alpha=1-0.38*pow(Re_L,-0.088)/F_Xtt+0.0361*pow(Re_L,-0.176)/pow(F_Xtt,2);

    if (Alpha<=0):
        Alpha=0.000000000001;
    if (Alpha>=1):
        Alpha=0.999999999999;
    
    return 1/(Alpha/P['vv']+(1-Alpha)/P['vl'])


def VolumeTP_Thom(TXP1,G,D,Ref):
    '''
    #B.S.----------------------------------------------------------
    #Thome
    '''

    VolParams = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
    P = VolParams

    TXP_prop={'T':0,'X':0,'P':0};
    
    #/********************************************
    #Rice&Thom model
    #S=6.9887-481.29*PI2+14797*PI2^2-179022*PI2^3+875059*PI2^4-1.6248E+06*PI2^5+914478*PI2^6    #this is for Thome&Rice
    #*************************************************/
    
    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vl'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    mul = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vv'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    muv = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]    
    
    
    a=[6.9887,-481.29,14797,-179022,875059,-1.6248e06,914478]; #a is an array of 7 elements 

    PI1=P['vl']/P['vv'];
    PI2=PI1*pow(mul/muv,0.2);    
    
    S=0
    for i in range(7):
        S += a[i]*pow(PI2,i)
    
    Alpha=1/(1+(1-TXP1['X'])/TXP1['X']*P['vl']/P['vv']*S);
    
    if (Alpha<=0):
        Alpha=0.000000000001;
    if (Alpha>=1):
        Alpha=0.999999999999;

    return 1/(Alpha/P['vv']+(1-Alpha)/P['vl'])


def VolumeTP_Premoli(TXP1,G,D,Ref):
    '''
    #B.S.----------------------------------------------
    #Premolli
    '''

    VolParams = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
    P = VolParams

    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    sigma = PropertyTXPtr('I',TXP_prop, Ref) #[N/m]
    P['vl'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    P['ul'] = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vv'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    P['uv'] = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]


    PI1=P['vl']/P['vv'];

    We=pow(G,2)*D/(sigma/P['vl']);

    Re_L=G*D/P['ul'];

    F_2=0.0273*We*pow(Re_L,-0.51)*pow(PI1,0.08);

    F_1=1.578*pow(Re_L,-0.19)*pow(PI1,-0.22);

    Y=(TXP1['X']/(1-TXP1['X']))*1/PI1;

    S=1+F_1*pow((Y/(1+F_2*Y)-F_2*Y),0.5);

    Alpha=1.0/(1.0+(1-TXP1['X'])/TXP1['X']*PI1*S);

    if (Alpha>=1): 
        Alpha=0.9999999999;
    if (Alpha<=0):
        Alpha=0.0000000001;

    return 1/(Alpha/P['vv']+(1-Alpha)/P['vl'])


def VolumeTP_Homo(TXP1,G,D,Ref):
    '''
    #B.S.------------------------------------------------
    #Homogeneous
    '''

    VolParams = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
    P = VolParams

    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vl'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vv'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]

    # homogeneous mixture specific volume
    return P['vv']*TXP1['X'] + P['vl']*(1-TXP1['X'])


def VolumeTP_Rouhani(TXP1,G,D,Ref):
    '''
    #B.S.------------------------------------------------
    #Rouhani
    '''

    VolParams = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
    P = VolParams

    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    sigma = PropertyTXPtr('I',TXP_prop, Ref) #[N/m]
    P['vl'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    P['vv'] = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    

    x=TXP1['X'];
    rho_v=1/P['vv'];
    rho_l=1/P['vl'];

    Alpha=x/rho_v*1.0/((1+0.12*(1-x))*(x/rho_v+(1-x)/rho_l)+1.18*(1-x)*pow((9.8*sigma*(rho_l-rho_v)),0.25)/(G*pow(rho_l,0.5)));
    
    if (Alpha>=1):
        Alpha=0.9999999999;
    if (Alpha<=0):
        Alpha=0.0000000001;

    return 1/(Alpha/P['vv']+(1-Alpha)/P['vl'])


def VolumeTP_Taitel(TXP1,G,D,Ref):

    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    vl = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    vv = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]


    PI1=vl/vv;
    C_sf=1.2;
    u_d=0.35*pow((9.8*D),0.5);
    Alpha=1.0/(C_sf+C_sf*((1-TXP1['X'])/TXP1['X'])*PI1+u_d/(TXP1['X']*G*vv));
    
    if (Alpha>=1):
        Alpha=0.9999999999;
    if (Alpha<=0):
        Alpha=0.0000000001;

    return 1/(Alpha/vv+(1-Alpha)/vl)


def VolumeTP_Thome(TXP1,G,D,q,Ref):
    
    FlowPattern = {'JudgPattern':0,'G_wavy':0.0,'G_strat':0.0,'G_mist':0.0,
                   'G_bub':0.0,'X_lA':0.0,'epsilon':0.0,'theta_dry':0.0,
                   'delta':0.0,'h_nb':0.0,'h_cb':0.0,'h_v':0.0,'h_wet':0.0,
                   'h_tp':0.0,'Pattern':0} #output struct variable from Kattan-Thome flow-pattern-dependent heat transfer model
    #initialize
    FlowPat = FlowPattern
    FlowPat['JudgPattern']=1;
    
    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    vl = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    vv = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    

    if (q<0):
        FlowPat = Cond_FlowPattern(TXP1,G,D,q,FlowPat, Ref);
    else:
        FlowPat = Eva_FlowPattern(TXP1,G,D,q,FlowPat, Ref);


    Alpha=FlowPat['epsilon'];

    if (Alpha>=1):
        Alpha=0.9999999999;
    if (Alpha<=0):
        Alpha=0.0000000001;

    return 1/(Alpha/vv+(1-Alpha)/vl)


def VolumeTP_FlowPattern(TXP1,G,D,q,Ref):
    
    FlowPattern = {'JudgPattern':0,'G_wavy':0.0,'G_strat':0.0,'G_mist':0.0,
                   'G_bub':0.0,'X_lA':0.0,'epsilon':0.0,'theta_dry':0.0,
                   'delta':0.0,'h_nb':0.0,'h_cb':0.0,'h_v':0.0,'h_wet':0.0,
                   'h_tp':0.0,'Pattern':0} #output struct variable from Kattan-Thome flow-pattern-dependent heat transfer model
    #initialize
    FlowPat = FlowPattern
    FlowPat['JudgPattern']=1;
    
    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    vl = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]
    
    #vapor specific volumes
    TXP_prop['P']=TXP1['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    vv = 1/PropertyTXPth('D',TXP_prop, Ref) #[m^3/kg]


    if (q<0):
        FlowPat = Cond_FlowPattern(TXP1,G,D,q,FlowPat, Ref);
    else:
        FlowPat = Eva_FlowPattern(TXP1,G,D,q,FlowPat, Ref);


    if (FlowPat['Pattern']==2):
        #V=VolumeTP_ACRC(TXP1,G,D);
        V = VolumeTP_Rouhani(TXP1,G,D, Ref);
    else:
        V = VolumeTP_Taitel(TXP1,G,D, Ref);
    

    return V

 
if __name__=='__main__':
    print('Hello world')
