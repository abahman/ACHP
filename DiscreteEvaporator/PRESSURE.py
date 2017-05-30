from __future__ import division, print_function, absolute_import
from math import log,pi,sqrt,exp,cos,sin,tan,log10
#from scipy.integrate import quad,quadrature,trapz,simps,fixed_quad
#from scipy.optimize import brentq,fsolve
#import numpy as np

#from CoolProp.CoolProp import PropsSI

from extra_functions import toTXP, PropertyTXPtr, PropertyTXPth, HPtoTXP, PreAcc
from VOLUME import VolumeALL
from CORR import Xtt

def dPdLfriction(TXPi,G,D,Ref):
    '''
    /********************************************************************
    Calculates friction pressure drop (Pa) for internal refrigerant flow.
    ********************************************************************/
    '''

    X1=0.1; X2=0.9;
    #X1=0.05; X2=0.95;       
    #X1=0.01; X2=0.99;
    
    if (TXPi['X']>=1.0 or TXPi['X']<=0):
        y=dPdLfricSP(TXPi,G,D, Ref) 
    
    elif (TXPi['X']<X1):
        TXP1=toTXP(TXPi['T'],0,TXPi['P']);
        TXP2=toTXP(TXPi['T'],X1,TXPi['P']);
        y1=dPdLfricSP(TXP1,G,D, Ref)
        y2=dPdLfricTP(TXP2,G,D, Ref)
        y=y1+TXPi['X']*(y2-y1)/X1;
    
    elif(TXPi['X']>X2):
        TXP1=toTXP(TXPi['T'],X2,TXPi['P']);
        TXP2=toTXP(TXPi['T'],1,TXPi['P']);
        y1=dPdLfricTP(TXP1,G,D, Ref)
        y2=dPdLfricSP(TXP2,G,D, Ref)
        y=y2-(1-TXPi['X'])*(y2-y1)/(1-X2);
    
    else:
        y=dPdLfricTP(TXPi,G,D, Ref)

    return y


def dPdLfricSP(TXP,G,D,Ref):
    '''
    /********************************************************************
    Calculates friction pressure drop (Pa) for single phase internal flow.
    ********************************************************************/
    '''

    mu = PropertyTXPtr('V',TXP, Ref) #[Pa-s]

    v = 1/PropertyTXPth('D',TXP, Ref) #[m^3/kg]

    Re = G*D/mu;
    
    if (Re<0.0): 
        print("dPdLfricSP","Re<0.0")

    f = 0.046*pow(Re,-0.2)
    
    if (G<0.0):
        print("dPdLfricSP","G<0.0")

    dPdL = -0.5*f/D*v*pow(G,2);
    
    return dPdL


def dPdLfricTP(TXP,G,D,Ref):
    '''
    /********************************************************************
    Calculates friction pressure drop for two-phase internal flow.
    ********************************************************************/
    '''   
     
    dPdLf0 = dPdLfricSP(toTXP(TXP['T'],0,TXP['P']),G,D, Ref);

    xtt = Xtt(TXP, Ref)

    arg = 1.0-TXP['X'];
    if (arg<0.0): 
        print("dPdLfricTP","1.0-TXP['X']<0.0")

    PhiSqf0 = pow(arg,1.8)*PhiSqfWALLIS(xtt)
    #PhiSqf0 = pow(arg,1.8)*PhiSqfHPSIM(xtt)

    dPdL = PhiSqf0*dPdLf0;

    return dPdL


def PhiSqfWALLIS(xtt):
    '''
    /********************************************************************
    Correlation factor for two phase frictional pressure drop.
    ********************************************************************/
    '''
    
    n=3.5;
    Phi = pow(1+pow(xtt,-2/n),n);

    return Phi


def dPmom(v1,v2,G):
    '''
    /********************************************************************
    Calculates the pressure drop (Pa) due to a change in specific volume
    usually resulting from a heat transfer process.  Assumptions include
    no body forces, no viscous forces, steady flow, and constant area.
    ********************************************************************/
    '''
    return pow(G,2)*(v2-v1)


def dPelbow(TXP,G,D,Ro,Ref):
    '''
    /********************************************************************
    Calculates the pressure drop for 180 return bend.
    Assumptions: (1) smooth wall, (2) incompressible flow.
    Reference:
        Idelchik, I.E.
        Handbook of Hydraulic Resistance, 2nd ed.
        Hemisphere Publishing Company, 1986.
        p.289
    ********************************************************************/
    '''

    delta=180

    dPdL = dPdLfriction(TXP,G,D,Ref)

    v = VolumeALL(TXP,G,D,0, Ref)

    A1 = 0.7+0.35*delta/90;
    B1 = 0.21/sqrt(Ro/D);
    C1 = 1;
    dP = 0.0175*Ro*delta*dPdL-(0.5)*pow(G,2)*v*A1*B1*C1;

    return dP


def dPdLfrictionEvap(TXPi,G,D,dx,L,Ref):
    '''
    /********************************************************************
    Calculates frictional pressure drop (Pa) for internal refrigerant flow
    in the evaporator.
    ********************************************************************/
    '''
    
    X1=0.1; X2=0.9;
    #X1=0.05; X2=0.95;       
    #X1=0.01; X2=0.99;;

    if (TXPi['X']>=1.0 or TXPi['X']<=0):
        y = dPdLfricSP(TXPi,G,D, Ref)
         
    elif (TXPi['X']<X1):
        TXP1=toTXP(TXPi['T'],0,TXPi['P']);
        TXP2=toTXP(TXPi['T'],X1,TXPi['P']);
        y1 = dPdLfricSP(TXP1,G,D, Ref)
        y2 = dPdLfricEvapTP(TXP2,G,D,dx,L, Ref)
        y=y1+TXPi['X']*(y2-y1)/X1;
     
    elif(TXPi['X']>X2):
        TXP1=toTXP(TXPi['T'],X2,TXPi['P']);
        TXP2=toTXP(TXPi['T'],1,TXPi['P']);
        y1 = dPdLfricEvapTP(TXP1,G,D,dx,L, Ref)
        y2 = dPdLfricSP(TXP2,G,D, Ref)
        y = y2-(1-TXPi['X'])*(y2-y1)/(1-X2);
    
    else:
        y = dPdLfricEvapTP(TXPi,G,D,dx,L, Ref)
    
    return y


def dPdLfricEvapTP(TXPm,G,D,dx,L,Ref):
    '''
    /********************************************************************
    Calculates friction pressure drop (Pa) for evaporation two phase
    internal flow.
    ********************************************************************/
    '''

    g=9.807;
    dx0=1e-2;
    TXP_prop={'T':0.0,'X':0.0,'P':0.0}

    TXP_prop['P']=TXPm['P'];
    TXP_prop['X']=0;
    Tsat = PropertyTXPth('T',TXP_prop, Ref) #[K]


    TXP_prop['P']=TXPm['P'];
    TXP_prop['T']=Tsat;
    TXP_prop['X']=0;
    mul = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]
    hl = PropertyTXPth('H',TXP_prop, Ref) #[J/kg]
    
    v = 1/PropertyTXPth('D',TXPm, Ref) #[m^3/kg]
    
    
    TXP_prop['P']=TXPm['P'];
    TXP_prop['X']=1;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    hv = PropertyTXPth('H',TXP_prop, Ref) #[J/kg]

    hlv = hv-hl;
    Re = G*D/mul;
    f0 = 0.0185*pow(hlv/(Re*L*g),0.25);

    if(dx>=dx0): 
        f = f0*pow(dx,0.25); 
    else:
        A = f0*pow(dx0,0.25)/exp(0.25);
        B = 1/(4*dx0);
        f = A*exp(B*dx);

    dPdL = -f*pow(G,2)*v/D;

    return dPdL


def FricFactorSP(TXP,G,D,Ref):
    '''
    /************************************************************************
    Calculates friction factor for single-phase internal refrigerant flow.
    ************************************************************************/
    '''

    # Option 1: recommended constant value
    #
    # return 0.02;

    # Option 2: smooth pipe
    #
    # char s[SL];
    # double mu = PropertyTXP(VISC,TXP,status);
    # if(*status) {sprintf(s,"FricFactorSP:\n"); bprintf(s); return -1;}
    # double Re=G*D/mu;
    # return pow(1.8*log10(Re)-1.64,-2);

    # Option 3: rough pipe
    #
    delta=0.01;
    return pow(2*log10(3.7/delta),-2)


def LossCoeffReturnBend(TXP,G,R,D,Ref):
    '''
    /***********************************************************************
    This function provides the local and frictional losses of a return bend.
    Reference:
        Idelchik, I.E.
        Handbook of Hydraulic Resistance, 2nd ed.
        Hemisphere Publishing Company, 1986.
        p.289
    Inputs:
        TXP = refrigerant state
        G = refrigerant mass flux (kg/s/m^2)
        R = radius of return bend (m)
        D = inside diameter of pipe (m)
    Outputs:
        (return value) return bend loss coefficient (K)
    ************************************************************************/
    '''

    delta=180;  # angle of return bend in degrees

    # local loss coefficient
    A1 = 0.7+0.35*delta/90;      # delta > 100 degrees
    B1 = 0.21/sqrt(R/D);          # R/D > 1.0
    C1 = 1.0;                   # circular or square x-section
    Kloc = A1*B1*C1;

    f = FricFactorSP(TXP,G,D,Ref) #(TXP, G, D and Ref) are passed, but actually do nothing
    
    Kfric = 0.0175*R/D*delta*f

    return Kfric+Kloc


def LossCoeffStraightPipe(TXP,G,L,D,Ref):
    '''
    /***********************************************************************
    This function provides the local and frictional losses of a straight
    piece of pipe.
    Inputs:
        TXP = refrigerant state
        G = refrigerant mass flux (kg/s/m^2)
        L = length of pipe (m)
        D = inside diameter of pipe (m)
    Outputs:
        (return value) return bend loss coefficient (K)
    ************************************************************************/
    '''
    
    f = FricFactorSP(TXP,G,D,Ref) #(TXP, G, D and Ref) are passed, but actually do nothing

    K = f*L/D

    return K


def FricFactorTP(TXP,G,D,Ref):
    ''' 
    /************************************************************************
    Calculates friction factor for two-phase internal refrigerant flow.
    ************************************************************************/
    '''

    fSP = FricFactorSP(toTXP(TXP['T'],0,TXP['P']),G,D,Ref)

    xtt = Xtt(TXP, Ref)

    PhiSqf0 = pow(1-TXP['X'],1.8)*PhiSqfWALLIS(xtt)
    #double PhiSqf0 = pow(1-TXP.X,1.8)*PhiSqfHPSIM(xtt);

    fTP = PhiSqf0*fSP;

    return fTP


def FricFactor(TXPi,G,D,Ref):
    '''
    /************************************************************************
    Calculates friction factor for internal refrigerant flow.
    ************************************************************************/
    '''
    
    X1=0.1; X2=0.9;
    #X1=0.05; X2=0.95;       
    #X1=0.01; X2=0.99;

    if (TXPi['X']>=1.0 or TXPi['X']<=0):
        y = FricFactorSP(TXPi,G,D, Ref)
        
    elif (TXPi['X']<X1):
        TXP1 = toTXP(TXPi['T'],0,TXPi['P']);
        TXP2 = toTXP(TXPi['T'],X1,TXPi['P']);
        y1 = FricFactorSP(TXP1,G,D, Ref)
        y2 = FricFactorTP(TXP2,G,D, Ref)
        y = y1+TXPi['X']*(y2-y1)/X1;

    elif (TXPi['X']>X2):
        TXP1=toTXP(TXPi['T'],X2,TXPi['P']);
        TXP2=toTXP(TXPi['T'],1,TXPi['P']);
        y1=FricFactorTP(TXP1,G,D, Ref)
        y2=FricFactorSP(TXP2,G,D, Ref)
        y = y2-(1-TXPi['X'])*(y2-y1)/(1-X2);

    else:
        y = FricFactorTP(TXPi,G,D, Ref)

    return y


def GET_PreAcc(DP_ACC, Ref, Params=None):
    '''
    /*******************************************************
    B.S. add for getting the acceleration pressure drop (Pa)
    ********************************************************/
    '''
    if (Params==None):
        Preacc = PreAcc()
    else:
        Preacc = Params

    DP_FR=Preacc['DP_FR'];
    G=Preacc['G'];
    H_OUT=Preacc['H_OUT'];
    P=Preacc['P_IN'];
    X=Preacc['X_IN'];
    TXP_prop={'T':0,'X':0,'P':0};

    TXP_prop['P']=P;
    TXP_prop['X']=1;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop,Ref)
    DV=PropertyTXPth('D',TXP_prop, Ref) #[kg/m^3]
    
    TXP_prop['P']=P;
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop,Ref)
    DL=PropertyTXPth('D',TXP_prop,Ref) #[kg/m^3]

    DP_OLD=DP_ACC;
    P_OUT=P-DP_FR-DP_ACC;

    HP_OUT={'H':0.0,'P':0.0};
    HP_OUT['H']=H_OUT;
    HP_OUT['P']=P_OUT;
    TXP_OUT=HPtoTXP(HP_OUT,Ref);

    X_OUT=TXP_OUT['X'];
    DP_ACC=pow(G,2)*(1/DV-1/DL)*(X_OUT-X);
    ERR_P=(DP_ACC-DP_OLD);

    return ERR_P

 
if __name__=='__main__':
    print('Hello world')
