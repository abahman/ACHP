from __future__ import division, print_function, absolute_import
#from math import pi,log,sqrt,exp,cos,sin,tan,log10
#from scipy.integrate import quad,quadrature,trapz,simps,fixed_quad
#from scipy.optimize import brentq,fsolve
#import numpy as np
#import CoolProp as CP
#import CoolProp
from CoolProp.CoolProp import PropsSI



def toTXP(T,X,P):
    '''
    /********************************************************************
    Convert T,X,P format to TXP format
    ********************************************************************/
    '''
    txp = {'T':0.0,'X':0.0,'P':0.0};
    txp['T']=T;
    txp['X']=X;
    txp['P']=P;
    
    return txp

def HPtoTXP(HP,Ref):
    '''Convert HP format to TXP format.'''

    TXP={'T':0,'X':0,'P':0};

    TXP['P'] = HP['P'];
    hl = PropsSI('H','P',HP['P'],'Q',0,Ref)
    hv = PropsSI('H','P',HP['P'],'Q',1,Ref)
    Tl = PropsSI('T','P',HP['P'],'Q',0,Ref)
    Tv = PropsSI('T','P',HP['P'],'Q',1,Ref)
    
    if (HP['H']>=hv):#superheated
        TXP['X']=1;
        TXP['T'] = PropsSI('T','P',HP['P'],'H',HP['H'],Ref)
    elif (HP['H']<=hl):#subcooled
        TXP['X']=0;
        TXP['T'] = PropsSI('T','P',HP['P'],'H',HP['H'],Ref) 
    else: #two-phase
        TXP['X'] = (HP['H']-hl)/(hv-hl)
        TXP['T'] = Tl + TXP['X']*(Tv-Tl) #use waieghted averaged to determin
    
    return TXP

def TXPtoHP(TXP,Ref):
    '''Convert TXP format to HP format.'''
    
    HP={'H':0.0,'P':0.0};

    HP['P'] = TXP['P'];
    HP['H'] = PropertyTXPth('H',TXP,Ref); #1 means that we wants to find ENTHALPY (i.e., ENTH) ##Now it takes string of 'H'

    return HP

def PropertyTXPth(prop,TXP,Ref):
    '''
    /********************************************************************
    Evaluates properties given TXP data.
    Refrigerant properties
    thermodynamic properties (th) - not transport properties
    pressure and temperature (PT) are independent variables

    prop = sting of 
        'T' Temperature, 0
        'H' Enthalpy, 1
        (No volume 2 >>> it should be 1/Density)
        'S' Entropy, 3
        'U' Internal Energy, 4
        'D' Density, 5
    ********************************************************************/
    '''

    if (TXP['X']>=1): #Superheated
        Tsat = PropsSI('T','P',TXP['P'],'Q',1,Ref)
        if (prop=='T'): #0 is TSAT             
            return Tsat

        if (TXP['T']>Tsat): #superheat vapor
            a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
        else: #saturated vapor
            a = PropsSI(prop,'P',TXP['P'],'Q',1,Ref)
            
    elif (TXP['X']<=0): #Subcooled
        Tsat = PropsSI('T','P',TXP['P'],'Q',0,Ref)
        if(prop=='T'): #0 is TSAT
            return Tsat

        if (TXP['T']<Tsat): #subcooled liquid
            a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
        else: #saturated liquid
            a = PropsSI(prop,'P',TXP['P'],'Q',0,Ref)
    
    else: #two-phase
        if (prop == 'T'): 
            al = PropsSI(prop,'P',TXP['P'],'Q',0,Ref)
            av = PropsSI(prop,'P',TXP['P'],'Q',1,Ref)
            a = al + TXP['X']*(av-al);
        else:
            a = PropsSI(prop,'P',TXP['P'],'Q',TXP['X'],Ref)

    return a

def PropertyTXPtr(prop,TXP,Ref):
    '''
    /********************************************************************
    Evaluates properties given TXP data.
    Refrigerant properties
    transport properties (tr) - not thermodynamic properties
    pressure and temperature (PT) are independent variables
    
    prop = sting of 
        'P' *pressure, 0
        'T' *Temperature, 1
        'L' conductivity, 2
        'V' viscosity, 3
        'C' specific heat, 4
        'I' surface tension, 5
    ********************************************************************/
    '''

    if (TXP['X']>=1): 
        Tsat = PropsSI('T','P',TXP['P'],'Q',1,Ref)

        if (TXP['T']>Tsat): #superheat vapor
            a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
        else: #saturated vapor
            a = PropsSI(prop,'P',TXP['P'],'Q',1,Ref)
    
    elif (TXP.X<=0):
        Tsat = PropsSI('T','P',TXP['P'],'Q',0,Ref)
        
        if (TXP['T']<Tsat): #subcooled liquid
            a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
        else: #saturated liquid
            a = PropsSI(prop,'P',TXP['P'],'Q',0,Ref)
        
    else: #two-phase
        a = PropsSI(prop,'P',TXP['P'],'Q',TXP['X'],Ref)
        
        #(commented) No need to average since PropsSI can find the property with given quality value
        #al = PropsSI(prop,'P',TXP['P'],'Q',0,Ref)
        #av = PropsSI(prop,'P',TXP['P'],'Q',1,Ref)
        #a = al + TXP['X']*(av-al);
           
    return a

if __name__=='__main__':
    print('Hello world')
