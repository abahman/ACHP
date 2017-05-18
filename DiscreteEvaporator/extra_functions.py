from __future__ import division, print_function, absolute_import
#from math import pi,log,sqrt,exp,cos,sin,tan,log10
#from scipy.integrate import quad,quadrature,trapz,simps,fixed_quad
#from scipy.optimize import brentq,fsolve
#import numpy as np
#import CoolProp as CP
#import CoolProp
from CoolProp.CoolProp import PropsSI


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

        if (TXP['T']>Tsat):
            a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
        else: 
            a = PropsSI(prop,'P',TXP['P'],'Q',1,Ref)
            
    elif (TXP['X']<=0): #Subcooled
        Tsat = PropsSI('T','P',TXP['P'],'Q',0,Ref)
        if(prop=='T'): #0 is TSAT
            return Tsat

        if (TXP['T']<Tsat):
            a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
        else:
            a = PropsSI(prop,'P',TXP['P'],'Q',0,Ref)
    
    else: #two-phase
        if (prop == 'T'): 
            al = PropsSI(prop,'P',TXP['P'],'Q',0,Ref)
            av = PropsSI(prop,'P',TXP['P'],'Q',1,Ref)
            a = al + TXP['X']*(av-al);
        else:
            a = PropsSI(prop,'P',TXP['P'],'Q',TXP['X'],Ref)

    return a
 
if __name__=='__main__':
    print('Hello world')
