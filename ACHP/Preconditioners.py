from __future__ import absolute_import
from scipy.optimize import fsolve, minimize

from math import pi
import numpy as np

import CoolProp as CP
from .convert_units import *
from CoolProp.CoolProp import HAPropsSI 

from .Solvers import MultiDimNewtRaph
from . import Correlations

def DXPreconditioner(Cycle,epsilon=0.96):
    
    #Assume the heat exchangers are highly effective
    #Condensing heat transfer rate from enthalpies
    rho_air=1.1
    Cp_air =1005 #[J/kg/K]
        
    #AbstractState
    if hasattr(Cycle,'Backend'): #check if backend is given
        AS = CP.AbstractState(Cycle.Backend, Cycle.Ref)
    else: #otherwise, use the defualt backend
        AS = CP.AbstractState('HEOS', Cycle.Ref)
        Cycle.Backend = 'HEOS'
    Cycle.AS = AS
    
    def OBJECTIVE(x):
        Tevap=x[0]
        Tcond=x[1]
        
        #Use fixed effectiveness to get a guess for the condenser capacity
        Qcond=epsilon*Cycle.Condenser.Fins.Air.Vdot_ha*rho_air*(Cycle.Condenser.Fins.Air.Tdb-Tcond)*Cp_air
        
        Cycle.AS.update(CP.QT_INPUTS,1.0,Tevap)
        pevap=Cycle.AS.p() #[pa]
        Cycle.AS.update(CP.QT_INPUTS,0.0,Tcond)
        P_l =Cycle.AS.p() #[pa]
        Cycle.AS.update(CP.QT_INPUTS,1.0,Tcond)
        P_v=Cycle.AS.p() #[pa]
        pcond= (P_l+P_v)/2.0
        Cycle.Compressor.pin_r=pevap
        Cycle.Compressor.pout_r=pcond
        Cycle.Compressor.Tin_r=Tevap+Cycle.Evaporator.DT_sh
        Cycle.Compressor.Ref=Cycle.Ref
        Cycle.Compressor.Backend=Cycle.Backend
        Cycle.Compressor.Calculate()
        W=Cycle.Compressor.W
        
        # Evaporator fully-dry analysis
        Qevap_dry=epsilon*Cycle.Evaporator.Fins.Air.Vdot_ha*rho_air*(Cycle.Evaporator.Fins.Air.Tdb-Tevap)*Cp_air
        
        #Air-side heat transfer UA
        Evap=Cycle.Evaporator
        Evap.mdot_r=Cycle.Compressor.mdot_r
        Evap.psat_r=pevap
        Evap.Ref=Cycle.Ref
        Evap.Backend=Cycle.Backend
        Evap.Initialize()
        UA_a=Evap.Fins.h_a*Evap.Fins.A_a*Evap.Fins.eta_a
        Tin_a=Evap.Fins.Air.Tdb
        Tout_a=Tin_a+Qevap_dry/(Evap.Fins.mdot_da*Evap.Fins.cp_da)
        #Refrigerant-side heat transfer UA
        UA_r=Evap.A_r_wetted*Correlations.ShahEvaporation_Average(0.5,0.5,Cycle.AS,Evap.G_r,Evap.ID,Evap.psat_r,Qevap_dry/Evap.A_r_wetted,Evap.Tbubble_r,Evap.Tdew_r)
        #Get wall temperatures at inlet and outlet from energy balance
        T_so_a=(UA_a*Evap.Tin_a+UA_r*Tevap)/(UA_a+UA_r)
        T_so_b=(UA_a*Tout_a+UA_r*Tevap)/(UA_a+UA_r)
        
        Tdewpoint=HAPropsSI('D','T',Cycle.Evaporator.Fins.Air.Tdb, 'P',101325, 'R',Evap.Fins.Air.RH)
        
        #Now calculate the fully-wet analysis
        #Evaporator is bounded by saturated air at the refrigerant temperature.
        h_ai=HAPropsSI('H','T',Cycle.Evaporator.Fins.Air.Tdb, 'P',101325, 'R', Cycle.Evaporator.Fins.Air.RH) #[J/kg_da]
        h_s_w_o=HAPropsSI('H','T',Tevap, 'P',101325, 'R', 1.0) #[J/kg_da]
        Qevap_wet=epsilon*Cycle.Evaporator.Fins.Air.Vdot_ha*rho_air*(h_ai-h_s_w_o)
        
        #Coil is either fully-wet, fully-dry or partially wet, partially dry
        if T_so_a>Tdewpoint and T_so_b>Tdewpoint:
            #Fully dry, use dry Q
            f_dry=1.0
        elif T_so_a<Tdewpoint and T_so_b<Tdewpoint:
            #Fully wet, use wet Q
            f_dry=0.0
        else:
            f_dry=1-(Tdewpoint-T_so_a)/(T_so_b-T_so_a)
        Qevap=f_dry*Qevap_dry+(1-f_dry)*Qevap_wet
        
        if Cycle.ImposedVariable == 'Subcooling': #if Subcooling impose
            Cycle.AS.update(CP.PT_INPUTS,pcond,Tcond-Cycle.DT_sc_target)
            h_target = Cycle.AS.hmass() #[J/kg]
            Qcond_enthalpy=Cycle.Compressor.mdot_r*(Cycle.Compressor.hout_r-h_target)
        else: #otherwise, if Charge impose
            Cycle.AS.update(CP.PT_INPUTS,pcond,Tcond-5)
            h_target = Cycle.AS.hmass() #[J/kg]
            Qcond_enthalpy=Cycle.Compressor.mdot_r*(Cycle.Compressor.hout_r-h_target)
        
        resids=[Qevap+W+Qcond,Qcond+Qcond_enthalpy]#,Qevap,f_dry]
        return resids
    
    Tevap_init=Cycle.Evaporator.Fins.Air.Tdb-15
    Tcond_init=Cycle.Condenser.Fins.Air.Tdb+8
    x=fsolve(OBJECTIVE,[Tevap_init,Tcond_init])
    DT_evap=Cycle.Evaporator.Fins.Air.Tdb-x[0]
    DT_cond=x[1]-Cycle.Condenser.Fins.Air.Tdb
    
    return DT_evap-4, DT_cond+5#sometimes DT_evap-2, DT_cond+2 can work with no error


def SecondaryLoopPreconditioner(Cycle,epsilon=0.9):
    rho_air= 1.1 #[kg/m^2]
    Cp_air = 1005 #[J/kg-K]
    
    def OBJECTIVE(x):
        Tevap=x[0]
        Tcond=x[1]
        Tin_CC=x[2]
        if Cycle.Mode=='AC':
            #Condenser heat transfer rate
            Qcond=epsilon*Cycle.Condenser.Fins.Air.Vdot_ha*rho_air*(Cycle.Condenser.Fins.Air.Tdb-Tcond)*Cp_air
            
            #Compressor power
            Cycle.AS.update(CP.QT_INPUTS,1.0,Tevap)
            pevap=Cycle.AS.p() #[Pa]
            Cycle.AS.update(CP.QT_INPUTS,1.0,Tcond)
            pcond=Cycle.AS.p() #[Pa]
            Cycle.Compressor.pin_r=pevap
            Cycle.Compressor.pout_r=pcond
            Cycle.Compressor.Tin_r=Tevap+Cycle.Compressor.DT_sh
            Cycle.Compressor.Ref=Cycle.Ref
            Cycle.Compressor.Calculate()
            W=Cycle.Compressor.W
            
            Qcoolingcoil_dry=epsilon*Cycle.CoolingCoil.Fins.Air.Vdot_ha*rho_air*(Cycle.CoolingCoil.Fins.Air.Tdb-Tin_CC)*Cp_air
            
            # Air-side heat transfer UA
            CC=Cycle.CoolingCoil
            CC.Initialize()
            UA_a=CC.Fins.h_a*CC.Fins.A_a*CC.Fins.eta_a
            #Air outlet temp from dry analysis
            Tout_a=CC.Tin_a-Qcoolingcoil_dry/(CC.Fins.mdot_a*CC.Fins.cp_a)
            
            # Refrigerant side UA
            f,h,Re=Correlations.f_h_1phase_Tube(Cycle.Pump.mdot_g/CC.Ncircuits, CC.ID, Tin_CC, CC.pin_g, CC.AS_g)
            UA_r=CC.A_g_wetted*h
            #Glycol specific heat
            Cycle.Pump.AS_g.update(CP.PT_INPUTS,Cycle.Pump.pin_g,Tin_CC)
            cp_g=Cycle.Pump.AS_g.cpmass() #[J/kg-K]
            #Refrigerant outlet temp
            Tout_CC=Tin_CC+Qcoolingcoil_dry/(Cycle.Pump.mdot_g*cp_g)
            
            #Get wall temperatures at inlet and outlet from energy balance
            T_so_a=(UA_a*CC.Tin_a+UA_r*Tout_CC)/(UA_a+UA_r)
            T_so_b=(UA_a*Tout_a+UA_r*Tin_CC)/(UA_a+UA_r)
            
            Tdewpoint=HAPropsSI('D','T',CC.Fins.Air.Tdb,'P',101325,'R',CC.Fins.Air.RH)
            #Now calculate the fully-wet analysis
            #Evaporator is bounded by saturated air at the refrigerant temperature.
            h_ai= HAPropsSI('H','T',CC.Fins.Air.Tdb,'P',101325,'R',CC.Fins.Air.RH)
            h_s_w_o=HAPropsSI('H','T',Tin_CC,'P',101325,'R',1.0)
            Qcoolingcoil_wet=epsilon*CC.Fins.Air.Vdot_ha*rho_air*(h_ai-h_s_w_o)
            
            #Coil is either fully-wet, fully-dry or partially wet, partially dry
            if T_so_a>Tdewpoint and T_so_b>Tdewpoint:
                #Fully dry, use dry Q
                f_dry=1.0
            elif T_so_a<Tdewpoint and T_so_b<Tdewpoint:
                #Fully wet, use wet Q
                f_dry=0.0
            else:
                f_dry=1-(Tdewpoint-T_so_a)/(T_so_b-T_so_a)
            Qcoolingcoil=f_dry*Qcoolingcoil_dry+(1-f_dry)*Qcoolingcoil_wet
        
            Tin_IHX=Tin_CC+Qcoolingcoil/(Cycle.Pump.mdot_g*cp_g)
            QIHX=epsilon*Cycle.Pump.mdot_g*cp_g*(Tin_IHX-Tevap)
            
            Cycle.AS.update(CP.PT_INPUTS,pcond,Tcond-Cycle.DT_sc_target)
            h_target = Cycle.AS.hmass() #[J/kg]
            Qcond_enthalpy=Cycle.Compressor.mdot_r*(Cycle.Compressor.hout_r-h_target)
            resids=[QIHX+W+Qcond,Qcond+Qcond_enthalpy,Qcoolingcoil-QIHX]
            return resids
        
        elif Cycle.Mode=='HP':
            #Evaporator heat transfer rate
            Qevap=epsilon*Cycle.Evaporator.Fins.Air.Vdot_ha*rho_air*(Cycle.Evaporator.Fins.Air.Tdb-Tevap)*Cp_air
            
            #Compressor power
            Cycle.AS.update(CP.QT_INPUTS,1.0,Tevap)
            pevap=Cycle.AS.p() #[pa]
            Cycle.AS.update(CP.QT_INPUTS,1.0,Tcond)
            pcond=Cycle.AS.p() #[pa]
            Cycle.Compressor.pin_r=pevap
            Cycle.Compressor.pout_r=pcond
            Cycle.Compressor.Tin_r=Tevap+Cycle.Evaporator.DT_sh
            Cycle.Compressor.Ref=Cycle.Ref
            Cycle.Compressor.Calculate()
            W=Cycle.Compressor.W
            
            #Evaporator will be dry
            Qcoolingcoil=epsilon*Cycle.CoolingCoil.Fins.Air.Vdot_ha*rho_air*(Tin_CC-Cycle.CoolingCoil.Fins.Air.Tdb)*Cp_air
            
            #Glycol specifi heat
            Cycle.Pump.AS_g.update(CP.PT_INPUTS,Cycle.Pump.pin_g,Tin_CC)
            cp_g=Cycle.Pump.AS_g.cpmass() #[J/kg/K]
            
            Tin_IHX=Tin_CC-Qcoolingcoil/(Cycle.Pump.mdot_g*cp_g)
            QIHX=epsilon*Cycle.Pump.mdot_g*cp_g*(Tin_IHX-Tcond)
            
            Cycle.AS.update(CP.PT_INPUTS,pcond,Tcond-Cycle.DT_sc_target)
            h_target = Cycle.AS.hmass() #[J/kg]
            QIHX_enthalpy=Cycle.Compressor.mdot_r*(Cycle.Compressor.hout_r-h_target)
            
            resids=[QIHX+W+Qevap,QIHX+QIHX_enthalpy,Qcoolingcoil+QIHX]
            return resids
          
    solverFunc=fsolve
    if Cycle.Mode=='AC':
        Tevap_init=Cycle.CoolingCoil.Fins.Air.Tdb-15
        Tcond_init=Cycle.Condenser.Fins.Air.Tdb+8
        Tin_CC=Tevap_init+1
        #First try using the fsolve algorithm
        try:
            x=fsolve(OBJECTIVE,[Tevap_init,Tcond_init,Tin_CC])
        except:
            #If that doesnt work, try the Mult-Dimensional Newton-raphson solver
            try:
                x=MultiDimNewtRaph(OBJECTIVE,[Tevap_init,Tcond_init,Tin_CC])
            except:
                x=[Tevap_init,Tcond_init,284]
        DT_evap=Cycle.CoolingCoil.Fins.Air.Tdb-x[0]
        DT_cond=x[1]-Cycle.Condenser.Fins.Air.Tdb
        Tin_CC=x[2]
    elif Cycle.Mode=='HP':
        Tevap_init=Cycle.Evaporator.Fins.Air.Tdb-8
        Tcond_init=Cycle.CoolingCoil.Fins.Air.Tdb+15
        Tin_CC=Tcond_init-1
        x=solverFunc(OBJECTIVE,[Tevap_init,Tcond_init,Tin_CC])
        DT_evap=Cycle.Evaporator.Fins.Air.Tdb-x[0]
        Tin_CC=x[2]
        DT_cond=x[1]-Tin_CC
    else:
        raise ValueError()
        
    return DT_evap-2,DT_cond+2,Tin_CC


def VICompPreconditioner(Cycle,epsilon=0.96):
    
    #Assume the heat exchangers are highly effective
    #Condensing heat transfer rate from enthalpies
    rho_air=1.1
    Cp_air =1005 #[J/kg/K]
        
    #AbstractState
    if hasattr(Cycle,'Backend'): #check if backend is given
        AS = CP.AbstractState(Cycle.Backend, Cycle.Ref)
    else: #otherwise, use the defualt backend
        AS = CP.AbstractState('HEOS', Cycle.Ref)
        Cycle.Backend = 'HEOS'
    Cycle.AS = AS
    
    def OBJECTIVE(x):
        Tevap=x[0]
        Tcond=x[1]
        Tdew_inj=x[2]
        
        #Use fixed effectiveness to get a guess for the condenser capacity
        Qcond=epsilon*Cycle.Condenser.Fins.Air.Vdot_ha*rho_air*(Cycle.Condenser.Fins.Air.Tdb-Tcond)*Cp_air
        
        Cycle.AS.update(CP.QT_INPUTS,1.0,Tevap)
        pevap=Cycle.AS.p() #[pa]
        Cycle.AS.update(CP.QT_INPUTS,1.0,Tcond)
        pcond=Cycle.AS.p() #[pa]
        Cycle.AS.update(CP.QT_INPUTS,1.0,Tdew_inj)
        pinj=Cycle.AS.p() #[pa]
        cp_dew=Cycle.AS.cpmass() #[J/kg-K]
        Cycle.AS.update(CP.PQ_INPUTS,pinj,0.0)
        Tbubble_inj=Cycle.AS.T() #[K]
        cp_bubble=Cycle.AS.cpmass() #[J/kg-K]
        cp_inj=(cp_dew+cp_bubble)/2 #[J/kg-K]
        Cycle.Compressor.pin_r=pevap
        Cycle.Compressor.pout_r=pcond
        Cycle.Compressor.Tinj_r=Tdew_inj+Cycle.PHEHX.DT_sh_target
        Cycle.Compressor.pinj_r=pinj
        Cycle.Compressor.Tin_r=Tevap+Cycle.Evaporator.DT_sh
        Cycle.Compressor.Ref=Cycle.Ref
        Cycle.Compressor.Backend=Cycle.Backend
        Cycle.Compressor.Calculate()
        W=Cycle.Compressor.W
        
        # Evaporator fully-dry analysis
        Qevap_dry=epsilon*Cycle.Evaporator.Fins.Air.Vdot_ha*rho_air*(Cycle.Evaporator.Fins.Air.Tdb-Tevap)*Cp_air
        
        #Air-side heat transfer UA
        Evap=Cycle.Evaporator
        Evap.mdot_r=Cycle.Compressor.mdot_r
        Evap.psat_r=pevap
        Evap.Ref=Cycle.Ref
        Evap.Backend=Cycle.Backend
        Evap.Initialize()
        UA_a=Evap.Fins.h_a*Evap.Fins.A_a*Evap.Fins.eta_a
        Tin_a=Evap.Fins.Air.Tdb
        Tout_a=Tin_a+Qevap_dry/(Evap.Fins.mdot_da*Evap.Fins.cp_da)
        #Refrigerant-side heat transfer UA
        UA_r=Evap.A_r_wetted*Correlations.ShahEvaporation_Average(0.5,0.5,Cycle.AS,Evap.G_r,Evap.ID,Evap.psat_r,Qevap_dry/Evap.A_r_wetted,Evap.Tbubble_r,Evap.Tdew_r)
        #Get wall temperatures at inlet and outlet from energy balance
        T_so_a=(UA_a*Evap.Tin_a+UA_r*Tevap)/(UA_a+UA_r)
        T_so_b=(UA_a*Tout_a+UA_r*Tevap)/(UA_a+UA_r)
        
        Tdewpoint=HAPropsSI('D','T',Cycle.Evaporator.Fins.Air.Tdb, 'P',101325, 'R',Evap.Fins.Air.RH)
        
        #Now calculate the fully-wet analysis
        #Evaporator is bounded by saturated air at the refrigerant temperature.
        h_ai=HAPropsSI('H','T',Cycle.Evaporator.Fins.Air.Tdb, 'P',101325, 'R', Cycle.Evaporator.Fins.Air.RH) #[J/kg_da]
        h_s_w_o=HAPropsSI('H','T',Tevap, 'P',101325, 'R', 1.0) #[J/kg_da]
        Qevap_wet=epsilon*Cycle.Evaporator.Fins.Air.Vdot_ha*rho_air*(h_ai-h_s_w_o)
        
        #Coil is either fully-wet, fully-dry or partially wet, partially dry
        if T_so_a>Tdewpoint and T_so_b>Tdewpoint:
            #Fully dry, use dry Q
            f_dry=1.0
        elif T_so_a<Tdewpoint and T_so_b<Tdewpoint:
            #Fully wet, use wet Q
            f_dry=0.0
        else:
            f_dry=1-(Tdewpoint-T_so_a)/(T_so_b-T_so_a)
        Qevap=f_dry*Qevap_dry+(1-f_dry)*Qevap_wet
        
        if Cycle.ImposedVariable == 'Subcooling': #if Subcooling impose
            Cycle.AS.update(CP.PT_INPUTS,pcond,Tcond-Cycle.DT_sc_target)
            h_target = Cycle.AS.hmass() #[J/kg]
            cp_cond = Cycle.AS.cpmass() #[J/kg-K]
            Cycle.AS.update(CP.PT_INPUTS,pcond,Tcond-Cycle.DT_sc_target-5)
            cp_econ = Cycle.AS.cpmass() #[J/kg-K] #assume economizer has 5K more subcooling
            Qcond_enthalpy=Cycle.Compressor.mdot_tot*(Cycle.Compressor.hout_r-h_target)
        else: #otherwise, if Charge impose
            Cycle.AS.update(CP.PT_INPUTS,pcond,Tcond-5)
            h_target = Cycle.AS.hmass() #[J/kg]
            cp_cond = Cycle.AS.cpmass() #[J/kg-K]
            Cycle.AS.update(CP.PT_INPUTS,pcond,Tcond-5-5)
            cp_econ = Cycle.AS.cpmass() #[J/kg-K] #assume economizer has 5K more subcooling
            Qcond_enthalpy=Cycle.Compressor.mdot_tot*(Cycle.Compressor.hout_r-h_target)
        
        #New analysis for the economizer
#         cp_c = cp_inj #[J/kg-K] taken at average injection state
#         cp_h = cp_cond #[J/kg-K] specific heat at the inlet of hot side
#         Cc = Cycle.Compressor.mdot_inj * cp_c
#         Ch = Cycle.Compressor.mdot_tot * cp_h
#         if Ch < Cc: Cmin = Ch
#         else: Cmin = Cc
#         #Temperature at the inlet of cold side (two phase)
#         Ty = (Cc * (Tdew_inj+Cycle.PHEHX.DT_sh_target) - epsilon * Cmin * (Tcond-Cycle.DT_sc_target)) /(Cc - epsilon*Cmin) 
#         #Enthalpy at the inlet of cold side (two phase), approximated qaulity becasue it should be 2-phase
#         xy=(Ty-Tbubble_inj)/(Tdew_inj-Tbubble_inj)
#         if xy <= 0:
#             xy = 0.00001
#         elif xy >= 1:
#             xy = 0.99999
#         Cycle.AS.update(CP.PQ_INPUTS,pinj,xy)
#         hy = Cycle.AS.hmass() #[J/kg]
        #Enthalpy at the inlet of evaporator
        hevap = Cycle.Compressor.hin_r - Qevap/Cycle.Compressor.mdot_r
        #Enthalpy at the exit of the hot side (single phase)
        #hx = (Cycle.Compressor.mdot_r * hevap + Cycle.Compressor.mdot_inj * hy) / Cycle.Compressor.mdot_tot
        Q_c = Cycle.Compressor.mdot_inj * (Cycle.Compressor.hinj_r - hevap)
        Q_h = Cycle.Compressor.mdot_tot * (h_target - hevap)
#         Cycle.AS.update(CP.HmassP_INPUTS,hevap,pcond)
#         Tz = Cycle.AS.T()
#         Tx = Ty = Tz
#         Cycle.AS.update(CP.PT_INPUTS,pcond,Tx)
#         hx=Cycle.AS.hmass()
#         hy =(hx*Cycle.Compressor.mdot_tot - hevap*Cycle.Compressor.mdot_r)/Cycle.Compressor.mdot_inj
# #         Cycle.AS.update(CP.HmassP_INPUTS,hy,pinj)
# #         print Cycle.AS.Q()
#         Q_c = Cycle.Compressor.mdot_inj * (Cycle.Compressor.hinj_r - hy)
#         Q_h = Cycle.Compressor.mdot_tot * (h_target - hx)
#         Cycle.AS.update(CP.PQ_INPUTS,pinj,0.5)
#         hy = Cycle.AS.hmass() #[J/kg]
#         Cycle.PHEHX.pin_c=pinj
#         Cycle.PHEHX.pin_h=pcond
#         Cycle.PHEHX.hin_c=hy
#         Cycle.PHEHX.hin_h=h_target
#         Cycle.PHEHX.mdot_c=Cycle.Compressor.mdot_inj
#         Cycle.PHEHX.mdot_h=Cycle.Compressor.mdot_tot
#         Cycle.PHEHX.Ref_c=Cycle.Ref
#         Cycle.PHEHX.Ref_h=Cycle.Ref
#         Cycle.PHEHX.Backend_c=Cycle.Backend
#         Cycle.PHEHX.Backend_h=Cycle.Backend
#         Cycle.PHEHX.Calculate()
#         Q_econ = Cycle.PHEHX.Q
        
        resids=[Qevap+W+Qcond-Q_h,Qcond+Qcond_enthalpy,Q_h-Q_c]#,Cycle.Compressor.mdot_inj*(Cycle.Compressor.hinj_r - Cycle.PHEHX.hout_c)]

        return resids
    
    Tevap_init=Cycle.Evaporator.Fins.Air.Tdb-15
    Tcond_init=Cycle.Condenser.Fins.Air.Tdb+15
    Tdew_inj_init=np.sqrt(Tevap_init*Tcond_init)

    #First try using the fsolve algorithm
#     try:
    x=fsolve(OBJECTIVE,[Tevap_init,Tcond_init,Tdew_inj_init])
#     cons = ()
#     bnds = ((None,Cycle.Evaporator.Fins.Air.Tdb), (Cycle.Condenser.Fins.Air.Tdb,None), (Cycle.Evaporator.Fins.Air.Tdb,Cycle.Condenser.Fins.Air.Tdb))
#     guess = (Tevap_init,Tcond_init,Tdew_inj_init)
#     res = minimize(OBJECTIVE, guess, method='SLSQP', bounds=bnds, constraints=cons, tol=1e-6)
#     except:
#         #If that doesnt work, try the Mult-Dimensional Newton-raphson solver
#         try:
#             print 'try using the Mult-Dimensional Newton-raphson solver'
#             x=MultiDimNewtRaph(OBJECTIVE,[Tevap_init,Tcond_init,Tdew_inj_init])
#         except:
#             #use the simplified equation propsed by Emerson (see Thomas thesis page 21)
#             print 'try using simplified equation propsed by Emerson'
#             Tdew_inj_init = 0.8 * K2F(Tevap_init) + 0.5 * K2F(Tcond_init) - 21 #[F]
#             x=[Tevap_init,Tcond_init,F2K(Tdew_inj_init)]
            
    DT_evap=Cycle.Evaporator.Fins.Air.Tdb-x[0]
    DT_cond=x[1]-Cycle.Condenser.Fins.Air.Tdb
    Tdew_inj=x[2]

    return DT_evap-1, DT_cond+5, Tdew_inj+20
