from __future__ import division, print_function, absolute_import
from math import pi, sqrt, exp, log, tanh, cos, log10
from scipy.optimize import brentq
from scipy.integrate import quad
import CoolProp as CP
from CoolProp.CoolProp import HAPropsSI


class DWSVals():
    """ 
    Empty Class for passing data with DryWetSegment
    """
    def __init__(self):
        #Don't do anything
        pass

def DryWetSegment(DWS):
    """
    Generic solver function for air-side and NTU heat rate
    """
    
    # Retrieve values from structures defined above
    Tin_a=DWS.Tin_a
    if DWS.h_a<0.000000001:
        print ("Warning: Dws.h_a was constrained to 0.001, original value: ", DWS.h_a)
        h_a=0.000000001
    else:
        h_a=DWS.h_a    
    cp_da=DWS.cp_da
    eta_a=DWS.eta_a  #from fin correlations, overall airside surface effectiveness
    A_a=DWS.A_a
    pin_a=DWS.pin_a
    RHin_a=DWS.RHin_a
    mdot_da=DWS.mdot_da

    Tin_r=DWS.Tin_r
    pin_r=DWS.pin_r
    if DWS.h_r<0.000000001:
        print ("Warning: Dws.h_r was constrained to 0.001, original value: ", DWS.h_r)
        h_r=0.000000001
    else:
        h_r=DWS.h_r
    cp_r=DWS.cp_r
    A_r=DWS.A_r
    mdot_r=DWS.mdot_r

    # Calculate the dewpoint (amongst others)
    omega_in=HAPropsSI('W','T',Tin_a,'P',pin_a,'R',RHin_a)
    Tdp=HAPropsSI('D','T',Tin_a,'P',pin_a,'W',omega_in)
    hin_a=HAPropsSI('H','T',Tin_a,'P',pin_a,'W',omega_in) #[J/kg_da]
    
    # Internal UA between fluid flow and outside surface (neglecting tube conduction)
    UA_i=h_r*A_r #[W/K], from Shah or f_h_1phase_Tube-fct -> Correlations.py
    # External UA between wall and free stream
    UA_o=eta_a*h_a*A_a #[W/K], from fin correlations
    # Internal Ntu
    Ntu_i=UA_i/(mdot_r*cp_r)   #[-]
    # External Ntu (multiplied by eta_a since surface is finned and has lower effectiveness)
    Ntu_o=eta_a*h_a*A_a/(mdot_da*cp_da) #[-]
    
    # Start calculation  
    # Overall UA
    UA = 1 / (1 / (UA_i) + 1 / (UA_o));
    # Min and max capacitance rates [W/K]
    Cmin = min([cp_r * mdot_r, cp_da * mdot_da])
    Cmax = max([cp_r * mdot_r, cp_da * mdot_da])
    # Capacitance rate ratio [-]
    C_star = Cmin / Cmax
    # Ntu overall [-]
    Ntu_dry = UA / Cmin
    
    if Ntu_dry<0.0000001:
        print("warning:  NTU_dry in dry wet segment was negative. forced it to positive value of 0.001!")
        Ntu_dry=0.0000001
    
    # Crossflow effectiveness
    if (cp_r * mdot_r)<(cp_da * mdot_da):
        epsilon_dry= 1-exp(-C_star**(-1)*(1-exp(-C_star*(Ntu_dry))))
        #Cross flow, single phase, cmax is airside, which is unmixed
    else:
        epsilon_dry=(1/C_star)*(1-exp(-C_star*(1-exp(-Ntu_dry))))
        #Cross flow, single phase, cmax is refrigerant side, which is mixed

    # Dry heat transfer [W]
    Q_dry = epsilon_dry*Cmin*(Tin_a-Tin_r)
    # Dry-analysis air outlet temp [K]
    Tout_a_dry=Tin_a-Q_dry/(mdot_da*cp_da)
    # Dry-analysis outlet temp [K]
    Tout_r=Tin_r+Q_dry/(mdot_r*cp_r)
    # Dry-analysis air outlet enthalpy from energy balance [J/kg]
    hout_a=hin_a-Q_dry/mdot_da
    # Dry-analysis surface outlet temp [K]
    Tout_s=(UA_o*Tout_a_dry+UA_i*Tin_r)/(UA_o+UA_i)
    # Dry-analysis surface inlet temp [K]
    Tin_s=(UA_o*Tin_a+UA_i*Tout_r)/(UA_o+UA_i)
    # Dry-analysis outlet refrigerant temp [K]
    Tout_r_dry=Tout_r
    # Dry fraction [-]
    f_dry=1.0
    # Air outlet humidity ratio [-]
    DWS.omega_out = omega_in

    # If inlet surface temp below dewpoint, whole surface is wetted 
    if Tin_s<Tdp:
        isFullyWet=True
    else:
        isFullyWet=False

    if Tout_s<Tdp or isFullyWet:
        # There is some wetting, either the coil is fully wetted or partially wetted 
        print('coil is fully or partially wet')
        if (Tin_s>Tdp and not isFullyWet):
            # Partially wet and partially dry with single-phase on refrigerant side
            print('coil is partially wet and partially dry')
            raise
    else:
        # Coil is fully dry
        Tout_a=Tout_a_dry
        Q=Q_dry
        Q_sensible=Q_dry
          

    DWS.f_dry=f_dry
    DWS.omega_out=HAPropsSI('W','T',Tout_a,'P',101325,'H',hout_a)
    DWS.RHout_a=HAPropsSI('R','T',Tout_a,'P',101325,'W',DWS.omega_out)
    DWS.Tout_a=Tout_a
    DWS.Q=Q
    DWS.Q_sensible=Q_sensible
    
    DWS.hout_a=hout_a
    DWS.hin_a=hin_a
    DWS.Tout_r=Tout_r
    DWS.Twall_s=Tout_r - Q/UA_i #inner wall temperature for gas cooler model

class FinsVals:
    pass    
class TubesVals:
    pass
class AirVals:
    pass

class FinInputs():
    """ 
    Empty Class for fin data
    """
    def __init__(self):
        #Don't do anything
        self.Tubes=TubesVals()
        self.Air=AirVals()
        self.Fins=FinsVals()
    def __repr__(self):
        string="Tubes::\n"
        for field in self.Tubes.__dict__.keys():
            string+=field+": "+repr(self.Tubes.__dict__[field])+"\n"
        string+="Fins::\n"
        for field in self.Fins.__dict__.keys():
            string+=field+": "+repr(self.Fins.__dict__[field])+"\n"
        string+="Air::\n"
        for field in self.Air.__dict__.keys():
            string+=field+": "+repr(self.Air.__dict__[field])+"\n"
        return string
            
def WavyLouveredFins(Inputs):
    """
    # Chi-Chuan Wang and Yu-Min Tsai and Ding-Chong Lu, 1998, "Comprehensive
    # Study of Convex-Louver and Wavy Fin-and-Tube Heat Exchangers", Journal 
    # of Thermophysics and Heat Transfer
    
    || -    xf    ->
    ^              ==                          ==
    |           ==  |  ==                  ==
    Pd       ==     |     ==            ==
    |     ==        |        ==     ==
    =  ==           s             ==  
                    |
                    |
                    |
                   ==                        ==
                ==     ==                 ==
             ==           ==           ==
          ==                 ==     ==
       ==                        ==  
    
     t: thickness of fin plate
     Pf: fin pitch (centerline-to-centerline distance between fins)
     Pd: indentation for waviness (not including fin thickness)
     s: fin spacing (free space between fins) = Pf-t
    
    
    
                 |--       Pl      -|
                ___                 |
              /     \               |
       =     |       |              |
       |      \ ___ /               |
       |                            |
       |                           ___      
       |                         /     \      |
      Pt                        |       |     D
       |                         \ ___ /      |
       |    
       |        ___
       |      /     \
       =     |       |
              \ ___ /
    """
    
    Ntubes_bank = Inputs.Tubes.NTubes_per_bank  #tubes per bank
    Nbank =       Inputs.Tubes.Nbank            #Number of banks
    Ltube =       Inputs.Tubes.Ltube            #length of a single tube
    D =           Inputs.Tubes.OD               #Outer diameter of tube
    Pl =          Inputs.Tubes.Pl               #Horizontal spacing between banks (center to center)
    Pt =          Inputs.Tubes.Pt               #Vertical spacing between tubes in a bank (center to center)

    FPI =         Inputs.Fins.FPI
    pd =          Inputs.Fins.Pd
    xf =          Inputs.Fins.xf
    t =           Inputs.Fins.t
    k_fin =       Inputs.Fins.k_fin

    Vdot_ha =     Inputs.Air.Vdot_ha
    p =           Inputs.Air.p
    
    if hasattr(Inputs,'Tmean'):
        Temp = Inputs.Air.Tmean
    else:
        Temp = Inputs.Air.Tdb
    
    if hasattr(Inputs,'RHmean'):
        RHin = Inputs.Air.RHmean
    else:
        RHin = Inputs.Air.RH

    # Check that cs_cp is defined, if so, set it to the value passed in
    if (hasattr(Inputs,'cs_cp') and Inputs.cs_cp>0) or (hasattr(Inputs,'WetDry') and Inputs.WetDry=='Wet'):
        isWet=True
        cs_cp=Inputs.Air.cs_cp
    else:
        isWet=False
        cs_cp=1.0

    # Film temperature [K]
    Tfilm = Temp
    # Fins per meter [1/m]
    FPM = FPI / 0.0254
    # Fin pitch (distance between centerlines of fins)
    pf = 1 / FPM
    # Spacing between fins
    s = 1 / FPM - t

    # Height of heat exchanger [m]
    Height = Pt * (Ntubes_bank+1)  #assuming that fin extends 1/2 pt above/below last tube in bundle
    # A_duct is the face area [m^2] equivalent to the duct cross-section
    A_duct = Height * Ltube  #neglecting the additional height of the fins above/below the last tubes in the bundle
    # Number of fins in the tube sheet [-]
    Nfin = Ltube * FPM
    # Secant of theta is the area enhancement factor [-]
    # It captures the increase in area due to the waviness of the fins 
    sec_theta = sqrt(xf*xf + pd*pd) / xf
    # Duct cross-sectional area that is not fin or tube [m^2]
    Ac = A_duct - t * Nfin * (Height-D*Ntubes_bank) - Ntubes_bank * D * Ltube
    # Total outer area of the tubes [m^2]
    Atube = Ntubes_bank * Nbank * pi * D * Ltube
    # Wetted Area of a single fin [m^2]
    A_1fin = 2.0 * (Height * Pl * (Nbank+1) * sec_theta  - Ntubes_bank*Nbank * pi*D*D/4) #assuming that fin extends 1/2 pt in front/after last tube in bundle
    # Total wetted area of the fins [m^2]
    Af = Nfin * A_1fin
    # Total area including tube and fins [m^2]
    A = Af + Ntubes_bank * Nbank * pi * D * (Ltube-Nfin*t)

    # Evaluate the mass flow rate based on inlet conditions
    # To convert a parameter from per kg_{dry air} to per kg_{humid air}, divide by (1+W)
    W=HAPropsSI('W','T',Inputs.Air.Tdb,'P',p,'R',Inputs.Air.RH)
    v_da=HAPropsSI('V','T',Inputs.Air.Tdb,'P',p,'W',W)
    h_da=HAPropsSI('H','T',Inputs.Air.Tdb,'P',p,'W',W) #[J/kg]
    rho_ha = 1 / v_da*(1+W) #[kg_ha/m^3]
    rho_da = 1 / v_da #[kg_da/m^3]
    mdot_ha = Vdot_ha * rho_ha #[kg_ha/s]
    mdot_da = Vdot_ha * rho_da #[kg_da/s]

    umax = mdot_ha / (rho_ha * Ac) #[m/s]

    # Use a forward difference to calculate cp from cp=dh/dT
    dT=0.0001 #[K]
    cp_da=(HAPropsSI('H','T',Inputs.Air.Tdb+dT,'P', p, 'W',W)-h_da)/dT #[J/kg_da/K]
    cp_ha=cp_da/(1+W) #[J/kg_ha/K]
    
    # Transport properties of humid air from CoolProp
    mu_ha=HAPropsSI('M','T',Inputs.Air.Tdb,'P',p,'W',W)
    k_ha=HAPropsSI('K','T',Inputs.Air.Tdb,'P',p,'W',W)
    
    # Dimensionless Groups
    Pr = cp_ha * mu_ha / k_ha
    Re_D = rho_ha * umax * D / mu_ha

    # Heat transfer
    j = 16.06 * pow(Re_D,-1.02 * (pf / D) - 0.256) * pow(A / Atube, -0.601) * pow(Nbank,-0.069) * pow(pf / D,0.84) #Colburn j-Factor
    h_a = j * rho_ha * umax * cp_ha / pow(Pr,2.0/3.0) #air side mean heat transfer coefficient
    
    # Air side pressure drop correlations
    if (Re_D<1e3):
        fa_total=0.264*(0.105+0.708*exp(-Re_D/225.0))*pow(Re_D,-0.637)*pow(A/Atube,0.263)*pow(pf/D,-0.317)
    else:
        fa_total=0.768*(0.0494+0.142*exp(-Re_D/1180.0))*pow(A/Atube,0.0195)*pow(pf/D,-0.121)
    
    # Calcs needed for specific fin types
    r = D / 2
    X_D = sqrt(Pl*Pl + Pt*Pt / 4) / 2
    X_T = Pt / 2
    rf_r = 1.27 * X_T / r * sqrt(X_D / X_T - 0.3)
    m = sqrt(2 * h_a * cs_cp / (k_fin * t)) #cs_cp is the correction for heat/mass transfer
    
    # Using the circular fin correlation of Schmidt
    phi = (rf_r - 1) * (1 + 0.35 * log(rf_r))
    eta_f = tanh(m * r * phi) / (m * r * phi)
    
    # Fin efficiency
    phi = (rf_r - 1) * (1 + (0.3+pow(m*(rf_r*r-r)/2.5,1.5-rf_r/12.0)*(0.26*pow(rf_r,0.3)-0.3)) * log(rf_r))
    
    # Finned surface efficiency
    eta_f = tanh(m * r * phi) / (m * r * phi)*cos(0.1 * m * r * phi)
    
    # Overall surface efficiency
    eta_o = 1 - Af / A * (1 - eta_f)

    G_c=mdot_ha/Ac #air mass flux
    DeltaP_air=A/Ac/rho_ha*G_c**2/2.0*fa_total #airside pressure drop
    
    # Write necessary values back into the given structure
    Inputs.A_a=A;
    Inputs.cp_da=cp_da
    Inputs.cp_ha=cp_ha
    if isWet==True:
        Inputs.eta_a_wet=eta_o
    else:
        Inputs.eta_a=eta_o
    Inputs.h_a=h_a
    Inputs.mdot_ha=mdot_ha
    Inputs.mdot_da=mdot_da
    Inputs.f_a=fa_total
    Inputs.dP_a=DeltaP_air
    Inputs.Re=Re_D

def Petterson_supercritical_average(Tout,Tin,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w):
    '''
    Petterson et al. (2000), Heat transfer and pressure drop for flow supercritical and subcritical CO2 in microchannel tubes
    '''

    def Petterson_supercritical(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w):
        AS.update(CP.PT_INPUTS,p,T_w)
        h_w = AS.hmass() #[J/kg]
        mu_w = AS.viscosity() #[Pa-s OR kg/m-s]
        cp_w = AS.cpmass() #[J/kg-K]
        k_w = AS.conductivity() #[W/m-K]
        rho_w = AS.rhomass() #[kg/m^3]
        Pr_w = cp_w * mu_w / k_w #[-]
        
        AS.update(CP.PT_INPUTS,p,T)
        h = AS.hmass() #[J/kg]
        mu = AS.viscosity() #[Pa-s OR kg/m-s]
        cp = AS.cpmass() #[J/kg-K]
        k = AS.conductivity() #[W/m-K]
        rho = AS.rhomass() #[kg/m^3]
        Pr = cp * mu / k #[-]
        
        Dh = OD - ID
        Area=pi*(OD**2-ID**2)/4.0
        u=mdot/(Area*rho)
        Re=rho*u*Dh/mu
        Re_w=Re#rho_w*u*Dh/mu_w
            
        if G > 350:
            e_D = 1e-6 #smooth pipe
            f = (-1.8*log10(6.9/Re + (1/3.7*e_D)**1.11))**(-2)/4
            Nu_m = (f/8)*(Re-1000)*Pr/(1+12.7*sqrt(f/8)*(Pr**(2/3)-1)) *(1+(D_l)**(2/3))
            Nu = Nu_m * (Pr/Pr_w)**0.11
         
        else: # G<350
             
            M = 0.001 #[kg/J]
            K = 0.00041 #[kg/J]
             
            cp_avg = (h-h_w)/(T-T_w)
             
            if cp_avg/cp_w <= 1:
                n = 0.66 - K*(q_flux_w/G)
            else: #cp_avg/cp_w >1
                n = 0.9 - K*(q_flux_w/G)
             
            f0 = (0.79*log(Re)-1.64)**(-2)
             
            g =9.81
            #coefficient of thermal expansion
            beta = AS.isobaric_expansion_coefficient() #[1/K]
            #Grashoff number
            Gr = g*beta*(T_w-T)*Dh**3/(mu/rho)**2
            if Gr/Re**2 < 5e-4:
                f = f0 * (mu_w/mu)**0.22
            elif  Gr/Re**2 >= 5e-4 and G/Re**2 < 0.3:
                f = 2.15 * f0 * (mu_w/mu)**0.22 * (Gr/Re)**0.1
            else: #use f0 for friction factor
                f = f0
                 
            Nu_w_ppk = (f0/8)*Re_w*Pr_w/(1.07+12.7*sqrt(f/8)*(Pr_w**(2/3)-1))
             
            Nu = Nu_w_ppk * (1-M*q_flux_w/G) * (cp_avg/cp_w)**n
             
        h = k*Nu/Dh #[W/m^2-K]
         
        return (h,f,cp,rho)

    def SuperCriticalCondensation_h(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w):
        '''return h value'''
        return Petterson_supercritical(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w)[0]
    def SuperCriticalCondensation_f(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w):
        '''return f value'''
        return Petterson_supercritical(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w)[1]
    def SuperCriticalCondensation_cp(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w):
        '''return cp value'''
        return Petterson_supercritical(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w)[2]
    def SuperCriticalCondensation_rho(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w):
        '''return rho value'''
        return Petterson_supercritical(T,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w)[3]
            
    if not Tout==Tin:
        #A proper range is given
        h = quad(SuperCriticalCondensation_h,Tin,Tout,args=(T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w))[0]/(Tout-Tin)
        f = quad(SuperCriticalCondensation_f,Tin,Tout,args=(T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w))[0]/(Tout-Tin)
        cp = quad(SuperCriticalCondensation_cp,Tin,Tout,args=(T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w))[0]/(Tout-Tin)
        rho = quad(SuperCriticalCondensation_rho,Tin,Tout,args=(T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w))[0]/(Tout-Tin)
        return (h,f,cp,rho)
    else:
        #A single value is given
        return Petterson_supercritical(Tout,T_w,AS,G,OD,ID,D_l,mdot,p,q_flux_w)
        
class GasCoolerClass():
    def __init__(self,**kwargs):
        # Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)      
    
    def Update(self,**kwargs):
        # Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    
    def OutputList(self):
        """
            Return a list of parameters for this component for further output
            
            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        return [
            ('Volumetric flow rate of humid air','m^3/s',self.Fins.Air.Vdot_ha),
            ('Inlet Dry bulb temp','K',self.Tin_a),
            ('Inlet Air pressure','Pa',self.Fins.Air.p),
            ('Inlet Air Relative Humidity','-',self.Fins.Air.RH),
            ('Tubes per bank','-',self.Fins.Tubes.NTubes_per_bank),
            ('Number of banks','-',self.Fins.Tubes.Nbank),
            ('Number circuits','-',self.Fins.Tubes.Ncircuits),
            ('Length of tube','m',self.Fins.Tubes.Ltube),
            ('Tube OD','m',self.OD),
            ('Tube ID','m',self.ID),
            ('Tube Long. Pitch','m',self.Fins.Tubes.Pl),
            ('Tube Transverse Pitch','m',self.Fins.Tubes.Pt),
            ('Fins per inch','1/in',self.Fins.Fins.FPI),
            ('Fin waviness pd','m',self.Fins.Fins.Pd),
            ('Fin waviness xf','m',self.Fins.Fins.xf),
            ('Fin thickness','m',self.Fins.Fins.t),
            ('Fin Conductivity','W/m-K',self.Fins.Fins.k_fin),
            ('Fins Type','-',self.FinsType),
            ('Q Total','W',self.Q),
            ('Q Supercritical','W',self.Q_supercritical),
            ('Q Supercritical_liquid','W',self.Q_subcool),
            ('Inlet Temp','K',self.Tin_r),
            ('Outlet Temp','K',self.Tout_r),
            ('Pressure Drop Total','Pa',self.DP_r),
            ('Pressure Drop Supercritical','Pa',self.DP_r_supercritical),
            ('Pressure Drop Supercritical_liquid','Pa',self.DP_r_subcool),
            ('Charge Total','kg',self.Charge),
            ('Charge Supercritical','kg',self.Charge_supercritical),
            ('Charge Supercritical_liquid','kg',self.Charge_subcool),
            ('Mean HTC Superheat','W/m^2-K',self.h_r_supercritical),
            ('Mean HTC Supercritical_liquid','W/m^2-K',self.h_r_subcool),
            ('Wetted Area Fraction Supercritical','-',self.w_supercritical),
            ('Wetted Area Fraction Supercritical_liquid','-',self.w_subcool),
            ('Mean Air HTC','W/m^2-K',self.Fins.h_a),
            ('Surface Effectiveness','-',self.Fins.eta_a),
            ('Air-side area (fin+tubes)','m^2',self.Fins.A_a),
            ('Mass Flow rate of Dry Air','kg/s',self.Fins.mdot_da),
            ('Mass Flow rate of Humid Air','kg/s',self.Fins.mdot_ha),
            ('Pressure Drop Air-side','Pa',self.Fins.dP_a),
            ('Approach temperature degree','K',self.DT_app),
        ]
        

    def Initialize(self):
        # AbstractState
        if hasattr(self,'Backend'): #check if backend is given
            AS = CP.AbstractState(self.Backend, self.Ref)
        else: #otherwise, use the defualt backend
            AS = CP.AbstractState('HEOS', self.Ref)
        self.AS = AS
    
        # Retrieve some parameters from nested structures 
        # for code compactness
        self.ID=self.Fins.Tubes.ID
        self.OD=self.Fins.Tubes.OD
        self.Ltube=self.Fins.Tubes.Ltube
        self.NTubes_per_bank=self.Fins.Tubes.NTubes_per_bank
        self.Nbank=self.Fins.Tubes.Nbank
        self.Ncircuits=self.Fins.Tubes.Ncircuits
        self.Tin_a=self.Fins.Air.Tdb
        
        # Calculate an effective length of circuit if circuits are 
        # not all the same length
        TotalLength=self.Ltube*self.NTubes_per_bank*self.Nbank
        self.Lcircuit=TotalLength/self.Ncircuits
        self.V_r = pi * self.ID**2 / 4.0 * self.Lcircuit * self.Ncircuits
        self.A_r_wetted = pi * self.ID * self.Ncircuits * self.Lcircuit
        self.G_r = self.mdot_r/(self.Ncircuits*pi*self.ID**2/4.0) 
         
        # Define known parameters
        AS.update(CP.PT_INPUTS, self.psat_r, self.Tin_r)
        self.hin_r=AS.hmass() #[J/kg]
        self.sin_r=AS.smass() #[J/kg-K]
        
        # Define critical pressure and temperature
        self.Pcr=AS.p_critical() #[Pa]
        self.Tcr=AS.T_critical() #[K]
        
        # Critical enthalpy at defined pressure
        AS.update(CP.PT_INPUTS, self.psat_r, self.Tcr)
        self.hcr = AS.hmass() #[J/kg]
        self.scr = AS.smass() #[J/kg]
        
        # Triple temperature
        self.Ttriple = AS.Ttriple()
        
        self.Fins.Air.RHmean=self.Fins.Air.RH
        
        # Update with user FinType
        if self.FinsType == 'WavyLouveredFins':
            WavyLouveredFins(self.Fins)

        
        self.mdot_ha=self.Fins.mdot_ha #[kg_ha/s]
        self.mdot_da=self.Fins.mdot_da #[kg_da/s]
        
        
    def Calculate(self):
        
        # Initialize
        self.Initialize()
        AS = self.AS
        
        # Assume we have all supercritical region
        self.Tout_r_cr=self.Tcr
        # Give an intial guess for the inner wall temperature
        self.T_w = (self.Tout_r_cr+self.Tin_a)/2
        # If we have already used too much of the HX (max possible sum of w is 1.0)
        if self._Supercritical_Forward(1.0)<0:
            self.existsSubcooled=False
            self.w_supercritical=1.0
            def OBJECTIVE(Tout_r_cr):
                self.Tout_r_cr=Tout_r_cr
                AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r_cr)
                hout = AS.hmass()
                Q_target=self.mdot_r*(hout-self.hin_r)
                self._Supercritical_Forward(self.w_supercritical)
                return self.Q_supercritical-Q_target
            brentq(OBJECTIVE,self.Tin_r,self.Tcr)
            # Zero out all the supercritical_liquid parameters
            self.Q_subcool=0.0
            self.DP_r_subcool=0.0
            self.Charge_subcool=0.0
            self.w_subcool=0.0
            self.h_r_subcool=0.0
            self.Re_r_subcool=0.0
            self.Tout_a_subcool=0.0
            self.fdry_subcool=0.0    
        else:
            # By definition then we have a supercritical_liquid portion, solve for it
            self.existsSubcooled=True 
            self.w_supercritical=brentq(self._Supercritical_Forward,0.00000000001,0.9999999999)
            def OBJECTIVE(Tout_r_sc):
                self.Tout_r_sc=Tout_r_sc
                AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r_sc)
                hout = AS.hmass()
                Q_target=self.mdot_r*(hout-self.hcr)
                self._Subcool_Forward(1-self.w_supercritical)
                return self.Q_subcool-Q_target
            brentq(OBJECTIVE,self.Tcr,self.Ttriple)
            
            
        # Overall calculations
        self.Q=self.Q_supercritical+self.Q_subcool
        self.DP_r=self.DP_r_supercritical+self.DP_r_subcool
        self.Charge=self.Charge_supercritical+self.Charge_subcool
        
        # Average air outlet temperature (area fraction weighted average) [K]
        self.Tout_a=self.w_supercritical*self.Tout_a_supercritical+self.w_subcool*self.Tout_a_subcool
        
        # Outlet enthalpy obtained from energy balance
        self.hout_r=self.hin_r+self.Q/self.mdot_r
        
        AS.update(CP.HmassP_INPUTS, self.hout_r, self.psat_r)
        self.Tout_r = AS.T() #[K]
        self.sout_r = AS.smass() #[J/kg-K]
        
        # Approach temperature
        self.DT_app = self.Tout_r - self.Tin_a
        
        self.hmean_r=self.w_supercritical*self.h_r_supercritical+self.w_subcool*self.h_r_subcool
        self.UA_r=self.hmean_r*self.A_r_wetted
        self.UA_a=self.Fins.h_a*self.Fins.A_a*self.Fins.eta_a

        
    def _Supercritical_Forward(self,w_supercritical):
        # **********************************************************************
        #                      SUPERCRITICAL PART 
        # **********************************************************************
        
        # AbstractState
        AS = self.AS
        
        DWS=DWSVals() #DryWetSegment structure (only dry-analysis, single phase is used)
    
        # Store temporary values to be passed to DryWetSegment
        DWS.Fins=self.Fins
        DWS.FinsType = self.FinsType                                            
        DWS.A_a=self.Fins.A_a*w_supercritical
        DWS.cp_da=self.Fins.cp_da
        DWS.eta_a=self.Fins.eta_a
        DWS.h_a=self.Fins.h_a  #Heat transfer coefficient, not enthalpy
        DWS.mdot_da=self.mdot_da*w_supercritical
        DWS.pin_a=self.Fins.Air.p

        DWS.Tin_a=self.Tin_a
        DWS.RHin_a=self.Fins.Air.RH
    
        DWS.Tin_r=self.Tin_r
        DWS.A_r=self.A_r_wetted*w_supercritical
        DWS.pin_r=self.psat_r
        DWS.mdot_r=self.mdot_r
        
        AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r_cr)
        hout = AS.hmass() #[J/kg]
        
        # Target heat transfer to go from inlet temperature to iterative outlet temperature
        Q_target=self.mdot_r*(hout-self.hin_r)
        print(-Q_target/DWS.A_r)
        if Q_target>0:
            raise ValueError('Q_target in Gas cooler must be negative')
        
        # This block calculates the average refrigerant heat transfer coefficient, average friction factor, average specific heat, and average density
        DWS.h_r, f_r_supercritical, DWS.cp_r, rho_supercritical = Petterson_supercritical_average(self.Tout_r_cr, self.Tin_r, self.T_w, self.AS, self.G_r, self.ID, 0, self.ID/self.Lcircuit, self.mdot_r / self.Ncircuits, self.psat_r, -Q_target/DWS.A_r);

        # Compute Fins Efficiency based on FinsType 
        DryWetSegment(DWS)
        
        self.T_w=DWS.Twall_s #inner surface wall temperature (refrigerant)
        self.Q_supercritical=DWS.Q
        self.h_r_supercritical=DWS.h_r
        self.fdry_supercritical=DWS.f_dry
        self.Tout_a_supercritical=DWS.Tout_a
        
        # Pressure drop calculations for supercritical refrigerant
        v_r=1./rho_supercritical
        # Pressure gradient using Darcy friction factor
        dpdz_r=-f_r_supercritical*v_r*self.G_r**2/(2*self.ID) #Pressure gradient
        self.DP_r_supercritical=dpdz_r*self.Lcircuit*w_supercritical
        # Charge for the supercritical portion
        self.Charge_supercritical = w_supercritical * self.V_r * rho_supercritical
        
        return Q_target-DWS.Q
    
    
    def _Subcool_Forward(self,w_subcool):
        # **********************************************************************
        #                      SUPERCRITICAL_LIQUID PART 
        # **********************************************************************
        self.w_subcool=w_subcool
        
        if self.w_subcool<0:
            raise ValueError('w_subcool in Gas cooler cannot be less than zero')
        
        # AbstractState
        AS = self.AS
        
        DWS=DWSVals() #DryWetSegment structure
        
        # Store temporary values to be passed to DryWetSegment
        DWS.A_a=self.Fins.A_a*w_subcool
        DWS.cp_da=self.Fins.cp_da
        DWS.eta_a=self.Fins.eta_a
        DWS.h_a=self.Fins.h_a  #Heat transfer coefficient
        DWS.mdot_da=self.mdot_da*w_subcool
        DWS.pin_a=self.Fins.Air.p
        DWS.Fins=self.Fins
        DWS.FinsType = self.FinsType           
    
        # Inputs on the air side to two phase region are inlet air again
        DWS.Tin_a=self.Tin_a
        DWS.RHin_a=self.Fins.Air.RH
    
        DWS.Tin_r=self.Tcr
        DWS.A_r=self.A_r_wetted*w_subcool
        
        DWS.pin_r=self.psat_r
        DWS.mdot_r=self.mdot_r
        
        AS.update(CP.PT_INPUTS, self.psat_r, self.Tout_r_sc)
        hout = AS.hmass() #[J/kg]
        # Target heat transfer to go from inlet temperature to iterative outlet temperature
        Q_target=self.mdot_r*(hout-self.hcr)
        
        if Q_target>0:
            raise ValueError('Q_target in Gas cooler must be negative')
        
        # Friction factor and HTC in the refrigerant portions.
        # Average fluid temps are used for the calculation of properties 
        # Average temp of refrigerant is average of sat. temp and outlet temp
        # Secondary fluid is air over the fins
        DWS.h_r, self.f_r_subcool, DWS.cp_r, rho_subcool = Petterson_supercritical_average(self.Tout_r_sc, self.Tcr, self.T_w, self.AS, self.G_r, self.ID, 0, self.ID/self.Lcircuit, self.mdot_r / self.Ncircuits, self.psat_r, -Q_target/DWS.A_r);
        
        # Run DryWetSegment
        DryWetSegment(DWS)
        
        # Positive Q is heat input to the refrigerant, negative Q is heat output from refrigerant. 
        # Heat is removed here from the refrigerant since it is condensing
        self.T_w = DWS.Twall_s #inner surface wall temperature (refrigerant)
        self.Q_subcool=DWS.Q
        self.h_r_subcool=DWS.h_r
        self.fdry_subcool=DWS.f_dry
        self.Tout_a_subcool=DWS.Tout_a
        self.Tout_r=DWS.Tout_r
        
        # Charge calculation
        self.Charge_subcool = self.w_subcool * self.V_r * rho_subcool
    
        # Pressure drop calculations for subcooled refrigerant
        v_r=1/rho_subcool
        # Pressure gradient using Darcy friction factor
        dpdz_r=-self.f_r_subcool*v_r*self.G_r**2/(2*self.ID)  #Pressure gradient
        self.DP_r_subcool=dpdz_r*self.Lcircuit*self.w_subcool
        
        return Q_target-DWS.Q
     
def SampleGasCooler():
    Fins=FinInputs()
    Fins.Tubes.NTubes_per_bank=18       #number of tubes per bank or row
    Fins.Tubes.Nbank=3                  #number of banks or rows
    Fins.Tubes.Ncircuits=1              #number of circuits
    Fins.Tubes.Ltube=0.61               #one tube length
    Fins.Tubes.OD=7.9/1000
    Fins.Tubes.ID=7.5/1000
    Fins.Tubes.Pl=12.5/1000             #distance between center of tubes in flow direction                                                
    Fins.Tubes.Pt=24.211/1000           #distance between center of tubes orthogonal to flow direction
    
    Fins.Fins.FPI=1/(1.5/1000/0.0254)   #number of fins per inch
    Fins.Fins.Pd=0.001                  #2* amplitude of wavy fin
    Fins.Fins.xf=0.001                  #1/2 period of fin
    Fins.Fins.t=0.13/1000               #thickness of fin material
    Fins.Fins.k_fin=237                 #thermal conductivity of fin material
    
    Fins.Air.Vdot_ha=1*0.281            #rated volumetric flowrate (m^3/s)
    Fins.Air.Tmean=29.4+273.15   
    Fins.Air.Tdb=29.4+273.15            #dry Bulb Temperature
    Fins.Air.p=101325                   #air pressure in Pa
    Fins.Air.RH=0.5                     #relative Humidity
    Fins.Air.RHmean=0.5
    Fins.Air.FanPower=160    
    
    params={
        'Ref': 'R744',
        'mdot_r': 0.038,
        'Tin_r': 118.1+273.15,
        'psat_r': 9*1000000,
        'Fins': Fins,
        'FinsType': 'WavyLouveredFins',
        'Verbosity':0,
        'Backend':'HEOS'
    }
    GasCooler=GasCoolerClass(**params)
    GasCooler.Calculate()
    return GasCooler
    
if __name__=='__main__':
    GasCooler=SampleGasCooler()
    print ('Refrigerant outlet temperature in gascooler is', GasCooler.Tout_r-273.15,'C')
    print ('Heat transfer rate in gascooler is', GasCooler.Q,'W')
    print ('Heat transfer rate in gascooler (supercritical section) is', GasCooler.Q_supercritical,'W')
    print ('Heat transfer rate in gascooler (supercritical_liquid section) is', GasCooler.Q_subcool,'W')
    print ('Fraction of circuit length in supercritical section is', GasCooler.w_supercritical)
    print ('Fraction of circuit length in supercritical_liquid section is', GasCooler.w_subcool)