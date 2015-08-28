from __future__ import division
from CoolProp.CoolProp import PropsSI
from Correlations import f_h_1phase_Tube,TrhoPhase_ph
from math import log,pi,exp

class SightGlassFilterDrierMicroMotionClass():
    def __init__(self,**kwargs):
        #Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    def Update(self,**kwargs):
        #Load the parameters passed in
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
            ('height of Sight Glass','m',self.h),
            ('Sight glass diameter','m',self.D),
            ('Charge','kg',self.Charge),
            ('Volume filter drier', 'm^3',self.V),
            ('Micromotion tube diameter', 'm', self.D_Micro),
            ('Micormotion tube length','m', self.L_Micro),
            ('Micormotion number of tubes','-', self.n_Micro),

         ]
    
    def Calculate(self):
        if not 'INCOMP' in self.Ref: #if not IsFluidType(self.Ref,'Brine'):
            #Figure out the inlet state
            self.Tbubble=PropsSI('T','P',self.pin,'Q',0.0,self.Ref)
            self.Tdew=PropsSI('T','P',self.pin,'Q',1.0,self.Ref)
        else:
            #It is a brine
            self.Tbubble = None
            self.Tdew = None
        
        self.Tin,self.rhoin,self.Phasein=TrhoPhase_ph(self.Ref,self.pin,self.hin,self.Tbubble,self.Tdew)
        ###Solver shows TwoPhase in the first iteration, the following if statement just to avoid ValueError with CoolProp for pseudo-pure refrigerants
        if self.Phasein =='TwoPhase':
            self.f_fluid, self.h_fluid, self.Re_fluid=f_h_1phase_Tube(self.mdot, self.ID, self.Tin-1, self.pin, self.Ref)
            # Specific heat capacity [J/kg-K]                        
            cp=PropsSI('C','T',self.Tin-1,'P',self.pin,self.Ref)
            # Density [kg/m^3]
            rho=PropsSI('D','T',self.Tin-1, 'P', self.pin, self.Ref)              
        else: #Single phase
            self.f_fluid, self.h_fluid, self.Re_fluid=f_h_1phase_Tube(self.mdot, self.D_Micro, self.Tin, self.pin, self.Ref)
            # Specific heat capacity [J/kg-K]                        
            cp=PropsSI('C','T',self.Tin,'P',self.pin,self.Ref) #*1000
            # Density [kg/m^3]
            rho=PropsSI('D','T',self.Tin, 'P', self.pin, self.Ref)

    
        #Pressure drop calculations for single phase refrigerant
        v=1./rho
        #G=self.mdot/(pi*self.ID**2/4.0)
        #Pressure gradient using Darcy friction factor
        #dpdz=-self.f_fluid*v*G**2/(2.*self.ID) #Pressure gradient
        #self.DP=dpdz*(2*self.B+self.E + self.h) #For total length of sight glass and micromotion only)
        
        #Charge in Sight Glass [kg]
        self.SightGlassCharge = pi*self.D**2/4.0*self.h*rho
        #Charge in FilterDrier [kg]
        self.FilterDrierCharge= self.V*rho
        #Charge in MicroMotion [kg]
        self.MicroMotionCharge= self.n_Micro*pi*self.D_Micro**2/4.0*self.L_Micro*rho
        #Total Chnarge [kg]
        self.Charge= 2*self.SightGlassCharge + self.FilterDrierCharge + self.MicroMotionCharge #multiply sight glass by 2 because we have 2 SG