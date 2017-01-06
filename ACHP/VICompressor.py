from __future__ import division
from CoolProp.CoolProp import PropsSI
import CoolProp as CP
from convert_units import *

class VICompressorClass():
    """
    Compressor Model based on 10-coefficient Model from `ANSI/AHRI standard 540 <http://www.ahrinet.org/App_Content/ahri/files/standards%20pdfs/ANSI%20standards%20pdfs/ANSI-ARI-540-2004%20latest.pdf>`_
    
    Required Parameters:
        
    ===========   ==========  ========================================================================
    Variable      Units       Description
    ===========   ==========  ========================================================================
    M             Ibm/hr      A numpy-like list of compressor map coefficients for mass flow
    P             Watts       A numpy-like list of compressor map coefficients for electrical power
    Ref           N/A         A string representing the refrigerant
    Tin_r         K           Refrigerant inlet temperature
    pin_r         Pa          Refrigerant suction pressure (absolute)
    pout_r        Pa          Refrigerant discharge pressure (absolute)
    fp            --          Fraction of electrical power lost as heat to ambient
    Vdot_ratio    --          Displacement Scale factor
    ===========   ==========  ========================================================================
    
    All variables are of double-type unless otherwise specified
        
    """
    def __init__(self,**kwargs):
        #Load up the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    def Update(self,**kwargs):
        #Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
        
    def OutputList(self):
        """
            Return a list of parameters for this component for further output
            
            It is a list of tuples, and each tuple is formed of items with indices:
                [0] Description of value
                
                [1] Units of value
                
                [2] The value itself
        """
        
        return [ 
            ('M1','-',self.M[0]),
            ('M2','-',self.M[1]),
            ('M3','-',self.M[2]),
            ('M4','-',self.M[3]),
            ('M5','-',self.M[4]),
            ('M6','-',self.M[5]),
            ('P1','-',self.P[0]),
            ('P2','-',self.P[1]),
            ('P3','-',self.P[2]),
            ('P4','-',self.P[3]),
            ('P5','-',self.P[4]),
            ('P6','-',self.P[5]),
            ('R1','-',self.R[0]),
            ('R2','-',self.R[1]),
            ('R3','-',self.R[2]),
            ('R4','-',self.R[3]),
            ('R5','-',self.R[4]),
            ('R6','-',self.R[5]),
            ('Te1','-',self.Te[0]),
            ('Te2','-',self.Te[1]),
            ('Te3','-',self.Te[2]),
            ('Te4','-',self.Te[3]),
            ('Te5','-',self.Te[4]),
            ('Te6','-',self.Te[5]),
            ('Heat Loss Fraction','-',self.fp),
            ('Displacement scale factor','-',self.Vdot_ratio),
            ('Power','W',self.W),
            ('Suction mass flow rate','kg/s',self.mdot_r),
            ('Injection mass flow rate','kg/s',self.mdot_inj),
            ('Total mass flow rate','kg/s',self.mdot_tot),
            ('Inlet Temperature','K',self.Tin_r),
            ('Injection Temperature','K',self.Tinj_r),
            ('Outlet Temperature','K',self.Tout_r),
            ('Inlet Enthalpy','J/kg',self.hin_r),
            ('Injection Enthalpy','J/kg',self.hinj_r),
            ('Outlet Enthalpy','J/kg',self.hout_r),
            ('Inlet Pressure','Pa',self.pin_r),
            ('Injection Pressure','Pa',self.pinj_r),
            ('Outlet Pressure','Pa',self.pout_r),
            ('Inlet Entropy','J/kg-K',self.sin_r),
            ('Injection Entropy','J/kg-K',self.sinj_r),
            ('Outlet Entropy','J/kg-K',self.sout_r),
            ('Overall isentropic efficiency','-',self.eta_oi),
            ('Pumped flow rate','m**3/s',self.Vdot_pumped),
            ('Ambient heat loss','W',self.Q_amb)
         ]
        
    def Calculate(self):
        #AbstractState
        if hasattr(self,'Backend'): #check if backend is given
            AS = CP.AbstractState(self.Backend, self.Ref)
        else: #otherwise, use the defualt backend
            AS = CP.AbstractState('HEOS', self.Ref)
        self.AS = AS
        
        #Local copies of coefficients
        P=self.P
        M=self.M
        R=self.R
        Te=self.Te
        
        #Calculate suction and injection dew temperatures
        AS.update(CP.PQ_INPUTS, self.pin_r, 0.0)
        h_bubble=AS.hmass() #[J/kg]
        
        AS.update(CP.PQ_INPUTS, self.pinj_r, 0.0)
        h_bubble_inj=AS.hmass() #[J/kg]
        
        AS.update(CP.PT_INPUTS, self.pinj_r, self.Tinj_r)
        self.hinj_r=AS.hmass() #[J/kg]
        self.sinj_r=AS.smass() #[J/kg-K]
        
        AS.update(CP.PT_INPUTS, self.pin_r, self.Tin_r)
        self.hin_r=AS.hmass() #[J/kg]
        self.sin_r=AS.smass() #[J/kg-K]
        self.vin_r = 1 / AS.rhomass() #[m**3/kg]
        
        #AS.update(CP.QT_INPUTS, 0.0, self.Tinj_r)
        AS.update(CP.PQ_INPUTS, self.pinj_r, 0.0)
        h_f_inj=AS.hmass()#[J/kg]
        #AS.update(CP.QT_INPUTS, 1.0, self.Tinj_r)
        AS.update(CP.PQ_INPUTS, self.pinj_r, 1.0)
        h_g_inj=AS.hmass()#[J/kg]
        
        #AS.update(CP.QT_INPUTS, 0.0, self.Tin_r)
        AS.update(CP.PQ_INPUTS, self.pin_r, 0.0)
        h_f=AS.hmass()#[J/kg]
        #AS.update(CP.QT_INPUTS, 1.0, self.Tin_r)
        AS.update(CP.PQ_INPUTS, self.pin_r, 1.0)
        h_g=AS.hmass()#[J/kg]
        
        #Injection:
        Dh_shinj_r=self.hinj_r-h_bubble_inj
        Dh_fg_inj=h_g_inj-h_f_inj
        #Suction
        Dh_sh=self.hin_r-h_bubble
        Dh_fg=h_g-h_f
        
        p_inj_norm = self.pinj_r/self.pin_r
        deltah_inj_norm = Dh_shinj_r/Dh_fg_inj
        p_dis_norm = self.pout_r/self.pin_r
        deltah_suc_norm = Dh_sh/Dh_fg
        
        #Power [Watts]
        #power = P[0] + P[1] * ((self.pinj_r/self.pin_r)**P[2]) * ((Dh_shinj_r/Dh_fg_inj)**P[3]) * ((self.pout_r/self.pin_r)**P[4]) * ((Dh_sh/Dh_fg)**P[5])
        #power=-221692.6142+226000.5146*p_dis_norm**(-0.02318188686)*p_inj_norm**(+0.02748851436)*deltah_inj_norm**(+0.1951974063)*deltah_suc_norm**(+0.03808133981)
        #power=+2255.855988+1.272237025*p_dis_norm**(+1.386084529)*p_inj_norm**(+0.4665091523)*deltah_suc_norm**(+108.6198468)*deltah_inj_norm**(-0.2734908995)
        #power=-486688.1319+490640.1614*p_dis_norm**(-0.01665460385)*p_inj_norm**(+0.01936780152)*deltah_suc_norm**(+0.03077386404)*deltah_inj_norm**(+0.1223111911)
        power=-319729.1256+321981.9721*p_dis_norm**(+0.003473718649)*p_inj_norm**(+0.001361941186)

        #suction mass flow rate [kg/s]
        #mdot = M[0] + M[1] * ((self.pinj_r/self.pin_r)**M[2]) * ((Dh_shinj_r/Dh_fg_inj)**M[3]) * ((self.pout_r/self.pin_r)**M[4]) * ((Dh_sh/Dh_fg)**M[5])
        #mdot=+32.03321447+1495.38647*p_dis_norm**(-2.323951449)*p_inj_norm**(+1.664814757)*deltah_inj_norm**(+11.87719284)*deltah_suc_norm**(+7.977699852)
        #mdot=+0.0006866353183+0.009703849695*p_dis_norm**(-0.3195095861)*p_inj_norm**(+0.06435925202)*deltah_suc_norm**(+49.42744463)*deltah_inj_norm**(+0.4097326195)
        #mdot=-0.001049098943+0.1684268596*p_dis_norm**(-2.803142126)*p_inj_norm**(+2.252610925)*deltah_suc_norm**(+9.139310451)*deltah_inj_norm**(+14.54322906)
        mdot=-49745.88214+50671.82075*p_dis_norm**(-0.005391951797)*p_inj_norm**(-0.0006888184216)

        #ratio of injection to suction mass flow rate [-]
        #ratio_mass = R[0] + R[1] * ((self.pinj_r/self.pin_r)**R[2]) * ((Dh_shinj_r/Dh_fg_inj)**R[3]) * ((self.pout_r/self.pin_r)**R[4]) * ((Dh_sh/Dh_fg)**R[5])
        #ratio_mass=-5.107639209+4.903278345*p_dis_norm**(+0.04064998309)*p_inj_norm**(+0.08813423426)*deltah_inj_norm**(-0.2250659425)*deltah_suc_norm**(+0.03886689634)
        #ratio_mass=-33.36701011+33.05046622*p_dis_norm**(+0.007937476525)*p_inj_norm**(+0.004692211409)*deltah_suc_norm**(+0.06654063672)*deltah_inj_norm**(-0.03190257993)
        #ratio_mass=-10.05302602+9.817948036*p_dis_norm**(-0.001821866536)*p_inj_norm**(+0.07187989147)*deltah_suc_norm**(+0.0681223862)*deltah_inj_norm**(+0.006131357799)
        ratio_mass=-0.2434501617+0.183719441*p_dis_norm**(+0.1700795496)*p_inj_norm**(+1.342340672)

        #injection mass flow rate [kg/s]
        mdot_inj = mdot * ratio_mass
        mdot = lbh2kgs(mdot)
        mdot_inj = lbh2kgs(mdot_inj)
        
        #Discharge tempearture [K]
        #T_dis = Te[0] + Te[1] * ((self.pinj_r/self.pin_r)**Te[2]) * ((Dh_shinj_r/Dh_fg_inj)**Te[3]) * ((self.pout_r/self.pin_r)**Te[4]) * ((Dh_sh/Dh_fg)**Te[5])
        #T_dis=+62.17255277+52.69594129*p_dis_norm**(-0.07656264047)*p_inj_norm**(+0.62997606)*deltah_inj_norm**(+4.96795673)*deltah_suc_norm**(+5.197706204)
        #T_dis=+321.7965732+0.6011743237*p_dis_norm**(+1.171991802)*p_inj_norm**(+0.3640605433)*deltah_suc_norm**(+45.74189473)*deltah_inj_norm**(+1.62244683)
        #T_dis=+294.5941967+24.21525236*p_dis_norm**(-0.3468596178)*p_inj_norm**(+0.953536384)*deltah_suc_norm**(+6.396931042)*deltah_inj_norm**(+6.811709684)
        T_dis=-12328.67021+12410.38063*p_dis_norm**(+0.005848322644)*p_inj_norm**(+0.0008722038865)

        self.Tout_r = F2K(T_dis) #[K]

        AS.update(CP.PT_INPUTS, self.pout_r, self.Tout_r)
        self.hout_r = AS.hmass() #[J/kg]
        self.sout_r = AS.smass() #[J/kg-K]        
        
        #define properites for isentropic efficency
        AS.update(CP.PSmass_INPUTS, self.pout_r, self.sin_r)
        h_2s=AS.hmass() #[J/kg]
        
        AS.update(CP.PSmass_INPUTS, self.pout_r, self.sinj_r)
        h_4s=AS.hmass() #[J/kg]
        
        #isentropic effeicincy defined by Groll
        self.eta_oi=(mdot*(h_2s-self.hin_r) + mdot_inj*(h_4s-self.hinj_r))/power
        
        self.mdot_r = mdot
        self.mdot_inj = mdot_inj
        self.mdot_tot = mdot + mdot_inj
        self.W=power
        self.CycleEnergyIn=power*(1-self.fp)
        self.Vdot_pumped= mdot*self.vin_r
        self.Q_amb=-self.fp*power
        
if __name__=='__main__':        
    import numpy as np
    Tin_dew = 273.15 #assume 5K superheat
    pin_r = PropsSI('P','T',Tin_dew,'Q',1,'R407C')
    pout_r = PropsSI('P','T',333.15,'Q',1,'R407C')
    pinj_r = np.sqrt(pin_r * pout_r)
    Tinj_dew = PropsSI('T','P',pinj_r,'Q',1,'R407C')
    
    kwds={
          'M':[-0.08135903196,0.0458453839,0.129432516,0.1325702995,-0.1485129565,9.144668476],
          'P':[2102.503527,0.6448139189,0.6561921643,0.1690042408,1.270243282,41.2967932],
          'R':[-22.45861341,22.05630275,0.007613286149,-0.01443163273,0.01145168107,0.057749393],
          'Te':[317.359575,0.6939083124,0.4131861045,0.6598603942,1.009091398,16.1910005],
          'Ref':'R407C',
          'Tin_r':Tin_dew+5,
          'pin_r':pin_r,
          'pout_r':pout_r,
          'pinj_r':pinj_r,
          'Tinj_r':Tinj_dew+5,
          'fp':0.15, #Fraction of electrical power lost as heat to ambient
          'Vdot_ratio': 1.0, #Displacement Scale factor
          'Backend':'HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
          }
    Comp=VICompressorClass(**kwds)
    Comp.Calculate()
    print Comp.W,'W'
    print Comp.Tout_r,'K'
    print Comp.mdot_r,'kg/s'
    print Comp.mdot_inj,'kg/s'
    print Comp.Q_amb, 'W'
    print Comp.eta_oi*100, '%'
    print pinj_r