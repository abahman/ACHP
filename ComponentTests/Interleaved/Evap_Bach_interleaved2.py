#3-circuit version of interleaved evaporator for RTU
from __future__ import division #Make integer 3/2 give 1.5 in python 2.x

from math import floor,ceil
from CoolProp.CoolProp import PropsSI
from FinCorrelations import FinInputs
from Christian_Evaporator import EvaporatorClass
import numpy as np, pylab
import copy
from CoolProp.Plots import PropertyPlot
from scipy.optimize import newton,fmin_slsqp
from ACHPTools import Write2CSV
from convert_units import *
import pylab as plt

#MultiCircuitEvaporator inherits things from the Evaporator base class
class MCE_N(EvaporatorClass):
    def __init__(self,**kwargs):
        self.__dict__.update(kwargs)
    def Update(self,**kwargs):
        self.__dict__.update(kwargs)

    def OutputList(self):
        """
            Return a list of parameters for this component for further output
            
            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        Output_list=[]
        num_evaps= self.num_evaps
        Output_list_i=[]
        for i in range(num_evaps):                                                                                                                              #Generate output list with all evaporators
            #Output_list_i.append([('Vol flow A'+' '+str(i),'m^3/s',float(self.EvapsA[i].Fins.Air.Vdot_ha)),('Vol flow B'+' '+str(i),'m^3/s',float(self.EvapsB[i].Fins.Air.Vdot_ha)),
            Output_list_i.append([('Vol flow A'+' '+str(i),'m^3/s',self.EvapsA[i].Fins.Air.Vdot_ha),('Vol flow B'+' '+str(i),'m^3/s',self.EvapsB[i].Fins.Air.Vdot_ha),
                                    ('Inlet DB_A'+' '+str(i),'K',self.EvapsA[i].Fins.Air.Tdb),('Inlet DB_B'+' '+str(i),'K',self.EvapsB[i].Fins.Air.Tdb),
                                        ('Inlet Air pressure_A'+' '+str(i),'kP_A',self.EvapsA[i].Fins.Air.p),    ('Inlet Air pressure_B'+' '+str(i),'kP_B',self.EvapsB[i].Fins.Air.p),
                                        ('Inlet Air Relative Humidity_A'+' '+str(i),'-',self.EvapsA[i].Fins.Air.RH),('Inlet Air Relative Humidity_B'+' '+str(i),'-',self.EvapsB[i].Fins.Air.RH),
                                        ('Tubes per bank_A'+' '+str(i),'-',self.EvapsA[i].Fins.Tubes.NTubes_per_bank), ('Tubes per bank_B'+' '+str(i),'-',self.EvapsB[i].Fins.Tubes.NTubes_per_bank),
                                        ('Number of banks_A'+' '+str(i),'-',self.EvapsA[i].Fins.Tubes.Nbank),('Number of banks_B'+' '+str(i),'-',self.EvapsB[i].Fins.Tubes.Nbank),
                                        ('Number circuits_A'+' '+str(i),'-',self.EvapsA[i].Fins.Tubes.Ncircuits),    ('Number circuits_B'+' '+str(i),'-',self.EvapsB[i].Fins.Tubes.Ncircuits),
                                        ('Length of tube_A'+' '+str(i),'m',self.EvapsA[i].Fins.Tubes.Ltube),    ('Length of tube_B'+' '+str(i),'m',self.EvapsB[i].Fins.Tubes.Ltube),
                                        ('Tube OD_A'+' '+str(i),'m',self.EvapsA[i].Fins.Tubes.OD),    ('Tube OD_B'+' '+str(i),'m',self.EvapsB[i].Fins.Tubes.OD),
                                        ('Tube ID_A'+' '+str(i),'m',self.EvapsA[i].Fins.Tubes.ID),    ('Tube ID_B'+' '+str(i),'m',self.EvapsB[i].Fins.Tubes.ID),
                                        ('Tube Long. Pitch_A'+' '+str(i),'m',self.EvapsA[i].Fins.Tubes.Pl),    ('Tube Long. Pitch_B'+' '+str(i),'m',self.EvapsB[i].Fins.Tubes.Pl),
                                        ('Tube Transverse Pitch_A'+' '+str(i),'m',self.EvapsA[i].Fins.Tubes.Pt),    ('Tube Transverse Pitch_B'+' '+str(i),'m',self.EvapsB[i].Fins.Tubes.Pt),
                                        ('Outlet superheat_A'+' '+str(i),'K',self.EvapsA[i].Tout_r-self.EvapsA[i].Tdew_r),    ('Outlet superheat_B'+' '+str(i),'K',self.EvapsB[i].Tout_r-self.EvapsB[i].Tdew_r),
                                        ('Fins per inch_A'+' '+str(i),'1/in',self.EvapsA[i].Fins.Fins.FPI),    ('Fins per inch_B'+' '+str(i),'1/in',self.EvapsB[i].Fins.Fins.FPI),
                                        ('Fin waviness pd_A'+' '+str(i),'m',self.EvapsA[i].Fins.Fins.Pd),    ('Fin waviness pd_B'+' '+str(i),'m',self.EvapsB[i].Fins.Fins.Pd),
                                        ('Fin waviness xf_A'+' '+str(i),'m',self.EvapsA[i].Fins.Fins.xf),    ('Fin waviness xf_B'+' '+str(i),'m',self.EvapsB[i].Fins.Fins.xf),
                                        ('Fin thickness_A'+' '+str(i),'m',self.EvapsA[i].Fins.Fins.t),    ('Fin thickness_B'+' '+str(i),'m',self.EvapsB[i].Fins.Fins.t),
                                        ('Fin Conductivity_A'+' '+str(i),'W/m-K',self.EvapsA[i].Fins.Fins.k_fin),    ('Fin Conductivity_B'+' '+str(i),'W/m-K',self.EvapsB[i].Fins.Fins.k_fin),
                                        ('Refrigerant flowrate_A'+' '+str(i),'kg/s',self.EvapsA[i].mdot_r),    ('Refrigerant flowrate_B'+' '+str(i),'kg/s',self.EvapsB[i].mdot_r),
                                        ('Q Total_A'+' '+str(i),'W',self.EvapsA[i].Q),    ('Q Total_B'+' '+str(i),'W',self.EvapsB[i].Q),
                                        ('Q Superheat_A'+' '+str(i),'W',self.EvapsA[i].Q_superheat),    ('Q Superheat_B'+' '+str(i),'W',self.EvapsB[i].Q_superheat),
                                        ('Q Two-Phase_A'+' '+str(i),'W',self.EvapsA[i].Q_2phase),    ('Q Two-Phase_B'+' '+str(i),'W',self.EvapsB[i].Q_2phase),
                                        ('Inlet ref. temp_A'+' '+str(i),'K',self.EvapsA[i].Tin_r),    ('Inlet ref. temp_B'+' '+str(i),'K',self.EvapsB[i].Tin_r),
                                        ('Outlet ref. temp_A'+' '+str(i),'K',self.EvapsA[i].Tout_r),    ('Outlet ref. temp_B'+' '+str(i),'K',self.EvapsB[i].Tout_r),
                                        ('Outlet air temp_A'+' '+str(i),'K',self.EvapsA[i].Tout_a),    ('Outlet air temp_B'+' '+str(i),'K',self.EvapsB[i].Tout_a),
                                        ('Pressure Drop Total_A'+' '+str(i),'P_A',self.EvapsA[i].DP_r),    ('Pressure Drop Total_B'+' '+str(i),'P_B',self.EvapsB[i].DP_r),
                                        ('Pressure Drop Superheat_A'+' '+str(i),'P_A',self.EvapsA[i].DP_r_superheat),    ('Pressure Drop Superheat_B'+' '+str(i),'P_B',self.EvapsB[i].DP_r_superheat),
                                        ('Pressure Drop Two-Phase_A'+' '+str(i),'P_A',self.EvapsA[i].DP_r_2phase),    ('Pressure Drop Two-Phase_B'+' '+str(i),'P_B',self.EvapsB[i].DP_r_2phase),
                                        ('Charge Total_A'+' '+str(i),'kg',self.EvapsA[i].Charge),    ('Charge Total_B'+' '+str(i),'kg',self.EvapsB[i].Charge),
                                        ('Charge Superheat_A'+' '+str(i),'kg',self.EvapsA[i].Charge_superheat),    ('Charge Superheat_B'+' '+str(i),'kg',self.EvapsB[i].Charge_superheat),
                                        ('Charge Two-Phase_A'+' '+str(i),'kg',self.EvapsA[i].Charge_2phase),    ('Charge Two-Phase_B'+' '+str(i),'kg',self.EvapsB[i].Charge_2phase),
                                        ('Mean HTC Superheat_A'+' '+str(i),'W/m^2-K',self.EvapsA[i].h_r_superheat),    ('Mean HTC Superheat_B'+' '+str(i),'W/m^2-K',self.EvapsB[i].h_r_superheat),
                                        ('Mean HTC Two-phase_A'+' '+str(i),'W/m^2-K',self.EvapsA[i].h_r_2phase),    ('Mean HTC Two-phase_B'+' '+str(i),'W/m^2-K',self.EvapsB[i].h_r_2phase),
                                        ('Wetted Area Fraction Superheat_A'+' '+str(i),'-',self.EvapsA[i].w_superheat),    ('Wetted Area Fraction Superheat_B'+' '+str(i),'-',self.EvapsB[i].w_superheat),
                                        ('Wetted Area Fraction Two-phase_A'+' '+str(i),'-',self.EvapsA[i].w_2phase),    ('Wetted Area Fraction Two-phase_B'+' '+str(i),'-',self.EvapsB[i].w_2phase),
                                        ('Mean Air HTC_A'+' '+str(i),'W/m^2-K',self.EvapsA[i].Fins.h_a),    ('Mean Air HTC_B'+' '+str(i),'W/m^2-K',self.EvapsB[i].Fins.h_a),
                                        ('Surface Effectiveness_A'+' '+str(i),'-',self.EvapsA[i].Fins.eta_a),    ('Surface Effectiveness_B'+' '+str(i),'-',self.EvapsB[i].Fins.eta_a),
                                        ('Air-side area (fin+tubes)_A'+' '+str(i),'m^2',self.EvapsA[i].Fins.A_a),    ('Air-side area (fin+tubes)_B'+' '+str(i),'m^2',self.EvapsB[i].Fins.A_a),
                                        ('Mass Flow rate dry Air_A'+' '+str(i),'kg/s',self.EvapsA[i].Fins.mdot_da),    ('Mass Flow rate dry Air_B'+' '+str(i),'kg/s',self.EvapsB[i].Fins.mdot_da),
                                        ('Mass Flow rate humid Air_A'+' '+str(i),'kg/s',self.EvapsA[i].Fins.mdot_ha),    ('Mass Flow rate humid Air_B'+' '+str(i),'kg/s',self.EvapsB[i].Fins.mdot_ha),
                                        ('Pressure Drop Air-side_A'+' '+str(i),'P_A',self.EvapsA[i].Fins.dP_a),('Pressure Drop Air-side_B'+' '+str(i),'P_B',self.EvapsB[i].Fins.dP_a),
                                        ('Refrigerant mass flow A'+' '+str(i),'kg/s',self.EvapsA[i].mdot_r),('Refrigerant mass flow B'+' '+str(i),'kg/s',self.EvapsB[i].mdot_r),
                                        ('Sensible Heat Ratio_A'+' '+str(i),'-',self.EvapsA[i].SHR),('Sensible Heat Ratio_B'+' '+str(i),'-',self.EvapsB[i].SHR)])
        for i in range(0,len(Output_list_i[0]),1): #sort output list, such that corresponding values are next to each other
            sumsi=0    #need sums and avgs
            for n in range(0,num_evaps):
                Output_list.append(Output_list_i[n][i])
                try:
                    if type(Output_list_i[n][i][2]) is not type("string"):
                        sumsi+=Output_list_i[n][i][2] #sum up values
                        if n == num_evaps-1:
                            Output_list.append((Output_list_i[n][i][0][:-2]+" sum",Output_list_i[n][i][1],sumsi))
                            Output_list.append((Output_list_i[n][i][0][:-2]+" avg",Output_list_i[n][i][1],sumsi/num_evaps))
                except:
                    print "something wrong in Output List, wrong size of element?"
                    print Output_list_i[n][i],Output_list_i[n][i][2],i,type(1.0)
                    print type(Output_list_i[n][i][2])
                    raise
        Output_List_tot=[]
        #append optional parameters, if applicable
        if hasattr(self,'TestName'):
            Output_List_tot.append(('Name','N/A',self.TestName)) 
        if hasattr(self,'TestDescription'):
            Output_List_tot.append(('Description','N/A',self.TestDescription))
        if hasattr(self,'Details'):
            Output_List_tot.append(('Details','N/A',self.Details))
        if hasattr(self,'md_severity'):
            Output_List_tot.append(('Maldistribution Severity','N/A',self.md_severity))
        if hasattr(self,'resids'):
            Output_List_tot.append(('Residuals','N/A',self.resids))
        if hasattr(self,'resid_eq_sh'):
            Output_List_tot.append(('Residuals_SH_equal','N/A',self.resid_eq_sh))
        else:
            Output_List_tot.append(('Residuals_SH_equal','N/A','0'))
        if hasattr(self,'maldistributed'):
            try: 
                len(self.maldistributed) #if it is a bool (=False), then we can't calculate anything
                m_dot_bar=np.average(self.maldistributed)
                m_dot_i_minus_m_dot_ave=self.maldistributed-m_dot_bar
                sum_squares=np.sum(np.multiply(m_dot_i_minus_m_dot_ave,m_dot_i_minus_m_dot_ave))
                self.md_Epsilon=np.sqrt(sum_squares/len(self.maldistributed))/m_dot_bar
            except:
                self.md_Epsilon=-999
        if hasattr(self,'md_Epsilon'):
            Output_List_tot.append(('Epsilon_md','-',self.md_Epsilon))
        else:
            Output_List_tot.append(('Epsilon_md','-','-999'))
        Output_List_tot=Output_List_tot+[('Q Total','W',self.Q),('m_dot Total','kg/s',self.mdot_r_tot),('h_out Total','kJ/kg',self.hout_r/1000.0)]
        for i in range(0,len(Output_list)):                             #append default parameters to output list
            Output_List_tot.append(Output_list[i])
        return Output_List_tot
    
    def Evaporator_LRCS_Fins(self):
        #define parameters for evaporator as used in LRCS
        Evaporator=EvaporatorClass()
        Evaporator.Fins=FinInputs()
         
        #--------------------------------------
        #--------------------------------------
        #           Evaporator
        #           -> see Condenser and GUI for explanations
        #--------------------------------------
        #--------------------------------------
        Evaporator.Fins.Tubes.NTubes_per_bank=1#8 (each cell 1 tube)
        Evaporator.Fins.Tubes.Nbank=2#4(half of actual number for a single cell)
        Evaporator.Fins.Tubes.Ncircuits=1#8 (each cell is part of 1 circuit)
        Evaporator.Fins.Tubes.Ltube=in2m(108) #measured fin pack length
        Evaporator.Fins.Tubes.OD=0.01307 #measured
        Evaporator.Fins.Tubes.ID=0.01307-0.001 #guess of 1 mm for wall thickness
        Evaporator.Fins.Tubes.Pl=0.03303       #distance between center of tubes in flow direction (measured)
        Evaporator.Fins.Tubes.Pt=0.02285
             
        Evaporator.Fins.Fins.FPI=8
        Evaporator.Fins.Fins.Pd=0.00065  #fins are basically flat; measured Pd in wrong direction (wavyness perpendicular to airflow direction)
        Evaporator.Fins.Fins.xf=0.0031 
        Evaporator.Fins.Fins.t=0.000198   #tuned; measurement with callipper, confirmed withmicrometer screw (0.0078inch=0.00019812m)
        Evaporator.Fins.Fins.k_fin=237 #Thermal conductivity of fin material, aluminum, from wikipedia (replace with other source)
         
        Evaporator.Fins.Air.Vdot_ha=(1/8)*cfm2cms(4440.0)*0.67 #4440rated cfm >set manually in liquid_receiver_cycle
        #Evaporator.Fins.Air.Tmean=C2K(2.0)  #this is not actually used
        Evaporator.Fins.Air.Tdb=C2K(2.0)
        Evaporator.Fins.Air.p=101325      #Air pressure
        #################################
        Evaporator.Fins.Air.RH=0.48  #0.48
        Evaporator.Fins.Air.RHmean=0.48 #0.48
        #################################
        Evaporator.Fins.Air.FanPower=0#378  #W, average from clean coil hybrid measurements
        return Evaporator.Fins
    
    def Evaporator_RAC_Fins(self):
        #define parameters for evaporator as used in RAC
        
        FinsTubes=FinInputs()

        #RAC has not as nice of a circuitry as the LRCS has,
        #therefore need to define 2 banks (instead of 3) with total
        #of 14 tubes, e.g. each cell 1 bank with 7 tubes

        FinsTubes.Tubes.NTubes_per_bank=7
        FinsTubes.Tubes.Ncircuits=1#7
        FinsTubes.Tubes.Nbank=1
        
        FinsTubes.Tubes.Ltube=in2m(33)
        FinsTubes.Tubes.OD=in2m(0.375)  #given by manufacturer
        FinsTubes.Tubes.ID=FinsTubes.Tubes.OD-2*in2m(0.0118)  #manufacturer gave wall thickness
        FinsTubes.Tubes.Pl=in2m(1.5/2.0)   #distance between tubes in air flow direction (measured distance between 3 tubes)
        FinsTubes.Tubes.Pt=in2m(1.0)        #distance between center of tubes orthogonal to flow direction
        
        FinsTubes.Fins.FPI=15     #given by manufacturer
        FinsTubes.Fins.Pd=0.001#not measurable without dismantling unit; assumed standard value
        FinsTubes.Fins.xf=0.001 #not measurable without dismantling unit; assumed standard value
        FinsTubes.Fins.t=0.0001 #measured
        FinsTubes.Fins.k_fin=237 #aluminum fins
        
        FinsTubes.Air.Vdot_ha=cfm2cms(1174)*(1/7) #125F-80F-41.2SF-50Damp, hybrid
        #FinsTubes.Air.Tmean=C2K(36.7)  #125F-80F-41.2SF-50Damp
        FinsTubes.Air.Tdb=C2K(36.6)#125F-80F-41.2SF-50Damp, hybrid
        FinsTubes.Air.p=101325
        FinsTubes.Air.RH=0.2  #keep dry fins
        FinsTubes.Air.RHmean=0.2
        FinsTubes.Air.FanPower=0#250 #manufacturer rating
        
        return FinsTubes

    def Calculate(self,evap_type='LRCS'):
        #common inputs; note: flowrates are "per circuit", not total!
        if not hasattr(self,'num_evaps'):
            self.num_evaps=2 #number of evaporators
        if evap_type=='LRCS':
            self.Ref='R404a'
            self.psat_r=445100
            if hasattr(self,'mdot_r'):
                self.mdot_r=self.mdot_r/float(self.num_evaps) #internally using individual circuit average flowrate
            else:
                self.mdot_r=(86.57/1000.0)/(8.0)*1.2   #later on add handling to automatically get back to flowrate of one circuit from total flowrate
            self.mdot_r_=self.mdot_r*1.0   #used as backup if first value in superheat iteration does not converge
            self.hin_r=PropsSI('H','P', 1429000,'Q',0,self.Ref)
            self.Verbosity=0
            self.cp_r_iter=False  #iterate for CP in evaporator
            self.FinsType = 'WavyLouveredFins'
            
        elif evap_type=='RAC':
            self.Ref='R410a'
            #Values from 125F-80F-41.2SF-50Damp hybrid test
            self.psat_r=1414000
            self.mdot_r=(68.0/1000.0)/(7.0)*1.0   #later on add handling to automatically get back to flowrate of one circuit from total flowrate
            self.mdot_r_=self.mdot_r*1.0   #used as backup if first value in superheat iteration does not converge
            self.hin_r=PropsSI('H','P', 3687000,'T',C2K(52.7),self.Ref)
            self.Verbosity=0
            self.cp_r_iter=False  #iterate for CP in evaporator
            self.FinsType = 'WavyLouveredFins'
        
        else:
            print "undefined evaporator type"
            raise()
        #self.interleaved=False  #this needs to be set by calling function!

        
        if not hasattr(self,'EvapsA'):   #if we don't already have run the calculate function once
            #generate dictionaries for evaporators
            EvapDict=self.__dict__
            ED=copy.deepcopy(EvapDict)
            self.EvapsA=[]  #first row at air inlet
            self.EvapsB=[]  #second row, air inlet is air outlet from first row
            #for both rows, [0] is on top, [1] is on bottom
            for i in range(self.num_evaps):
                #Create new evaporator instanciated with new deep copied dictionary
                ####first row
                ED2=copy.deepcopy(ED)
                E=EvaporatorClass(**ED2)
                if evap_type=='LRCS':
                    dict_tmp=copy.deepcopy(self.Evaporator_LRCS_Fins())
                elif evap_type=='RAC':
                    dict_tmp=copy.deepcopy(self.Evaporator_RAC_Fins())
                    E.h_a_tuning=1.0#1.4 was used for other RAC simulations
                else:
                    print "undefined evaporator type"
                    raise()
                E.Fins=dict_tmp
                self.Vdot_ha=dict_tmp.Air.Vdot_ha*1.0  #save air flowrate in main structure
                self.Tdb=dict_tmp.Air.Tdb*1.0  #save air temperature in main structure
                #Add to list of evaporators
                self.EvapsA.append(E)
                ####second row
                ED2=copy.deepcopy(ED)
                E=EvaporatorClass(**ED2)
                if evap_type=='LRCS':
                    dict_tmp=copy.deepcopy(self.Evaporator_LRCS_Fins())
                elif evap_type=='RAC':
                    dict_tmp=copy.deepcopy(self.Evaporator_RAC_Fins())
                    E.h_a_tuning=1.0
                else:
                    print "undefined evaporator type"
                    raise()
                E.Fins=dict_tmp
                #Add to list of evaporators
                self.EvapsB.append(E)

        if hasattr(self,'maldistributed'):    
            #apply airside FLOW maldistribution
            try: 
                float(self.maldistributed[i])
                air_flow_rat=self.maldistributed
            except: #no airside maldistribution
                air_flow_rat=np.linspace(1.0,1.0,self.num_evaps)
                print "invalid vector for aiside maldistribution, proceeding without maldistribution"
            print "using air flow maldistribution, volumetric version",air_flow_rat,"air_flow_rat"
            for i in range(self.num_evaps):
                self.EvapsA[i].Fins.Air.Vdot_ha=self.Vdot_ha*air_flow_rat[i]
                self.EvapsB[i].Fins.Air.Vdot_ha=self.Vdot_ha*air_flow_rat[i]
                if hasattr(self,'interleaved'):
                    #adjust air flowrates to second row (is interleaved)
                    if self.interleaved==True:
                        self.EvapsB[i].Fins.Air.Vdot_ha=self.Vdot_ha*air_flow_rat[self.num_evaps-1-i]
                print "circuit",i,air_flow_rat[i],"circ Ai",self.EvapsB[i].Fins.Air.Vdot_ha,"circ Bi",self.EvapsA[i].Fins.Air.Vdot_ha
            print " " #newline
        if hasattr(self,'air_temp_maldistributed'):    
            #apply airside temperature maldistribution
            #50% from nominal airflowrate as default value flowrate on top circuit
            try: 
                float(self.air_temp_maldistributed[i])
                air_temp_change=self.air_temp_maldistributed
            except:
                air_temp_change=np.linspace(-0.5,+0.5,self.num_evaps)  #K
                print "applying 0.5K  maldistribution on airside inlet temperature "
            for i in range(self.num_evaps):
                self.EvapsB[i].Fins.Air.Tdb=self.Tdb + air_temp_change[i] #second row actually gets inlet temperature during solving process...
                self.EvapsA[i].Fins.Air.Tdb=self.Tdb +air_temp_change[i] 
            print "###applied airside temperatuere maldistribution of",air_temp_change
        
        def adjust_flowrate_EXV(): 
            #apply refrigerant side maldistribution
            #function to adjust flowrates for EXV control
            #adjust refrigerant flowrates according to equal flow assumption or maldistribution
            #parallel flow for normal EXV
            try: 
                float(self.ref_maldistributed[self.num_evaps-1])
                ref_flow_rat=self.ref_maldistributed
                print "applying maldistribution on refrigerant side"
            except:
                ref_flow_rat=np.ones(self.num_evaps)
            for i in range(self.num_evaps):
                self.EvapsA[i].mdot_r=self.mdot_r*ref_flow_rat[i]
                self.EvapsB[i].mdot_r=self.mdot_r*ref_flow_rat[i]
            
        def adjust_flowrates():  #function to adjust flowrates for hybrid control
            if hasattr(self,'Hybrid'):
                if self.Hybrid=='equal_flow':
                    v_dot_avg=0.0
                    for i in range(self.num_evaps):
                        v_dot_avg+=self.EvapsA[i].Fins.Air.Vdot_ha
                    v_dot_avg/=(1.0*self.num_evaps)  #average circuit inlet flowrate
                    #adjust refrigerant flowrates according to air flowrates
                    for i in range(self.num_evaps):
                        self.EvapsA[i].mdot_r=self.mdot_r*self.EvapsA[i].Fins.Air.Vdot_ha/v_dot_avg
                        self.EvapsB[i].mdot_r=self.mdot_r*self.EvapsB[i].Fins.Air.Vdot_ha/v_dot_avg
                elif self.Hybrid=='equal_flow_DT_in':
                    Tsat_r=PropsSI('T','P',self.psat_r,'Q',1.0,self.Ref)  #neglect temperature dependency
                    DT_avg=0.0
                    for i in range(self.num_evaps):
                        DT_avg+=self.EvapsA[i].Fins.Air.Tdb
                    DT_avg=(DT_avg)/(1.0*self.num_evaps)-Tsat_r
                    #adjust refrigerant flowrates according to air inlet temperature difference
                    for i in range(self.num_evaps):
                        self.EvapsA[i].mdot_r=self.mdot_r*(self.EvapsA[i].Fins.Air.Tdb-Tsat_r)/DT_avg
                        self.EvapsB[i].mdot_r=self.mdot_r*(self.EvapsA[i].Fins.Air.Tdb-Tsat_r)/DT_avg  #parallel circuitry
                elif self.Hybrid=='adjust_superheat' or self.Hybrid=='adjust_area_fraction' or 'adjust_superheat_iter':
                        if not hasattr(self,'Hybrid_ref_distribution'):
                            raise()
                        for i in range(self.num_evaps):
                            self.EvapsA[i].mdot_r=self.mdot_r*self.Hybrid_ref_distribution[i]
                            self.EvapsB[i].mdot_r=self.mdot_r*self.Hybrid_ref_distribution[i]  #parallel circuitry
                else:
                    adjust_flowrate_EXV()  #already done at first call of Calculate, but needs to be repeated for calculation of mass flowrate with given target SH
                    print "Wrong input for self.Hybrid - temporarily using EXV adjustment for hybrid", self.Hybrid
            else:
                adjust_flowrate_EXV()  #already done at first call of Calculate, but needs to be repeated for calculation of mass flowrate with given target SH
                
        adjust_flowrates()  #run once at startup, needed for calculation of guess values
        
        #self.iterationloopnum=0

        def residual(hin_rA):
            #objective function to calculate residual of evaporators
            #self.iterationloopnum+=1
            #print self.iterationloopnum, "self.iterationloopnum"
            for i in range(self.num_evaps):  #update and calculate first row
                #print "A - evaps"
                if self.Verbosity: print "residual evap calcs first coil sheet",i
                if hin_rA[i]<PropsSI('H','Q',0.0,'P',self.psat_r,self.Ref):
                    print "something wrong with inlet enthalpy, too small- is",hin_rA[i]/1000.0,i,"but saturated liquid would be",PropsSI('H','Q',0.0,'P',self.psat_r,self.Ref),
                    print "fixing it to slight amount of quality to allow solver to proceed"
                    hin_rA[i]=PropsSI('H','Q',0.01,'P',self.psat_r,self.Ref)
                if hin_rA[i]>PropsSI('H','T',self.EvapsA[i].Fins.Air.Tdb,'P',self.psat_r,self.Ref):
                    print "something wrong with inlet enthalpy, too large"
                    print "hin_rA[i]",hin_rA[i],i
                    hin_rA[i]=PropsSI('H','T',self.EvapsA[i].Fins.Air.Tdb-0.1,'P',self.psat_r,self.Ref)
                    print "limited to 0.1K les than the air inlet temperature"
                self.EvapsA[i].hin_r=hin_rA[i]
                self.EvapsA[i].Calculate()
                if self.Verbosity: print "first coil sheet",i
            
            hB_out_for_residue=np.zeros(self.num_evaps)
            for i in range(self.num_evaps):   #update and calculate second row
                #print "B - evaps"
                #if self.iterationloopnum==7:
                #    self.EvapsB[0].plotit=True
                if self.Verbosity: print "residual evap calcs second coil sheet",i
                if self.interleaved:
                    self.EvapsB[i].Fins.Air.RH= self.EvapsA[self.num_evaps-1-i].Fins.Air.RH_out
                    self.EvapsB[i].Fins.Air.Tdb= self.EvapsA[self.num_evaps-1-i].Tout_a
                else:
                    self.EvapsB[i].Fins.Air.RH= self.EvapsA[i].Fins.Air.RH_out
                    self.EvapsB[i].Fins.Air.Tdb= self.EvapsA[i].Tout_a
                self.EvapsB[i].Calculate()
                hB_out_for_residue[i]=self.EvapsB[i].hout_r
            #calculate the error between estimated and actual inlet enthalpy to first row
            residue=hin_rA-hB_out_for_residue
            self.resids=residue
            return residue
        
        print self.EvapsA[0].Fins.Air.Tdb
        h_guess_max=PropsSI('H','P',self.psat_r,'T',self.EvapsA[0].Fins.Air.Tdb,self.Ref)-5.0
        guess_value=300000*np.ones(self.num_evaps)
        guess_value=1000.0*h_guess_max**np.ones(self.num_evaps)
        print self.EvapsA[0].Fins.Air.Tdb,h_guess_max,guess_value
        print residual(guess_value)
        
        def solve_for_exit_sh(self):
            "solve for the mass flow rate for a given target super-heat"
            Target_SH=np.float(self.Target_SH)  #check if it is a float
            #import solver and solve
            from scipy.optimize import fsolve
            def objective_SH_out(mdot_guess):
                print "mdot_guess",mdot_guess
                if mdot_guess<0.001:
                    print 'warning - mass flowrate has a negative value during solving process - constrained to 0.002'
                    mdot_guess=[0.001]
                self.mdot_r=mdot_guess[0]
                adjust_flowrates()  #adjust flowrates for hybrid or control
                fsolve(residual, guess_value)  #actually solve evaporator
                self.mdot_r_tot=0.0
                self.hout_r=0.0
                for i in range(self.num_evaps):
                    self.mdot_r_tot+=self.EvapsA[i].mdot_r
                    self.hout_r+=self.EvapsA[i].mdot_r*self.EvapsA[i].hout_r
                self.hout_r/=self.mdot_r_tot
                self.resids=self.hout_r-self.hout_r_target #store nested for csv output
                return self.hout_r-self.hout_r_target
            T_sat=PropsSI('T','P',self.psat_r,'Q',1.0,self.Ref)
            self.hout_r_target=PropsSI('H','T',self.Target_SH+T_sat,'P',self.psat_r,self.Ref)
            fsolve(objective_SH_out, self.mdot_r*1.0)  #note - flowrate is based on single evaporator up to here; mdot_r used as initial guess
            print "residual for superheat is",self.resids
            if self.resids>0.5:
                fsolve(objective_SH_out, self.mdot_r_*0.8)  #backup massflow
            if self.resids<-0.5:
                fsolve(objective_SH_out, self.mdot_r_*1.2)  #backup massflow
            print "result of fsolve is that target superheat is ",np.float(self.Target_SH), "actual value is ", self.EvapsA[0].Tout_r-T_sat,self.EvapsA[1].Tout_r-T_sat
            
        if hasattr(self,'Target_SH'):
            if hasattr(self,'Hybrid'):
                if not self.Hybrid=='adjust_superheat':
                    if not self.Hybrid=='adjust_area_fraction':
                        print "calculating exit superheat, since neither ajust_sh option nor adjust area_fraction"
                        solve_for_exit_sh(self) #otherwise we can skip this calculation
            else:
                solve_for_exit_sh(self)
        else:  #just calculate the normal output values
            #import solver and solve
            from scipy.optimize import fsolve
            actual_hinA=fsolve(residual, guess_value)
        
        #equalize exit superheat for hybrid control, if applicable
        def EQ_SH_OBJECTIVE(x):
            #objective is equalization of exit superheats
            self.Hybrid_ref_distribution=x #update distribution
            print x
            solve_for_exit_sh(self)
            #calculate residual is the difference of the individual circuit exit enthalpies
            resid_eq_sh=np.zeros(self.num_evaps)
            for i in range(self.num_evaps):
                resid_eq_sh[i]=self.EvapsA[i].hout_r/10000.0 #bring on similar scale as inputs
            resid_eq_sh=np.sum(resid_eq_sh*resid_eq_sh)
            self. resid_eq_sh=resid_eq_sh
            return resid_eq_sh

        if hasattr(self,'Hybrid'):
            if self.Hybrid=='adjust_superheat':
                if not hasattr(self,'Hybrid_ref_distribution'):
                    #use equal flow distribution on refrigerant and airside
                    self.Hybrid_ref_distribution=np.zeros(self.num_evaps)
                    v_dot_avg=0.0
                    for i in range(self.num_evaps):
                        v_dot_avg+=self.EvapsA[i].Fins.Air.Vdot_ha
                    v_dot_avg/=(1.0*self.num_evaps)  #average circuit inlet flowrate
                    #adjust refrigerant flowrates according to air flowrates
                    for i in range(self.num_evaps):
                        self.Hybrid_ref_distribution[i]=self.EvapsA[i].Fins.Air.Vdot_ha/v_dot_avg
                #Starting guess for mass flow weighting parameters are air volume flow weighting parameters
                print "warning, the Hybrid_ref_distribution does not seem to work properly"
                x0=self.Hybrid_ref_distribution
                print "x0=self.Hybrid_ref_distribution",self.Hybrid_ref_distribution
                #Actually solve constrained optimization problem to obtain equal exit superheats
                Bounds=[]; 
                B_min=np.min(x0*0.5)
                B_max=np.max(x0*1.5)
                for n in range(0,len(x0)): Bounds.append((B_min,B_max)) #create the solver-boundaries for the flowrates
                x=fmin_slsqp(EQ_SH_OBJECTIVE,x0,eqcons=[lambda x: np.sum(x)-len(x0)],bounds=Bounds)
                print "residual for fining equal exit superheats is",self. resid_eq_sh
            elif self.Hybrid=='adjust_area_fraction':
                print "self.Hybrid  ia adjust_area_fraction"
                #Iteratively find flowrates for same area fraction
                #use equal flow distribution on refrigerant and airside
                self.Hybrid_ref_distribution=np.zeros(self.num_evaps)
                v_dot_avg=0.0
                for i in range(self.num_evaps):
                    v_dot_avg+=self.EvapsA[i].Fins.Air.Vdot_ha
                v_dot_avg/=(1.0*self.num_evaps)  #average circuit inlet flowrate
                #adjust refrigerant flowrates according to air flowrates
                for i in range(self.num_evaps):
                    self.Hybrid_ref_distribution[i]=self.EvapsA[i].Fins.Air.Vdot_ha/v_dot_avg
                #solve for exit superheat
                solve_for_exit_sh(self)
                #iterate
                if not hasattr(self,'adjust_area_fraction_iternum'):
                    self.adjust_area_fraction_iternum=1 #do one iteration
                for i in range(self.adjust_area_fraction_iternum):
                    print "in pseudo I-control loop, iteration",i
                    w_evaps=np.zeros(self.num_evaps)
                    for i in range(self.num_evaps): w_evaps[i]=self.EvapsA[i].w_2phase+self.EvapsB[i].w_2phase
                    w_evaps_ave=np.average(w_evaps)
                    for i in range(self.num_evaps): self.Hybrid_ref_distribution[i]=(0.5+0.5*(w_evaps_ave-w_evaps[i])/w_evaps_ave)*self.Hybrid_ref_distribution[i]  #using 50% "integration" value as pseudo I-control
                    self.Hybrid_ref_distribution=self.Hybrid_ref_distribution/np.sum(self.Hybrid_ref_distribution)*float(self.num_evaps)
                    solve_for_exit_sh(self) #recalculate
            elif self.Hybrid=='adjust_superheat_iter':
                print "self.Hybrid  ia adjust_area_fraction"
                #Iteratively find flowrates for same exit superheat
                #use equal flow distribution on refrigerant and airside
                self.Hybrid_ref_distribution=np.zeros(self.num_evaps)
                v_dot_avg=0.0
                for i in range(self.num_evaps):
                    v_dot_avg+=self.EvapsA[i].Fins.Air.Vdot_ha
                v_dot_avg/=(1.0*self.num_evaps)  #average circuit inlet flowrate
                #adjust refrigerant flowrates according to air flowrates
                for i in range(self.num_evaps):
                    self.Hybrid_ref_distribution[i]=self.EvapsA[i].Fins.Air.Vdot_ha/v_dot_avg
                #solve for exit superheat
                solve_for_exit_sh(self)
                #iterate
                if not hasattr(self,'adjust_area_fraction_iternum'):
                    self.adjust_area_fraction_iternum=1 #do one iteration
                for i in range(self.adjust_area_fraction_iternum):
                    print "in pseudo I-control loop, iteration",i
                    SH_evaps=np.zeros(self.num_evaps)
                    for i in range(self.num_evaps): SH_evaps[i]=self.EvapsA[i].Tout_r-self.EvapsA[i].Tdew_r
                    SH_evaps_ave=np.average(SH_evaps)
                    for i in range(self.num_evaps): 
                        tmp=self.Hybrid_ref_distribution[i] #check if it  is an ovberwriting issue
                        self.Hybrid_ref_distribution[i]=tmp+0.07*(SH_evaps[i]-SH_evaps_ave)/SH_evaps_ave*tmp  #using 7% "integration" value as pseudo I-control
                    self.Hybrid_ref_distribution=self.Hybrid_ref_distribution/np.sum(self.Hybrid_ref_distribution)*float(self.num_evaps)
                    solve_for_exit_sh(self) #recalculate

            else:
                #already covered by previous calculations for equal flow and DT in
                print "maybe wrong input for Hybrid?"

        ##Calculate overall outputs
        self.Q=0.0
        self.mdot_r_tot=0.0
        self.mdot_r_totB=0.0  #second coil sheet, as a check
        self.hout_r=0.0
        for i in range(self.num_evaps):
            self.Q+=self.EvapsA[i].Q+self.EvapsB[i].Q
            self.mdot_r_tot+=self.EvapsA[i].mdot_r
            self.mdot_r_totB+=self.EvapsB[i].mdot_r
            self.hout_r+=self.EvapsA[i].mdot_r*self.EvapsA[i].hout_r  #flow weighted average
        self.hout_r/=self.mdot_r_tot
        
        print "flowrate at first and second coil sheet is ",self.mdot_r_tot,self.mdot_r_totB,"relative error is", (self.mdot_r_tot-self.mdot_r_totB)/self.mdot_r_tot

        if False:  #plot results
            from CoolProp.Plots import PropertyPlot
            fig = plt.figure()
            PropertyPlot(self.Ref,'PH')
            ax = fig.add_subplot(111)
            if self.interleaved==True:
                subplottitleis="interleaved"
            else:
                try:
                    self.Hybrid
                    subplottitleis="non-interleaved"+' '+self.Hybrid
                except:
                    subplottitleis="non-interleaved"
            try: 
                subplottitleis=subplottitleis+' air '+str(self.maldistributed)
            except:
                try:
                    subplottitleis=subplottitleis+' ref '+str(self.ref_maldistributed)
                except:
                    subplottitleis=subplottitleis
            ax.set_title(subplottitleis)
            ax.set_yscale('log')
            ax.plot(guess_value/1000.0,self.psat_r*np.ones(self.num_evaps),'x',label='guess hin_rA')
            actualH=[]
            linestiles=['-',':','--','-.','-',':','--','-.']
            markerstiles=['D','o','S','p','h','*','D','o','S','p','h','*']
            colors=['r','b','g','y','b','g','y','r']
            #linestyle=linestiles[i]
            for i in range(self.num_evaps):
                actualH.append(np.array([self.EvapsB[i].hin_r,self.EvapsA[i].hin_r,self.EvapsA[i].hout_r]))
                ax.scatter(actualH[i]/1000.0,self.psat_r*np.ones(3),s=np.linspace(5,50,3),label='actual enthalpy branch'+str(i),linestyle='-',c=colors[i],alpha=0.8)
            ax.plot(self.hout_r/1000.0*np.ones(1),self.psat_r*np.ones(1),'x',markerfacecolor='None',color='black',markersize=15,label='outlet enthalpy')
            print self.psat_r*np.ones(1)
            print np.array([PropsSI('H','P',self.psat_r,'T',self.EvapsA[0].Fins.Air.Tdb,self.Ref)]),self.EvapsA[0].hout_r/1000.0,"temps",self.EvapsA[0].Tout_r,self.EvapsA[0].Tin_a
            ax.plot(np.array([PropsSI('H','P',self.psat_r,'T',self.EvapsA[0].Fins.Air.Tdb,self.Ref)]),self.psat_r*np.ones(1),'x',label='theoretical maximum at T_in_air A[0]')
            ax.legend()
            #plt.show()


def refside_maldistribution_study(evap_type='LRCS'):
    #run different refrigerant side maldistributions to check effect on performance
    refside_maldistributions=[0,0.05,0.1,0.2,0.3,0.4,0.5]
    MD_severity=refside_maldistributions
    if evap_type=='RAC':
        num_evaps=7
    elif evap_type=='LRCS':
        num_evaps=8
    else:
        num_evaps=2 #number of evaporators
        
    for i in range(len(refside_maldistributions)):
        refside_maldistributions[i]=np.linspace(1.+refside_maldistributions[i],1.-refside_maldistributions[i],num_evaps)
    print refside_maldistributions

    Target_SH=5.0
    filenameMDref =evap_type+'-NCircuit_refMD_linear.csv'
    
    evap=MCE_N()
    evap.num_evaps=num_evaps #update evaporator
    evap.Target_SH=Target_SH
    evap.interleaved=False
    evap.maldistributed=False
    evap.Calculate(evap_type)
    evap.TestDescription='Standard' #to use for plotting in Excel Details
    evap.md_severity=str(0) #to use for plotting in Excel 
    evap.Details='Standard without MD' 
    Write2CSV(evap,open(filenameMDref,'w'),append=False)
    Q_base=evap.Q

    for i in range(len(refside_maldistributions)):
        evap=MCE_N()
        evap.Target_SH=Target_SH
        evap.interleaved=False
        evap.num_evaps=num_evaps #update evaporator
        evap.ref_maldistributed=refside_maldistributions[i]
        evap.Calculate(evap_type)
        evap.TestDescription='Standard' #to use for plotting in Excel Details
        evap.md_severity=str(MD_severity[i]) #to use for plotting in Excel 
        evap.Details=make_name('Standard ',str(np.round(refside_maldistributions[i],2)),'ref flow MD') 
        Write2CSV(evap,open(filenameMDref,'a'),append=True)
        Q_noninterleaved=evap.Q

    for i in range(len(refside_maldistributions)):
        evap=MCE_N()
        evap.Target_SH=Target_SH
        evap.interleaved=True
        evap.num_evaps=num_evaps #update evaporator
        evap.ref_maldistributed=refside_maldistributions[i]
        evap.Calculate(evap_type)
        evap.TestDescription='Interleaved' #to use for plotting in Excel Details
        evap.md_severity=str(MD_severity[i]) #to use for plotting in Excel 
        evap.Details=make_name('Interleaved ',str(np.round(refside_maldistributions[i],2)),'ref flow MD') 
        Write2CSV(evap,open(filenameMDref,'a'),append=True)
        Q_interleaved=evap.Q
        
    print "capacity-non-interleaved",Q_noninterleaved,"capacity, interleaved",Q_interleaved,"ratio",(Q_interleaved/Q_noninterleaved),"performance improvement over non-interleaved",((Q_interleaved-Q_noninterleaved)/Q_noninterleaved)*100,'%'
    print "Capacity of basecase without maldistribution",Q_base,"performance degradation caused by maldistribution",((Q_base-Q_noninterleaved)/Q_base)*100
    #plt.show()

def flow_maldistribution_profiles(num_evaps,type,severity=[0,0.05,0.1,0.2,0.3,0.4,0.5],parametric_study=False,custom=False,profile=np.array(range(5))):
    """
    generates a scaled maldistribution profile, returns an array with the distribution factors
    to be used on air or refrigerant side for flowdistribution
    
    num_evaps -> number of circuits
    type -> maldistribution type (linear, pyramid, halflinear A/B, etc.
    severety -> list with severities (0 = no maldistribution, 1=original profile) or
    severety -> float that gives severety
    
    Note: distribution is given as multiplicator to nominal circuit flowrate
    
    Distribution profiles see also:
    D:\Purdue\Thesis\Tex-document\source files and links\interleaved circuitry\maldistribution profiles.xlsx
    """

    if custom==True:
        #use a custom generated profile, note: max_dev is used as severity!
        profile=profile #just to show what is going on
    elif type=='linear':
        #use linear profile
        profile=np.linspace(2,0,num_evaps)  #maximum possible maldistribution (0 flowrate at one circuit)
    elif type=='pyramid':
        #use pyramidial profile
        num_steps=np.ceil(float(num_evaps)/2.0) #make such that have an overlap for odd numbers of circuits
        profile_up=np.linspace(0,2,num_steps)
        profile_down=np.flipud(profile_up)
        if np.mod(num_evaps,2):
            profile_down=profile_down[1:] #remove overlap for odd number of evaps
        profile=np.hstack((profile_up,profile_down)) #add the increasing and decreasing arrays 
    elif type=='Halflinear A':
        ###
        num_steps=np.floor_divide(num_evaps,2) #make such that odd number circuit evaps have one less maldistributed circuit
        profile_up=np.linspace(0,1,num_steps)
        profile_cnst=np.ones(num_evaps-num_steps)*1.0
        profile=np.hstack((profile_up,profile_cnst))  
    elif type=='Halflinear B':
        ###
        num_steps=np.floor_divide(num_evaps,2) #make such that odd number circuit evaps have one less maldistributed circuit
        profile_cnst=np.ones(num_evaps-num_steps)*1.0
        profile_up=np.linspace(1,2,num_steps)
        profile=np.hstack((profile_cnst,profile_up))
    profile=(profile/np.sum(profile))*num_evaps #normalize
    #print num_evaps,type,len(profile),severity,np.round(profile,2)
    maldistribution=maldistribution_scaler(profile,severity=severity,parametric_study=parametric_study)
    
    try:
        len(maldistribution[0])
        dim_md=np.zeros(len(maldistribution))
        for i in range(len(maldistribution)):
            dim_md[i]=dim_md_calc_array(maldistribution[i])
    except:
        dim_md=[dim_md_calc_array(maldistribution)]    
    return maldistribution,dim_md

def flow_maldistribution_profiles_tester():
    #tester)
    from matplotlib import rc
    font = {'family' : 'sans-serif'}
        #'weight' : 'bold',
        #'size'   : 'larger'}

    rc('font', **font)  # pass in the font dict as kwargs
    
    for severity in [0.8,[0,0.05,0.1,0.2,0.3,0.4,0.5]]:
        for num_evaps in [7,8]:
            for type in ['linear','pyramid','Halflinear A','Halflinear B']:
                fig=plt.figure()
                ax = fig.add_subplot(111)
                plt.rc('text', usetex=True)
                use_fontsize=11
                plt.rcParams['font.size']=use_fontsize 
                plt.rc('legend',**{'fontsize':use_fontsize})
                try:
                    len(severity)
                    parametric_study=True
                except:
                    parametric_study=False
                profiles,dim_md=flow_maldistribution_profiles(num_evaps,type,severity=severity,parametric_study=parametric_study)
                if parametric_study==True:
                    #profiles=np.rot90(profiles)
                    for i in range(len(profiles)):
                        print i,profiles[i],dim_md[i],(range(num_evaps)),type,num_evaps
                        #plt.plot(range(num_evaps),profiles[i],'o',ls='-',label=type+r', $\varepsilon$ = '+str(np.round(dim_md[i],2)))#+" sev="+str(severity))  #x=circuit
                        plt.plot(profiles[i],range(num_evaps),'o',ls='-',label=type+r', $\varepsilon$ = '+str(np.round(dim_md[i],2)))#+" sev="+str(severity))     #y=circuit
                else:
                    if type=='pyramid':type='Pyramid' #capitalize..
                    if type=='linear':type='Linear' #capitalize..
                    #plt.plot(range(num_evaps),profiles,'o',ls='-',label=type+r', $\varepsilon$ = '+str(np.round(dim_md[0],2)))#+" sev= "+str(severity)) #x=circuit
                    plt.plot(profiles,range(num_evaps),'o',ls='-',label=type+r', $\varepsilon$ = '+str(np.round(dim_md[0],2)))#+" sev= "+str(severity))  #y=circuit
                plt.legend(loc='best',fancybox=True)#,title='LEGEND')
                plt.xlim(0,2)
                plt.ylim(num_evaps-0.8,-0.2)
                for label in ax.xaxis.get_ticklabels():label.set_fontsize(use_fontsize)# label is a Text instance #label.set_color('red')#label.set_rotation(45)
                for label in ax.yaxis.get_ticklabels(): label.set_fontsize(use_fontsize)
                ax.set_xlabel('Normalized Flowrate',fontsize=use_fontsize)
                ax.set_ylabel('Circuit Number',fontsize=use_fontsize)
                #modify labels
                labels = [str(int(item+1)) for item in ax.yaxis.get_majorticklocs()]
                ax.set_yticklabels(labels)
                labels = [str(item) for item in ax.xaxis.get_majorticklocs()]
                #labels[0]='0'
                ax.set_xticklabels(labels)
                plt.savefig(str(num_evaps)+type+str(np.round(dim_md[0],2))+'.pdf',bbox_inches='tight')
        plt.close('all')
        #plt.show()

def air_temp_maldistribution_profiles(num_evaps,type,severity=[0,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.8,1.0,1.1],parametric_study=False,custom=False):
    #test airside temperature maldistribution
    if type=='RAC':
        Original_Profile=np.array([8.963222053,8.493439286,4.582012002,-0.840190032,-5.842297047,-8.493439286,-6.862746977])  
        maldistribution=maldistribution_scaler(Original_Profile,severity=severity,parametric_study=parametric_study)
    try:
        len(maldistribution[0])
        dim_md=np.zeros(len(maldistribution))
        for i in range(len(maldistribution)):
            dim_md[i]=dim_md_calc_array(maldistribution[i])
    except:
        dim_md=[dim_md_calc_array(maldistribution)]
    return maldistribution,dim_md

def air_temp_maldistribution_profiles_tester():
    for severity in [1.0,[0,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.8,1.0,1.1]]:
        for num_evaps in [7]:
            plt.figure()
            plt.rc('text', usetex=True)
            plt.rc('legend',**{'fontsize':16})
            plt.rcParams['font.size']=20 
            plt.rcParams['axes.labelsize']=18
            plt.rcParams['axes.titlesize']=18
            for type in ['RAC']:
                try:
                    len(severity)
                    parametric_study=True
                except:
                    parametric_study=False
                print severity,parametric_study
                profiles,dim_md=air_temp_maldistribution_profiles(num_evaps,type,severity=severity,parametric_study=parametric_study)
                if parametric_study==True:
                    #profiles=np.rot90(profiles)
                    for i in range(len(profiles)):
                        print i,profiles[i],dim_md[i],(range(num_evaps)),type,num_evaps
                        plt.plot(range(num_evaps),profiles[i],label=type+r'severity'+str(np.round(severity[i],2)))
                else:
                    plt.plot(range(num_evaps),profiles,label=type+r'$severity$'+str(np.round(severity,2)))
                plt.legend(loc='best',fancybox=True,title='LEGEND')
                plt.show()
        
def dim_md_calc_array(m_dot):
    #equation (4.10 in prelim)
    #is recalculated in output list for backwards compatebility
    #input is a single array-like m_dot
    m_dot_bar=np.average(m_dot)
    m_dot_i_minus_m_dot_ave=m_dot-m_dot_bar
    sum_squares=np.sum(np.multiply(m_dot_i_minus_m_dot_ave,m_dot_i_minus_m_dot_ave))
    #print 'm_dot',m_dot
    epsilon=np.sqrt(sum_squares/len(m_dot))/m_dot_bar
    return epsilon 

def maldistribution_scaler(maldistribution_profile,severity=1,parametric_study=False):
    """
    scales the np.array(maldistribution_profile):
    severity = 0 --> equal distribution
    severty= 1--> given profile
    severity>0 --> extrapolate (need to take care if gets smaller than 0....)
    """

    def scaler(maldistribution_profile,severity):
        equal_distribution=np.ones(len(maldistribution_profile))*np.sum(maldistribution_profile)/np.size(maldistribution_profile)
        maldistrib_tmp=equal_distribution+(maldistribution_profile-equal_distribution)*severity
        return (maldistrib_tmp/np.sum(maldistrib_tmp))*np.sum(maldistribution_profile) #normalize to original sum
    
    if parametric_study:
        maldistributions=[]
        for severity in severity:  #severity is iterable
            maldistributions.append(scaler(maldistribution_profile,severity))
        return maldistributions
    else:
        return scaler(maldistribution_profile,severity)  #severity=scalar

def make_name(startstring,maldistribution,endstring):
    #helperfunction to set testdescription while preventing linebreaks
    try: 
        if len(maldistribution)>5:
            TestDescription=startstring+" "+str(maldistribution[0:5])[:-1]+" "+str(maldistribution[5:])[1:]+" "+endstring
        else:
            TestDescription=startstring+str(maldistribution)+endstring
        return TestDescription
    except:  #maldistribution is something else but not an array!
        return startstring+" "+endstring  #assume no maldistribution

def airside_maldistribution_study(evap_type='LRCS',MD_Type=None,MD_severity=None,airside_maldistributions=None,num_evaps=2,filenameMDair='debug.csv',Hybrid='equal_flow',adjust_area_fraction_iternum=10):
    #run different airside flow maldistributions to check effect on performance
    if MD_Type==None:
        airside_maldistributions=[0,0.05,0.1,0.2,0.3,0.4,0.5]
        MD_severity=airside_maldistributions
        for i in range(len(airside_maldistributions)):
            airside_maldistributions[i]=np.linspace(1.+airside_maldistributions[i],1.-airside_maldistributions[i],num_evaps)
        filenameMDair =evap_type+'-NCircuit_airMD_linear.csv'
    elif MD_Type=="LRCS_Type_A":  #see D:\Purdue\Thesis\Tex-document\source files and links\interleaved circuitry\LRCS\maldistribution profiles.xlsx
        Original_Profile=np.array([0.0135415976822403,0.0221506896994024,0.0369272399580833,0.111895731459975,0.106096750782192,0.265750418904745,0.196007841404425,0.247629730108938])*8.0  #different definition compared to normal ACHP MCE
        MD_severity=[0,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.8,0.9,1.0]
        #MD_severity=[1.0]
        airside_maldistributions=maldistribution_scaler(Original_Profile,severity=MD_severity,parametric_study=True)
        num_evaps=8 #number of evaporators
        filenameMDair =evap_type+'-NCircuit_airMD_typeA.csv'
    elif MD_Type=="RAC_avg":  #see D:\Purdue\Thesis\Tex-document\source files and links\interleaved circuitry\LRCS\maldistribution profiles.xlsx
        Original_Profile=np.array([0.12741462,0.154399656,0.132419274,0.128817928,0.178563757,0.15644018,0.121944584])*7.0  #different definition compared to normal ACHP MCE, average between hybrid and PEXV
        MD_severity=[0,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.8,1.0,1.1]
        airside_maldistributions=maldistribution_scaler(Original_Profile,severity=MD_severity,parametric_study=True)
        num_evaps=7 #number of evaporators
        filenameMDair =evap_type+'-NCircuit_airMD_avg.csv' 
    elif MD_Type=='custom':
        print "using custum maldistribution as passed in"
    Target_SH=5.0

    evap=MCE_N()
    evap.Target_SH=Target_SH
    evap.interleaved=False
    evap.maldistributed=False
    evap.num_evaps=num_evaps #update evaporator
    evap.Calculate(evap_type)
    evap.TestDescription='Standard' #to use for plotting in Excel Details
    evap.md_severity=str(0) #to use for plotting in Excel 
    evap.Details="without maldistribution"
    Write2CSV(evap,open(filenameMDair,'w'),append=False)
    Q_base=evap.Q

    for i in range(len(airside_maldistributions)):
        evap=MCE_N()
        evap.Target_SH=Target_SH
        evap.interleaved=False
        evap.num_evaps=num_evaps #update evaporator
        evap.maldistributed=airside_maldistributions[i]
        evap.Calculate(evap_type)
        evap.TestDescription='Standard' #to use for plotting in Excel Details
        evap.md_severity=str(MD_severity[i]) #to use for plotting in Excel 
        evap.Details=make_name('Standard ',str(np.round(airside_maldistributions[i],2)),'air flow MD') 
        Write2CSV(evap,open(filenameMDair,'a'),append=True)
        Q_noninterleaved=evap.Q
    
    for i in range(len(airside_maldistributions)):
        evap=MCE_N()
        evap.Target_SH=Target_SH
        evap.interleaved=True
        evap.num_evaps=num_evaps #update evaporator
        evap.maldistributed=airside_maldistributions[i]
        evap.Calculate(evap_type)
        evap.TestDescription='Interleaved' #to use for plotting in Excel Details
        evap.md_severity=str(MD_severity[i]) #to use for plotting in Excel 
        evap.Details=make_name('Interleaved ',str(np.round(airside_maldistributions[i],2)),'air flow MD') 
        Write2CSV(evap,open(filenameMDair,'a'),append=True)
        Q_interleaved=evap.Q
    
    for i in range(len(airside_maldistributions)):
        evap=MCE_N()
        evap.Target_SH=Target_SH
        evap.Hybrid=Hybrid
        if evap.Hybrid=='adjust_superheat_iter':
            evap.Hybrid_ref_distribution=airside_maldistributions[i]
            evap.adjust_area_fraction_iternum=adjust_area_fraction_iternum
        evap.interleaved=False
        evap.num_evaps=num_evaps #update evaporator
        evap.maldistributed=airside_maldistributions[i]
        evap.Calculate(evap_type)
        evap.TestDescription='Equal flow' #to use for plotting in Excel Details
        evap.md_severity=str(MD_severity[i]) #to use for plotting in Excel 
        evap.Details=make_name('Equal flow ',str(np.round(airside_maldistributions[i],2)),'air flow MD') 
        Write2CSV(evap,open(filenameMDair,'a'),append=True)
        Q_hybrid=evap.Q
        
    print "capacity-non-interleaved",Q_noninterleaved,"capacity, interleaved",Q_interleaved,"ratio",(Q_interleaved/Q_noninterleaved),"performance improvement over non-interleaved",((Q_interleaved-Q_noninterleaved)/Q_noninterleaved)*100,'%'
    print "capacity-non-interleaved",Q_noninterleaved,"capacity, hybrid",Q_hybrid,"ratio",(Q_hybrid/Q_noninterleaved),"performance improvement over non-interleaved",((Q_hybrid-Q_noninterleaved)/Q_noninterleaved)*100,'%'
    print "Capacity of basecase without maldistribution",Q_base,"performance degradation caused by maldistribution",((Q_base-Q_noninterleaved)/Q_base)*100
    #plt.show()
    
def airside_temp_maldistribution_study(evap_type='RAC',MD_Type=None):
    #run different airside maldistributions to check effect on performance
    if MD_Type==None:
        airside_temp_maldistributions=[0.0,0.2,0.3,0.4,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0]
        MD_severity=airside_temp_maldistributions #needed for plotting with pandas later on
        num_evaps=2 #number of evaporators
        for i in range(len(airside_temp_maldistributions)):
            airside_temp_maldistributions[i]=np.linspace(airside_temp_maldistributions[i],-airside_temp_maldistributions[i],num_evaps)
        filenameMDair =evap_type+'-NCircuit_airTMD_default.csv'
    elif MD_Type=="RAC_Temp":  #see D:\Purdue\Thesis\Tex-document\source files and links\interleaved circuitry\LRCS\maldistribution profiles.xlsx
        Original_Profile=np.array([8.963222053,8.493439286,4.582012002,-0.840190032,-5.842297047,-8.493439286,-6.862746977])  
        MD_severity=[0,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.8,1.0,1.1]
        airside_temp_maldistributions=maldistribution_scaler(Original_Profile,severity=MD_severity,parametric_study=True)
        num_evaps=7 #number of evaporators
        filenameMDair =evap_type+'-NCircuit_airTMD_typeRAC_Temp.csv'
    else:
        print 'something wrong in evap_N_circ_cross_counter_flow'

    Target_SH=5.0
    
    evap=MCE_N()
    evap.num_evaps=num_evaps #update evaporator
    evap.Target_SH=Target_SH
    evap.interleaved=False
    evap.maldistributed=False
    evap.Calculate(evap_type)
    evap.TestDescription='Standard'
    evap.md_severity=0
    evap.Details='without maldistribution'
    Write2CSV(evap,open(filenameMDair,'w'),append=False)
    Q_base=evap.Q

    for i in range(len(airside_temp_maldistributions)):
        evap=MCE_N()
        evap.num_evaps=num_evaps #update evaporator
        evap.Target_SH=Target_SH
        evap.interleaved=False
        evap.air_temp_maldistributed=airside_temp_maldistributions[i]
        evap.Calculate(evap_type)
        evap.TestDescription='Standard' #to use for plotting in Excel Details
        evap.md_severity=str(MD_severity[i]) #to use for plotting in Excel 
        evap.Details=make_name('Standard ',str(np.round(airside_temp_maldistributions[i],2)),'air flow temp MD') #to use for plotting in Excel 
        Write2CSV(evap,open(filenameMDair,'a'),append=True)
        Q_noninterleaved=evap.Q

    for i in range(len(airside_temp_maldistributions)):
        evap=MCE_N()
        evap.num_evaps=num_evaps #update evaporator
        evap.Target_SH=Target_SH
        evap.interleaved=True
        evap.air_temp_maldistributed=airside_temp_maldistributions[i]
        evap.Calculate(evap_type)
        evap.TestDescription='Interleaved' #to use for plotting in Excel Details
        evap.md_severity=str(MD_severity[i]) #to use for plotting in Excel 
        evap.Details=make_name('Interleaved ',str(np.round(airside_temp_maldistributions[i],2)),'air flow temp MD') #to use for plotting in Excel 
        Write2CSV(evap,open(filenameMDair,'a'),append=True)
        Q_interleaved=evap.Q
        

    for i in range(len(airside_temp_maldistributions)):
        evap=MCE_N()
        evap.num_evaps=num_evaps #update evaporator
        evap.Target_SH=Target_SH
        evap.Hybrid='equal_flow_DT_in'
        evap.interleaved=False
        evap.air_temp_maldistributed=airside_temp_maldistributions[i]
        evap.Calculate(evap_type)
        evap.TestDescription='Equal flow' #to use for plotting in Excel Details
        evap.md_severity=str(MD_severity[i]) #to use for plotting in Excel 
        evap.Details=make_name('equal_flow_DT_in ',str(np.round(airside_temp_maldistributions[i],2)),'air flow temp MD') #to use for plotting in Excel 
        Write2CSV(evap,open(filenameMDair,'a'),append=True)
        Q_hybrid=evap.Q
        
    print "capacity-non-interleaved",Q_noninterleaved,"capacity, interleaved",Q_interleaved,"ratio",(Q_interleaved/Q_noninterleaved),"performance improvement over non-interleaved",((Q_interleaved-Q_noninterleaved)/Q_noninterleaved)*100,'%'
    print "capacity-non-interleaved",Q_noninterleaved,"capacity, hybrid",Q_hybrid,"ratio",(Q_hybrid/Q_noninterleaved),"performance improvement over non-interleaved",((Q_hybrid-Q_noninterleaved)/Q_noninterleaved)*100,'%'
    print "Capacity of basecase without maldistribution",Q_base,"performance degradation caused by maldistribution",((Q_base-Q_noninterleaved)/Q_base)*100
    #plt.show()

def combined_maldistribution_study(evap_type='RAC',MD_Type='RAC_combined'):
    #study combination of airside and refrigerant side maldistribution
    if MD_Type!='RAC_combined':
        print "something wrong in evap_N_circ_cross_counter_flow"
    elif MD_Type=='RAC_combined':
        MD_severity=[0,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.8,1.0,1.1]
        Original_Profile=np.array([8.963222053,8.493439286,4.582012002,-0.840190032,-5.842297047,-8.493439286,-6.862746977])  
        airside_temp_maldistributions=maldistribution_scaler(Original_Profile,severity=MD_severity,parametric_study=True) #generate array for airside temperature maldst
        Original_Profile=np.array([0.12741462,0.154399656,0.132419274,0.128817928,0.178563757,0.15644018,0.121944584])*7.0  #different definition compared to normal ACHP MCE, average between hybrid and PEXV
        airside_flow_maldistributions=maldistribution_scaler(Original_Profile,severity=MD_severity,parametric_study=True)
        num_evaps=7 #number of evaporators
        filenameMDair =evap_type+'-NCircuit_combined_RTU_MD.csv'
        Target_SH=5.0
        
    for i in range(len(airside_temp_maldistributions)):
        #standard evaporator
        evap=MCE_N()
        evap.num_evaps=num_evaps #update evaporator
        evap.Target_SH=Target_SH
        #evap.Hybrid='equal_flow_DT_in' do not set attribute
        evap.interleaved=False
        evap.air_temp_maldistributed=airside_temp_maldistributions[i] #apply airside maldistribution
        evap.maldistributed=airside_flow_maldistributions[i] #apply refrigerant side maldistribution
        evap.Calculate(evap_type)
        evap.TestDescription='Standard' #to use for plotting in Excel Details
        evap.md_severity=str(MD_severity[i]) #to use for plotting in Excel 
        evap.Details=make_name('TempMD ',str(np.round(airside_temp_maldistributions[i],2),' '))+make_name('FlowMD ',np.round(airside_flow_maldistributions[i],2),'') #to use for plotting in Excel 
        if i==0:
            Write2CSV(evap,open(filenameMDair,'a'),append=False)
        else:
            Write2CSV(evap,open(filenameMDair,'a'),append=True)
        Q_hybrid=evap.Q
        
    for i in range(len(airside_temp_maldistributions)):
        #interleaved evaporator
        evap=MCE_N()
        evap.num_evaps=num_evaps #update evaporator
        evap.Target_SH=Target_SH
        #evap.Hybrid='equal_flow_DT_in' do not set attribute
        evap.interleaved=True
        evap.air_temp_maldistributed=airside_temp_maldistributions[i] #apply airside maldistribution
        evap.maldistributed=airside_flow_maldistributions[i] #apply refrigerant side maldistribution
        evap.Calculate(evap_type)
        evap.TestDescription='Interleaved' #to use for plotting in Excel Details
        evap.md_severity=str(MD_severity[i]) #to use for plotting in Excel 
        evap.Details=make_name('TempMD ',str(np.round(airside_temp_maldistributions[i],2)),' ')+make_name('FlowMD ',np.round(airside_flow_maldistributions[i],2),'') #to use for plotting in Excel 
        if i==0:
            Write2CSV(evap,open(filenameMDair,'a'),append=True)
        else:
            Write2CSV(evap,open(filenameMDair,'a'),append=True)
        Q_hybrid=evap.Q
    
    for i in range(len(airside_temp_maldistributions)):
        # pseudo hybrid scheme
        evap=MCE_N()
        evap.num_evaps=num_evaps #update evaporator
        evap.Target_SH=Target_SH
        evap.Hybrid='equal_flow_DT_in'
        evap.interleaved=False
        evap.air_temp_maldistributed=airside_temp_maldistributions[i] #apply airside maldistribution
        evap.maldistributed=airside_flow_maldistributions[i] #apply refrigerant side maldistribution
        evap.Calculate(evap_type)
        evap.TestDescription='Equal flow' #to use for plotting in Excel Details
        evap.md_severity=str(MD_severity[i]) #to use for plotting in Excel 
        evap.Details=make_name('TempMD ',str(np.round(airside_temp_maldistributions[i],2)),' ')+make_name('FlowMD ',np.round(airside_flow_maldistributions[i],2),'') #to use for plotting in Excel 
        if i==0:
            Write2CSV(evap,open(filenameMDair,'a'),append=True)
        else:
            Write2CSV(evap,open(filenameMDair,'a'),append=True)
        Q_hybrid=evap.Q

def sh_equalizer_tester(evap_type='LRCS',num_evaps=8,md_type='linear',Target_SH=5.0,MD_severity=[0.2]):
    #test if superheat equalizer works
    if len(MD_severity)>1:
        parametric_study=True
    else:
        parametric_study=False
    maldistributions=flow_maldistribution_profiles(num_evaps,md_type,severity=MD_severity,parametric_study=parametric_study,custom=False,profile=np.array(range(5)))
    Append=False
    for i_iternum in [10]:
        for i_MD in range(len(MD_severity)):
            evap=MCE_N()
            evap.Target_SH=Target_SH
            evap.interleaved=False
            evap.num_evaps=num_evaps #update evaporator
            evap.maldistributed=maldistributions[0][i_MD]
            #evap.Hybrid='adjust_area_fraction' #works fine
            #evap.Hybrid='adjust_superheat' #doesn't work correctly
            evap.Hybrid='adjust_superheat_iter' #10 iterations is good for 5% integration
            filename_sh_eq='sh_equalizer_tester_'+evap.Hybrid+'_'+md_type+'.csv'
            evap.adjust_area_fraction_iternum=i_iternum
            if not parametric_study:
                evap.Hybrid_ref_distribution=maldistributions[0][i_MD]
            else:
                evap.Hybrid_ref_distribution=maldistributions[0][i_MD]
                evap.md_Epsilon=maldistributions[1][i_MD] #to use for plotting in Excel 
            evap.md_severity=str(MD_severity[i_MD]) #to use for plotting in Excel 
            #evap.mdot_r=0.114424706615
            evap.Calculate(evap_type)
            evap.TestDescription='Standard' #to use for plotting in Excel Details
            evap.Details=make_name('Standard '+str(i_iternum),str(np.round(maldistributions[0][i_MD],2)),'air flow MD') 
            Write2CSV(evap,open(filename_sh_eq,'a'),append=Append)
            if Append==False: Append=True #do append after first run
            Q_noninterleaved=evap.Q


if __name__=='__main__':
    if 0:
        maldistribution_profile=np.array([0.0135415976822403,0.0221506896994024,0.0369272399580833,0.111895731459975,0.106096750782192,0.265750418904745,0.196007841404425,0.247629730108938])
        tmp=maldistribution_scaler(np.array(maldistribution_profile),severity=[0.5,1.0,1.1],parametric_study=True)
        print tmp
        print make_name('bla',maldistribution_profile,'endstring')#.TestDescription
        print len(maldistribution_profile), np.sum(tmp[0]),np.sum(tmp[1]),np.sum(tmp[2]),"and original size was",np.sum(maldistribution_profile)
    if 0: #test profiles
        flow_maldistribution_profiles_tester()
        air_temp_maldistribution_profiles_tester()
    if 1: #run parametric study for 2-circuit cases only
        airside_maldistribution_study(evap_type='LRCS',MD_Type=None,Hybrid='adjust_superheat_iter',adjust_area_fraction_iternum=30)  #this runs the 2-circuit case with the only possible maldistribution for that case (code is ugly...)
    if 0: #run parametric studies
        airside_maldistribution_study(evap_type='LRCS',MD_Type="LRCS_Type_A")
        refside_maldistribution_study(evap_type='LRCS')
        airside_temp_maldistribution_study(evap_type='RAC',MD_Type="RAC_Temp")
        refside_maldistribution_study(evap_type='RAC')
        airside_maldistribution_study(evap_type='RAC',MD_Type="RAC_avg")
        combined_maldistribution_study(evap_type='RAC',MD_Type='RAC_combined')
    if 0: #test superheat equalizer
        sh_equalizer_tester()
    if 0: #run different flow distribution profiles for LRCS
        MD_severity=[0,0.05,0.1,0.2,0.3,0.4,0.5]
        #MD_severity=[0.5]
        for md_type in ['linear','pyramid','Halflinear A','Halflinear B']:
            maldistributions=flow_maldistribution_profiles(8,md_type,severity=MD_severity,parametric_study=True,custom=False,profile=np.array(range(5)))
            #print maldistributions[1][1],md_type
            if 0:
                sh_equalizer_tester(evap_type='LRCS',num_evaps=8,md_type=md_type,Target_SH=5.0,MD_severity=MD_severity)
            if 1:
                airside_maldistribution_study(evap_type='LRCS',MD_Type=md_type,MD_severity=MD_severity,airside_maldistributions=maldistributions[0],num_evaps=8,filenameMDair=md_type+'LRCS'+'.csv')
        
    #plt.show() #show plots, if any
