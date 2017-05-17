from __future__ import division, print_function, absolute_import
from math import pi,log,exp

from scipy.optimize import brentq #solver to find roots (zero points) of functions
import numpy as np
import pandas as pd

import CoolProp as CP
from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI

from ACHP.convert_units import *

#from .Correlations import f_h_1phase_Tube,ShahEvaporation_Average, LockhartMartinelli,LMPressureGradientAvg,AccelPressureDrop,TwoPhaseDensity
#from .FinCorrelations import WavyLouveredFins,FinInputs,IsFinsClass, HerringboneFins, PlainFins
#from .DryWetSegment import DWSVals, DryWetSegment
#from .ACHPTools import ValidateFields


# class EvaporatorClass():
#     
#     def __init__(self,**kwargs):
#         self.__dict__.update(kwargs)
#         
#     def Update(self,**kwargs):
#         self.__dict__.update(kwargs)
#         
#     def OutputList(self):
#         """
#             Return a list of parameters for this component for further output
#             
#             It is a list of tuples, and each tuple is formed of items:
#                 [0] Description of value
#                 [1] Units of value
#                 [2] The value itself
#         """
#         Output_List=[]
#         #append optional parameters, if applicable
#         if hasattr(self,'TestName'):
#             Output_List.append(('Name','N/A',self.TestName)) 
#         if hasattr(self,'TestDescription'):
#             Output_List.append(('Description','N/A',self.TestDescription))
#         if hasattr(self,'TestDetails'):
#             Output_List.append(('Details','N/A',self.TestDetails))
#         Output_List_default=[                                                   #default output list
#             ('Volumetric flow rate','m^3/s',self.Fins.Air.Vdot_ha),
#             ('Inlet Dry bulb temp','K',self.Tin_a),
#             ('Inlet Air pressure','Pa',self.Fins.Air.p),
#             ('Inlet Air Relative Humidity','-',self.Fins.Air.RH),
#             ('Tubes per bank','-',self.Fins.Tubes.NTubes_per_bank),
#             ('Number of banks','-',self.Fins.Tubes.Nbank),
#             ('Number circuits','-',self.Fins.Tubes.Ncircuits),
#             ('Length of tube','m',self.Fins.Tubes.Ltube),
#             ('Tube OD','m',self.Fins.Tubes.OD),
#             ('Tube ID','m',self.Fins.Tubes.ID),
#             ('Tube Long. Pitch','m',self.Fins.Tubes.Pl),
#             ('Tube Transverse Pitch','m',self.Fins.Tubes.Pt),
#             ('Outlet superheat','K',self.Tout_r-self.Tdew_r),
#             ('Fins per inch','1/in',self.Fins.Fins.FPI),
#             ('Fin waviness pd','m',self.Fins.Fins.Pd),
#             ('Fin waviness xf','m',self.Fins.Fins.xf),
#             ('Fin thickness','m',self.Fins.Fins.t),
#             ('Fin Conductivity','W/m-K',self.Fins.Fins.k_fin),
#             ('Fins Type','-',self.FinsType),
#             ('Q Total','W',self.Q),
#             ('Q Superheat','W',self.Q_superheat),
#             ('Q Two-Phase','W',self.Q_2phase),
#             ('Inlet ref. temp','K',self.Tin_r),
#             ('Outlet ref. temp','K',self.Tout_r),
#             ('Outlet air temp','K',self.Tout_a),
#             ('Inlet Entropy','J/kg-K',self.sin_r),
#             ('Outlet Entropy','J/kg-K',self.sout_r),
#             ('Inlet Pressure (sat)','Pa',self.psat_r),
#             ('Outlet pressure (sat+DP)','Pa',self.psat_r+self.DP_r),
#             ('Evaporator inlet quality','-',self.xin_r),
#             ('Evaporator ref. flowrate','kg/s',self.mdot_r),
#             ('Pressure Drop Total','Pa',self.DP_r),
#             ('Pressure Drop Superheat','Pa',self.DP_r_superheat),
#             ('Pressure Drop Two-Phase','Pa',self.DP_r_2phase),
#             ('Charge Total','kg',self.Charge),
#             ('Charge Superheat','kg',self.Charge_superheat),
#             ('Charge Two-Phase','kg',self.Charge_2phase),
#             ('Mean HTC Superheat','W/m^2-K',self.h_r_superheat),
#             ('Mean HTC Two-phase','W/m^2-K',self.h_r_2phase),
#             ('Wetted Area Fraction Superheat','-',self.w_superheat),
#             ('Wetted Area Fraction Two-phase','-',self.w_2phase),
#             ('Mean Air HTC','W/m^2-K',self.Fins.h_a),
#             ('Surface Effectiveness','-',self.Fins.eta_a),
#             ('Air-side area (fin+tubes)','m^2',self.Fins.A_a),
#             ('Mass Flow rate dry Air','kg/s',self.Fins.mdot_da),
#             ('Mass Flow rate humid Air','kg/s',self.Fins.mdot_ha),
#             ('Pressure Drop Air-side','Pa',self.Fins.dP_a),
#             ('Superheat','K',self.DT_sh_calc),
#             ('Sensible Heat Ratio','-',self.SHR),
#         ]
#         for i in range(0,len(Output_List_default)):         #append default parameters to output list
#             Output_List.append(Output_List_default[i])
#         return Output_List
#         
#     def AirSideCalcs(self):
#         
#         pass
#     
#     def Initialize(self):
#         
#         pass
#         
#     def Calculate(self):
#         
#         #Initialize
#         self.Initialize()
#         
#         pass



class StructEvap():
    
    def __init__(self,filename):
        self.REV = 0 #reverse calculation
   
    #def StructEvap(self,filename):      #reading the data file in
        #initialize
        EvapNode={
            'NodNo':int(0),#node number
            'InNum':int(0),#number of tubes flowing into the node
            'OutNum':int(0),#number of the tubes flowing out of the node
            'BranIN':[],#index of the tube branch flowing in
            'BranOUT':[],#index of the tube branch flowing out
            }
        EvapBranch = {
            'BranNo':int(0),     #branch number
            'EqulNo':int(0),     #the equivalent branch number
            'Ini':int(0),        #signal variable
            'GrFac':0.0,    #mass flow distribution factor of the branch
            'Gr':0.0,       #mass flux of the branch
            'TubNum':int(0),     #total tube numbers in the branch
            'TubNo':[],      #index of the tubes in the branch
            'HPi':{'H':0.0,'P':0.0},#inlet enthalpy and pressure 
            'HPo':{'H':0.0,'P':0.0},#outlet enthalpy and pressure
            'm':{'m':0.0,'V':0.0},  #mass and volume
            'Para_Struc':[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]#output parameters of this condenser branch
            }
        TubeEvap = {
            'TubNo':int(0),             #tube No
            'RowNo':int(0),             #the row no where the tube located
            'Refdownstream':int(0),     #downstream tube No at the refrigerant side (flow direction)
            'AirUpstreamUpper':int(0),  #Upstream tube No at the air upper side
            'AirUpstreamLower':int(0),  #Upstream tube No at the air lower side
            'even':int(0),              #tube flow direction
            'GaFac':0.0,                #air flow distribution factor
            'Ga':0.0,                   #maximum air flow rate in the segment
            'HPi':{'H':0.0,'P':0.0},    #refrigerant inlet state, flowing direction
            'HPo':{'H':0.0,'P':0.0},    #refrigerant outlet state, flowing direction
            'm':{'m':0.0,'V':0.0},      #mass and volume of the tube
            'Seg':[]
            }
        
        TubEvpSeg = {'TPi':{'T':0.0,'P':0.0}, 'WHo':{'W':0.0,'H':0.0}}
        
        #endSign=int(0);
        #error=0;
           
        df = pd.read_excel(filename,sheetname=0,header = 0)
        df_node = pd.read_excel(filename,sheetname='node',header = 0)
        df_branch = pd.read_excel(filename,sheetname='branch',header = 0)
        df_tube = pd.read_excel(filename,sheetname='tube',header = 0)
           
        #fgets(scanstr,2048,fp);
        self.NodNum = int(df.EvapStruc[0]);          #Node Number
        self.BranNum = int(df.EvapStruc[1]);         #Branch Number
        self.RowNum = int(df.EvapStruc[2]);          #Row Number
        self.TubeNum = int(df.EvapStruc[3]);          #overall tube Number
        self.SegNum = int(df.EvapStruc[4]);          #Segment Number of per tube
        self.AreaFront = int(df.EvapStruc[5]);       #frontal area
        self.Volum = int(df.EvapStruc[6]);           #total nominal air mass flow rate [m^3/s]
       
        self.Nod = [EvapNode.copy() for k in range(self.NodNum)]
        self.Bra = [EvapBranch.copy() for k in range(self.BranNum)]
        self.GaRow = [0.0 for k in range(self.RowNum)]
        self.Tub = [TubeEvap.copy() for k in range(self.TubeNum)]
       
        #fgets(scanstr,2048,fp);
        for i in range(self.NodNum):   
            self.Nod[i]['NodNo'] = df_node.NodNo[i]
            self.Nod[i]['OutNum'] = df_node.OutNum[i]    #get outlet branches at first
            self.Nod[i]['InNum'] = df_node.InNum[i]      #get inlet branches second
            self.Nod[i]['BranOUT'] = list(df_node.ix[i][3 : 3+self.Nod[i]['OutNum']]) #get outlet branches first
            self.Nod[i]['BranIN'] = list(df_node.ix[i][3+self.Nod[i]['OutNum'] : 3+self.Nod[i]['OutNum']+self.Nod[i]['InNum']]) #get inlet branches second
       
        #fgets(scanstr,2048,fp);
        for i in range(self.BranNum):
            self.Bra[i]['BranNo'] = df_branch.BranNo[i]
            self.Bra[i]['EqulNo'] = df_branch.EqulNo[i]
            self.Bra[i]['GrFac'] = df_branch.GrFac[i]
            self.Bra[i]['TubNum'] = df_branch.TubNum[i]
            self.Bra[i]['TubNo'] = list(df_branch.ix[i][4 : 4+self.Bra[i]['TubNum']])#important, the tube number in branch, inputted first from the compressor suction
       
        #fgets(scanstr,2048,fp);
        #fgets(scanstr,2048,fp);
        for i in range(self.TubeNum):
            self.Tub[i]['Seg'] = [TubEvpSeg.copy() for k in range(self.SegNum)]
            self.Tub[i]['TubNo'] = df_tube.TubNo[i]
            self.Tub[i]['RowNo'] = df_tube.RowNo[i]
            self.Tub[i]['Refdownstream'] = df_tube.Refdownstream[i]
            self.Tub[i]['AirUpstreamUpper'] = df_tube.AirUpstreamUpper[i]
            self.Tub[i]['AirUpstreamLower'] = df_tube.AirUpstreamLower[i]
            self.Tub[i]['GaFac'] = df_tube.GaFac[i]
            self.Tub[i]['even'] = df_tube.even[i]
       
        return 0
       
       
    def DeStructEvap(self):
        ''''delete the memory'''
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
         
        return 0
   
    def _EvapCircuit_Fwd(self,mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,P):

        Gr=0.0;#mass flux
        H_in=0.0;
        
        if (self.Rev): #reversed calculation
            self._EvapCircuit_Rev(mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,P)
        return 0
    
        TPo=TPi;#inlet air state
    
        H_in=HPi['H'];#inlet refrigerant enthalpy
        
        #air side initializing
        rho_air=1/HAPropsSI("V", "T", TPi['T'], "P", 101325, "R", TPi['P']) #[kg dry air/ m^3]
        #if(errorLog.IsError()) {
        #errorLog.Add("StructEvap::_EvapCircuit","rho_air");
        #return -1;}
    
        #const double Ma =Ga*rho_air*4.719e-4; 
        Ga=Ga/self.AreaFront*P['vsp']*P['Ls']/((P['vsp']-P['Do'])*(P['Ls']-P['N']*P['th']));
    
    
        for i in range(self.RowNum):
            self.GaRow[i]=0 #this step can be deleted because GaRow[i] already initialized with zero
            
        for j in range(self.RowNum):
            N=0;
            for i in range(self.TubeNum):
                
                if(self.Tub[i]['RowNo']==j):
                
                    for k in range(self.SegNum):
                        self.Tub[i]['Seg'][k]['TPi']=TPi;
                        self.Tub[i]['Seg'][k]['WHo']['W'] = HAPropsSI("W", "T", TPi['T'], "P", 101325, "R", TPi['P']) #[kg water/kg dry air] ####wair.HumidityRatio(TPi.T,TPi.P)
                        self.Tub[i]['Seg'][k]['WHo']['H'] = HAPropsSI("H", "T", TPi['T'], "P", 101325, "R", TPi['P']) #[J/kg dry air] ####wair.h(TPi.T,TPi.P);
                        #if(errorLog.IsError()) {
                        #    errorLog.Add("StructEvap::_EvapCircuit","WHo");
                        #    return -1;}
        
                    if(self.Tub[i]['RowNo']>0):#not the first row
                        Upper = self.Tub[i]['AirUpstreamUpper']
                        Lower = self.Tub[i]['AirUpstreamLower']
                        if (Upper>=0 and Lower>=0):
                            self.Tub[i]['GaFac']=(self.Tub[Upper]['GaFac']+self.Tub[Lower]['GaFac'])/2
                        elif (Upper>=0):
                            self.Tub[i]['GaFac']=self.Tub[Upper]['GaFac']
                        else:
                            self.Tub[i]['GaFac']=self.Tub[Lower]['GaFac']
                    
                self.GaRow[j]=self.GaRow[j]+self.Tub[i]['GaFac']
                N=N+1;
                #ifend    
            #i loop end
            self.GaRow[j]=self.GaRow[j]/N
        #j loop end
        
        for i in range(self.TubeNum):
            RowN=self.Tub[i]['RowNo'];
            self.Tub[i]['Ga']=self.Tub[i]['GaFac']/self.GaRow[RowN]*Ga
    
        H_out=0.0;  #intermidiate evaporator exit enthalpy
        Res=0;
        DP = 20000.0;
        P_in=0.0;
        P_out=0.0;
        HPi['P']=HPo['P']+DP;   #temperary inlet refrigerant pressure
        HPi['H']=H_in;          #evaporator inlet refrigerant enthalpy
        #*HPi is the intermediate variable
        P_in = HPo['P']+DP;       #temperary inlet refrigerant pressure
        P_out = HPo['P'];         #evaporator outlet refrigerant pressure
        IterN=0;

    
        #refrigerant state intialize first time
        for i in range(self.NodNum):#inlet nodes
            if (self.Nod[i]['BranIN'][0]<0):#no inlet branch, only from the distributor
                Gr = mr/(P['Ax']*self.Nod[i]['OutNum']);    
                #if(errorLog.IsError()) {
                #errorLog.Add("StructCond::_CondCircuit","CondMan0");
                #return -1;}
                for j in range(self.Nod[i]['OutNum']):#states flowing out from the node
                    jj = self.Nod[i]['BranOUT'][j];#index of the outlet branches
                    self.Bra[jj]['HPi']=HPi;
                    self.Bra[jj]['HPo']=HPo;
                    self.Bra[jj]['Gr']=Gr*self.Bra[jj]['GrFac'];
                    self.Bra[jj]['Ini']=1;
                #end j loop
            #endif
        #end i loop
    
        while (abs(Res)>1e-7 and IterN<20):
            
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
            HPi['H']=H_in;#evaporator inlet refrigerant enthalpy
            HPi['P']=P_in;
    
            for i in range(self.BranNum):
                self.Bra[i]['Ini']=0
        
            for i in range(self.NodNum): #inlet nodes
                if (self.Nod[i]['BranIN'][0]<0):#from the distributor
                    Gr = mr/(P['Ax']*self.Nod[i]['OutNum']);    
                    for j in range(self.Nod[i]['OutNum']): #states flowing out from the node
                        jj = self.Nod[i]['BranOUT'][j];     #index of the outlet branches
                        self.Bra[jj]['HPi']['H']=HPi['H'];
            
                        if (self.Nod[i]['OutNum']==self.BranNum):#no middle nodes
                            DDP = (self.Bra[jj]['HPi']['P']-self.Bra[jj]['HPo']['P']);
                            self.Bra[jj]['HPi']['P']=HPo['P']+DDP;
                            P['DPr'][jj]=DDP;#output the pressure drop of each branch
                        else:
                            self.Bra[jj]['HPi']['P']=HPi['P'];
                        self.Bra[jj]['Gr']=Gr*self.Bra[jj]['GrFac'];
                        self.Bra[jj]['Ini']=1;
                    #end j circle
                #endif
            #end i loop
        
            self.Cal_HP(P,0,HPi);    #heat transfer and pressure drop calculation
        
            for i in range(self.BranNum):
                self.Bra[i]['Ini']=0;
        
            Gr=0;
            HPi['H']=0;
            HPi['P']=0;
               
            #nodes in the middle
            for i in range(self.NodNum):
                if (self.Nod[i]['BranIN'][0]>=0 and self.Nod[i]['BranOUT'][0]>=0): #nodes in the middle
                    for j in range(self.Nod[i]['InNum']):   #node inlet state
                        jj= self.Nod[i]['BranIN'][j];
                        Gr=Gr+self.Bra[jj]['Gr'];
                        HPi['H']=self.Bra[jj]['HPo']['H']*self.Bra[jj]['Gr']+HPi['H'];
                        HPi['P']=self.Bra[jj]['HPo']['P']*self.Bra[jj]['Gr']+HPi['P'];
        
                    HPi['H']=HPi['H']/Gr;
                    HPi['P']=HPi['P']/Gr;
                    Gr=Gr/self.Nod[i]['InNum'];
                    Gr=Gr*self.Nod[i]['InNum']/self.Nod[i]['OutNum'];
        
                    for j in range(self.Nod[i]['OutNum']):
                        jj = self.Nod[i]['BranOUT'][j];   #index of outlet branches
                        self.Bra[jj]['HPi']=HPi;
                        self.Bra[jj]['Gr']=Gr*self.Bra[jj]['GrFac'];
                        self.Bra[jj]['Ini']=1;
                    #end j circle
                #endif
            #end i circle
        
        
            self.Cal_HP(P,1,HPi);#heat transfer and pressure drop calculation
            
            Gr=0;
            HPi['H']=0;
            HPi['P']=0;
        
            #exit nodes
            for i in range(self.NodNum):
                if (self.Nod[i]['BranOUT'][0]<0):#no outlet branch except the suction line
                    for j in range(self.Nod[i]['InNum']):
                        jj=self.Nod[i]['BranIN'][j];
                        Gr=Gr+self.Bra[jj]['Gr'];
                        HPi['H']=self.Bra[jj]['HPo']['H']*self.Bra[jj]['Gr']+HPi['H'];
                        HPi['P']=self.Bra[jj]['HPo']['P']*self.Bra[jj]['Gr']+HPi['P'];
                        P['Hout8'][jj] = self.Bra[jj]['HPo']['H'];  #output the exit enthalpy of each branch
                    
                    HPi['H']=HPi['H']/Gr;
                    HPi['P']=HPi['P']/Gr;
                    Gr=Gr/self.Nod[i]['InNum'];
                    Gr=Gr*self.Nod[i]['InNum']/self.Nod[i]['OutNum'];
                #endif
            #end i circle
            
            if (self.RowNum==1):
                Res=0;
            else:
                Res=2*(HPi['H']-H_out)/(HPi['H']+H_out);
            H_out=HPi['H'];
            P_out=HPi['P'];
            P_in=HPo['P']+(P_in-P_out);
        
            
        if (abs(Res)>1e-5):
            print('_EvapCircuit_Fwd, Can NOT reach the required tolerance, Res= '+str(Res))
            raise
    
        Sm['m']=0;
        Sm['V']=0;
    
        for i in range(self.BranNum):
            Sm['m']=Sm['m']+self.Bra[i]['m']['m'];
            Sm['V']=Sm['V']+self.Bra[i]['m']['V'];
    
        WH_out={'W':0,'H':0};
        Ma_out=0.0;
    
        #air outputs
        for i in range(self.BranNum):
            if (self.Bra[i]['EqulNo']<0): #no equivalent circuits
                for j in range(self.Bra[i]['TubNum']):
                    Tubj=0;
                    if (self.Rev):   #counted from the point connected to compressor suction
                        Tubj=j;
                    else:
                        Tubj=self.Bra[i]['TubNum']-1-j;#counted from the point of evaporator entrance
                    TubeN=self.Bra[i]['TubNo'][Tubj];
    
                    if (self.Tub[TubeN]['RowNo']==self.RowNum-1)#outlet row tube at air side
                        for k in range(self.SegNum):
                            WH_out['W']=self.Tub[TubeN]['Ga']*self.Tub[TubeN]['Seg'][k]['WHo']['W']+WH_out['W'];
                            WH_out['H']=self.Tub[TubeN]['Ga']*self.Tub[TubeN]['Seg'][k]['WHo']['H']+WH_out['H'];
                            Ma_out=Ma_out+self.Tub[TubeN]['Ga'];
    
        WH_out['W'] = WH_out['W']/Ma_out;
        WH_out['H'] = WH_out['H']/Ma_out;
    
        TPo['T']=TPi['T']-5;
        TPo['P']=0.8;
        try:
            TPo= {'T': HAPropsSI('T','P',101325,'H',WH_out['H'],'W',WH_out['W']), 'P': HAPropsSI('R','P',101325,'H',WH_out['H'],'W',WH_out['W'])} #WHtoTP(WH_out,*TPo);
        except:
            print('check TPo that need to be numerically solved!')
            raise
        #if(errorLog.IsError()) {
        #errorLog.ClearError("StructEvap::_EvapCircuit","WH_out to TPo");
        #return -1;}
    
        #refrigerant outputs
        HPo['H']=H_out;
        HPo['P']=HPo['P'];
        HPi['H']=H_in;
        HPi['P']=P_in;
    
        #Sm->m=Sm->m + (2.87-2.766)+(3.10-3.05602)/(2.0792-7.955)*(P->VapL-7.955);#for three points charge tuning
    
        return 0
    
    
    def _EvapCircuit_Rev(self, mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,P):
        '''reversed evaporator solver'''
        pass
    
    def Cal_HP(self, P, Pos, HPi):
        ''''function to calculate the heat transfer and pressure drop'''
        
        #heat transfer and pressure drop of the two-phase region
        Logic=0;
        WHo = {'W':0,'H':0};
        Sm = {'m':0,'V':0};
        mi = {'m':0,'V':0};
    
        Bak=P.copy();#keeping the information of the evaporator
        
        if(Pos<0||Pos>1):
            print("Cal_HP, Wrong position")
        
    
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
        
            if Pos == 0: 
                if (self.Bra[i]['Ini']==1):#this branch has been initialized
                    Logic=1;
                else:
                    Logic=0;  
            elif Pos == 1:
                if(self.Bra[i]['Ini']==1):#this branch has been initialized
                    Logic=1;
                else:
                    Logic=0;
            
            if (Logic):
                if(Rev):
                    HPi=self.Bra[i]['HPo']; #opposite to flow direction
                else:
                    HPi=self.Bra[i]['HPi']; #paralell to flow direction
            
                if(self.Bra[i]['EqulNo']<0):    #no equivalent branch

                    for j in range(self.Bra[i]['TubNum']):
                        Tubj=0;
                        if (Rev):   #counted from the point connected to compressor suction
                            Tubj=j;
                        else:
                            Tubj=self.Bra[i]['TubNum']-1-j;#counted from the point of evaporator entrance
            
                        TubeN=self.Bra[i]['TubNo'][Tubj];
            
                        self.Tub[TubeN]['HPo']=HPi;
                        self.Tub[TubeN]['m']['m']=0;
                        self.Tub[TubeN]['m']['V']=0;
            
                        if (self.Tub[TubeN]['RowNo']>0): #not the first row, to get the air state  
                            Upper = self.Tub[TubeN]['AirUpstreamUpper'];
                            Lower = self.Tub[TubeN]['AirUpstreamLower'];
                            for k in range(self.SegNum):
                                if(Upper>=0 and Lower>=0):
                                    WHo['W']=(self.Tub[Upper]['Seg'][k]['WHo']['W']*self.Tub[Upper]['Ga']+self.Tub[Lower]['Seg'][k]['WHo']['W']*self.Tub[Lower]['Ga'])/(self.Tub[Upper]['Ga']+self.Tub[Lower]['Ga']);
                                    WHo['H']=(self.Tub[Upper]['Seg'][k]['WHo']['H']*self.Tub[Upper]['Ga']+self.Tub[Lower]['Seg'][k]['WHo']['H']*self.Tub[Lower]['Ga'])/(self.Tub[Upper]['Ga']+self.Tub[Lower]['Ga']);
                                elif (Upper>=0):
                                    WHo=self.Tub[Upper]['Seg'][k]['WHo'];
                                else:
                                    WHo=self.Tub[Lower]['Seg'][k]['WHo'];
            
                                try:
                                    self.Tub[TubeN]['Seg'][k]['TPi'] = {'T': HAPropsSI('T','P',101325,'H',WHo['H'],'W',WHo['W']), 'P': HAPropsSI('R','P',101325,'H',WHo['H'],'W',WHo['W'])} #WHtoTP(WHo,self.Tub[TubeN]['Seg'][k]['TPi']);
                                except:
                                    print('Ammar: I guess the error is becasue the humidity is so high and can not converge, try to write the numerical solver that Boshen used')
                                    raise
                            #end k circle
                        #endif
                 
                        for k in range(self.SegNum):
                            #for debugging
                            #if(k==9 and TubeN==13):    
                            #    shenb=0;
                
                            realk=0;
                            if (self.Tub[TubeN]['even']):
                                realk=(self.SegNum-1-k);
                            else:
                                realk=k;
                        
                            if (self.Rev):
                                EvapTubeL_Rev(self.Bra[i]['Gr'],HPi,self.Tub[TubeN]['Ga'],self.Tub[TubeN]['Seg'][realk]['TPi'],WHo,mi,P)
                            else:
                                EvapTubeL_Fwd(elf.Bra[i]['Gr'],HPi,self.Tub[TubeN]['Ga'],self.Tub[TubeN]['Seg'][realk]['TPi'],WHo,mi,P)
                
                            #if(errorLog.IsError()) {
                            #errorLog.Add("StructEvap::Cal_HP","EvapTube");
                            #return -1;}
                            self.Tub[TubeN]['Seg'][realk]['WHo']=WHo;
                            self.Tub[TubeN]['m']['m']=self.Tub[TubeN]['m']['m']+mi['m'];
                            self.Tub[TubeN]['m']['V']=self.Tub[TubeN]['m']['V']+mi['V'];
                            Sm['m']=Sm['m']+mi['m'];
                            Sm['V']=Sm['V']+mi['V'];
                        #end k circle
                        
                        
                        if (self.Rev):
                            EvapTubeBend(self.Bra[i]['Gr'],HPi,mi,P);
                        else:
                            EvapTubeBend_Fwd(self.Bra[i]['Gr'],HPi,mi,P);
                        #if(errorLog.IsError()) {
                        #errorLog.Add("StructCond::Cal_HP","CondTubeBend");
                        #return -1;}
                
                        self.Tub[TubeN]['m']['m']=self.Tub[TubeN]['m']['m']+mi['m'];
                        self.Tub[TubeN]['m']['V']=self.Tub[TubeN]['m']['V']+mi['V'];
                        Sm['m']=Sm['m']+mi['m'];
                        Sm['V']=Sm['V']+mi['V'];
                        if (self.Rev):
                            self.Tub[TubeN]['HPi']=HPi;
                        else:
                            self.Tub[TubeN]['HPo']=HPi;
                    #end j circle
            
                    #output of this branch
                    if(Rev):
                        self.Bra[i]['HPi']=HPi;
                    else:
                        self.Bra[i]['HPo']=HPi;
                    self.Bra[i]['m']=Sm;
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
                    if (self.Rev):
                        self.Bra[i]['HPi']=self.Bra[NoBra]['HPi'];
                    else:
                        self.Bra[i]['HPo']=self.Bra[NoBra]['HPo'];
                    self.Bra[i]['m']=self.Bra[NoBra]['m'];
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
                    #end NN loop
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
            
        #end i loop
    
        P=Bak;
        
        return 0;
    
    
#     char filename[128];      #data file name
#     int error;      #for keeping track of the errorLog
#     int TubeNum;          #number of all tubes
#     int BranNum;       # number of branch
#     int NodNum;        # number of node
#     int RowNum;           #number of Row
#     int SegNum;       #number of segment
#     int Rev;
#     double AreaFront;      #front area
#     double Volum;      #air flow volume
#     double *GaRow;      #average air flow distribution factor of each row
#     EvapBranch *Bra;      #array of branch objects
#     EvapNode *Nod;      #array of node objects
#     TubeEvap *Tub;      #array of tube objects
     

def EvapCircuit(type,mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D):

    EvapInfo = StructEvap("../InputDoc/EvapStruc.xlsx");
    EvapInfo.Rev=D['REV'];

    if type == 301: #user defined
        EvapInfo._EvapCircuit_Fwd(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D);
#     elif type == 302: #Carrier RTU at Purdue
#         EvapCircuit301(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D); 
#     elif type == 303: #Single finned tube
#         EvapCircuit302(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D);
#     elif type == 304: #this function is for CS30TR
#         EvapCircuit303(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D); 
    else:
        print("Evaporator model "+str(type)+" is not found");

    return 0


def Evaporator(filename,
               mr,#refrigerant mass flow rate
               HPo,#refrigerant outlet state
               GaI,#air mass flux
               TPi,#inlet air state
               HPi,#inlet refrigerant state of the segment
               TPo,#outlet air state
               Sm,#mass and volume in the evaporator
               Aflow,#air flow area
               Evap_struc,
               Prms, #adjustment factors
               NSeg,evapInit):
   
    '''    
    Model of a direct expansion evaporator  coil. The model is comprised 
    of small pieces of single finned tubes put together to form rows.  
    Rows of finned tubes are put together using elbows and manifolds to 
    make an entire coil. This function defines the coil's circuit pattern.
    Model parameters are stored in the file "evap.dat".
       
    INPUTS
    mr - refrigerant mass flow rate (kg/s).
    HPo - outlet refrigerant thermodynamic state (J/kg,kPa). 
    Ga - air mass flux (kg/s*m^2).
    TPi - inlet air thermodynamic state (K,-).
    Prms - array of three multiplicative correlation correction factors.
        Prms[0] - corrects air side convection coefficient (dry coil).
        Prms[1] - corrects refrigerant side convection coefficient.
        Prms[2] - corrects refrigerant side pressure drop correlation.
        Prms[3] - corrects air side convection coefficient (wet coil).
        Prms[4] - reduction in air flow rate (0-1) to simulate fouling.
        Set these parameters to 1 (Prms={1,1,1,1}) for no corrections.
       
    OUTPUTS
    HPi - inlet refrigerant thermodynamic state (J/kg,kPa).
    TPo - outlet air thermodynamic state (K,-).
    Sm - mass of charge in evaporator (kg) and volume of charge (m^3).
    Aflow = air flow area (ma=Ga*Aflow) (m^2).
    status - error flag.  Set to 1 to indicate error has occurred.
    '''
    D = Evap_struc.copy()
     
 
    if (evapInit > 0 or NSeg>0):
         
        df = pd.read_excel(filename,header = 0) #open the file and all the data are under column 'Evaporator'
         
        D['type'] = df.Evaporator[0]
        D['Di'] = df.Evaporator[1]             # tube inside diameter
        D['L'] = df.Evaporator[2]              # tube length 
        #printf("evap tube length = %lf\n",D.L);
        D['xp'] = df.Evaporator[3]             # tube thickness
        D['Df'] = df.Evaporator[4]              #outside fin diameter (m)
        D['z'] = df.Evaporator[5]            #it is used as the fin pitch
        D['th'] = df.Evaporator[6]          #spacing between tubes in bank (m)
        D['vsp'] = df.Evaporator[7]       #actually it is P_t
        D['P_l'] = df.Evaporator[8]         #spacing between tubs in the longitudual direction (m)
        D['NSeg'] = int(df.Evaporator[9])   #number of segments per tube in finite difference model
        if(NSeg>0): 
            D['NSeg'] = NSeg
        D['Brad']=df.Evaporator[10]     #radius of return bend (m)
         
        GaNom = df.Evaporator[11]      #nominal air mass flux (kg/s/m^2), correponding to the maximum air flux
         
        D['NBranchs'] = int(df.Evaporator[12]) # number of equivalent branches
        D['NBraTube'] = int(df.Evaporator[13]) # number of tube in each branch
        D['Nrows'] = int(df.Evaporator[14])        # number of rows high
        D['Ndeep'] = int(df.Evaporator[15])        # number of rows deep
        D['Frontal_A'] = df.Evaporator[16]    # frontal area
         
        #correction factors
        hAirAdj = df.Evaporator[17]
        hRefAdj = df.Evaporator[18]
        PRefAdj = df.Evaporator[19]
        WAirAdj = df.Evaporator[20]
        FoulFac = df.Evaporator[21]
   
        #B.S. ------------------------
        #new input parameter for evaporator
        D['microfin'] = int(df.Evaporator[22])     #microfin type, 0=smooth tube, 1=helical, 2=cross-grooved, 3=herringbone
        D['finN'] = df.Evaporator[23]              #fin number in a micro-fin tube
        D['gama'] = df.Evaporator[24]              #fin apex angle in a micro-fin tube
        D['beta'] = df.Evaporator[25]              #fin helix angle in a micro-fin tube
        D['finH'] = df.Evaporator[26]              #fin height in a micro-fin tube
        D['w_b'] = df.Evaporator[27]                #base width of a single fin
        D['w_e'] = df.Evaporator[28]                #top width of a single fin
        D['w_z'] = df.Evaporator[29]                #base distance between two neighboring fins
        # airside fin type
        D['airFin'] = int(df.Evaporator[30])       #1-plain, 2-corrugated, 3-slit, 4-louvered, 5-convex louvered, 6-smooth wavy, 7-spine
        D['sub1'] = df.Evaporator[31]              #fin substructures, for the possible use later, such as the waffle height of the wavy fin, louver number of the louverd fin
        D['sub2'] = df.Evaporator[32]               #fin substructures
        D['sub3'] = df.Evaporator[33]                #fin substructures
        D['sub4'] = df.Evaporator[34]                #fin substructures
        D['sub5'] = df.Evaporator[35]                #fin substructures
        D['K_T'] = df.Evaporator[36]                #400-- tube wall conductance which is the conductance factor of copper, for calculating tube resistance
        D['K_F'] = df.Evaporator[37]                    #237-- fin conductance
        #------------------------------B.S.
       
        endSign=int(df.Evaporator[38])      #-256
        if(endSign != -256):
            print ("Evaporator parameter input wrong")
            raise
   
        evapInit = 0
         
        D['Ls']=D['L']/float(D['NSeg'])
           
        #B.S.-----------------------------------
        # tube inside geometry
        if (D['microfin'] > 0):
            D['D_b'] = D['Di'] + 2*D['finH']        #tube diameter at the base of the fin
            D['Do'] = D['D_b'] + 2*D['xp']           #Pipe outside diameter.
            D['D_m'] = (D['Di'] + D['D_b'])/2       #mean diameter of the micro-fin tube
            D['P_H'] = D['finN']*(D['w_z']+2*pow((pow(D['finH'],2.0)+pow((0.5*(D['w_b']-D['w_e'])),2.0)),0.5)+D['w_e'])# the hydraulical circumference
            D['Acs'] = pi*pow((D['D_b']/2.0),2.0)-D['finN']*0.5*D['finH']*(D['w_b']+D['w_e'])#cross area of the micro-fin tube, this is the actual cross-section area
            D['Dh_i'] = 4*D['Acs']/D['P_H']               #Inside hydraulical diameter 
            D['Ax'] = pi/4.0*D['Di']*D['Di']              #Inside pipe cross sectional area, based on fin tips
            D['Api'] = pi*D['Di']*D['Ls']                 #Inside pipe surface area (one tube segment), based on fin tips
        else:
            D['finH']=0                     #tube without fins
            D['D_b'] = D['Di']              #tube diameter at the base of the fin
            D['Do'] = D['D_b']+2*D['xp']          #Pipe outside diameter.
            D['D_m'] = D['Di']                 #averge diameter of the micro-fin tube
            D['Ax'] = pi/4.0*D['Di']*D['Di']      # Inside pipe cross sectional area
            D['Api'] = pi*D['Di']*D['Ls']         # Inside pipe surface area (one tube segment)
            D['P_H'] = pi*D['D_b']             # the hydraulical circumference
            D['Acs'] = D['Ax']                 #cross area of the micro-fin tube, this is the actual cross-section area
            D['Dh_i'] = D['Di']                #inside hydraulical diameter
        #-----------------------------------------B.S.
        D['Vs'] = D['Acs']*D['Ls']
        D['Blen'] = pi*D['Brad']
        D['BVs'] = D['Blen']*D['Ax']
   
        ##*External Tube Parameters*##
        #B.S.-----------------------------
        #get the fin reference length, according to Schimidt equation
        M = D['vsp']/2e0
        R_O = D['Do']/2
        L = 0.5e0*pow((pow(M,2e0)+pow(D['P_l'],2e0)),0.5)
        BETA = L/M
        PHI = M/R_O
        R_EFF = R_O*1.27e0*PHI*pow((BETA-0.3e0),0.5e0)
        D['Df'] = 2.0*R_EFF                             #confirm the outside reference fin diameter
        D['L_F'] = (R_EFF-R_O)*(1e0+0.35e0*log(R_EFF/R_O))   #get the fin reference length for calculating the fin efficiency, according to Schimidt equation
        #--------------------------------------B.S.
   
        D['y'] = (D['Df']-D['Do'])/2        #Distance from outside of pipe to edge of fin. 
        D['N'] = int(D['Ls']/D['z'])        #B.S., change the original, (D.Ls/(D.z+D.th))    
        D['Apo'] = pi*D['Do']*(D['Ls']-D['N']*D['th'])       #B.S., correct pi*D.Do*(D.Ls)
   
        x = D['Df'] + D['th']/2.0
        D['Af'] = D['N']*pi/2.0*(x*x-D['Do']*D['Do'])                   #fin has faces
        D['Aflow'] = (D['vsp']-D['Do'])*(D['Ls']-D['N']*D['th'])        #minimum airflow area per segment
        # printf("Minimum free flow area = %lf (m^2)\n",D.Aflow)
         
        #* hydraulic diameter *#
        D['Dh']=4*D['Aflow']/(D['Af']+D['Apo'])*D['vsp']     
     
   
    D['hAirAdj'] = Prms[0]*hAirAdj
    D['hRefAdj'] = Prms[1]*hRefAdj
    D['PRefAdj'] = Prms[2]*PRefAdj
    D['WAirAdj'] = Prms[3]*WAirAdj
   
    if (GaI <= 0):
        Ga = FoulFac*GaNom*Prms[4]
    else:
        Ga = GaI
       
    # Calculate air cfm after Ga is established.
    A_total = D['Aflow']*D['NSeg']*float(D['Nrows']) #total area air flows throguh for entire coil (m^2)
    v_air = HAPropsSI("V", "T", TPi['T'], "P", 101325, "R", TPi['P']) # specific volume of air [m^3/kg dry air]
    vfra = Ga*A_total*v_air # nominal volumetric flow rate air
    D['cfma'] = vfra*(60.0/((12.0*0.0254)*(12.0*0.0254)*(12.0*0.0254))) # nominal cfm of air flow
   
    #B.S., prepare to set up the simple model (A.B., this might not be necessary because they are already initialized with 0 values in Evap_Sturc dictionary
    D['V_TOT']=0; D['V_TP']=0; D['V_Liq']=0; D['V_Vap']=0;
    D['m_TOT']=0; D['m_TP']=0; D['m_Liq']=0; D['m_Vap']=0;
    D['mr']=0; D['ma_TOT']=0;
    D['UA_TOT']=0; D['UA_Liq']=0; D['UA_Vap']=0; D['UA_TP']=0;
    D['UAw_TOT']=0; D['UAw_Liq']=0; D['UAw_Vap]']=0; D['UAw_TP']=0;
    D['LiqL']=0; D['TPL']=0; D['VapL']=0;
    D['L_dry']=0; D['L_wet']=0;
    D['count1']=0; D['count2']=0;
    D['HP_TP1']['H']=0; D['HP_TP1']['P']=0;
    D['HP_TP2']['H']=0; D['HP_TP2']['P']=0;
    D['Qtp_dry']=0; D['Qtp_wet']=0;
    #-----------------------------------------B.S.
       
    D['GetP'] = 0     #don't calculate the airside pressure drop in the following procedure
       
    D['REV'] = Evap_struc['REV']    #iteration loop direction
     
    EvapCircuit(D['type'],mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,D)
   
    DP_circuit = Circuit_DP(Ga,TPi,TPo,D)#B.S., calculate the airside pressure drop across the circuit
    #gOutData.Eva_AirDP = DP_circuit#for output
   
    #B.S.---------------------------------------
    D['HP_out'] = HPo      #suciton site of the evaporator
    D['HP_in'] = HPi       #entrance of the evaporator
    D['mr'] = mr
    D['Ga'] = Ga
    D['TPi'] = TPi
    D['ma_TOT'] = Ga*D['Aflow']*D['NSeg']*D['Nrows']
    #Build_Eva(&D)        #correlate the moving boundary and lumped evaporator model
    Evap_struc = D          #output the evaporator construct
    #---------------------------------------B.S.
  
    return Evap_struc

    
def ExerciseEvaporator():
    
    Ref = 'R410A'
    #design conditions
    SH = 23.004 #superheat (F)
    #SC = 10.9512 #subcooling F (***no effect)
    ET_SL = 48.3264-23 # F, evaporation temperature at the suction line
    ET_DT = 45+4 # F, add 4F due to evaporator coild pressure drop
    #CTOA_DL = 29.511+4 # F, add 4F due to pressure drop in cond coil, 18F is at the liquid line (** no changes)
    #CTOA_LL = 29.511 # F, at the liquid line  (** no changes)
    AMB = 95 # F
    RA = 83.048#80 # indoor room dry bulb temperature (F)
    RARH = 0.51 # indoor room humidity (-)
    RAWB = 67.856 #  indoor room wetbulb temoperature (F)
    MR = 955.572 # lbm/hr, refrigerant mass flow rate
    #GOA = 3.82 # kg/s/m^2, outdoor air mass flux (this is for condenser component)
    GIA = 5.6 # kg/s/m^2, indoor air mass flux (0.9237/0.3623)
    #DSH = 67.13 # F, discharge superheat (***nothing changes)
    
    
    #initialize the values
    hpi = {'H':0.0,'P':0.0}
    m = {'m':0.0,'V':0.0}
    TPo = {'T':0.0,'P':0.0}
    Aflow = 0.0
    Evap_prms = [1,1,1,1,1] #correction factor (see Evaporator model for more details)
    NSeg = -1
    Evap_struc={'Di':0.0,'L':0.0,'xp':0.0,'Df':0.0,'z':0.0,'th':0.0,'vsp':0.0,'Brad':0.0,
                'NSeg':int(0),
                'Dh':0.0,'Do':0.0,'y':0.0,'Ls':0.0,'N':0.0,'Ax':0.0,'Api':0.0,'Apo':0.0,'Af':0.0,'Aflow':0.0,'Vs':0.0,'Blen':0.0,'BVs':0.0,
                'Ro:':0.0,
                'Gr':0.0,'Ga':0.0,
                'HPo':{'H':0.0,'P':0.0},
                'TXPo':{'T':0.0,'X':0.0,'P':0.0},
                'TPi':{'T':0.0,'P':0.0},
                'hAirAdj':0.0,'hRefAdj':0.0,'PRefAdj':0.0,'WAirAdj':0.0,
                'type':int(0),
                'Nrows':int(0),'Ndeep':int(0),
                'NBranchs':int(0), 'NBraTube':int(0),
                'P_l':0.0, #spacing between tubs in the longitudual direction (m)
                'microfin':int(0), #new parameters for micro-fin tubes, and specially configured fins #decide whether is micro-fin tubes microfin=1
                'w_b':0.0, #base length for single micro-fin
                'w_e':0.0, #top length for single micro-fin
                'w_z':0.0,    #width between the neighboring micro-fins
                'finH':0.0,    #micro-fin height
                'beta':0.0,#micro-fin helical angle
                'gama':0.0,    #micro-fin apex angle
                'finN':0.0,    #total micro-fin number
                'Acs':0.0,    # micro-fin cross-sectional area
                'P_H':0.0,    #micro-fin hydraulic circumference
                'Dh_i':0.0,    #micro-fin tube hydraulic diameter
                'K_T':0.0,    #conductance factor of tube wall
                'K_F':0.0,    #conductance factor of fin material
                'L_F':0.0,    #reference fin length for Schimidt fin efficiency calculation
                'D_b':0.0,'D_m':0.0,    #fin base diameter, and mean diameter
                'airFin':int(0),
                'sub1':0.0,'sub2':0.0,'sub3':0.0,'sub4':0.0,'sub5':0.0,#for inputing sub-structures of fin surface
                'ho':0.0, 'wetadj':0.0,#airside heat transfer coefficient
                'Frontal_A':0.0,#frontal area
                'GetP':0.0,#calculate the airside pressure drop
                #variables for generating the simple evaporator function
                'V_TOT':0.0, 'V_TP':0.0, 'V_Liq':0.0, 'V_Vap':0.0,#inner volume of different phase
                'L_TOT':0.0,'LiqL':0.0,'VapL':0.0,'TPL':0.0,#tube length of different phase
                'L_dry':0.0, 'L_wet':0.0,#tube length of dry and wet heat transfer
                'A_TOT':0.0,'A_Liq':0.0,'A_Vap':0.0,'A_TP':0.0,#heat transfer surface area of different phase
                'm_TOT':0.0, 'm_TP':0.0, 'm_Liq':0.0, 'm_Vap':0.0,#airflow rate across different phase
                'rho_TOT':0.0, 'rho_TP':0.0, 'rho_Liq':0.0, 'rho_Vap':0.0,#average density of different phase
                'U_TP':0.0, 'U_Liq':0.0, 'U_Vap':0.0,#averge dry heat conductance of different phase 
                'Uw_TP':0.0, 'Uw_Liq':0.0, 'Uw_Vap':0.0,#average wet heat conductance of different phase
                'DP_TOT':0.0, 'DP_TP':0.0, 'DP_Liq':0.0, 'DP_Vap':0.0,#average pressure gradient of different phase
                'UA_TOT':0.0, 'UA_Liq':0.0,'UA_Vap':0.0,'UA_TP':0.0,#overall dry heat conductance of different phase
                'UAw_TOT':0.0, 'UAw_Liq':0.0,'UAw_Vap':0.0,'UAw_TP':0.0,#overall wet heat conductance of different phase
                'mr':0.0, 'ma_TOT':0.0,#overall refrigerant and air mass flow rate
                'Ga_meanL':0.0,#average air flow rate per tube length
                'HP_out':{'H':0.0,'P':0.0}, 'HP_dry':{'H':0.0,'P':0.0}, 'HP_wet':{'H':0.0,'P':0.0}, 'HP_in':{'H':0.0,'P':0.0},'HP_TP1':{'H':0.0,'P':0.0},'HP_TP2':{'H':0.0,'P':0.0},#state parameters at important locations
                'count1':0.0,'count2':0.0,#count the state points of two-phase flow begin point and end point 
                'Qtp_dry':0.0, 'Qtp_wet':0.0,#two-phase dry heat transfer and wet heat transfer amount
                'r_v':0.0, 'r_l':0.0,'r_tp':0.0,#parameters for adjusting the theoretical heat transfer effectiveness of each phase
                'H2_residual':0.0,#for the consistency of the moving boundary model analysis
                #------------------------------B.S.
                'q_flux':0.0, 'T_w':0.0,#segment heat flux and inside tube wall temperature
                'wet':int(0),#wet=0 to calculate dry heat transfer, wet=1 to calculate wet heat transfer
                'REV':int(0),
                'cfma':0.0,
                'Hout8':[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],'DPr':[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],#superheat and pressure drop at each evaporator branch (array of 10 elements)
                }

    #inlet air state
    TPi = {'T':F2K(RA), 'P':RARH} #T: air inlet (return air) dry bulb temperature [K], P: air inlet (return air) relative humidity [-]
    #air mass flux (kg/s/m^2)
    Ga = GIA
    #refrigerant mass flow rate [kg/s]
    mr = lbh2kgs(MR)
    #outlet refrigerant state
    ETBRA = DeltaF2K(RA-ET_SL) #evaporating temperature below return air [K]
    Tsh = DeltaF2K(SH) #suction line superheat [K]
    
    Tsat_sl = TPi['T']-ETBRA #Temperature of saturated liquid [K]
    P_sl = PropsSI("P", "T", Tsat_sl, "Q", 1, Ref)
    txpo = {'T':Tsat_sl+Tsh,'X':1,'P':P_sl} #initially assume outlet temperature equals to saturated liquid pressure (txpo_P)
    hpo = {'H': PropsSI("H", "T", txpo['T'], "P", txpo['P'], Ref),'P':txpo['P']}  #outlet refregerant enthalpy
    Twbi = HAPropsSI("Twb", "T", TPi['T'], "P", 101325, "R", TPi['P']) #inlet air wet bulb temperature 

    #*Run the evaporator model*#
    Evaporator('../InputDoc/Evap_Data.xlsx',mr,hpo,Ga,TPi,hpi,TPo,m,Aflow,Evap_struc,Evap_prms,NSeg,evapInit=1)
    
    #inlet refrigerant state and   
    #saturation temperature computed from the inlet refrigerant state
    try:
        txpi = {'T':PropsSI("T", "H", hpi['H'], "P", hpi['P'], Ref),'X':PropsSI("Q", "H", hpi['H'], "P", hpi['P'], Ref),'P':hpi['P']}
        Tsat_dt =  PropsSI("T", "Q", txpi['X'], "P", txpi['P'], Ref) #[K]
    except:
        print('this is two-phase, need to find quality first, Ammar you need to change the previous line to calculate temperature from quality')
        raise
     
    Twbo = HAPropsSI("Twb", "T", TPo['T'], "P", 101325, "R", TPo['P']) #outlet air wet bulb temperature 
     
    #Calculate the total heat transfer and then the UA = Qdot/DT
    Qdot = W2BTUh(mr*(hpo['H']-hpi['H'])) #[BTU/hr] - for all 4 sections of the heat exchanger
    DT = DeltaK2F(TPi['T'] - Tsat_sl) #[F]
    UA = Qdot/DT #[BTU/hr/F] - for all 4 sections of the heat exchanger 
         
    #print inputs
    print("Tai (F), RHai (%)", C2F(TPi['T']),TPi['P']*100.0)
    print("CFM (cfm)",Evap_struc['cfma'])
    print("mr (lbm/hr)", MR)
    print("ETBRA_sl (F)",DeltaK2F(TPi['T']-Tsat_sl))
    print("Tsh_sl (F)",DeltaK2F(txpo['T']-Tsat_sl))
    print("X_sl (0-1)",txpo['X'])
    print("NSeg",NSeg)
     
    #print outputs
    print("ETBRA_dt (F)",DeltaK2F(TPi['T']-Tsat_dt)) #ETBRA_dt
    print("X_dt (0-1)",txpi['X'])
    print("dP (psi)",kPa2psi(txpi['P']-txpo['P'])) ###BEWARE that txpo_P was set to saturatewd liquid pressure, check if the code modifies it
    print("DTa (F)",DeltaK2F(TPi['T']-TPo['T']))
    print("DTwb (F)",DeltaK2F(Twbi-Twbo))
    print("mass (kg)",m['m']) #charge
    print("Qdot (BTU/hr)",Qdot)
    print("UA (BTU/hr/F)",UA)
    print("Hin (kJ/kg)",hpi['H']/1000.0)
    print("Hout (kJ/kg)",hpo['H']/1000.0)
    print("Tsat_sl (C)",K2C(Tsat_sl))
    print("Tsat_dt (C)",K2C(Tsat_dt))

    return EvapD
    
def ExerciseComponents():
    ExerciseEvaporator()

if __name__=='__main__':    
#     ExerciseComponents()
    df = pd.read_excel("../InputDoc/EvapStruc.xlsx",sheetname='tube', header = 0)
    print (df)
