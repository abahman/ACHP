from __future__ import division, print_function, absolute_import
from math import pi,log,exp

#from scipy.optimize import brentq #solver to find roots (zero points) of functions
#import numpy as np
import pandas as pd

#import CoolProp as CP
#from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI

from EVAP import EvapTubeBend, EvapTubeBend_Fwd

def EvapCircuit(Ref,type,mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D):

    EvapInfo = StructEvap("../InputDoc/EvapStruc.xlsx", Ref);
    EvapInfo.Rev=D['REV'];

    if type == 301: #user defined
        EvapInfo._EvapCircuit_Fwd(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D)
#     elif type == 302: #Carrier RTU at Purdue
#         EvapCircuit301(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D); 
#     elif type == 303: #Single finned tube
#         EvapCircuit302(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D);
#     elif type == 304: #this function is for CS30TR
#         EvapCircuit303(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D); 
    else:
        print("Evaporator model "+str(type)+" is not found")

    return (EvapInfo.HPo, EvapInfo.HPi, EvapInfo.TPo, EvapInfo.Sm, EvapInfo.Aflow, EvapInfo.D)

class StructEvap():
    
    def __init__(self,filename,Ref):
        self.REV = 0 #reverse calculation
        self.Ref = Ref
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
        self.TubeNum = int(df.EvapStruc[3]);         #overall tube Number
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
          
    def _EvapCircuit_Fwd(self,mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,P):
        
        #initialize as those will be in the output 
        self.HPo = HPo
        self.HPi = HPi
        self.TPo = TPo
        self.Sm = Sm
        self.Aflow = Aflow
        self.P = P
        
        Gr=0.0;#mass flux
        
        if (self.Rev): #reversed calculation
            self._EvapCircuit_Rev(mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,P)
        return 0
    
        self.TPo=TPi;#inlet air state
    
        H_in=self.HPi['H'];#inlet refrigerant enthalpy
        
        #air side initializing
        rho_air=1/HAPropsSI("V", "T", TPi['T'], "P", 101325, "R", TPi['P']) #[kg dry air/ m^3]
    
        #const double Ma =Ga*rho_air*4.719e-4; 
        Ga=Ga/self.AreaFront*self.P['vsp']*self.P['Ls']/((self.P['vsp']-self.P['Do'])*(self.P['Ls']-self.P['N']*self.P['th']));
    
    
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
        DP = 20.0;
        P_in=0.0;
        P_out=0.0;
        self.HPi['P']=self.HPo['P']+DP;   #temperary inlet refrigerant pressure
        self.HPi['H']=H_in;          #evaporator inlet refrigerant enthalpy
        #*HPi is the intermediate variable
        P_in = self.HPo['P']+DP;       #temperary inlet refrigerant pressure
        P_out = self.HPo['P'];         #evaporator outlet refrigerant pressure
        IterN=0;

    
        #refrigerant state intialize first time
        for i in range(self.NodNum):#inlet nodes
            if (self.Nod[i]['BranIN'][0]<0):#no inlet branch, only from the distributor
                Gr = mr/(self.P['Ax']*self.Nod[i]['OutNum']);    
                #if(errorLog.IsError()) {
                #errorLog.Add("StructCond::_CondCircuit","CondMan0");
                #return -1;}
                for j in range(self.Nod[i]['OutNum']):#states flowing out from the node
                    jj = self.Nod[i]['BranOUT'][j];#index of the outlet branches
                    self.Bra[jj]['HPi']=self.HPi;
                    self.Bra[jj]['HPo']=self.HPo;
                    self.Bra[jj]['Gr']=Gr*self.Bra[jj]['GrFac'];
                    self.Bra[jj]['Ini']=1;
                #end j loop
            #endif
        #end i loop
    
        while True:
            
            self.P['VapL'] =0;
            self.P['TPL'] = 0;
            self.P['LiqL'] = 0;
            self.P['V_Vap'] =0;
            self.P['V_TP'] =0;
            self.P['V_Liq'] = 0;
            self.P['m_Vap'] = 0;
            self.P['m_TP'] = 0;
            self.P['m_Liq'] = 0;    
            self.P['UA_Vap'] = 0;
            self.P['UA_TP'] = 0;
            self.P['UA_Liq'] = 0;
    
            IterN=IterN+1;
            self.HPi['H']=H_in; #evaporator inlet refrigerant enthalpy
            self.HPi['P']=P_in;
    
            for i in range(self.BranNum):
                self.Bra[i]['Ini']=0
        
            for i in range(self.NodNum): #inlet nodes
                if (self.Nod[i]['BranIN'][0]<0):#from the distributor
                    Gr = mr/(self.P['Ax']*self.Nod[i]['OutNum']);    
                    for j in range(self.Nod[i]['OutNum']): #states flowing out from the node
                        jj = self.Nod[i]['BranOUT'][j];     #index of the outlet branches
                        self.Bra[jj]['HPi']['H']=self.HPi['H'];
            
                        if (self.Nod[i]['OutNum']==self.BranNum):#no middle nodes
                            DDP = (self.Bra[jj]['HPi']['P']-self.Bra[jj]['HPo']['P']);
                            self.Bra[jj]['HPi']['P']=self.HPo['P']+DDP;
                            P['DPr'][jj]=DDP;#output the pressure drop of each branch
                        else:
                            self.Bra[jj]['HPi']['P']=self.HPi['P'];
                        self.Bra[jj]['Gr']=Gr*self.Bra[jj]['GrFac'];
                        self.Bra[jj]['Ini']=1;
                    #end j circle
                #endif
            #end i loop
        
            self.Cal_HP(0);    #heat transfer and pressure drop calculation
        
            for i in range(self.BranNum):
                self.Bra[i]['Ini']=0;
        
            Gr=0;
            self.HPi['H']=0;
            self.HPi['P']=0;
               
            #nodes in the middle
            for i in range(self.NodNum):
                if (self.Nod[i]['BranIN'][0]>=0 and self.Nod[i]['BranOUT'][0]>=0): #nodes in the middle
                    for j in range(self.Nod[i]['InNum']):   #node inlet state
                        jj= self.Nod[i]['BranIN'][j];
                        Gr=Gr+self.Bra[jj]['Gr'];
                        HPi['H']=self.Bra[jj]['HPo']['H']*self.Bra[jj]['Gr']+HPi['H'];
                        HPi['P']=self.Bra[jj]['HPo']['P']*self.Bra[jj]['Gr']+HPi['P'];
        
                    self.HPi['H']=self.HPi['H']/Gr;
                    self.HPi['P']=self.HPi['P']/Gr;
                    Gr=Gr/self.Nod[i]['InNum'];
                    Gr=Gr*self.Nod[i]['InNum']/self.Nod[i]['OutNum'];
        
                    for j in range(self.Nod[i]['OutNum']):
                        jj = self.Nod[i]['BranOUT'][j];   #index of outlet branches
                        self.Bra[jj]['HPi']=self.HPi;
                        self.Bra[jj]['Gr']=Gr*self.Bra[jj]['GrFac'];
                        self.Bra[jj]['Ini']=1;
                    #end j circle
                #endif
            #end i circle
        
        
            self.Cal_HP(1);#heat transfer and pressure drop calculation
            
            Gr=0;
            self.HPi['H']=0;
            self.HPi['P']=0;
        
            #exit nodes
            for i in range(self.NodNum):
                if (self.Nod[i]['BranOUT'][0]<0):#no outlet branch except the suction line
                    for j in range(self.Nod[i]['InNum']):
                        jj=self.Nod[i]['BranIN'][j];
                        Gr=Gr+self.Bra[jj]['Gr'];
                        HPi['H']=self.Bra[jj]['HPo']['H']*self.Bra[jj]['Gr']+HPi['H'];
                        HPi['P']=self.Bra[jj]['HPo']['P']*self.Bra[jj]['Gr']+HPi['P'];
                        P['Hout8'][jj] = self.Bra[jj]['HPo']['H'];  #output the exit enthalpy of each branch
                    
                    self.HPi['H']=self.HPi['H']/Gr;
                    self.HPi['P']=self.HPi['P']/Gr;
                    Gr=Gr/self.Nod[i]['InNum'];
                    Gr=Gr*self.Nod[i]['InNum']/self.Nod[i]['OutNum'];
                #endif
            #end i circle
            
            if (self.RowNum==1):
                Res=0;
            else:
                Res=2*(self.HPi['H']-H_out)/(self.HPi['H']+H_out);
            H_out=self.HPi['H'];
            P_out=self.HPi['P'];
            P_in=self.HPo['P']+(P_in-P_out);
            
            if (abs(Res)<1e-7 and IterN>20): #condition to break the while loop
                break
        #end while loop
            
        if (abs(Res)>1e-5):
            print('_EvapCircuit_Fwd, Can NOT reach the required tolerance, Res= '+str(Res))
            raise
    
        self.Sm['m']=0;
        self.Sm['V']=0;
    
        for i in range(self.BranNum):
            self.Sm['m']=self.Sm['m']+self.Bra[i]['m']['m'];
            self.Sm['V']=self.Sm['V']+self.Bra[i]['m']['V'];
    
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
    
                    if (self.Tub[TubeN]['RowNo']==self.RowNum-1):#outlet row tube at air side
                        for k in range(self.SegNum):
                            WH_out['W']=self.Tub[TubeN]['Ga']*self.Tub[TubeN]['Seg'][k]['WHo']['W']+WH_out['W'];
                            WH_out['H']=self.Tub[TubeN]['Ga']*self.Tub[TubeN]['Seg'][k]['WHo']['H']+WH_out['H'];
                            Ma_out=Ma_out+self.Tub[TubeN]['Ga'];
    
        WH_out['W'] = WH_out['W']/Ma_out;
        WH_out['H'] = WH_out['H']/Ma_out;
    
        self.TPo['T']=TPi['T']-5;
        self.TPo['P']=0.8;
        try:
            #self.TPo = WHtoTP(WH_out,self.TPo)
            self.TPo = {'T': HAPropsSI('T','P',101325,'H',WH_out['H'],'W',WH_out['W']), 'P': HAPropsSI('R','P',101325,'H',WH_out['H'],'W',WH_out['W'])} 
        except:
            print('_EvapCircuit_Fwd :: check TPo that need to be numerically solved!')
            print('Other source of exception might be due to relative humidity ~100%, try to solve TPo with relative humidity of 95%')
            self.TPo = {'T': HAPropsSI('T','P',101325,'H',WH_out['H'],'R',0.995), 'P': 0.995}
        #if(errorLog.IsError()) {
        #errorLog.ClearError("StructEvap::_EvapCircuit","WH_out to TPo");
        #return -1;}
    
        #refrigerant outputs
        self.HPo['H']=H_out;
        self.HPo['P']=self.HPo['P'];
        self.HPi['H']=H_in;
        self.HPi['P']=P_in;
    
        #Sm->m=Sm->m + (2.87-2.766)+(3.10-3.05602)/(2.0792-7.955)*(P->VapL-7.955);#for three points charge tuning
    
        return 0
    
    def _EvapCircuit_Rev(self, mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,P):
        '''reversed evaporator solver'''
        pass
    
    def Cal_HP(self, Pos):
        ''''function to calculate the heat transfer and pressure drop'''
        
        #heat transfer and pressure drop of the two-phase region
        Logic=0;
        WHo = {'W':0,'H':0};
        Sm = {'m':0,'V':0};
        mi = {'m':0,'V':0};
        
        Bak=self.P.copy();#keeping the information of the evaporator
        
        if(Pos<0 or Pos>1):
            print("Cal_HP, Wrong position")
        
    
        for i in range(self.BranNum):
            #branch parameter clear zero
            self.P['VapL'] =0;
            self.P['TPL'] = 0;
            self.P['LiqL'] = 0;
            self.P['V_Vap'] =0;
            self.P['V_TP'] =0;
            self.P['V_Liq'] = 0;
            self.P['m_Vap'] = 0;
            self.P['m_TP'] = 0;
            self.P['m_Liq'] = 0;    
            self.P['UA_Vap'] = 0;
            self.P['UA_TP'] = 0;
            self.P['UA_Liq'] = 0;
        
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
                if(self.Rev):
                    self.HPi=self.Bra[i]['HPo']; #opposite to flow direction
                else:
                    self.HPi=self.Bra[i]['HPi']; #paralell to flow direction
            
                if(self.Bra[i]['EqulNo']<0):    #no equivalent branch

                    for j in range(self.Bra[i]['TubNum']):
                        Tubj=0;
                        if (self.Rev):   #counted from the point connected to compressor suction
                            Tubj=j;
                        else:
                            Tubj=self.Bra[i]['TubNum']-1-j;#counted from the point of evaporator entrance
            
                        TubeN=self.Bra[i]['TubNo'][Tubj];
            
                        self.Tub[TubeN]['HPo']=self.HPi;
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
                                    #self.Tub[TubeN]['Seg'][k]['TPi'] = WHtoTP(WHo,self.Tub[TubeN]['Seg'][k]['TPi']);
                                    self.Tub[TubeN]['Seg'][k]['TPi'] = {'T': HAPropsSI('T','P',101325,'H',WHo['H'],'W',WHo['W']), 'P': HAPropsSI('R','P',101325,'H',WHo['H'],'W',WHo['W'])}
                                except:
                                    print('Cal_HP :: check self.Tub that need to be numerically solved!')
                                    print('Other source of exception might be due to relative humidity ~100%, try to solve TPo with relative humidity of 95%')
                                    self.Tub[TubeN]['Seg'][k]['TPi'] = {'T': HAPropsSI('T','P',101325,'H',WHo['H'],'R',0.995), 'P': 0.995}
                                    
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
                                EvapTubeL_Rev(self.Bra[i]['Gr'],self.HPi,self.Tub[TubeN]['Ga'],self.Tub[TubeN]['Seg'][realk]['TPi'],WHo,mi,self.P)
                            else:
                                EvapTubeL_Fwd(self.Bra[i]['Gr'],self.HPi,self.Tub[TubeN]['Ga'],self.Tub[TubeN]['Seg'][realk]['TPi'],WHo,mi,self.P)
                
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
                            self.HPi, mi, self.P = EvapTubeBend(self.Ref,self.Bra[i]['Gr'],self.HPi,mi,self.P);
                        else:
                            self.HPi, mi, self.P = EvapTubeBend_Fwd(self.Ref,self.Bra[i]['Gr'],self.HPi,mi,self.P);
                        #if(errorLog.IsError()) {
                        #errorLog.Add("StructCond::Cal_HP","CondTubeBend");
                        #return -1;}
                
                        self.Tub[TubeN]['m']['m']=self.Tub[TubeN]['m']['m']+mi['m'];
                        self.Tub[TubeN]['m']['V']=self.Tub[TubeN]['m']['V']+mi['V'];
                        Sm['m']=Sm['m']+mi['m'];
                        Sm['V']=Sm['V']+mi['V'];
                        if (self.Rev):
                            self.Tub[TubeN]['HPi']=self.HPi;
                        else:
                            self.Tub[TubeN]['HPo']=self.HPi;
                    #end j circle
            
                    #output of this branch
                    if(self.Rev):
                        self.Bra[i]['HPi']=self.HPi;
                    else:
                        self.Bra[i]['HPo']=self.HPi;
                    self.Bra[i]['m']=Sm;
                    Sm['m']=0;
                    Sm['V']=0;
                    self.Bra[i]['Para_Struc'][0]=self.P['VapL'];
                    self.Bra[i]['Para_Struc'][1]=self.P['TPL'];
                    self.Bra[i]['Para_Struc'][2]=self.P['LiqL'];
                    self.Bra[i]['Para_Struc'][3]=self.P['V_Vap'];
                    self.Bra[i]['Para_Struc'][4]=self.P['V_TP'];
                    self.Bra[i]['Para_Struc'][5]=self.P['V_Liq'];
                    self.Bra[i]['Para_Struc'][6]=self.P['m_Vap'];
                    self.Bra[i]['Para_Struc'][7]=self.P['m_TP'];
                    self.Bra[i]['Para_Struc'][8]=self.P['m_Liq'];    
                    self.Bra[i]['Para_Struc'][9]=self.P['UA_Vap'];
                    self.Bra[i]['Para_Struc'][10]=self.P['UA_TP'];
                    self.Bra[i]['Para_Struc'][11]=self.P['UA_Liq'];
                
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
    
        self.P=Bak;
        
        return 0
