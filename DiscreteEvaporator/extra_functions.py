from __future__ import division, print_function, absolute_import
from math import pi,log,sqrt,exp,cos,sin,tan,log10

from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.CoolProp import PropsSI

def EvapNode():
    evap_node={
        'NodNo':int(0),#node number
        'InNum':int(0),#number of tubes flowing into the node
        'OutNum':int(0),#number of the tubes flowing out of the node
        'BranIN':[],#index of the tube branch flowing in
        'BranOUT':[],#index of the tube branch flowing out
        }
    
    return evap_node

def EvapBranch():
    evap_branch = {
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
    
    return evap_branch

def TubeEvap():
    tube_evap = {
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
    
    return tube_evap

def TubEvpSeg():
    
    tub_evp_seg = {'TPi':{'T':0.0,'P':0.0}, 'WHo':{'W':0.0,'H':0.0}}
    
    return tub_evp_seg
    
        
def FlowPattern():
    
    flowpattern = {'JudgPattern':0,'G_wavy':0.0,'G_strat':0.0,'G_mist':0.0,
                   'G_bub':0.0,'X_lA':0.0,'epsilon':0.0,'theta_dry':0.0,
                   'delta':0.0,'h_nb':0.0,'h_cb':0.0,'h_v':0.0,'h_wet':0.0,
                   'h_tp':0.0,'Pattern':0} #output struct variable from Kattan-Thome flow-pattern-dependent heat transfer model
    
    return flowpattern

def PreAcc():
    
    pre_acc = {'DP_FR':0.0,'P_IN':0.0,'H_OUT':0.0,'G':0.0,'X_IN':0.0}

    return pre_acc

def InsideTube_Dim():
    '''
    returns dictionary for storing the parameters for micro-fin tube
    '''
    inside_tube_dim = {'Microfin':0, #microfin type, 0=smooth tube, 1=helical, 2=cross-grooved, 3=herringbone
                      'finN':0.0,   #fin number in a micro-fin tube
                      'Di':0.0,     #inside diameter at the fin tip
                      'gama':0.0,   #fin apex angle in a micro-fin tube
                      'beta':0.0,   #fin helix angle in a micro-fin tube
                      'finH':0.0,   #fin height in a micro-fin tube
                      'w_b':0.0,    #base width of a single fin
                      'w_e':0.0,    #top width of a single fin
                      'w_z':0.0,    #base distance between two neighboring fins
                      'K_T':0.0,    #400, this is the conductance factor of copper
                      'Ls':0.0,     #tube unit length
                      'D_b':0.0,    #tube diameter at the base of the fin
                      'Do':0.0,     #Pipe outside diameter.
                      'D_m':0.0,    #mean diameter of the micro-fin tube
                      'P_H':0.0,    #the hydraulical circumference
                      'Acs':0.0,    #cross area of the micro-fin tube, this is the actual cross-section area
                      'Dh_i':0.0,   #Inside hydraulical diameter 
                      'Ax':0.0,     #Inside pipe cross sectional area, based on fin tips
                      'Api':0.0}    #Inside pipe surface area (one tube segment), based on fin tips
    
    return inside_tube_dim

def Airside_Dim():
    '''
    returns dictionary for storing the parameters for airside fins
    '''
    airside_dim = {'evap':0,    #1 represent evaporator
                   'T':0.0,     #air temperature
                   'P':0.0,     #air humidity ratio
                   'Do':0.0,    #tube outside diameter
                   'Df':0.0,    #fin outside diameter
                   'z':0.0,     #space between fins (m)
                   'th':0.0,    #fin thickness (m)
                   'vsp':0.0,   #vertical tube spacing (m)
                   'P_l':0.0,   #B.S. tube spacing along the airflow direction (m)
                   'Apo':0.0,   #nominal air mass flux (kg/s/m^2)
                   'Ndeep':0.0, #B.S. number of high rows
                   'airFin':0,  #fin type
                   'K_F':0.0,   #conductance factor of the fin material         
                   'y':0.0,     #Distance from outside of pipe to edge of fin.
                   'N':0.0,     #Number of fins along one tube segment.
                   'Af':0.0,    #Fin wetted area (one tube segment).
                   'Aflow':0.0, #Air flow area (one tube segment).
                   'Dh':0.0,    #Hydrolic diameter
                   'sub1':0.0,'sub2':0.0,'sub3':0.0,'sub4':0.0,'sub5':0.0,  #possible sub-structures for different fin surface
                   'ho':0.0,'wetadj':0.0,   #for output the calculation result
                   'Ls':0.0}
    
    return airside_dim

def ETdim():
    '''
    This function return an initialized dictionary (with zeros) for all evaporator structure
    '''
    evap_struc={'Di':0.0,'L':0.0,'xp':0.0,'Df':0.0,'z':0.0,'th':0.0,'vsp':0.0,'Brad':0.0,
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

    return evap_struc

def EVA_Get_Q_dic():
    '''
    This function return an initialized dictionary (with zeros)
    #add for evaporative heat transfer iteration, since most of the evporating heat transfer calculation contains the nucleate boiling component, 
    #which needs to know the tube wall temperature at first
    '''
    eva_get_q = {'ma':0.0,  #air mass flow rate
                 'mr':0.0, #refrigerant mass flow rate
                 'TXPo':{'T':0.0,'X':0.0,'P':0.0},   #refrigerant state in the segment
                 'TPi':{'T':0.0,'P':0.0},     #air inlet state of the segment
                 'Gr':0.0,  #mass flow flux
                 'q':0.0,   #heat flux
                 'T_w':0.0, #inside tube wall temperature
                 'P':ETdim(), #evaporator struct
                 'W':0.0,   #outlet air humidity 
                 'Cmin':0.0,#Cmine
                 'hi':0.0,  #inside refrigerant heat transfer coefficent
                 'HPo':{'T':0.0,'P':0.0},     #inlet state of the refrigerant
                 'T_S_O':0.0}#effective tube surface temperature
    
    
    return eva_get_q

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

def HPtoTP(H,P):
    '''
    /********************************************************************
    Converts format of thermodynamic state representation for moist
    air from enthalpy (H) and relative humidity (P) to the
    standard format for ACMODEL: temperature (T) and relative
    humidity (P).
    ********************************************************************/
    '''
    TP = {'T':0.0,'P':0.0}
    HP = {'H':0.0,'P':0.0}
    
    HP['H'] = H;
    HP['P'] = P;
    
    TP['T'] = HAPropsSI('T','P',101325,'H',H,'R',P) #[K]
    TP['P'] = P;

    return TP

def WPtoTP(W,P):
    '''
    /********************************************************************
    Converts format of thermodynamic state representation for moist
    air from humidity ratio (W) and relative humidity (P) to the
    standard format for ACMODEL: temperature (T) and relative
    humidity (P).
    ********************************************************************/
    '''

    TP = {'T':0.0,'P':0.0}
    WP = {'W':0.0,'P':0.0}
    
    WP['W'] = W; 
    WP['P'] = P;
    
    TP['T'] = HAPropsSI('T','P',101325,'W',W,'R',P) #[K]
    TP['P'] = P;

    return TP

def THtoTP(T,H):
    '''
    /********************************************************************
    Converts format of thermodynamic state representation for moist
    air from temperature (T) and enthalpy (H) to the standard format
    for ACMODEL: temperature (T) and relative humidity (P).
    ********************************************************************/
    '''

    TP = {'T':0.0,'P':0.0}
    TH = {'T':0.0,'H':0.0};

    TH['T'] = T;
    TH['H'] = H;
    
    TP['P'] = HAPropsSI('R','P',101325,'T',T,'Hha',H) #[-]
    TP['T'] = T;

    return TP

def PropertyTXPth(prop,TXP,Ref):
    '''
    /********************************************************************
    Evaluates properties given TXP data.
    Refrigerant properties
    thermodynamic properties (th) - not transport properties
    pressure and temperature (PT) are independent variables

    prop = sting of 
        'T' Temperature [K], 0
        'H' Enthalpy [J/kg], 1
        (No volume 2 >>> it should be 1/Density)
        'S' Entropy [J/kg/K], 3
        'U' Internal Energy [J/kg], 4
        'D' Density [kg/m^3], 5
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
#         if (prop == 'T'): 
#             al = PropsSI(prop,'P',TXP['P'],'Q',0,Ref)
#             av = PropsSI(prop,'P',TXP['P'],'Q',1,Ref)
#             a = al + TXP['X']*(av-al);
#         else:
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
        'P' *pressure [Pa], 0
        'T' *Temperature [K], 1
        'L' conductivity [W/m/K], 2
        'V' viscosity [Pa-s], 3
        'C' specific heat [J/kg/K], 4
        'I' surface tension [N/m], 5
    ********************************************************************/
    '''

    if (TXP['X']>=1): 
        Tsat = PropsSI('T','P',TXP['P'],'Q',1,Ref)

        if (TXP['T']>Tsat): #superheat vapor
            a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
        else: #saturated vapor
            a = PropsSI(prop,'P',TXP['P'],'Q',1,Ref)
    
    elif (TXP['X']<=0):
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
