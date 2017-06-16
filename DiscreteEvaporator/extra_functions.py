from __future__ import division, print_function, absolute_import
from math import pi,log,sqrt,exp,cos,sin,tan,log10

from scipy.optimize import brentq #solver to find roots (zero points) of functions

from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.CoolProp import PropsSI

def ACRP():
    acrp = {'HPi':{'H':0.0,'P':0.0},
            'Go':0.0,'Gi':0.0,'Di':0.0,'vi':0.0,'si':0.0}
    return acrp

def CLRP():
    clrp = {'HPi':{'H':0.0,'P':0.0},
            'q':0.0,'K':0.0,'G':0.0,'D':0.0,'vi':0.0}
    return clrp

def VolParams():
    
    vol_params = {'vv':0,'vl':0,'ul':0,'uv':0,'Beta':0,'G':0,'D':0,'x':0}
    
    return vol_params

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
        'Para_Struc':[0.0 for k in range(12)]#output parameters of this condenser branch
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
    
def CondNode():
    cond_node={
        'NodNo':int(0),#node number
        'InNum':int(0),#number of tubes flowing into the node
        'OutNum':int(0),#number of the tubes flowing out of the node
        'BranIN':[],#index of the tube branch flowing in
        'BranOUT':[],#index of the tube branch flowing out
        }
    
    return cond_node

def CondBranch():
    cond_branch = {
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
        'Para_Struc':[0.0 for k in range(12)]#output parameters of this condenser branch
        }
    
    return cond_branch


def TubeCond():
    tube_cond = {
        'TubNo':int(0),             #tube No
        'RowNo':int(0),             #the row no where the tube located
        'RefUpstream':int(0),       #Upstream tube No at the refrigerant side (flow direction)
        'AirUpstreamUpper':int(0),  #Upstream tube No at the air upper side
        'AirUpstreamLower':int(0),  #Upstream tube No at the air lower side
        'GaFac':0.0,                #air flow distribution factor
        'Ga':0.0,                   #maximum air flow rate in the segment
        'even':int(0),              #tube flow direction
        'Seg':[],
        'HPi':{'H':0.0,'P':0.0},    #refrigerant inlet state
        'HPo':{'H':0.0,'P':0.0},    #refrigerant outlet state
        'm':{'m':0.0,'V':0.0},      #mass and volume of the tube
        }
    
    return tube_cond

def TubCndSeg():
    
    tub_cnd_seg = {'Tai':{'T':0.0,'P':0.0}, 'hao':{'W':0.0,'H':0.0}}
    
    return tub_cnd_seg

def FlowPattern():
    
    flowpattern = {'JudgPattern':0,'G_wavy':0.0,'G_strat':0.0,'G_mist':0.0,
                   'G_bub':0.0,'X_lA':0.0,'epsilon':0.0,'theta_dry':0.0,
                   'delta':0.0,'h_nb':0.0,'h_cb':0.0,'h_v':0.0,'h_wet':0.0,
                   'h_tp':0.0,'Pattern':0} #output struct variable from Kattan-Thome flow-pattern-dependent heat transfer model
    
    return flowpattern

def PreAcc():
    
    pre_acc = {'DP_FR':0.0,'P_IN':0.0,'H_OUT':0.0,'G':0.0,'X_IN':0.0,'rho_IN':0.0}

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
                'TPo':{'T':0.0,'P':0.0},    #A.B.
                'Eva_AirDP':0.0,    #A.B.
                'hAirAdj':0.0,'hRefAdj':0.0,'PRefAdj':0.0,'WAirAdj':0.0,
                'type':int(0),
                'Nrows':int(0),'Ndeep':int(0),
                'NBranchs':int(0), 'NBraTube':int(0),
                'P_l':0.0, #spacing between tubs in the longitudual direction (m)
                'Microfin':int(0), #new parameters for micro-fin tubes, and specially configured fins #decide whether is micro-fin tubes microfin=1
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
                'Hout8':[0.0 for k in range(10)],'DPr':[0.0 for k in range(10)],#superheat and pressure drop at each evaporator branch (array of 10 elements)
                }

    return evap_struc

def CGP():
    '''
    This function return an initialized dictionary (with zeros) for all evaporator structure
    '''
    cond_struc={'Di':0.0,'L':0.0,'xp':0.0,'Df':0.0,'z':0.0,'th':0.0,'vsp':0.0,'Brad':0.0,
                'NSeg':int(0.0),
                'Dh':0.0,'Do':0.0,'y':0.0,'Ls':0.0,'N':0.0,'Ax':0.0,'Api':0.0,'Apo':0.0,'Af':0.0,'Aflow':0.0,'Vs':0.0,'Blen':0.0,
                'Ro':0.0,'R':0.0, #A.B. overall heat resistance
                'hAirAdj':0.0,'hRefAdj':0.0,'PRefAdj':0.0,
                'hRefAdj_Sub':0.0,#B.S. parameter for adjusting the heat transfer ratio between the two-phase heat transfer and single-phase heat transfer
                'Nbranchs':int(0.0),'Nmaintubes':int(0.0),'Nsubtubes':int(0.0),
                'Ndeep':int(0.0),#B.S. accounting for the row number in the airflow direction
                'cfmA':0.0,'cfmB':0.0,'cfmC':0.0,
                'Microfin':int(0.0),# B.S. for judging if it is micro-fin tube
                'finN':0.0, 'gama':0.0, 'beta':0.0, 'finH':0.0,#B.S. micro-fin geometry
                'w_b':0.0, 'w_e':0.0, 'w_z':0.0, #B.S. micro-fin geometry
                'D_b':0.0, 'D_m':0.0, #B.S. micro-fin geometry
                'Acs':0.0,'P_H':0.0, 'Dh_i':0.0,#B.S. micro-fin geometry
                'airFin':int(0.0), #B.S. airside fin type
                'sub1':0.0,'sub2':0.0,'sub3':0.0,'sub4':0.0,'sub5':0.0,#for inputing sub-structures of fin surface
                'P_l':0.0, #airside fin geometry
                'K_T':0.0, 'K_F':0.0,# conductance factor
                'Ga':0.0,#airside velocity
                'airT':0.0,#air temperature
                'Frontal_A':0.0, #frontal area
                'GetP':0.0,#calculate the airside pressure drop
                
                #variables for generating the simple condenser function
                'V_TOT':0.0, 'V_TP':0.0, 'V_Liq':0.0, 'V_Vap':0.0,#B.S., inner volume of different phase
                'L_TOT':0.0,'LiqL':0.0,'VapL':0.0,'TPL':0.0,#B.S., tube length of different phase
                'A_TOT':0.0,'A_Liq':0.0,'A_Vap':0.0,'A_TP':0.0,#heat tranfer surface area of different phase
        
                'm_TOT':0.0, 'm_TP':0.0, 'm_Liq':0.0, 'm_Vap':0.0,#air mass flow rate across different phase region
                'rho_TOT':0.0, 'rho_TP':0.0, 'rho_Liq':0.0, 'rho_Vap':0.0,#average density of different phases
        
                'U_TP':0.0, 'U_Liq':0.0, 'U_Vap':0.0,#average heat transfer conductance per tube length
                'DP_TOT':0.0, 'DP_TP':0.0, 'DP_Liq':0.0, 'DP_Vap':0.0,#average pressure drop gradient of different phase
                'UA_TOT':0.0, 'UA_Liq':0.0,'UA_Vap':0.0,'UA_TP':0.0,#heat transfer conductance of diffferent phase
        
                'Tai':{'T':0.0,'P':0.0},#inlet air state A.B.
                'Tao':{'T':0.0,'P':0.0}, #outlet air state A.B.
                'Cond_AirDP':0.0,    #A.B.
                'mr':0.0, 'ma_TOT':0.0,#refrigerant and air mass flow rate
        
                'HP_in':{'H':0.0,'P':0.0}, 'HP_TP1':{'H':0.0,'P':0.0}, 'HP_TP2':{'H':0.0,'P':0.0}, 'HP_out':{'H':0.0,'P':0.0},#state parameters at important locations
                'Ga_meanL':0.0,#averager mass flow rate per tube length
    
                'r_v':0.0, 'r_l':0.0,'r_tp':0.0,#parameters for adjusting the theoritical heat transfer effectivenss of different phase
                'count1':0.0,'count2':0.0,#count the state points of two-phase flow end and beginning
                'H1_residual':0.0, 'H2_residual':0.0,#for consistency of the moving boundary model analysis
                #------------------------------B.S.
                'T_w':0.0, #A.B. wall temperature
                'q_w':0.0, #A.B. heat flux
                'fi':0.0, #A.B. frcition factor
                'rho':0.0, #A.B. density to be passed for friction pressure drop
                #Lumped model included above---------------------------B.S. 
                'cfma':0.0} # cfm air

    return cond_struc

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
    if (Ref=='R744'):
        TXP['X']=1;
        TXP['T'] = PropsSI('T','P',HP['P'],'H',HP['H'],Ref)
        return TXP
    
    hl = PropsSI('H','P',HP['P'],'Q',0,Ref)
    hv = PropsSI('H','P',HP['P'],'Q',1,Ref)
    Tl = PropsSI('T','P',HP['P'],'Q',0,Ref)
    Tv = PropsSI('T','P',HP['P'],'Q',1,Ref)
    
    if (HP['H']>hv):#superheated
        TXP['X']=1;
        TXP['T'] = PropsSI('T','P',HP['P'],'H',HP['H'],Ref)
    elif (HP['H']<hl):#subcooled
        TXP['X']=0;
        TXP['T'] = PropsSI('T','P',HP['P'],'H',HP['H'],Ref) 
    else: #two-phase
        try:
            TXP['X'] = PropsSI('Q','P',HP['P'],'H',HP['H'],Ref)
            TXP['T'] = PropsSI('T','P',HP['P'],'H',HP['H'],Ref)
        except:
            TXP['X'] = (HP['H']-hl)/(hv-hl)
            TXP['T'] = Tl + TXP['X']*(Tv-Tl)
    
    return TXP

def TXPtoHP(TXP,Ref):
    '''Convert TXP format to HP format.'''
    
    HP={'H':0.0,'P':0.0};

    HP['P'] = TXP['P']; #[Pa]
    HP['H'] = PropertyTXPth('H',TXP,Ref); #[J/kg]

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
    def HPFunc(T,Params):
        '''
        /********************************************************************
        Used by {\it HPtoTP} to convert thermodynamic state representations.
        ********************************************************************/
        '''
        P = Params;
        H = HAPropsSI('H','P',101325,'T',T,'R',P['P'])
        Ref=P['H']; #B.S.
        if(abs(Ref)<1e4):
            Ref=1e4; #B.S.
        return (H-P['H'])/Ref #B.S.


    TP = {'T':0.0,'P':0.0}
    HP = {'H':0.0,'P':0.0}
    HP['H']=H;
    HP['P']=P
    TMIN = -20+273.15 #lower bound [K]
    TMAX = 50+273.15 #upper bound [K]
    
    TP['T'] = brentq(HPFunc,TMIN,TMAX,args=(HP),xtol=1e-7,rtol=6e-8,maxiter=40)
    #TP['T'] = HAPropsSI('T','P',101325,'H',H,'R',P) #[K]
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
    def WPFunc(T,Params):
        '''
        /********************************************************************
        Used by {WPtoTP} to convert thermodynamic state representations.
        ********************************************************************/
        '''
        P = Params
        W = HAPropsSI('W','P',101325,'T',T,'R',P['P'])
        return (W-P['W'])/P['W'];

    TP = {'T':0.0,'P':0.0}
    WP = {'W':0.0,'P':0.0}
    WP['W']=W; 
    WP['P']=P
    TMIN = -20+273.15 #lower bound [K]
    TMAX = 50+273.15 #upper bound [K]
    
    TP['T'] = brentq(WPFunc,TMIN,TMAX,args=(WP),xtol=1e-7,rtol=6e-8,maxiter=40)
    #TP['T'] = HAPropsSI('T','P',101325,'W',W,'R',P) #[K]
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
    def THFunc(P,Params):
        '''
        /********************************************************************
        Used by {\it THtoTP} to convert thermodynamic state representations.
        ********************************************************************/
        '''
        Q = Params;
        H = HAPropsSI('H','P',101325,'T',Q['T'],'R',P)
        Ref=Q['H']; #B.S.
        if(abs(Ref)<1e4):
            Ref=1e4;#B.S.
        return (H-Q['H'])/Ref #B.S.

    TP = {'T':0.0,'P':0.0}
    TH = {'T':0.0,'H':0.0}
    TH['T'] = T; 
    TH['H'] = H 
    
    TP['P'] = brentq(THFunc,0.05,1.0,args=(TH),xtol=1e-7,rtol=6e-8,maxiter=40)
    #TP['P'] = HAPropsSI('R','P',101325,'T',T,'H',H) #[-]
    TP['T'] = T;

    return TP


def WHtoTP(WH,TPi):
    '''
    /********************************************************************
    Converts format of thermodynamic state representation for moist
    air from humidity ratio (W) and enthalpy (H) to the standard format
    for ACMODEL: temperature (T) and relative humidity (P).
    ********************************************************************/
    '''
#     def WHFunc(X,Params):
#         '''
#         /********************************************************************
#         Used by {WHtoTP} to convert thermodynamic state representations.
#         ********************************************************************/
#         '''
#         P = Params;
#         W = HAPropsSI('W','P',101325,'T',X[1],'R',X[0])
#         H = HAPropsSI('H','P',101325,'T',X[1],'R',X[0])
#         F = np.zeros(2)
#         F[0]=(W-P['W'])/P['W'];
#         Ref=P['H']; #B.S.
#         if(abs(Ref)<1e4):
#             Ref=1e4; #B.S.
#         F[1]=(H-P['H'])/Ref;#B.S.
#         return np.dot(F,F)

#     def WHFuncConst(X):
#         '''
#         /********************************************************************
#         Tests constraint on relative humidity when solving WHFunc.
#         Returns 1 if constraint is violated.
#         ********************************************************************/
#         '''
#         TWAIRMAX = 50+273.15 #upper bound [K]
#         TWAIRMIN = -20+273.15 #lower bound [K]
#         F = np.zeros(4)
#         F[0] = X[0]-1
#         F[1] = 0.05-X[0]
#         F[2] = X[1]-TWAIRMAX
#         F[3] = TWAIRMIN-X[1]
#         return F
    
    
    TPo = {'T':0.0,'P':0.0}
    WHo = {'W':WH['W'],'H':WH['H']}
#     TWAIRMAX = 50+273.15 #upper bound [K]
#     TWAIRMIN = -20+273.15 #lower bound [K]
    
    #initial guess
    X = (TPi['P'], TPi['T'])
    try:
#         bnds = ((0.05, 1.0), (-20+273.15, 50+273.15))
#         cons = {'type':'eq', 'fun':WHFuncConst}
#         cons = ({'type': 'eq', 'fun': lambda x: x[0]-1},
#                 {'type': 'eq', 'fun': lambda x: 0.05-x[0]},
#                 {'type': 'eq', 'fun': lambda x: x[1]-TWAIRMAX},
#                 {'type': 'eq', 'fun': lambda x: TWAIRMIN-x[1]})
#         res = minimize(WHFunc, X, args=(WHo), method='SLSQP', bounds=bnds, constraints=(), options={'eps':1e-4, 'maxiter': 100, 'ftol':1e-6})
#         TPo['P']=res.x[0]
#         TPo['T']=res.x[1]
        TPo ={'P': HAPropsSI('R','P',101325,'W',WHo['W'],'H',WHo['H']), 'T':HAPropsSI('T','P',101325,'W',WHo['W'],'H',WHo['H'])}
    except:
        #print ("Error in WHtoTP() conversion TPo= {}".format(TPo['T'],TPo['P']))
        if(TPo['P']>0.995):
            TPo = HPtoTP(WH['H'],0.995)
        else:
            print ("Error in WHtoTP() conversion TPo= {}".format(TPo['T'],TPo['P']))
            raise
        
    return TPo

def HatoTa(ha):
    '''
    /********************************************************************
    Converts the enthalpy of dry air to temperature of dry air.
    ********************************************************************/
    '''
    def HaTaFunc(T,Params):
        '''
        /********************************************************************
        Function called by HatoTa().  It returns a residual that is zero
        when the temperature corresponding to known enthalpy is selected.
        ********************************************************************/
        '''
        ha = Params;
        h = HAPropsSI('H','P',101325,'T',T,'R',0)
        dh = (ha['H']-h)/(ha['H']);
        return dh;
    
    Ta = {'T':0.0,'P':0.0}
    
#     TAIRMAX = 80+273.15
#     TAIRMIN = -20+273.15
#     dh = HaTaFunc(TAIRMAX,ha);
#     if(dh>0):
#         return 1e20;
# 
#     Ta['T'] = brentq(HaTaFunc,TAIRMIN,TAIRMAX,args=(ha),xtol=1e-8,rtol=6e-8,maxiter=40)
    Ta['T'] = HAPropsSI('T','P',101325,'H',ha['H'],'R',0)
    
    return Ta


def PropertyTXPth(prop,TXP,Ref):
    '''
    /********************************************************************
    Evaluates properties given TXP data.
    Refrigerant properties
    thermodynamic properties (th) - not transport properties
    pressure and temperature (PT) are independent variables

    prop = sting of 
        'T' Temperature [K],
        'H' Enthalpy [J/kg],
        'S' Entropy [J/kg/K],
        'U' Internal Energy [J/kg],
        'D' Density [kg/m^3],
    ********************************************************************/
    '''

    if (TXP['X']>=1): #Superheated
        try:
            Tsat = PropsSI('T','P',TXP['P'],'Q',1,Ref)
        except:
            if (Ref =='R744'):
                a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
                return a
            else:
                raise
        
        if (prop=='T'): #0 is TSAT             
            return Tsat

        if (TXP['T']>Tsat): #superheat vapor
            a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
        else: #saturated vapor
            try:
                a = PropsSI(prop,'P',TXP['P'],'Q',1,Ref)
            except:
                if (Ref =='R744'):
                    a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
                    return a
                else:
                    raise
            
            
    elif (TXP['X']<=0): #Subcooled
        try:
            Tsat = PropsSI('T','P',TXP['P'],'Q',0,Ref)
        except:
            if (Ref =='R744'):
                a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
                return a
            else:
                raise
            
        if(prop=='T'): #0 is TSAT
            return Tsat

        if (TXP['T']<Tsat): #subcooled liquid
            a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
        else: #saturated liquid
            try:
                a = PropsSI(prop,'P',TXP['P'],'Q',0,Ref)
            except:
                if (Ref =='R744'):
                    a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
                    return a
                else:
                    raise
    
    else: #two-phase
        try:
            a = PropsSI(prop,'P',TXP['P'],'Q',TXP['X'],Ref)
        except:
            al = PropsSI(prop,'P',TXP['P'],'Q',0,Ref)
            av = PropsSI(prop,'P',TXP['P'],'Q',1,Ref)
            a = al + TXP['X']*(av-al);

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
        try:
            Tsat = PropsSI('T','P',TXP['P'],'Q',1,Ref)
        except:
            if (Ref=='R744'):
                a = PropsSI(prop,'P',TXP['P'],'T',TXP['T'],Ref)
                return a
            else:
                raise
                
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
        if (prop=='I'):
            try:
                a = PropsSI(prop,'P',TXP['P'],'Q',TXP['X'],Ref)
            except:
                al = PropsSI(prop,'P',TXP['P'],'Q',0,Ref)
                av = PropsSI(prop,'P',TXP['P'],'Q',1,Ref)
                a = al + TXP['X']*(av-al);
        else: #B.S only the surface tension can be checked in two-phase state
            print('PropertyTXPtr :: only the surface tension can be checked in two-phase state') 
            raise  
    return a

if __name__=='__main__':
    print('Hello world')
