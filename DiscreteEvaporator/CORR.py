from __future__ import division, print_function, absolute_import
from math import log,pi,sqrt,exp,cos,sin,tan,log10,tanh
#from scipy.integrate import quad,quadrature,trapz,simps,fixed_quad
from scipy.optimize import brentq

from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI

from extra_functions import PropertyTXPth, PropertyTXPtr, Airside_Dim, InsideTube_Dim, FlowPattern, toTXP
from VOLUME import Xtt
#from MYBESSEL import bessk0, bessk1


def ConvCoeffAir_EVA(TPi,#air state
                     G,#air flux
                     P):#evaporator struct
    '''
    #B.S., ------------------------------------------------
    #interface function for calling the airside heat transfer and friction factor for evaporator
    '''
    
    D = Airside_Dim(); #Dictionary for storing the parameters for airside fins
    D['evap'] =1;#evaporating conditions
    D['T'] = TPi['T'];#air temperature
    D['P'] = TPi['P'];#air humidity ratio
    D['Do'] = P['Do'];#tube outside diameter
    D['Df'] = P['Df'];#fin outside diameter
    D['z'] = P['z'];#space between fins (m)
    D['th'] = P['th'];# fin thickness (m)
    D['vsp'] = P['vsp'];# vertical tube spacing (m)
    D['P_l'] = P['P_l'];# B.S. tube spacing along the airflow direction (m)
    D['Apo'] = P['Apo']; # tube outside area
    D['Ndeep'] = P['Ndeep'];# B.S. number of high rows
    D['airFin'] = P['airFin'];#fin type
    D['K_F'] = P['K_F'];#conductance factor of the fin material         
    D['y'] = P['y'];# Distance from outside of pipe to edge of fin.
    D['N'] = P['N']; # Number of fins along one tube segment.
    D['Af']= P['Af']; #Fin wetted area (one tube segment).
    D['Aflow'] = P['Aflow'];# Air flow area (one tube segment).
    D['Dh'] = P['Dh'];     # Hydrolic diameter
    D['sub1'] = P['sub1'];#possible sub-structures for different fin surface
    D['sub2'] = P['sub2'];#possible sub-structures for different fin surface
    D['sub3'] = P['sub3'];#possible sub-structures for different fin surface
    D['sub4'] = P['sub4'];#possible sub-structures for different fin surface
    D['sub5'] = P['sub5'];#possible sub-structures for different fin surface
    D['Ls'] = P['Ls'];#length of the segment
    
    if (P['wet']):
        D['wetadj']=0;
    else:
        D['wetadj']=1;
     
    if (P['GetP']<1):#calculate the heat transfer
        try:
            if D['airFin'] == 1:#plain 
                h, D = ConvCoeffAir_Plain(G, D);
            elif D['airFin'] == 2:#wavy plate
                h, D = ConvCoeffAir_Corrugated(G, D);
            elif D['airFin'] == 3:#slit
                h, D = ConvCoeffAir_Slit(G, D);
            elif D['airFin'] == 4:#louverd
                h, D = ConvCoeffAir_Louvered(G, D);
            elif D['airFin'] == 5:#convex louvered
                h, D = ConvCoeffAir_ConvexLouvered(G, D);#empty now
            elif D['airFin'] == 6:#smooth wavy
                h, D = ConvCoeffAir_SmoothWavy(G, D);
            elif D['airFin'] == 7:#spine
                h, D = ConvCoeffAir_Corrugated(G, D);#for temporary, since the correlations is not available
                #h, D = ConvCoeffAir_Spine(G, D);#empty now   
        except:
            h= 67.4
                     
        if (P['wet']):
            h=D['wetadj']*h;
            D['ho'] = h;#67.4;
        else:
            D['ho'] = h;#67.4;
     
        return h, P
 
    else: #calculate the airside pressure drop 
        try:
            if D['airFin'] == 1:#plain 
                f, D = FricAir_Plain(G, D);
            elif D['airFin'] == 2:#wavy plate
                f, D = FricAir_Corrugated(G, D);
            elif D['airFin'] == 3:#slit
                f, D = FricAir_Slit(G, D);
            elif D['airFin'] == 4:#louverd
                f, D = FricAir_Louvered(G, D);
            elif D['airFin'] == 5:#convex louvered
                f, D = FricAir_ConvexLouvered(G, D);#empty now
            elif D['airFin'] == 6:#smooth wavy
                f, D = FricAir_SmoothWavy(G, D);
            elif D['airFin'] == 7:#spine
                f, D = FricAir_Corrugated(G, D);#for temporary, since the correlations is not available
                #f, D = FricAir_Spine(G, D);#empty now   
        except:
            f= 1e-10
            
        return f, P

def ConvCoeffAir_CON(T,#air temperature
                     G,#air flux
                     P):#condenser struct
    '''
    #B.S., ------------------------------------------------
    #interface function for calling the airside heat transfer and friction factor for condenser
     
    '''
    
    C = Airside_Dim(); #dictionary for storing the parameters for airside fins
    C['evap'] =0;#evaporating conditions
    C['T'] = T;#air temperature
    C['P'] = 0.1;#air humidity ratio
    C['Do'] = P['Do'];#tube outside diameter
    C['Df'] = P['Df'];#fin outside diameter
    C['z'] = P['z'];#space between fins (m)
    C['th'] = P['th'];# fin thickness (m)
    C['vsp'] = P['vsp'];# vertical tube spacing (m)
    C['P_l'] = P['P_l'];# B.S. tube spacing along the airflow direction (m)
    C['Apo'] = P['Apo']; # tube outside area
    C['Ndeep'] = P['Ndeep'];# B.S. number of high rows
    C['airFin'] = P['airFin'];#fin type
    C['K_F'] = P['K_F'];#conductance factor of the fin material         
    C['y'] = P['y'];# Distance from outside of pipe to edge of fin.
    C['N'] = P['N']; # Number of fins along one tube segment.
    C['Af']= P['Af']; #Fin wetted area (one tube segment).
    C['Aflow'] = P['Aflow'];# Air flow area (one tube segment).
    C['Dh'] = P['Dh'];     # Hydrolic diameter
    C['sub1'] = P['sub1'];#possible sub-structures for different fin surface
    C['sub2'] = P['sub2'];#possible sub-structures for different fin surface
    C['sub3'] = P['sub3'];#possible sub-structures for different fin surface
    C['sub4'] = P['sub4'];#possible sub-structures for different fin surface
    C['sub5'] = P['sub5'];#possible sub-structures for different fin surface
    C['Ls'] = P['Ls'];#length of the segment
    
     
    if (P['GetP']<1):#calculate the heat transfer
        try:
            if C['airFin'] == 1:#plain 
                h, C = ConvCoeffAir_Plain(G, C);
            elif C['airFin'] == 2:#wavy plate
                h, C = ConvCoeffAir_Corrugated(G, C);
            elif C['airFin'] == 3:#slit
                h, C = ConvCoeffAir_Slit(G, C);
            elif C['airFin'] == 4:#louverd
                h, C = ConvCoeffAir_Louvered(G, C);
            elif C['airFin'] == 5:#convex louvered
                h, C = ConvCoeffAir_ConvexLouvered(G, C);#empty now
            elif C['airFin'] == 6:#smooth wavy
                h, C = ConvCoeffAir_SmoothWavy(G, C);
            elif C['airFin'] == 7:#spine
                h, C = ConvCoeffAir_Corrugated(G, C);#for temporary, since the correlations is not available
                #h, C = ConvCoeffAir_Spine(G, C);#empty now   
        except:
            h= 86
                     
        C['ho'] = h;#86;
        C['wetadj'] = 1;#for output the calculation result
     
        return h, P
  
    else: #calculate the airside pressure drop 
        try:
            if C['airFin'] == 1:#plain 
                f, C = FricAir_Plain(G, C);
            elif C['airFin'] == 2:#wavy plate
                f, C = FricAir_Corrugated(G, C);
            elif C['airFin'] == 3:#slit
                f, C = FricAir_Slit(G, C);
            elif C['airFin'] == 4:#louverd
                f, C = FricAir_Louvered(G, C);
            elif C['airFin'] == 5:#convex louvered
                f, C = FricAir_ConvexLouvered(G, C);#empty now
            elif C['airFin'] == 6:#smooth wavy
                f, C = FricAir_SmoothWavy(G, C);
            elif C['airFin'] == 7:#spine
                f, C = FricAir_Corrugated(G, C);#for temporary, since the correlations is not available
                #f, C = FricAir_Spine(G, C);#empty now   
        except:
            f= 1e-10
            
        return f, P


#===============================================================================
# the followings are the heat transfer calculations for different fin types
#===============================================================================

def ConvCoeffAir_Plain(G,#air flux 
                       P):#parameters of the airside fin
    '''
    /**************************************************
    Plain fin air heat transfer, B.S.
    <<Wang, C. C., et al, 2000. Heat Transfer And Friction Characteristics Of Plain Fin-and-Tube Heat Exchangers, 
    Part II: Correlation, International Journal Of Heat And Mass Transfer, Volume: 43, Issue: 15 August 1, 
    2000, pp. 2693-2700>>
    wet heat transfer adjustment from
    "empirical airside correlations of fin-and-tube heat exchangers under dehumidifying conditions"
    chi-chuan wang, wei-song lee, wen-jenn sheu and jane-sunn Liaw, 2001, INTERNATIONAL JOURNAL OF HEAT EXCHANGERS 
    ******************************************************/
    '''
 
    th = P['th'];#fin thickness
    D = P['Do']+2*th;#collar diameter
    D_h= P['Dh'];#hydraulic diameter
    F_p = P['z'];#fin pitch
    P_t=P['vsp'];#tube pitch perpendicular to the airflow
    P_l=P['P_l'];#tube pitch in the airflow direction
    N=P['Ndeep'];#tube number in the airflow direction
    F_s=F_p-th;#fin space
    T = P['T'];#air temperature
     
    # calc Prantl and Reynold number
    mu = HAPropsSI('mu','T',P['T'],'P',101325,'R',P['P']) #air.mu(T);#air viscosity [Pa-s]
    
    #air specific heat
    if (P['evap']): #during the conditions of evaporation
        cp = HAPropsSI('cp_ha','T',P['T'],'P',101325,'R',P['P']) #wair.Cp(P->T,P->P); #[J/kg humid air/K]
    else: #during the conditions of condensation
        cp = HAPropsSI('cp','T',P['T'],'P',101325,'R',P['P']) #air.Cp(T); #[J/kg dry air/K]
    
    k = HAPropsSI('K','T',P['T'],'P',101325,'R',P['P']) #air.k(T);#air heat conductance [W/m/K]
    
    Pr = mu*cp/k;#prandtl number
    if (Pr<0.0):
        print("ConvCoeffAir_Plain::ConvCoeffAir_Plain","Pr<0.0")
        
    Re = G*D/mu;#Reynold's number based on the collar diameter
    if (Re<0.0):
        print("ConvCoeffAir_Plain::ConvCoeffAir_Plain","Re<0.0")
    
    #jbraun number
    if (N<=1):#tube row less than 2
        P1 = 1.9-0.23*log(Re);
        P2 = -0.236+0.126*log(Re);
        j = 0.108*pow(Re,-0.29)*pow(P_t/P_l,P1)*pow(F_p/D,-1.084)*pow(F_p/D_h,-0.786)*pow(F_p/P_t,P2);
    else: #tube row larger than 1
        P3 = -0.361-0.042*N/(log(Re))+0.1581*log(N*pow(F_p/D,0.41));
        P4 = -1.224-0.076*pow(P_l/D_h,1.42)/log(Re);
        P5 = -0.083 + 0.058*N/log(Re);
        P6 = -5.735+1.211*log(Re/N);
        j = 0.086*pow(Re,P3)*pow(N,P4)*pow(F_p/D,P5)*pow(F_p/D_h,P6)*pow(F_p/P_t,-0.93);
    
    #airside heat transfer coefficient
    h = j*G*cp/pow(Pr,0.666667);
 
    if (not P['wetadj']):
        W4 = -0.213+0.0176*log(Re)+0.0172*(D/F_p);
        W5 = 0.551-2.63*(F_p/D)-0.012*N;
        P['wetadj']=0.968*pow(Re,W4)*pow((F_p/D),W5);
    
    return h, P

def ConvCoeffAir_Corrugated(G,P):
    '''
    /*****************************************************
    Corrugated fin, B.S.
    For tube diameters >= 1/2", 
    Wang, C. C., Y. T. Lin, C. J. Lee, and Y. J. Chang, "Investigation of 
    Wavy Fin-and-Tube Heat Exchangers: A Contribution to Databank", Experimental 
    Heat Transfer, 12:73-89(1999) 
    For tube diameters < 1/2", 
    Wang, C. C., J. Y. Jang, and N. F. Chiou, "A Heat Transfer and Friction Correlation 
    For Wavy Fin-and-tube Heat Exchangers", International Journal of Heat and Mass 
    Transfer, Vol 42(1999) pp.1919-1924. 
    wet heat transfer adjustment from
    "empirical airside correlations of fin-and-tube heat exchangers under dehumidifying conditions"
    chi-chuan wang, wei-song lee, wen-jenn sheu and jane-sunn Liaw, 2001, INTERNATIONAL JOURNAL OF HEAT EXCHANGERS
    *********************************************************/
    '''
     
    th = P['th'];
    D = P['Do']+2*th;
    D_h= P['Dh'];
    F_p = P['z'];
    P_t=P['vsp'];
    P_l=P['P_l'];
    N=P['Ndeep'];
    T = P['T'];
    F_s=F_p-th;
    X_f=P_l/4;#projected fin pattern length for one-half wavy length, typical structure from the paper
    P_d=1.32*0.001;#waffle height, from the paper
    
    # sec of the corrugation angle
    sec_angle=pow((pow(X_f,2.0)+pow(P_d,2.0)),0.5)/X_f;# calculate the sec of the corrugation angle
    # tan of the corrugation angle
    tg_angle=pow((pow(sec_angle,2.0)-1.0),0.5);#calculate the tan of the corrugation angle
 
    #beta=3.1415926*pow(D,2.0)/(4*P_l*P_t);#possible method for calculating hydraulic diameter for wavy fin, but not used here
    #D_h=2*F_p*(1-beta)/((1-beta)*sec_angle+2*F_p*beta/D);
    # calc Prantl and Reynold number
    
    # calc Prantl and Reynold number
    mu = HAPropsSI('mu','T',P['T'],'P',101325,'R',P['P']) #air.mu(T);#air viscosity [Pa-s]
    
    #air specific heat
    if (P['evap']): #during the conditions of evaporation
        cp = HAPropsSI('cp_ha','T',P['T'],'P',101325,'R',P['P']) #wair.Cp(P->T,P->P); #[J/kg humid air/K]
    else: #during the conditions of condensation
        cp = HAPropsSI('cp','T',P['T'],'P',101325,'R',P['P']) #air.Cp(T); #[J/kg dry air/K]
    
    k = HAPropsSI('K','T',P['T'],'P',101325,'R',P['P']) #air.k(T);#air heat conductance [W/m/K]
     
    Pr = mu*cp/k;#prandtl number
    if (Pr<0.0):
        print("ConvCoeffAir_Corrugated::ConvCoeffAir_Corrugated","Pr<0.0")
        
    Re = G*D/mu;#Reynold's number based on the collar diameter
    if (Re<0.0):
        print("ConvCoeffAir_Corrugated::ConvCoeffAir_Corrugated","Re<0.0")
    
    if (D<=25.4e-3/2):
        J_1=-0.229+0.115*pow((F_p/D),0.6)*pow((P_l/D_h),0.54)*pow(N,(-0.284))*log(0.5*tg_angle);
        J_2=-0.251+0.232*pow(N,1.37)/(log(Re)-2.303);
        J_3=-0.439*pow((F_p/D_h),0.09)*pow((P_l/P_t),(-1.75))*pow(N,(-0.93));
        J_4=0.502*(log(Re)-2.54);
        # the 'j-factor' correlation
        j =0.324*pow(Re,(J_1))*pow((F_p/P_l),(J_2))*pow((tg_angle),(J_3))*pow((P_l/P_t),J_4)*pow(N,0.428);
    else:
        J_5=-0.1707-1.374*pow(P_l/th,-0.493)*pow(F_p/D,-0.886)*pow(N,-0.143)*pow(P_d/X_f,-0.0296);
        # the 'j-factor' correlation
        j = 1.79097*pow(Re,J_5)*pow(P_l/th,-0.456)*pow(N,-0.27)*pow(F_p/D,-1.343)*pow(P_d/X_f,0.317);
 
    # calc conv coeff (h) from 'j-factor'(j)
    h = j*G*cp/pow(Pr,0.666667);
     
    if (not P['wetadj']):
        W5=0.0374-0.0018*log(Re)-0.00685*(D/F_p);
        P['wetadj']=0.794*pow(Re,W5)*pow((P_t/P_l),0.308)*pow((P_d/X_f),-0.119);
 
    return h, P

def ConvCoeffAir_Slit(G,P):
    '''
    /*****************************************************
    Slit fin, B.S.
    Wang, C. C., 2001, "A Comparative Study of Compact Enhanced Fin-and-Tube 
    Heat Exchangers", Int. J. Heat and Mass Transfer, Vol. 44, pp. 3565-3573. 
    wet heat transfer adjustment from
    "empirical airside correlations of fin-and-tube heat exchangers under dehumidifying conditions"
    chi-chuan wang, wei-song lee, wen-jenn sheu and jane-sunn Liaw, 2001, INTERNATIONAL JOURNAL OF HEAT EXCHANGERS
    *********************************************************/
    '''
    
    th = P['th'];
    D = P['Do']+2*th;
    D_h= P['Dh'];
    F_p = P['z'];
    P_t=P['vsp'];
    P_l=P['P_l'];
    N=P['Ndeep'];
    T = P['T'];
    F_s=F_p-th;
    S_h=0.99*0.001;#height of slit, typical structure from the paper
    S_s=2.2*0.001;#breadth of a slit in the direction of the airflow
    S_w=11*0.001;#width of slit
    S_n=4;#number of slits in an enhanced zone
    
    # calc Prantl and Reynold number
    mu = HAPropsSI('mu','T',P['T'],'P',101325,'R',P['P']) #air.mu(T);#air viscosity [Pa-s]
    
    #air specific heat
    if (P['evap']): #during the conditions of evaporation
        cp = HAPropsSI('cp_ha','T',P['T'],'P',101325,'R',P['P']) #wair.Cp(P->T,P->P); #[J/kg humid air/K]
    else: #during the conditions of condensation
        cp = HAPropsSI('cp','T',P['T'],'P',101325,'R',P['P']) #air.Cp(T); #[J/kg dry air/K]
 
    k = HAPropsSI('K','T',P['T'],'P',101325,'R',P['P']) #air.k(T);#air heat conductance [W/m/K]
     
    Pr = mu*cp/k;#prandtl number
    if (Pr<0.0):
        print("ConvCoeffAir_Slit::ConvCoeffAir_Slit","Pr<0.0")
        
    Re = G*D/mu;#Reynold's number based on the collar diameter
    if (Re<0.0):
        print("ConvCoeffAir_Slit::ConvCoeffAir_Slit","Re<0.0")
         
 
    J_4=-0.535+0.017*P_t/P_l-0.0107*N;
    J_5=0.4115+5.5756*pow(double(N)/(Re),0.5)*log(N/Re)+24.2028*pow(N/(Re),0.5);
    J_6=0.2646+1.0491*(S_s/S_h)*log(S_s/S_h)-0.216*pow((S_s/S_h),3.0);
    J_7=0.3749+0.0046*pow(Re,0.5)*log(Re)-0.0433*pow(Re,0.5);
 
    J_1=-0.255-0.0312/(F_s/D)-0.0487*N;
    J_2=0.9703 - 0.0455*pow(Re,0.5)-0.4986*pow(log(P_t/P_l),2.0);
    J_3=0.2405-0.003*Re+5.5349*(F_s/D);
     
    if (Re<700):
        j=0.9047*pow(Re,J_1)*pow(F_s/D,J_2)*pow(P_t/P_l,J_3)*pow((S_s/S_h),-0.0305)*pow(N,0.0782);
    else: # the below is the most important calculation
        j=1.0691*pow(Re,J_4)*pow(F_s/D,J_5)*pow(S_s/S_h,J_6)*pow(N,J_7);
 
    # calc conv coeff (h) from 'j-factor'(j)
    h = j*G*cp/pow(Pr,0.666667);
 
    if (not P['wetadj']):
        W5=-0.143+0.013*log(Re)+0.166*(F_p/D);
        P['wetadj']= 0.937*pow(Re,W5)*pow((S_s/S_h),0.1344)*pow(N,0.0657);
 
    return h, P

def ConvCoeffAir_Louvered(G,P):
    '''
    /*****************************************************
    Louvered fin, B.S.
    Wang, C. C., C. J. Lee, C. T. Chang, and S. P. Lin, 
    "Heat Transfer and Friction Correlation for Compact Louvered Fin-and-tube
     Heat Exchangers", International Journal of Heat and Mass Transfer, 
     Vol 42 (1999) pp.1945-1956. 
     wet heat transfer adjustment from
    "empirical airside correlations of fin-and-tube heat exchangers under dehumidifying conditions"
    chi-chuan wang, wei-song lee, wen-jenn sheu and jane-sunn Liaw, 2001, INTERNATIONAL JOURNAL OF HEAT EXCHANGERS
    *********************************************************/
    '''
 
    th = P['th'];
    D = P['Do']+2*th;
    D_h= P['Dh'];
    F_p = P['z'];
    P_t=P['vsp'];
    P_l=P['P_l'];
    N=P['Ndeep'];
    T = P['T'];
    F_s=F_p-th;
    L_h=1.07*0.001;#Louver height, typical structure from the paper
    L_p=2.4*0.001;#major louver pitch
 
    # calc Prantl and Reynold number
    mu = HAPropsSI('mu','T',P['T'],'P',101325,'R',P['P']) #air.mu(T);#air viscosity [Pa-s]
    
    #air specific heat
    if (P['evap']): #during the conditions of evaporation
        cp = HAPropsSI('cp_ha','T',P['T'],'P',101325,'R',P['P']) #wair.Cp(P->T,P->P); #[J/kg humid air/K]
    else: #during the conditions of condensation
        cp = HAPropsSI('cp','T',P['T'],'P',101325,'R',P['P']) #air.Cp(T); #[J/kg dry air/K]
 
    k = HAPropsSI('K','T',P['T'],'P',101325,'R',P['P']) #air.k(T);#air heat conductance [W/m/K]
     
    Pr = mu*cp/k;#prandtl number
    if (Pr<0.0):
        print("ConvCoeffAir_Louvered::ConvCoeffAir_Louvered","Pr<0.0")
        
    Re = G*D/mu;#Reynold's number based on the collar diameter
    if (Re<0.0):
        print("ConvCoeffAir_Louvered::ConvCoeffAir_Louvered","Re<0.0")
        
 
    j1 = -0.991-0.1055*pow(P_l/P_t,3.1)*log(L_h/L_p);
    j2 = -0.7344+2.1059*(pow(N,0.55)/(log(Re)-3.2));
    j3 = 0.08485*pow(P_l/P_t,-4.4)*pow(N,-0.68);
    j4 = -0.1741*log(N);
    j5 = -0.6027 + 0.02593*pow(P_l/D_h,0.52)*pow(N,-0.5)*log(L_h/L_p);
    j6 = -0.4776+0.40774*pow(N,0.7)/(log(Re)-4.4);
    j7 = -0.58655*pow(F_p/D_h,2.3)*pow(P_l/P_t,-1.6)*pow(N,-0.65);
    j8 = 0.0814*(log(Re)-3.0);
 
    if (Re<1000):
        j = 14.3117*pow(Re,j1)*pow(F_p/D,j2)*pow(L_h/L_p,j3)*pow(F_p/P_l,j4)*pow(P_l/P_t,-1.724);
    else:
        j = 1.1373*pow(Re,j5)*pow(F_p/P_l,j6)*pow(L_h/L_p,j7)*pow(P_l/P_t,j8)*pow(N,0.3545);
     
    # calc conv coeff (h) from 'j-factor'(j)
    h = j*G*cp/pow(Pr,0.666667);
 
    if (not P['wetadj']):
        W4=-0.0746+0.00115*pow(Re,0.5)+0.00028*pow(F_p/D,-2.0);
        W5=0.303-0.726*(L_h/L_p)+0.041*(L_p/F_p);
        P['wetadj'] = 0.263*pow(Re,W4)*pow((L_h/L_p),W5)*pow((F_p/D),-0.72)*pow((P_t/P_l),1.11)*pow((L_p/F_p),-0.742);
 
    return h, P

def ConvCoeffAir_ConvexLouvered(G,P):
    '''
    /*****************************************************
    Convexlouvered fin, B.S.
    Wang, C. C., Y. M. Tsai, and D. C. Lu, 
    "Comprehensive Study of Convex-louver and Wavy Fin-and-tube Heat Exchangers", 
    Journal of Thermophysics and Heat Transfer, Vol 12, No. 3, July-September 1998,
    pp.423-430. 
    wet heat transfer adjustment from
    "empirical airside correlations of fin-and-tube heat exchangers under dehumidifying conditions"
    chi-chuan wang, wei-song lee, wen-jenn sheu and jane-sunn Liaw, 2001, INTERNATIONAL JOURNAL OF HEAT EXCHANGERS
    *********************************************************/
    '''
 
    th = P['th'];
    D = P['Do']+2*th;
    D_h= P['Dh'];
    F_p = P['z'];
    P_t=P['vsp'];
    P_l=P['P_l'];
    N=P['Ndeep'];
    T = P['T'];
    F_s=F_p-th;
    
    # calc Prantl and Reynold number
    mu = HAPropsSI('mu','T',P['T'],'P',101325,'R',P['P']) #air.mu(T);#air viscosity [Pa-s]
    
    #air specific heat
    if (P['evap']): #during the conditions of evaporation
        cp = HAPropsSI('cp_ha','T',P['T'],'P',101325,'R',P['P']) #wair.Cp(P->T,P->P); #[J/kg humid air/K]
    else: #during the conditions of condensation
        cp = HAPropsSI('cp','T',P['T'],'P',101325,'R',P['P']) #air.Cp(T); #[J/kg dry air/K]
 
    k = HAPropsSI('K','T',P['T'],'P',101325,'R',P['P']) #air.k(T);#air heat conductance [W/m/K]
     
    Pr = mu*cp/k;#prandtl number
    if (Pr<0.0):
        print("ConvCoeffAir_ConvexLouvered::ConvCoeffAir_ConvexLouvered","Pr<0.0")
        
    Re = G*D/mu;#Reynold's number based on the collar diameter
    if (Re<0.0):
        print("ConvCoeffAir_ConvexLouvered::ConvCoeffAir_ConvexLouvered","Re<0.0")
        
        
    Coe_fin = ((P['Apo']+P['Af'])/(P['Ls']*P['Do']*pi));
    j1 = -1.02*(F_p/D)-0.256;
    j = 16.06*pow(Re,j1)*pow(Coe_fin,-0.601)*pow(N,-0.069)*pow(F_p/D,0.84);
     
    # calc conv coeff (h) from 'j-factor'(j)
    h = j*G*cp/pow(Pr,0.666667);
 
    if (not P['wetadj']):#borrow the factor of plain fin
        W4 = -0.213+0.0176*log(Re)+0.0172*(D/F_p);
        W5=0.551-2.63*(F_p/D)-0.012*N;
        P['wetadj']=0.968*pow(Re,W4)*pow((F_p/D),W5);
   
 
    return h, P

def ConvCoeffAir_SmoothWavy(G,P):
    '''
    /*****************************************************
    Smoothwavy fin, B.S.
    Mirth, D. R. and S. Ramadhyani, 1994,"Correlations for Predicting the Airside 
    Nusselt Numbers and Friction Factors in Chilled-water Cooling Coils", 
    Experimental Heat Transfer, Vol 7: 143-162.
    ---the author suggested predict both the dry heat transfer and wet heat transfer with the dry correlations
    wet heat transfer adjustment from
    "empirical airside correlations of fin-and-tube heat exchangers under dehumidifying conditions"
    chi-chuan wang, wei-song lee, wen-jenn sheu and jane-sunn Liaw, 2001, INTERNATIONAL JOURNAL OF HEAT EXCHANGERS
    *********************************************************/
    '''
 
    th = P['th'];
    D = P['Do']+2*th;
    D_h= P['Dh'];
    F_p = P['z'];
    P_t=P['vsp'];
    P_l=P['P_l'];
    N=P['Ndeep'];
    T = P['T'];
    F_s=F_p-th;
    
    # calc Prantl and Reynold number
    mu = HAPropsSI('mu','T',P['T'],'P',101325,'R',P['P']) #air.mu(T);#air viscosity [Pa-s]
    
    #air specific heat
    if (P['evap']): #during the conditions of evaporation
        cp = HAPropsSI('cp_ha','T',P['T'],'P',101325,'R',P['P']) #wair.Cp(P->T,P->P); #[J/kg humid air/K]
    else: #during the conditions of condensation
        cp = HAPropsSI('cp','T',P['T'],'P',101325,'R',P['P']) #air.Cp(T); #[J/kg dry air/K]
 
    k = HAPropsSI('K','T',P['T'],'P',101325,'R',P['P']) #air.k(T);#air heat conductance [W/m/K]
     
    Pr = mu*cp/k;#prandtl number
    if (Pr<0.0):
        print("ConvCoeffAir_SmoothWavy::ConvCoeffAir_SmoothWavy","Pr<0.0")
        
    Re = G*D/mu;#Reynold's number based on the collar diameter
    if (Re<0.0):
        print("ConvCoeffAir_SmoothWavy::ConvCoeffAir_SmoothWavy","Re<0.0")
        
     
    Nu = 0.0197*pow(Re, 0.94)*pow((P_t-D)/(2*F_s), -0.3)*(1+111.9/pow(Re*N*P_l/(2*F_s),1.2))*pow(Pr,1.0/3.0);
    h=Nu*k/(2*F_s);#it appears too small, around 34
 
    if (not P['wetadj']):#borrow the factor of plain fin
        W4 = -0.213+0.0176*log(Re)+0.0172*(D/F_p);
        W5=0.551-2.63*(F_p/D)-0.012*N;
        P['wetadj']=0.968*pow(Re,W4)*pow((F_p/D),W5);
      
    P['wetadj']=1.0;
     
    return h, P


def ConvCoeffAir_Spine(G,P):
    '''
    /*****************************************************
    Spine fin, empty now
    *********************************************************/
    '''
     
    h=0.0;
     
    return h, P


def Circuit_DP_EVAP(Ga,TPi,TPo,P):
    '''
    #B.S.--------------------------------------------------------
    ########the followings are the friction factors for different fin types
    '''

    P['GetP'] =1;#set the request for calculating the friction factor.
    
    rho1 = 1/HAPropsSI('Vda','T',TPi['T'],'P',101325,'R',TPi['P']) #air.v(TPi.T)#inlet air density
    
    rho2 = 1/HAPropsSI('Vda','T',TPo['T'],'P',101325,'R',TPo['P']) #air.v(TPo.T);#outlet air density
    
    rho_m=(rho1+rho2)/2;#mean air density
    
    f, P = ConvCoeffAir_EVA(TPi, Ga, P);#calculate the friction factor
    
    mu = HAPropsSI('mu','T',TPi['T'],'P',101325,'R',TPi['P']) #air.mu(TPi.T);#air viscosity
    
    Re = Ga*P['Do']/mu;#Reynolds number based on tube outside diameter
    if (Re<0.0):
        print("Circuit_DP_EVAP::Circuit_DP_E","Re<0.0")
        
    D1 = P['Aflow']/(P['Apo']+P['Af'])*rho1/rho_m;
    sigma = P['Ls']*P['vsp']/P['Aflow'];#frontal area contraction ratio
    D2 = (1+pow(sigma,2.0))*(rho1/rho2-1);
    
    if (P['airFin'] == 6):#the correlation for smooth wavy fin is different
        f_tubes = 2.17/pow(P['vsp']/P['Do'], 1.08) - 0.174*log(Re)/pow(P['vsp']/P['Do'],1.24);
        DP_tubes = f_tubes*P['Ndeep']*pow(Ga,2.0)/(2*rho1);
        DP_fins = f*rho1/rho_m*P['Af']/P['Aflow']*pow(Ga,2.0)/(2*rho1);
        DP_acc = (1+pow(sigma,2.0))*(rho1/rho2-1)*(pow(Ga,2.0)/(2*rho2));
        DP = DP_tubes + DP_fins + DP_acc;
    else:
        DP = (f/D1+D2)*pow(Ga,2.0)*rho1/2;#pressure drop calculations for other fin types
     
    return (DP, TPo, P)


def Circuit_DP_COND(Ga,Tai,Tao,P):
    
    P['GetP'] =1;#set the request for calculating the friction factor.
    
    rho1 = 1/HAPropsSI('Vda','T',Tai['T'],'P',101325,'R',Tai['P']) #air.v(Tai.T)#inlet air density
    
    rho2 = 1/HAPropsSI('Vda','T',Tao['T'],'P',101325,'R',Tao['P']) #air.v(Tao.T);#outlet air density
    
    rho_m=(rho1+rho2)/2;#mean air density
    
    f, P = ConvCoeffAir_CON(Tai, Ga, P);#calculate the friction factor
    
    mu = HAPropsSI('mu','T',Tai['T'],'P',101325,'R',Tai['P']) #air.mu(TPi.T);#air viscosity
    
    Re = Ga*P['Do']/mu;#Reynolds number based on tube outside diameter
    if (Re<0.0):
        print("Circuit_DP_COND::Circuit_DP_C","Re<0.0")
        
    D1 = P['Aflow']/(P['Apo']+P['Af'])*rho1/rho_m;
    sigma = P['Ls']*P['vsp']/P['Aflow'];#frontal area contraction ratio
    D2 = (1+pow(sigma,2.0))*(rho1/rho2-1);
    
    if (P['airFin'] == 6):#the correlation for smooth wavy fin is different
        f_tubes = 2.17/pow(P['vsp']/P['Do'], 1.08) - 0.174*log(Re)/pow(P['vsp']/P['Do'],1.24);
        DP_tubes = f_tubes*P['Ndeep']*pow(Ga,2.0)/(2*rho1);
        DP_fins = f*rho1/rho_m*P['Af']/P['Aflow']*pow(Ga,2.0)/(2*rho1);
        DP_acc = (1+pow(sigma,2.0))*(rho1/rho2-1)*(pow(Ga,2.0)/(2*rho2));
        DP = DP_tubes + DP_fins + DP_acc;
    else:
        DP = (f/D1+D2)*pow(Ga,2.0)*rho1/2;#pressure drop calculations for other fin types
        
    return (DP, Tao, P)

def FricAir_Plain(G,#air flux 
                  P):#parameters of the airside fin
    '''
    /**************************************************
    Plain fin air heat transfer, B.S.
    <<Wang, C. C., et al, 2000. Heat Transfer And Friction Characteristics Of Plain Fin-and-Tube Heat Exchangers, 
    Part II: Correlation, International Journal Of Heat And Mass Transfer, Volume: 43, Issue: 15 August 1, 
    2000, pp. 2693-2700>>
    ******************************************************/
    '''
    th = P['th'];#fin thichness
    D = P['Do']+2*th;#collar diameter
    D_h= P['Dh'];#hydraulic diameter
    F_p = P['z'];#fin pitch
    P_t=P['vsp'];#tube pitch perpendicular to the airflow
    P_l=P['P_l'];#tube pitch in the airflow direction
    N=P['Ndeep'];#tube number in the airflow direction
    F_s=F_p-th;#fin space
    T = P['T'];#air temperature
    
    mu = HAPropsSI('mu','T',P['T'],'P',101325,'R',P['P']) #air.mu(T);#air viscosity [Pa-s]
            
    Re = G*D/mu;#Reynold's number based on the collar diameter
    if (Re<0.0):
        print("FricAir_Plain::FricAir_Plain","Re<0.0")
        
     
    F1=-0.764+0.739*(P_t/P_l)+0.177*F_p/D-0.00758/N;
    F2=-15.689+64.021/log(Re);
    F3=1.696-15.695/log(Re);
    
    #airside heat transfer coefficient
    f = 0.0267*pow(Re,F1)*pow(P_t/P_l, F2)*pow(F_p/D,F3);
 
    return f, P


def FricAir_Corrugated(G,P):
    '''
    /*****************************************************
    Corrugated fin, B.S.
    For tube diameters >= 1/2", 
    Wang, C. C., Y. T. Lin, C. J. Lee, and Y. J. Chang, "Investigation of 
    Wavy Fin-and-Tube Heat Exchangers: A Contribution to Databank", Experimental 
    Heat Transfer, 12:73-89(1999) 
    For tube diameters < 1/2", 
    Wang, C. C., J. Y. Jang, and N. F. Chiou, "A Heat Transfer and Friction Correlation 
    For Wavy Fin-and-tube Heat Exchangers", International Journal of Heat and Mass 
    Transfer, Vol 42(1999) pp.1919-1924. 
    *********************************************************/
    '''
 
    th = P['th'];#fin thichness
    D = P['Do']+2*th;#collar diameter
    D_h= P['Dh'];#hydraulic diameter
    F_p = P['z'];#fin pitch
    P_t=P['vsp'];#tube pitch perpendicular to the airflow
    P_l=P['P_l'];#tube pitch in the airflow direction
    N=P['Ndeep'];#tube number in the airflow direction
    F_s=F_p-th;#fin space
    T = P['T'];#air temperature
    X_f=P_l/4;#projected fin pattern length for one-half wavy length, typical structure from the paper
    P_d=1.32*0.001;#waffle height, from the paper
    sec_angle=pow((pow(X_f,2.0)+pow(P_d,2.0)),0.5)/X_f;# calculate the sec of the corrugation angle
    tg_angle=pow((pow(sec_angle,2.0)-1.0),0.5);#calculate the tan of the corrugation angle
 
    #beta=3.1415926*pow(D,2.0)/(4*P_l*P_t);#possible method for calculating hydraulic diameter for wavy fin, but not used here
    #D_h=2*F_p*(1-beta)/((1-beta)*sec_angle+2*F_p*beta/D);
    
    mu = HAPropsSI('mu','T',P['T'],'P',101325,'R',P['P']) #air.mu(T);#air viscosity [Pa-s]
            
    Re = G*D/mu;#Reynold's number based on the collar diameter
    if (Re<0.0):
        print("FricAir_Corrugated::FricAir_Corrugated","Re<0.0")
     
     
    Coe_fin = ((P['Apo']+P['Af'])/(P['Ls']*P['Do']*pi));
 
    if (D<=25.4e-3/2):
        F1=0.4604-0.01336*pow(F_p/P_l, 0.58)*log(Coe_fin)*pow(tg_angle,-1.5);
        F2 = 3.247*pow(F_p/P_t,1.4)*log(Coe_fin);
        F3= -20.113/log(Re);
        f = 0.01915*pow(Re,F1)*pow(tg_angle,F2)*pow(F_p/P_l,F3)*pow(log(Coe_fin),-5.35)*pow(P['Dh']/D,1.3796)*pow(N,-0.0916);
    else:
        f1=0.1714-0.07372*pow(F_p/P_l,0.25)*log(Coe_fin)*pow(P_d/X_f, -0.2);
        f2 = 0.426*pow(F_p/P_t,0.3)*log(Coe_fin);
        f3 = -10.2192/(log(Re));
        f = 0.05273*pow(Re,f1)*pow(P_d/X_f, f2)*pow(F_p/P_t,f3)*pow(log(Coe_fin),-2.726)*pow(P['Dh']/D,0.1325)*pow(N,0.02305);
 
    return f, P


def FricAir_Slit(G,P):
    '''
    /*****************************************************
    Slit fin, B.S.
    Wang, C. C., 2001, "A Comparative Study of Compact Enhanced Fin-and-Tube 
    Heat Exchangers", Int. J. Heat and Mass Transfer, Vol. 44, pp. 3565-3573. 
    *********************************************************/
    '''
 
    th = P['th'];#fin thichness
    D = P['Do']+2*th;#collar diameter
    D_h= P['Dh'];#hydraulic diameter
    F_p = P['z'];#fin pitch
    P_t=P['vsp'];#tube pitch perpendicular to the airflow
    P_l=P['P_l'];#tube pitch in the airflow direction
    N=P['Ndeep'];#tube number in the airflow direction
    F_s=F_p-th;#fin space
    T = P['T'];#air temperature
    S_h=0.99*0.001;#height of slit, typical structure from the paper
    S_s=2.2*0.001;#breadth of a slit in the direction of the airflow
    S_w=11*0.001;#width of slit
    S_n=4;#number of slits in an enhanced zone
     
    mu = HAPropsSI('mu','T',P['T'],'P',101325,'R',P['P']) #air.mu(T);#air viscosity [Pa-s]
            
    Re = G*D/mu;#Reynold's number based on the collar diameter
    if (Re<0.0):
        print("FricAir_Slit::FricAir_Slit","Re<0.0")

 
    f1 = -0.1401+0.2567*log(F_s/D)+4.399*exp(-1*S_n);
    f2 = -0.383+0.7998*log(F_s/D)+5.1772/S_n;
    f3 = -1.7266 - 0.1102*log(Re)-1.4501*(F_s/D);
    f4 = 0.4034-0.199*(S_s/S_h)/log(S_s/S_h)+0.4208*log(S_s/S_h)/pow(S_s/S_h,2.0);
    N_ref=N;
    if (N==1):
        N_ref = 2;
    f5 = -9.0566+0.6199*log(Re)+32.8057/log(Re)-0.2881/log(N_ref)+0.9583/pow(N,1.5);
    f6 = -1.4994+1.209*(P_t/P_l)+1.4601/S_n;
     
    f = 1.201*pow(Re,f1)*pow(F_s/D,f2)*pow(P_t/P_l,f3)*pow(S_s/S_h,f4)*pow(N,f5)*pow(S_n,f6);
 
    return f, P


def FricAir_Louvered(G,P):
    '''
    /*****************************************************
    Louvered fin, B.S.
    Wang, C. C., C. J. Lee, C. T. Chang, and S. P. Lin, 
    "Heat Transfer and Friction Correlation for Compact Louvered Fin-and-tube
     Heat Exchangers", International Journal of Heat and Mass Transfer, 
     Vol 42 (1999) pp.1945-1956. 
    *********************************************************/
    '''
 
    th = P['th'];#fin thichness
    D = P['Do']+2*th;#collar diameter
    D_h= P['Dh'];#hydraulic diameter
    F_p = P['z'];#fin pitch
    P_t=P['vsp'];#tube pitch perpendicular to the airflow
    P_l=P['P_l'];#tube pitch in the airflow direction
    N=P['Ndeep'];#tube number in the airflow direction
    F_s=F_p-th;#fin space
    T = P['T'];#air temperature
    L_h=1.07*0.001;#Louver height, typical structure from the paper
    L_p=2.4*0.001;#major louver pitch
 
     
    mu = HAPropsSI('mu','T',P['T'],'P',101325,'R',P['P']) #air.mu(T);#air viscosity [Pa-s]
            
    Re = G*D/mu;#Reynold's number based on the collar diameter
    if (Re<0.0):
        print("FricAir_Louvered::FricAir_Louvered","Re<0.0")
 
    
    Coe_fin = ((P['Apo']+P['Af'])/(P['Ls']*P['Do']*pi));
 
    if(N>1):
        F5 = 0.1395-0.0101*pow(F_p/P_l, 0.58)*pow(L_h/L_p,-2.0)*log(Coe_fin)*pow(P_l/P_t, 1.9);
        F6 = -6.4367/log(Re);
        F7 = 0.07191*log(Re);
        F8 = -2.0585*pow(F_p/P_t, 1.67)*log(Re);
        F9 = 0.1036*log(P_l/P_t);
        f = 0.06393*pow(Re, F5)*pow(F_p/D, F6)*pow(P['Dh']/D, F7)*pow(L_h/L_p, F8)*pow(N, F9)*pow(log(Re)-4.0, -1.093);
    else:
        F1 = 0.1691+4.4118*pow(F_p/P_l,-0.3)*pow(L_h/L_p,-2.0)*log(P_l/P_t)*pow(F_p/P_t,3.0);
        F2 = -2.6642-14.3809/log(Re);
        F3 = -0.6816*log(F_p/P_l);
        F4 = 6.4668*pow(F_p/P_t,1.7)*log(Coe_fin);
        f = 0.00317*pow(Re,F1)*pow(F_p/P_l,F2)*pow(P['Dh']/D,F3)*pow(L_h/L_p,F4)*pow(log(Coe_fin), -6.0483);
 
    return f, P


def FricAir_ConvexLouvered(G,P):
    '''
    /*****************************************************
    Convexlouvered fin, B.S.
    Wang, C. C., Y. M. Tsai, and D. C. Lu, 
    "Comprehensive Study of Convex-louver and Wavy Fin-and-tube Heat Exchangers", 
    Journal of Thermophysics and Heat Transfer, Vol 12, No. 3, July-September 1998,
    pp.423-430. 
    *********************************************************/
    '''
 
    th = P['th'];#fin thichness
    D = P['Do']+2*th;#collar diameter
    D_h= P['Dh'];#hydraulic diameter
    F_p = P['z'];#fin pitch
    P_t=P['vsp'];#tube pitch perpendicular to the airflow
    P_l=P['P_l'];#tube pitch in the airflow direction
    N=P['Ndeep'];#tube number in the airflow direction
    F_s=F_p-th;#fin space
    T = P['T'];#air temperature
 
    mu = HAPropsSI('mu','T',P['T'],'P',101325,'R',P['P']) #air.mu(T);#air viscosity [Pa-s]
            
    Re = G*D/mu;#Reynold's number based on the collar diameter
    if (Re<0.0):
        print("FricAir_ConvexLouvered::FricAir_ConvexLouvered","Re<0.0")
        
        
    Coe_fin = ((P['Apo']+P['Af'])/(P['Ls']*P['Do']*pi));
    
    if (Re<1e3):
        f = 0.264*(0.105+0.708*exp(-1*Re/225))*pow(Re,-0.637)*pow(Coe_fin,0.263)*pow(F_p/D, -0.317);
    else:
        f = 0.768*(0.0494+0.142*exp(-1*Re/1180))*pow(Coe_fin,0.0195)*pow(F_p/D,-0.121);

     
    return f, P


def FricAir_SmoothWavy(G,P):
    '''
    /*****************************************************
    Smoothwavy fin, B.S.
    Mirth, D. R. and S. Ramadhyani, 1994,"Correlations for Predicting the Airside 
    Nusselt Numbers and Friction Factors in Chilled-water Cooling Coils", 
    Experimental Heat Transfer, Vol 7: 143-162.
    ---the author suggested predict both the dry heat transfer and wet heat transfer with the dry correlations
    *********************************************************/
    '''
 
    th = P['th'];#fin thichness
    D = P['Do']+2*th;#collar diameter
    D_h= P['Dh'];#hydraulic diameter
    F_p = P['z'];#fin pitch
    P_t=P['vsp'];#tube pitch perpendicular to the airflow
    P_l=P['P_l'];#tube pitch in the airflow direction
    N=P['Ndeep'];#tube number in the airflow direction
    F_s=F_p-th;#fin space
    T = P['T'];#air temperature
    wh = 2.38e-3;#typical geometry from the paper, with respect to first three coils in the paper
     
    mu = HAPropsSI('mu','T',P['T'],'P',101325,'R',P['P']) #air.mu(T);#air viscosity [Pa-s]
            
    Re = G*D/mu;#Reynold's number based on the collar diameter
    if (Re<0.0):
        print("FricAir_ConvexLouvered::FricAir_ConvexLouvered","Re<0.0")
        
 
    f = 8.64/pow(Re,0.457)*pow(2*F_s/wh,0.473)*pow(N*P_l/wh,-0.545);
     
    if (P['P']>0.4):#wet condition
        f = 2.71/pow(Re,0.737) +f;
    else: #dry condition
        f = f;
 
    return f, P


def FricAir_Spine(G,P):
    '''
    /*****************************************************
    Spine fin, empty now
    *********************************************************/
    '''
 
    f=0;
     
    return f, P


#===========================================================================
# single-phase heat transfer
#===========================================================================

def ConvCoeffSP(TXP,G,P, Ref):
    '''
    #interface function for calling the single-phase heat transfer calculation for the EVAPORATOR & CONDENSER
    '''
 
    D = InsideTube_Dim(); #initiate dictionarty for storing the parameters for micro-fin tube
    D['Microfin'] = P['microfin']; #microfin type, 0=smooth tube, 1=helical, 2=cross-grooved, 3=herringbone
    D['finN'] = P['finN']; #fin number in a micro-fin tube
    D['Di'] = P['Di'];#inside diameter at the fin tip
    D['gama'] =P['gama'] ;#fin apex angle in a micro-fin tube
    D['beta'] = P['beta'] ;    #fin helix angle in a micro-fin tube
    D['finH'] = P['finH']; #fin height in a micro-fin tube
    D['w_b'] = P['w_b'] ; #base width of a single fin
    D['w_e'] = P['w_e']; #top width of a single fin
    D['w_z'] = P['w_z']; #base distance between two neighboring fins
    D['K_T'] = P['K_T'];#400, this is the conductance factor of copper
    D['Ls'] = P['Ls'];#tube unit length
    D['D_b'] = P['D_b'];#tube diameter at the base of the fin
    D['Do'] = P['Do']; #Pipe outside diameter.
    D['D_m'] = P['D_m']; #mean diameter of the micro-fin tube
    D['P_H'] = P['P_H'];# the hydraulical circumference
    D['Acs'] = P['Acs'] ;#cross area of the micro-fin tube, this is the actual cross-section area
    D['Dh_i'] = P['Dh_i'];#inside hydraulical diameter 
    D['Ax'] = P['Ax'];# Inside pipe cross sectional area, based on fin tips
    D['Api'] =P['Api'];# Inside pipe surface area (one tube segment), based on fin tips
    h_sp, D = ConvCoeffSP_Microfin(TXP,G,D, Ref);

    return h_sp, P

#The following function is the same as the previous, therefroe it is commented out
# def ConvCoeffSP(TXP TXP,double G,CGP* P):
#     '''
#     #interface function for calling the single-phase heat transfer calculation for the condenser
#     '''
# 
#     InsideTube_Dim D;
#     D.Microfin = P->Microfin; #microfin type, 0=smooth tube, 1=helical, 2=cross-grooved, 3=herringbone
#     D.finN = P->finN; #fin number in a micro-fin tube
#     D.Di = P->Di;#inside diameter at the fin tip
#     D.gama =P->gama ;#fin apex angle in a micro-fin tube
#     D.beta = P->beta ;    #fin helix angle in a micro-fin tube
#     D.finH = P->finH; #fin height in a micro-fin tube
#     D.w_b = P->w_b ; #base width of a single fin
#     D.w_e = P->w_e; #top width of a single fin
#     D.w_z = P->w_z; #base distance between two neighboring fins
#     D.K_T = P->K_T;#400, this is the conductance factor of copper
#     D.Ls = P->Ls;#tube unit length
#     D.D_b = P->D_b;#tube diameter at the base of the fin
#     D.Do = P->Do; #Pipe outside diameter.
#     D.D_m = P->D_m; #mean diameter of the micro-fin tube
#     D.P_H = P->P_H;# the hydraulical circumference
#     D.Acs = P->Acs ;#cross area of the micro-fin tube, this is the actual cross-section area
#     D.Dh_i = P->Dh_i;#inside hydraulical diameter 
#     D.Ax = P->Ax;# Inside pipe cross sectional area, based on fin tips
#     D.Api =P->Api;# Inside pipe surface area (one tube segment), based on fin tips
#     const double h_sp = ConvCoeffSP_Microfin(TXP,G, &D);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffSP","h_sp_con");
#         return -1;
#     }
#     return h_sp


def ConvCoeffSP_Microfin(TXP,#refrigerant state
                         G, #refrigerant mass flux
                         P, Ref): #micro-fin tube struct
    '''
    /*****************************************************
    Single-phase heat transfer for the micro-fin tubes, B.S.,
    Tang, Liangyou. Ohadi, Michael M. Johnson, Arthur T. Flow condensation in smooth 
    and micro-fin tubes with HCFC-22, HFC-134a and HFC-410 refrigerants. 
    Part II: design equations Journal of Enhanced Heat Transfer. v 7 n 5 2000. 
    p 311-325
    *********************************************************/
    '''
     
    mu = PropertyTXPtr('V',TXP, Ref) #[Pa-s]
    Cp = PropertyTXPtr('C',TXP, Ref) #[J/kg/K]
    k = PropertyTXPtr('L',TXP, Ref) #[W/m/K]
    
    #double d_bottom,d_top,d_m,c,a;
    #double h;
    d_top = P['Di'];#diameter at the fin tip
    d_bottom = P['D_b'];#diameter at the fin bottom
    d_m = P['D_m'];#mean diameter
     
    
    if (TXP['X']>0.99 or P['Microfin']==0): #Dittus-Boelter
        G=G*pow(d_top/d_bottom,2.0);
        h=ConvCoeffSP_Smooth(TXP,G,d_bottom)
        h=h*d_bottom/d_top;
        return h, P
 
    Pr = mu*Cp/k;
    G_adj=G*pow(d_top/d_m,2.0);#the author mentioned the mean diamter was used to calculate the flow cross-sectional area
    Re = G_adj*d_bottom/mu;#the author mentioned the inside diamter at the fin bottom was used to calculate inner tube heat transfer surface area
 
    if (P['Microfin']==2):#cross-grooved tube
        c=0.012;
        a=0.95;
    elif (P['beta']<2): #axial fin tube
        c=0.0136;
        a=0.91;
    else: #helical fin tube
        a=0.8;
        c=0.0479; 
     
    h=c*pow(Re,a)*pow(Pr,0.4)*k/(d_bottom);
    h=d_bottom/d_top*h;#normalized the heat transfer coefficient to the inner surface area at the fin tip
     
    return h, P


def ConvCoeffSP_Smooth(TXP,G,D, Ref):
    '''
    /********************************************************************
    Dittus-Boelter correlation for convective heat transfer
    from flow inside of circular pipes.
    Refernece:
        1. Incopera & Dewitt, 3rd ed., p.496
        2. Dittus, F.W. L.M.K. Boelter
            University of California Publications on Engineering
            Vol. 2, p.443, 1930
    TXP = thermodynamic state of refrigerant (K/-/Pa)
    G = refrigerant mass fux (kg/s/m^2)
    D = inside diameter of pipe. (m)
    return = the convection coefficient (h) in (W/K/m^2)
    ********************************************************************/
    '''
 
    mu = PropertyTXPtr('V',TXP, Ref) #[Pa-s]
    Cp = PropertyTXPtr('C',TXP, Ref) #[J/kg/K]
    k = PropertyTXPtr('L',TXP, Ref) #[W/m/K]
    
    Pr = mu*Cp/k;#Boelter correlation
    Re = G*D/mu;
    Nu = 0.023*pow(Re,0.8)*pow(Pr,0.333);
    
    if (TXP['X']>0.9):
        h=Nu*k/D;
    else:
        h=Nu*k/D;    
     
    return h


#===============================================================================
# condensation two-phase heat transfer
#===============================================================================

# def ConvCoeffInside(TXP TXPi,#refrigerant inlet state
#                 double G,#refrigerant mass flux
#                 double D,#inside diameter
#                 CGP *P):#condenser struct
#     '''
#     /********************************************************************
#     Inside convection coefficient for the refrigerant inside the
#     condenser.
#     ********************************************************************/
#     '''
# 
#     double y;
#     const double X1=0.1,X2=0.9;
# 
#     if(TXPi.X>=1.0 || TXPi.X<=0) {
#         y = ConvCoeffSP(TXPi,G,P);
#     } else if(TXPi.X<X1) {
#         TXP TXP1 = toTXP(TXPi.T,0,TXPi.P);
#         TXP TXP2 = toTXP(TXPi.T,X1,TXPi.P);
#         double y1 = ConvCoeffSP(TXP1,G,P);
#         if(errorLog.IsError()) {
#             errorLog.Add("ConvCoeffInside");
#             return -1;
#         }
#         double y2 = ConvCoeffCondTP_Microfin(TXP2,G,P);
#         y = y1+TXPi.X*(y2-y1)/X1 ;
#         if(y2<y1) y=y1;#B.S.
# 
#     } else if(TXPi.X>X2) {
#         TXP TXP1 = toTXP(TXPi.T,X2,TXPi.P);
#         TXP TXP2 = toTXP(TXPi.T,1,TXPi.P);
#         double y1 = ConvCoeffCondTP_Microfin(TXP1,G,P);
#         if(errorLog.IsError()) {
#             errorLog.Add("ConvCoeffInside");
#             return -1;
#         }
#         double y2 = ConvCoeffSP(TXP2,G,P);
#         y = y2-(1-TXPi.X)*(y2-y1)/(1-X2);
#         if(y1<y2) y=y2;#B.S.
#     } else {
#         y = ConvCoeffCondTP_Microfin(TXPi,G,P);#two-phase heat transfer
#     }
# 
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffInside");
#         return -1;
#     }
# 
#     return y


# def ConvCoeffCondTP_Microfin(TXP TXPm,#refrigerant inlet state
#                      double G,#mass flux
#                      CGP* P):#condenser struct
#     '''
#     /***************************************************************************
#     condensing heat transfer correlation for micro-fin tube, B.S.
#     A. Cavallini, D. Del Col, L. Doretti, G. A. Longo and L. Rossetto 
#     "Heat transfer and pressure drop during condensation of refrigerants inside horizontal 
#     enhanced tubes;" Transfert de chaleur et chute de pression lors de la condensation de 
#     frigorignes l'intrieur de tubes horizontaux surface augme
#     nte, International Journal of Refrigeration, Volume 23, Issue 1, January 2000, Pages 4
#     ************************************************************************/
#     '''
# 
#     const double D = P->Di;
#     TXP TXP_prop={0,0,0};
# 
#     if(P->Microfin<1)#smooth tube
#     {
#     const double h_smooth = ConvCoeffCondTP_Smooth(TXPm,G,D);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffTP_Smooth","h_smooth");
#         return -1;
#     }
#     return h_smooth;
#     }
#     if(P->Microfin==3)#herringbone
#     {
#     const double h_Herringbone = ConvCoeffCondTP_Herringbone(TXPm,G,P);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffTP_Herringbone","h_Herringbone");
#         return -1;
#     }
#     return h_Herringbone;
#     }
# 
#     #Liquid refrigerant properties
#     TXP_prop.P= TXPm.P;
#     TXP_prop.X=0;
#     TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
# 
#     const double Tsat_l=PropertyTXPth(TSAT,TXP_prop);
# 
#     const double h_l=PropertyTXPth(ENTH,TXP_prop);
# 
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffTP_Microfin","h_l");
#         return -1;
#     }
# 
#     const double Tension=PropertyTXPtr(TENSION,TXP_prop);#reftpltrP.Tension(TXPm.P);#refrigerant surface tension
# 
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffTP_Microfin","Tension");
#         return -1;
#     }
# 
# 
#     const double rho_l=1.0/PropertyTXPth(VOL,TXP_prop);#1.0/reftplthP.v(TXPm.P);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffTP_Microfin","vl");
#         return -1;
#     }
# 
#     const double mu_l=PropertyTXPtr(VISC,TXP_prop);#refsctrPT.mu(TXPm.P,Tsat);#refrigerant liquid viscosity
#     
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffTP_Microfin","mul");
#         return -1;
#     }
# 
#     
#     const double Cp_l= PropertyTXPtr(SPEC,TXP_prop);#refsctrPT.Cp(TXPm.P,Tsat);#refrigerant liquid specific heat
# 
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffTP_Microfin","Cpl");
#         return -1;
#     }
# 
#     const double k_l=PropertyTXPtr(COND,TXP_prop);#refsctrPT.k(TXPm.P,Tsat);#refrigerant liquid heat conductance
# 
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffTP_Microfin","kl");
#         return -1;
#     }
# 
#     #vapor refrigerant properties
#     TXP_prop.P=TXPm.P;
#     TXP_prop.X=1.0;
#     TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
# 
#     const double Tsat_v=PropertyTXPth(TSAT,TXP_prop);
# 
#     const double h_g=PropertyTXPth(ENTH,TXP_prop);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffTP_Microfin","h_g");
#         return -1;
#     }
# 
#     const double rho_g=1.0/PropertyTXPth(VOL,TXP_prop);#1.0/reftpvthP.v(TXPm.P);#refrigerant gas density
# 
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffTP_Microfin","vv");
#         return -1;
#     }
# 
# 
#     const double mu_g=PropertyTXPtr(VISC,TXP_prop);#refshtrPT.mu(TXPm.P,Tsat);#refrigerant gas viscosity
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffTP_Microfin","muv");
#         return -1;
#     }
# 
#     const double Cp_g= PropertyTXPtr(SPEC,TXP_prop);#refsctrPT.Cp(TXPm.P,Tsat);#refrigerant liquid specific heat
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffTP_Microfin","Cpg");
#         return -1;
#     }
# 
#     const double k_g=PropertyTXPtr(COND,TXP_prop);#refsctrPT.k(TXPm.P,Tsat);#refrigerant liquid heat conductance
# 
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffTP_Microfin","kg");
#         return -1;
#     }
# 
#     const double h_fg=h_g-h_l;
#     #parameters for the correlation
# 
#     const double x=TXPm.X;#vapor quality
# 
#     const double pi=3.1415926;
#     const double gama=P->gama/180.0*pi;#fin apex angle, be attention to the unit of the angle
#     const double e=P->finH;#fin height
#     const double d_e=P->Di;#inner diameter at the fin tip
#     const double d_i=P->D_b;#d_e+2*e, innder diameter at the fin bottom
#     const double d_M=P->D_m;#mean inner diameter of the micro-fin tube
#     const double n_g=P->finN;#fin number
#     const double beta=P->beta/180.0*pi;#fin helical angle
# 
#     double s=2.00;
#     double t=-0.26;
#     
#     #empirical parameters corresponding to different tube type
#     if(P->Microfin==2)#cross-grooved
#     {s = 2.1; t=-0.26;}
#     else if(e/d_e>0.045)#low-fin tube, the original standard is e/d_e>0.04, here enlarge the range a bit.
#     {s = 1.4; t=-0.08;}
#     else#micro-fin tube
#     {s = 2.0; t=-0.26;}    
# 
# 
#     const double sigma=Tension;
# 
#     const double u_go=G/rho_g;#gas velocity based the whole mass flow rate
# 
#     const double Fr=pow(u_go,2.0)/(9.8*d_e);
# 
#     const double Bo=9.8*rho_l*e*pi*d_e/(8*sigma*n_g);
# 
#     const double Rx=((2*e*n_g*(1-sin(gama/2))/(pi*d_e*cos(gama/2))+1))/cos(beta);#geometry enhancement factor
#     const double Re_eq=4*G*pi*(pow(d_e,2.0)/4.0)*((1-x)+x*pow((rho_l/rho_g),0.5))/(pi*d_e*mu_l);#Reynolds number
# 
#     const double Pr_l=mu_l*Cp_l/k_l;#Prontal number
#     const double Pr_g=mu_g*Cp_g/k_g;#Prontal number
# 
#     const double Nusselt=0.05*pow(Re_eq,0.8)*pow(Pr_l,(0.33333333333))*pow(Rx,double(s))*pow((Bo*Fr),double(t));#Nusselts number
# 
#     const double h_cava=Nusselt*k_l/d_e;#heat transfer coefficient based on the surface area at the fin tip
# 
#     
#     double h=h_cava;
# 
# #ifdef _RefMix
#     #if zerotropic refrigerant, correct the overall condensation coefficient by considering the mass transfer resistance between vapora phase and liquid phase
#     {
#     const double T_delta = Tsat_v-Tsat_l;
#     double Corr_FlowBoiling=1.0;#parameter for the mass transfer resistance between the vapor phase and liquid phase
#     Corr_FlowBoiling = Correct_FLOW_Boiling(TXPm.X,Cp_g,T_delta,h_fg);#correction parameter for the mass transfer resistance between the liquid phase and vapor phase 
#     const double Re_vaporphase = G*TXPm.X*(d_e+2*e)/(mu_g);#Reynolds number, asuming the vapor only flowing in the tube
#     const double Nu = 0.023*pow(Re_vaporphase,0.8)*pow(Pr_g,0.333);#Dittus-Boelter equation to calculate the vapor phase coeffcient
#     const double h_vaporphase=Nu*k_g/(d_e+2*e);#vapor phase heat transfer coefficent
#     h=1.0/(1.0/h+Corr_FlowBoiling/h_vaporphase);#correct the overall flow boiling coefficient
#     }
# #endif
# 
#     return h


# def ConvCoeffCondTP_Smooth(TXP TXPm,#refrigerant inlet state
#                    double G,#refrigerant mass flux
#                    double D):#inside diameter
#     '''
#     /********************************************************************
#     Condensation in horizontal tubes, part 1: two-phase flow pattern map  ? ARTICLE
#     International Journal of Heat and Mass Transfer, Volume 46, Issue 18, August 2003, Pages 3349-3363 
#     J. El Hajal, J. R. Thome and A. Cavallini
#     Condensation in horizontal tubes, part 2: new heat transfer model based on flow regimes  ? ARTICLE
#     International Journal of Heat and Mass Transfer, Volume 46, Issue 18, August 2003, Pages 3365-3387 
#     J. R. Thome, J. El Hajal and A. Cavallini
#     ********************************************************************/
#     '''
# 
#     FlowPattern Cd;
#     const double q=0;    
#     Cd.JudgPattern=0;
#     Cond_FlowPattern(TXPm,G,D,q,&Cd);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffCondTP_Smooth","h");
#         return -1;}
#         
#     const double h=Cd.h_tp;
#     
#     return h


# def ConvCoeffCondTP_Herringbone(TXP TXPm,#refrigerant state
#                             double G,#refrigerant mass flux
#                             CGP *P):#condenser struct
#     '''
#     #B.S.-------------------------------------------------------
#     #Herringbone
#     /***********************************************************
#     Miyara, Akio Kengo Nonaka and Mitsunori Taniguchi 
#     "Condensation heat transfer and flow pattern inside a herringbone-type 
#     micro-fin tube;" Transfert de chaleur lors de la condensation et configuration de
#      l'coulement l'intrieur d'un tube microailettes che
#      vrons, International Journal of Refrigeration, Volume 23, Issue 2, March 2000, 
#      Pages 141-152
#     **********************************************************/
#     '''
# 
#     double H_I;    
#     const double d_M = P->D_m;
#     double delta_T=3;
#     double delta_T1=0;
#     const double PI = 4*atan(1.0);
#     TXP TXP_prop={0,0,0};
# 
#     G=G*pow(P->Di/P->D_m,2.0);
# 
#     #liquid refrigerant properties
#     TXP_prop.P=TXPm.P;
#     TXP_prop.X=0.0;
#     TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
#     const double H_F = PropertyTXPth(ENTH,TXP_prop);#reftplthP.h(TXPm.P);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffCond_Herringbone","hl");
#         return -1;
#     }
# 
# 
#     const double DL=1.0/PropertyTXPth(VOL,TXP_prop);#1.0/reftplthP.v(TXPm.P);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffCondTP_Herringbone","vl");
#         return -1;
#     }
# 
# 
#     const double MU_F=PropertyTXPtr(VISC,TXP_prop);#refsctrPT.mu(TXPm.P,TXPm.T);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffCondTP_Herringbone","ul");
#         return -1;
#     }
# 
#     const double CP_F= PropertyTXPtr(SPEC,TXP_prop);#refsctrPT.Cp(TXPm.P,T_sat);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffCondTP_Herringbone","Cpl");
#         return -1;
#     }
# 
#     const double TC_F=PropertyTXPtr(COND,TXP_prop);#refsctrPT.k(TXPm.P,T_sat);
# 
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffCondTP_Herringbone","kl");
#         return -1;
#     }
# 
#     #vapor refrigerant properties
#     
#     TXP_prop.P=TXPm.P;
#     TXP_prop.X=1.0;
#     TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
#     const double H_G = PropertyTXPth(ENTH,TXP_prop);#reftpvthP.h(TXPm.P);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffCondTP_Herringbone","hv");
#         return -1;
#     }
#     
#     const double DV=1.0/PropertyTXPth(VOL,TXP_prop);#1.0/reftpvthP.v(TXPm.P);
#     
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffCondTP_Herringbone","vv");
#         return -1;
#     }
# 
# 
# 
#     const double CP_A = air.Cp(P->airT);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffCondTP_Herringbone","CP_A");
#         return -1;
#     }
#     const double C = P->Ga*P->Aflow*CP_A;
# 
#     const double h_fg=H_G-H_F;
#     double xx=TXPm.X;
#     const double alfa=pow((1+((1-xx)*DV/(xx*DL))*(0.4+0.6*pow((xx*(DL/DV)+0.4*(1-xx)),0.5)/pow((xx+0.4*(1-xx)),0.5))),(-1));
#     const double Hfuc=alfa+(10*pow((1-alfa),0.1)-8.0)*pow((alfa),0.5)*(1-pow(alfa,0.5));
#     const double Gaa=9.8*pow(DL,2.0)*pow(d_M,3.0)/pow(MU_F,2.0);
#     const double eta_A=P->P_H/(PI*d_M);# "this is the heat transfer surface area increasing ratio" "pi*d_M*10"
#     
# 
#     const double PR_F=MU_F*CP_F/TC_F;
# 
#     while(fabs(delta_T1-delta_T)>=1e-7)#simple iteration for temporary use
#     {
#     delta_T1=delta_T;
#     const double Ja=h_fg/(CP_F*delta_T);
#     const double Nusselt_B=0.725*Hfuc*pow((Gaa*PR_F*Ja/(eta_A)),0.25);
#     const double Re_l_yu=G*(1-xx)*d_M/MU_F;
#     const double X_tt=Xtt(TXPm);
#     const double phi_g=1.2+1.65*(pow(G,0.35)*pow(X_tt,0.35)/pow((9.8*d_M*DV*(DL-DV)),0.175));
#     const double Nusselt_F=0.152*(phi_g/X_tt)*pow(Re_l_yu,0.68)*(0.3+0.1*pow(PR_F,1.1));
#     const double Nusselt_yu=pow((pow(Nusselt_B,2.0)+pow(Nusselt_F,2.0)),0.5);
#     const double h_yu=Nusselt_yu*TC_F/d_M;
#     H_I=h_yu*d_M/P->Di;
#     const double R_CONV=1e0/(H_I*P->P_H*P->Ls);
#     const double UA=1e0/(R_CONV+P->Ro);
#     const double NTU=1*UA/C;
#     const double EPSILON=1e0-exp(-NTU);
#     const double Q=EPSILON*C*(P->airT-TXPm.T);
#     delta_T=Q*R_CONV;
#     if(delta_T<0.01) delta_T=0.01;
#     }
#     return H_I

#===============================================================================
# evaporation two-phase heat transfer
#===============================================================================

def ConvCoeffEvapTP_microfin(TXPm,#refrigerant state
                             G,#refrigerant mass flux
                             P, Ref):#evaporator struct
    '''
    /***********************************************
    flow boiling correlation for micro-fin tube
    Cavallini A.. Del Col D., Doretti L. Longo G. A. Rossetto L. 
    "Refrigerant vaporization inside enhanced tubes, a heat transfer model" 
    Heat and technology, Vol. 17. n.2 1999
    *********************************************/
    '''
 
    Di = P['Di'];#inside diameter at the fin tip
    TXP_prop={'T':0.0,'X':0.0,'P':0.0};
 
    q=P['q_flux'];#for smooth tube evaporation model iteration
    T_w = P['T_w'];#for micro-fin tube evaporation model iteration
 
    if (not P['microfin']):#smooth tube
        h_smooth = ConvCoeffEvapTP_Smooth(TXPm,G,Di,q, Ref)
        return h_smooth, P 
     
    gama =P['gama'] ;#apex angle of the fin
    e = P['finH'];#fin height
    n_g = P['finN'];#fin number
    beta = P['beta'];#fin helical angle
    g=9.807;#gravational constant
    d_e=Di;#tube inside diameter at the fin tip
 
    #parameters
    A=1.36;
    B=0.36;
    C=0.38; 
    SS=2.14;
    if (G<500):
        T=-0.15;
    else:
        T=-0.21;
    V=0.59;
    Z=0.36;
    G_0=100;
 
    gama=gama*pi/180.0;#apex angle
    if (beta>30.0):
        beta = 30;#the author suggest the applicable range for this correlation
    beta=beta*pi/180;#helical angle
    d_0=0.01;
 
    if (G<100):
        A=1.36*sin(beta);
        B=0.36*pow((G/100.0),4.0);
        Z=-3.0;
 
    P_cr = PropsSI('PCRIT', Ref)     #critical pressure [Pa]
    T_cr = PropsSI('TCRIT', Ref)     #critical temperature [K]
    M = PropsSI('TCRIT', Ref)*1000   #molecular mass [kg/kmol]

    P_sat=TXPm['P'];#saturation pressure
 
    delta_T=T_w-TXPm['T'];#temperture difference between the tube wall and the refrigerant
    
    X_tt=Xtt(TXPm, Ref);#Martinelli parameter
    if (X_tt>1.0):
        X_tt=1.0;# the author suggested this restriction
     
    S=A*pow(X_tt,B);
    F=pow((d_0/d_e),C);
    P_R=P_sat/P_cr;
    
    
    #liquid refrigerant properties
    TXP_prop['P']=TXPm['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    sigma=PropertyTXPtr('I',TXP_prop, Ref) #[N/m]
    rho_l = PropertyTXPth('D',TXP_prop, Ref) #[kg/m^3]
    mu_l = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]
    Cp_l = PropertyTXPtr('C',TXP_prop, Ref) #[J/kg/K]
    k_l = PropertyTXPtr('L',TXP_prop, Ref) #[W/m/K]
    hl = PropertyTXPth('H',TXP_prop, Ref) #[J/kg/K]
    Tsat_l = PropertyTXPth('T',TXP_prop, Ref) #[K]
    
    #vapor refrigerant properties
    TXP_prop['P']=TXPm['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    rho_g = PropertyTXPth('D',TXP_prop, Ref) #[kg/m^3]
    mu_g = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]
    Cp_g = PropertyTXPtr('C',TXP_prop, Ref) #[J/kg/K]
    k_g = PropertyTXPtr('L',TXP_prop, Ref) #[W/m/K]
    hv = PropertyTXPth('H',TXP_prop, Ref) #[J/kg/K]
    Tsat_v = PropertyTXPth('T',TXP_prop, Ref) #[K]
    
         
    h_fg = hv-hl;#latent heat
    T_delta=Tsat_v-Tsat_l;#gliding temperature difference for zeotropic refrigerant
    Pr_l=mu_l*Cp_l/k_l;#prontal number
    Pr_g=mu_g*Cp_g/k_g;#prontal number
 
    if (delta_T<=0.0):
        delta_T=1e-20;# to remove the wrong result
    h_nb=55*pow(P_R,0.12)*pow(M,(-0.5))*pow((-log10(P_R)),(-0.55))*pow(q,0.67)*S*F;#calculate the nucleate boiling coefficient
    q_nb=q;#nucleate boiling heat flux
 
    
    if (abs(T_delta)>=0.1):#if zerotropic refrigerant, correct the nucleate boiling coefficent
        Corr_NUC=1.0;#parameter for the effect of mass transfer resistance on nucleate boiling
        B_0=1.0;#scaling factor
        beta_l = 3e-4;#m/s, mass transfer coefficient
        Corr_NUC = Correct_NUC_Boiling(h_nb,q_nb,T_delta,B_0,beta_l,rho_l,h_fg);#nucleate boiling correction
        h_nb=h_nb*Corr_NUC;#corrected
 
    u_go=G/rho_g;#all gas phase velocity
    Fr=pow(u_go,2)/(9.8*d_e);
    Bo=9.8*rho_l*e*pi*d_e/(8*sigma*n_g);
    Rx=((2*e*n_g*(1-sin(gama/2))/(pi*d_e*cos(gama/2))+1))/cos(beta);#geometrical parameter of the microfin tube
    F2=pow((d_0/d_e),V);
    F3=pow((G_0/G),Z);
    x=TXPm.X;
    Nusselt_cvsmooth=(0.023*pow((G*d_e/mu_l),0.8)*pow(Pr_l,(0.333333333)))*(pow(((1-x)+2.63*x*pow((rho_l/rho_g),0.5)),0.8));
    h_cv=k_l/d_e*Nusselt_cvsmooth*pow(Rx,SS)*pow((Bo*Fr),T)*F2*F3;#convective heat transfer coefficient
 
    h_tp_cavallini=h_cv+h_nb;#superposition form of the flow boiling
    h=h_tp_cavallini;
 
    if(G<100):
        h_cap = 0.332*k_l/e*pow(G*h_fg*sin(beta)/q,0.4326)*(1-pow(G/G_0,3.0));
        h=h+h_cap;
 
    if (abs(T_delta)>=0.1):#if zerotropic refrigerant, correct the overall flow boiling coefficient by considering the mass transfer resistance between vapora phase and liquid phase
        Corr_FlowBoiling=1.0;#parameter for the mass transfer resistance between the vapor phase and liquid phase
        Corr_FlowBoiling = Correct_FLOW_Boiling(TXPm['X'],Cp_g,T_delta,h_fg);#correction parameter for the mass transfer resistance between the liquid phase and vapor phase 
        Re_vaporphase = G*TXPm.X*(d_e+2*e)/(mu_g);#Reynolds number, asuming the vapor only flowing in the tube
        Nu = 0.023*pow(Re_vaporphase,0.8)*pow(Pr_g,0.333);#Dittus-Boelter equation to calculate the vapor phase coeffcient
        h_vaporphase=Nu*k_g/(d_e+2*e);#vapor phase heat transfer coefficent
        h=1.0/(1.0/h+Corr_FlowBoiling/h_vaporphase);#correct the overall flow boiling coefficient

 
    return h, P


def ConvCoeffEvapTP_Smooth(TXPm,#refrigerant state
                           G,#mass flux
                           D,#tube inside diameter
                           q, Ref):#tube wall
    '''
    #B.S.------------------------------------------------------
    /********************************************************************
    TXPm = mean refrigerant thermodynamic state defined by
        temperature (C), quality (-), and pressure (kPa).
    G = refrigerant mass flux (kg/m^2/s)
    D = pipe diameter (m)
    ********************************************************************/
    '''
 
    Ev = FlowPattern()
    Ev['JudgPattern']=0;
    Ev = Eva_FlowPattern(TXPm,G,D,q,Ev, Ref);
    
    h = Ev['h_tp'];

    return h


def Eva_FlowPattern(TXPm,#refrigerant state
                    G,#mass flux
                    d,#tube inside diameter
                    q,#heat flux (overall)
                    Ev,Ref): #for returning all the flow-pattern-dependent parameters
    '''
    /**************************************************
    This is a two-phase evaporating refrigerant heat transfer
    correlation.  
    from Thome-kattern flow-pattern-dependent evaporation model for smooth tubes 
    
    Thome J.R. and Jean Ei Hajal, "On recent advances in modelling of two-phase flow and heat transfer",
    1st International Conference on Heat Transfer, Fluid mechanics, and Thermodynamics, 8-10 April 2002, Kruger Park,
    south Africa TJ1
    
    Thome J.R. and Jean Ei Hajal, "Two-phase flow pattern map for evaporation in horizontal tubes: lastest version",
    1st International Conference on Heat Transfer, Fluid mechanics, and Thermodynamics, 8-10 April 2002, Kruger Park,
    south Africa TJ2
    
    TXPm = mean refrigerant thermodynamic state defined by
        temperature (C), quality (-), and pressure (kPa).
    G = refrigerant mass flux (kg/m^2/s)
    d = pipe diameter (m)
    q = heat flux (W/m^2)
    
    return 
    Ev = struct for returning the flow-pattern-dependent results
    **************************************************/
    '''
                       
    P=TXPm['P'];#corresponding pressure
    x=TXPm['X'];#local quality
    
    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid refrigerant properties
    TXP_prop['P']=TXPm['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    sigma=PropertyTXPtr('I',TXP_prop, Ref) #[N/m]
    rho_l = PropertyTXPth('D',TXP_prop, Ref) #[kg/m^3]
    mu_l = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]
    Cp_l = PropertyTXPtr('C',TXP_prop, Ref) #[J/kg/K]
    k_l = PropertyTXPtr('L',TXP_prop, Ref) #[W/m/K]
    hl = PropertyTXPth('H',TXP_prop, Ref) #[J/kg/K]
    Tsat_l = PropertyTXPth('T',TXP_prop, Ref) #[K]
    
    #vapor refrigerant properties
    TXP_prop['P']=TXPm['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    rho_v = PropertyTXPth('D',TXP_prop, Ref) #[kg/m^3]
    mu_v = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]
    Cp_v = PropertyTXPtr('C',TXP_prop, Ref) #[J/kg/K]
    k_v = PropertyTXPtr('L',TXP_prop, Ref) #[W/m/K]
    hv = PropertyTXPth('H',TXP_prop, Ref) #[J/kg/K]
    Tsat_g = PropertyTXPth('T',TXP_prop, Ref) #[K]
    

    P_cr = PropsSI('PCRIT', Ref)     #critical pressure [Pa]
    T_cr = PropsSI('TCRIT', Ref)     #critical temperature [K]
    M = PropsSI('TCRIT', Ref)*1000   #molecular mass [kg/kmol]

    P_r=P/P_cr;             #reduced pressure
    Pr_v=mu_v/(k_v/Cp_v);   #saturated vapor Prontal number
    Pr_l=mu_l/(k_l/Cp_l);   #saturated liquid Prontal number
    h_LV= hv-hl;            #latent heat

    delta_T = Tsat_g-Tsat_l;    #temperature glide of the refrigerant, used for zeotropic refrigerant

    
    #void fraction
    epsilon=x/rho_v*1.0/((1.0+0.12*(1-x))*(x/rho_v+(1.0-x)/rho_l)+1.18*(1-x)*pow((9.8*sigma*(rho_l-rho_v)),0.25)/(G*pow(rho_l,0.5)));
    Ev['epsilon']=epsilon;#keep the void fraction result

    #cross-sectional area
    A=pi/4.0*d*d;#inside cross-sectional area
    A_L=A*(1.0-epsilon);#liquid cross-sectional area
    A_V=A*epsilon;#vapor cross-sectional area
    A_Ld=A_L/(d*d);#dimensionless liquid cross-sectional area
    A_Vd=A_V/(d*d);#dimensionless vapor cross-sectional area

    
    theta_iter = {'A_Ld':0.0};#struct for iteration to get the stratified angle
    theta_iter['A_Ld']=A_Ld;#known dimensionless liquid cross-sectional area
    
    theta_max=2.0*pi;
    theta_min=0;
    theta_strat = brentq(theta_ev,theta_max,theta_min,args=(theta_iter),xtol=1e-7,rtol=6e-8,maxiter=40) #xtol, rtol and maxiter are changed to match "Zbrent" solver in ACMODEL
    
    
    h_Ld=0.5*(1.0-cos((2.0*pi-theta_strat)/2.0));#dimensionless liquid height
    P_id=sin((2.0*pi-theta_strat)/2.0);#dimensionless interface length between the vapor phase and liquid phase
    WeFL=9.8*d*d*rho_l/sigma;#ratio of liquid weber number to the liquid Froude number
    Xi=pow((1.138+2.0*log10(pi/(1.5*A_Ld))),(-2.0));#internal friction factor
    q_crit=0.131*pow(rho_v,0.5)*h_LV*pow((9.8*(rho_l-rho_v)*sigma),0.25);#critial heat flux
    Fq_1=646.0*pow((q/q_crit),2.0)+64.8*(q/q_crit);#parameter for considering the effect of heat flux
    Fq_2=18.8*(q/q_crit)+1.023;#parameter for considering the effect of heat flux

    Ev['G_wavy']=pow(((16.0*pow(A_Vd,3.0)*9.8*d*rho_l*rho_v)/(pow(x,2.0)*pow(pi,2.0)*pow((1.0-pow((2.0*h_Ld-1.0),2.0)),0.5))*(pow(pi,2.0)/(25.0*pow(h_Ld,2.0))*pow((1.0-x),(-1.0*Fq_1))*pow((WeFL),(-1.0*Fq_2))+1.0)),0.5)+50-75*exp(-pow((pow(x,2.0)-0.97),2.0)/(x*(1-x)));#transition curve between stratifed-wavy flow and annular flow

    Ev['G_strat']=pow(((pow(226.3,2.0)*A_Ld*pow(A_Vd,2.0)*rho_v*(rho_l-rho_v)*mu_l*9.8)/(pow(x,2.0)*(1.0-x)*pow(pi,3.0))),(1.0/3.0))+20.0*x;#transition curve between stratifed flow and stratified-wavy flow
    Ev['X_lA']=1.0/((0.2914*pow((rho_v/rho_l),(-1.0/1.75))*pow((mu_l/mu_v),(-1.0/7.0)))+1.0);#transition curve between the intermittent flow and the annular flow

    Ev['G_mist']=pow(((7680.0*pow(A_Vd,2.0)*9.8*d*rho_l*rho_v)/(pow(pi,2.0)*pow(x,2.0)*Xi)*1.0/WeFL),0.5);#transition curve between the annular flow and mist flow
    Ev['G_bub']=pow(((256.0*A_Vd*pow(A_Ld,2.0)*pow(d,1.25)*rho_l*(rho_l-rho_v)*9.8)/(0.3164*pow((1.0-x),1.75)*pow(pi,2.0)*P_id*pow(mu_l,0.25))),(1.0/1.75));#transition curve between the annular flow and bubbly flow
    
    if (G>Ev['G_mist']): #mist flow, dry angle is zero
        Ev['Pattern']=4;
        Ev['theta_dry']=0;
    elif (G<Ev['G_strat']):#stratified flow, dry angle is the stratified angle
        Ev['Pattern']=0;
        Ev['theta_dry']=theta_strat;
    elif (G<Ev['G_wavy'] and G>Ev['G_strat']): #stratified wavy flow, dry angle is the interpolation between the annular flow and the stratified flow
        Ev['Pattern']=1;
        G_high=Ev['G_wavy'];#high mass flux is from annular flow transition curve
        G_low=Ev['G_strat'];#low mass flux is from stratified flow transition curve.
        Ev['theta_dry']=theta_strat*(G_high-G)/(G_high-G_low);
    elif(G>Ev['G_wavy'] and x<Ev['X_lA']): #intermittent flow, dry angle is zero
        Ev['Pattern']=2;
        Ev['theta_dry']=0;
#     elif(G>Ev['G_wavy'] and G<Ev['G_bub'] and x<Ev['X_lA']): #intermittent flow, dry angle is zero
#     Ev['Pattern']=2;
#     Ev['theta_dry']=0;
    elif (G>Ev['G_wavy'] and G<Ev['G_mist'] and x>Ev['X_lA']): #annular flow, dry angle is zero
        Ev['Pattern']= 3;
        Ev['theta_dry']=0;
    else: #all others are supposed to be annular flow 
        Ev['Pattern']= 3;
        Ev['theta_dry']=0;

    if (Ev['JudgPattern']):
        return Ev

    delta=pi*d*(1.0-epsilon)/(2.0*(2.0*pi-Ev['theta_dry']));#equivalent liquid thickness
    Re_v=G*x*d/(epsilon*mu_v);#vapor Reynolds number, based on the void fraction, for calculating flow boiling coefficent
    Re_l=4.0*G*(1.0-x)*delta/((1.0-epsilon)*mu_l);#liquid phase Reynolds number, based on the liquid film thickness
    Ev['delta']=delta;
    Ev['h_nb']=55.0*pow(P_r,0.12)*pow((-log10(P_r)),(-0.55))*pow(M,(-0.5))*pow(q,0.67);#nucleate boiling coefficient from Copper equition
    
    if (abs(delta_T)>0.1): #if zerotropic refrigerant, correct the nucleate boiling coefficent
        Corr_NUC=1.0;#parameter for the effect of mass transfer resistance on nucleate boiling
        B_0=1.0;#scaling factor
        beta_l = 3e-4;#m/s, mass transfer coefficient
        Corr_NUC = Correct_NUC_Boiling(Ev['h_nb'],q,delta_T,B_0,beta_l,rho_l,h_LV);#nucleate boiling correction
        Ev['h_nb']=Ev['h_nb']*Corr_NUC;#corrected
    
    Ev['h_v']=0.023*pow(Re_v,0.8)*pow(Pr_v,0.4)*k_v/d;#gas phase heat transfer coeffcient
    Ev['h_cb']=0.0133*pow(Re_l,0.69)*pow(Pr_l,0.4)*k_l/delta;#convective boiling coefficient
    Ev['h_wet']=pow((pow(Ev['h_nb'],3.0)+pow(Ev['h_cb'],3.0)),(1.0/3.0));#asympotic equation that combines the nucleate boiling coeffcient and the convective boiling coefficent
    Ev['h_tp']=(Ev['theta_dry']*Ev['h_v']+(2.0*pi-Ev['theta_dry'])*Ev['h_wet'])/(2.0*pi);#overall two-phase flow boiling coefficient
    
    if (abs(delta_T)>0.1): #if zerotropic refrigerant, correct the overall flow boiling coefficient by considering the mass transfer resistance between vapora phase and liquid phase
        Corr_FlowBoiling=1.0;#parameter for the mass transfer resistance between the vapor phase and liquid phase
        Corr_FlowBoiling = Correct_FLOW_Boiling(TXPm['X'],Cp_v,delta_T,h_LV);#correction parameter for the mass transfer resistance between the liquid phase and vapor phase 
        Re_vaporphase = G*x*d/(mu_v);#Reynolds number, asuming the vapor only flowing in the tube
        Nu = 0.023*pow(Re_vaporphase,0.8)*pow(Pr_v,0.333);#Dittus-Boelter equation to calculate the vapor phase coeffcient
        h_vaporphase=Nu*k_v/d;#vapor phase heat transfer coefficent
        Ev['h_tp']=1.0/(1.0/Ev['h_tp']+Corr_FlowBoiling/h_vaporphase);#correct the overall flow boiling coefficient


    return Ev


def Cond_FlowPattern(TXPm,#refrigerant state
                     G,#mass flux
                     d,#tube inside diameter
                     q,#heat flux (overall)
                     Cd, Ref):#for returning all the flow-pattern-dependent parameters
    '''
    #not totally finished, need iterating by tube wall temperature
    '''

    P=TXPm['P'];#corresponding pressure
    x=TXPm['X'];#local quality
    
    TXP_prop={'T':0,'X':0,'P':0}
    
    #liquid refrigerant properties
    TXP_prop['P']=TXPm['P'];
    TXP_prop['X']=0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    sigma=PropertyTXPtr('I',TXP_prop, Ref) #[N/m]
    rho_l = PropertyTXPth('D',TXP_prop, Ref) #[kg/m^3]
    mu_l = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]
    Cp_l = PropertyTXPtr('C',TXP_prop, Ref) #[J/kg/K]
    k_l = PropertyTXPtr('L',TXP_prop, Ref) #[W/m/K]
    hl = PropertyTXPth('H',TXP_prop, Ref) #[J/kg/K]
    Tsat_l = PropertyTXPth('T',TXP_prop, Ref) #[K]
    
    #vapor refrigerant properties
    TXP_prop['P']=TXPm['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    rho_v = PropertyTXPth('D',TXP_prop, Ref) #[kg/m^3]
    mu_v = PropertyTXPtr('V',TXP_prop, Ref) #[Pa-s]
    Cp_v = PropertyTXPtr('C',TXP_prop, Ref) #[J/kg/K]
    k_v = PropertyTXPtr('L',TXP_prop, Ref) #[W/m/K]
    hv = PropertyTXPth('H',TXP_prop, Ref) #[J/kg/K]
    Tsat_g = PropertyTXPth('T',TXP_prop, Ref) #[K]
    
    
    P_cr = PropsSI('PCRIT', Ref)     #critical pressure [Pa]
    T_cr = PropsSI('TCRIT', Ref)     #critical temperature [K]
    M = PropsSI('TCRIT', Ref)*1000   #molecular mass [kg/kmol]

    P_r=P/P_cr;#reduced pressure
    Pr_v=mu_v/(k_v/Cp_v);#saturated vapor Prontal number
    Pr_l=mu_l/(k_l/Cp_l);#saturated liquid Prontal number
    h_LV= hv-hl;#latent heat

    delta_T = Tsat_g-Tsat_l;#reftpvthP.Tsat(TXPm.P)-reftplthP.Tsat(TXPm.P);#temperature glide of the refrigerant, used for zeotropic refrigerant
    
    A=pi/4.0*d*d;#inside cross-sectional area

    #determine the minimum mass flow rate and quality for wavy flow
    x_min=0; 
    G_min=1000000;
    P_static=0;

    if (abs(P_static-TXPm['P'])<0.1):
        ini=0;
    else:
        ini=1;
        P_static=TXPm['P'];

    if (ini):
        for i in range(100):
            x1=0.01*i;
            epsilon_h1=pow((1+(1-x1)/x1*(rho_v/rho_l)),(-1.0));
            epsilon_ra1=x1/rho_v*pow(((1+0.12*(1-x1))*(x1/rho_v+(1-x1)/rho_l)+1.18*(1-x1)*pow((9.8*sigma*(rho_l-rho_v)),0.25)/(G*pow(rho_l,0.5))),(-1.0));
            epsilon1=(epsilon_h1-epsilon_ra1)/log(epsilon_h1/epsilon_ra1);
            A_L1=A*(1-epsilon1);
            A_V1=A*epsilon1;
            A_Ld1=A_L1/pow(d,2.0);
            A_Vd1=A_V1/pow(d,2.0);
            theta_strat1=2.0*pi-2.0*(pi*(1.0-epsilon1)+pow((3.0/2.0*pi),(1.0/3.0))*(1-2.0*(1-epsilon1)+pow((1-epsilon1),(1.0/3.0))-pow(epsilon1,(1.0/3.0)))-1.0/200.0*(1-epsilon1)*epsilon1*(1-2*(1-epsilon1))*(1+4*(pow((1.0-epsilon1),2.0)+pow(epsilon1,2.0))));
            h_Ld1=0.5*(1.0-cos((2*pi-theta_strat1)/2.0));
            G_bak =pow(((16*pow(A_Vd1,3.0)*9.8*d*rho_l*rho_v)/(pow(x1,2.0)*pow(pi,2.0)*pow((1-pow((2*h_Ld1-1),2.0)),0.5))*(pow(pi,2.0)/(25*pow(h_Ld1,2.0))*pow((9.8*pow(d,2.0)*rho_l/sigma),(-1.023))+1)),0.5)+50-75.0*exp(-pow((pow(x1,2.0)-0.97),2.0)/(x1*(1-x1)));
            if (G_min>G_bak):
                G_min=G_bak;
                x_min=x1;

    #void fraction
    epsilon_h=pow((1+(1-x)/x*(rho_v/rho_l)),(-1.0));
    epsilon_ra=x/rho_v*pow(((1+0.12*(1-x))*(x/rho_v+(1-x)/rho_l)+1.18*(1-x)*pow((9.8*sigma*(rho_l-rho_v)),0.25)/(G*pow(rho_l,0.5))),(-1.0));
    Cd['epsilon']=(epsilon_h-epsilon_ra)/log(epsilon_h/epsilon_ra);
    #cross-sectional area
    A_L=A*(1-Cd['epsilon']);
    A_V=A*Cd['epsilon'];
    A_Ld=A_L/pow(d,2.0);
    A_Vd=A_V/pow(d,2.0);
    theta_strat=2.0*pi-2.0*(pi*(1.0-Cd['epsilon'])+pow((3.0/2.0*pi),(1.0/3.0))*(1-2.0*(1-Cd['epsilon'])+pow((1-Cd['epsilon']),(1.0/3.0))-pow(Cd['epsilon'],(1.0/3.0)))-1.0/200.0*(1-Cd['epsilon'])*Cd['epsilon']*(1-2*(1-Cd['epsilon']))*(1+4*(pow((1.0-Cd['epsilon']),2.0)+pow(Cd['epsilon'],2.0))));
    h_Ld=0.5*(1.0-cos((2*pi-theta_strat)/2.0));
    P_id=sin((2.0*pi-theta_strat)/2.0);

    #Flow pattern determination
    if (x>x_min):
        Cd['G_wavy']=G_min;
    else:
        Cd['G_wavy'] = pow(((16*pow(A_Vd,3)*9.8*d*rho_l*rho_v)/(pow(x,2.0)*pow(pi,2.0)*pow((1-pow((2*h_Ld-1),2.0)),0.5))*(pow(pi,2.0)/(25*pow(h_Ld,2.0))*pow((9.8*pow(d,2.0)*rho_l/sigma),(-1.023))+1)),0.5)+50-75.0*exp(-pow((pow(x,2.0)-0.97),2.0)/(x*(1-x)));


    Cd['G_strat']=pow(((pow(226.3,2.0)*A_Ld*pow(A_Vd,2.0)*rho_v*(rho_l-rho_v)*mu_l*9.8)/(pow(x,2.0)*(1-x)*pow(pi,3.0))),(1.0/3.0))+20*x;
    Cd['X_lA']=pow(((0.2914*pow((rho_v/rho_l),(-1.0/1.75))*pow((mu_l/mu_v),(-1.0/7.0)))+1),(-1.0));
    Xi=pow((1.138+2*log10(pi/(1.5*A_Ld))),(-2.0));
    WeFL=9.8*pow(d,2.0)*rho_l/sigma;
    Cd['G_mist']=pow(((7680*pow(A_Vd,2.0)*9.8*d*rho_l*rho_v)/(pow(pi,2.0)*pow(x,2.0)*Xi)*1/WeFL),0.5);
    Cd['G_bub']=pow(((256*A_Vd*pow(A_Ld,2.0)*pow(d,1.25)*rho_l*(rho_l-rho_v)*9.8)/(0.3164*pow((1-x),1.75)*pow(pi,2.0)*P_id*pow(mu_l,0.25))),(1.0/1.75));

    u_l=G*(1-x)/(rho_l*(1-Cd['epsilon']));
    u_v=G*x/(rho_v*Cd['epsilon']);

    if (G>Cd['G_mist']): #mist flow, dry angle is zero
        Cd['Pattern']=4;
        Cd['theta_dry']=0;
    elif (G<Cd['G_strat']):#stratified flow, dry angle is the stratified angle
        Cd['Pattern']=0;
        Cd['theta_dry']=theta_strat;
    elif (G<Cd['G_wavy'] and G>Cd['G_strat']): #stratified wavy flow, dry angle is the interpolation between the annular flow and the stratified flow
        Cd['Pattern']=1;
        Cd['theta_dry']= theta_strat*pow(((Cd['G_wavy']-G)/(Cd['G_wavy']-Cd['G_strat'])),0.5);
    elif (G>Cd['G_wavy'] and x<Cd['X_lA']):#intermittent flow, dry angle is zero
        Cd['Pattern']=2;
        Cd['theta_dry']=0;
#     elif (G>Cd['G_wavy'] and G<Cd['G_bub'] and x<Cd['X_lA']): #intermittent flow, dry angle is zero
#     Cd['Pattern']=2;
#     Cd['theta_dry']=0;
    elif (G>Cd['G_wavy'] and x>Cd['X_lA']): #annular flow, dry angle is zero
        Cd['Pattern']= 3;
        Cd['theta_dry']=0;
    else: #all others are supposed to be annular flow 
        Cd['Pattern']= 3;
        Cd['theta_dry']=0;

    if (Cd['JudgPattern']):
        return Cd

    Cd['delta']=(d-pow((pow(d,2.0)-A_L*8/(2*pi-Cd['theta_dry'])),0.5))/2.0;
    
    if (G<Cd['G_strat']):
        f_i=1+pow((u_v/u_l),0.5)*pow(((rho_l-rho_v)*9.8*pow(Cd['delta'],2.0)/sigma),0.25)*(G/Cd['G_strat']);
    else:
        f_i=1+pow((u_v/u_l),0.5)*pow(((rho_l-rho_v)*9.8*pow(Cd['delta'],2.0)/sigma),0.25);

    
    Re_l=4.0*G*(1-x)*Cd['delta']/(mu_l*(1-Cd['epsilon']));
    h_c=0.003*pow(Re_l,0.74)*pow(Pr_l,0.5)*k_l/Cd['delta']*f_i;
    #T_w=TXPm['T']-2;#need to input by iteration later
    #h_f=0.728*pow((rho_l*(rho_l-rho_v)*9.8*h_LV*pow(k_l,3.0)/(mu_l*d*(TXPm['T']-T_w))),0.25);
    h_f=0.655*pow((rho_l*(rho_l-rho_v)*9.8*h_LV*pow(k_l,3.0)/(mu_l*d*q)),1.0/3.0);
    Cd['h_tp']=(h_f*Cd['theta_dry']+(2*pi-Cd['theta_dry'])*h_c)/(2*pi);
    
    return Cd


def theta_ev(theta_strat, Params=None):
    '''
    #this function is for getting dry circumferencial angle
    /***********************************
    this subroutine work together with zbrent algorithm to calculate the stratifed angle by iteration
    *************************************/
    '''
    if (Params == None):
        theta_strat_ev = {'A_Ld': 0.0}#dictionary variable for getting the stratified angle in Kattan-Thome flow map by iteration
        theta = theta_strat_ev
    else:
        theta = Params 
        
    A_Ld1=1.0/8.0*((2.0*pi-theta_strat)-sin(2.0*pi-theta_strat));#dimensionless liquid cross-sectional area at the given stratified angle
    residual = A_Ld1-theta['A_Ld'];#deviation between the dimensionless liquid cross-sectional area at the given stratified angle and the real dimensionless liquid cross-sectional area
    
    return residual



def Correct_NUC_Boiling(h_nb,#nucleate boiling coefficient
                        q, #flow boiling heat flux
                        delta_T,#temperture glide of the zeotropic refrigerant
                        B_0,#scaling factor
                        beta_l,#mass transfer coefficient
                        rho_l,#liquid density
                        h_fg):#latent heat of the zeotropic refrigerant
    '''
    /************************************************
    For zerotropic refrigerant, accounting for the mass transfer resistance on nucleate boiling
    Kattan, N. Thome, J R. Favrat, D. Flow boiling in horizontal tubes: Part 3 - development of a new heat
    transfer model based on flow pattern Journal of Heat Transfer-Transactions of the ASME. v 120 n 1 Feb
    1998. p 156-165
    ***************************************************/
    '''

    F_c=1/(1+h_nb/q*delta_T*(1-exp(-1*B_0*q/(rho_l*h_fg*beta_l))));#correction parameter
    
    return F_c


def Correct_FLOW_Boiling(x,#quality
                         Cp_g,#specific heat of gas phase
                         delta_T,#temperature glide of the zeotropic refrigerant
                         h_fg):#latent heat of the zeotropic refrigerant
    '''
    /**************************************************************
    For zerotropic refrigerant, accounting for the mass transfer resistance between the liquid phase and vapor phase
    Bell and Ghaly correction for mixture thermal resistance, direct reference
    Cavallini A.. Del Col D., Doretti L. Longo G. A. Rossetto L. 
    "Refrigerant vaporization inside enhanced tubes, a heat transfer model" 
    Heat and technology, Vol. 17. n.2 1999
    **************************************************************/
    '''
    
    Corr=x*Cp_g*delta_T/h_fg;
    
    return Corr


#===============================================================================
# frictional pressure drop of refirgerant
#===============================================================================

def FricDP(TXPi, #refrigerant state
           Gr, #refrigerant mass flux
           q,#heat flux
           P, Ref):#evaporator struct
    '''
    #the following function is the interface to call the pressure drop calculation in the EVAPORATOR and CONDENSER
    '''
 
    D = InsideTube_Dim();
    
    TXP1={'T':0,'X':0,'P':0};
    TXP2={'T':0,'X':0,'P':0};
    
    X1=0.1; X2=0.95;
 
    D['Microfin'] = P['microfin']; #microfin type, 0=smooth tube, 1=helical, 2=cross-grooved, 3=herringbone
    D['finN'] = P['finN']; #fin number in a micro-fin tube
    D['Di'] = P['Di'];#inside diamete at the fin tip
    D['gama'] = P['gama'] ;#fin apex angle in a micro-fin tube
    D['beta'] = P['beta'] ;    #fin helix angle in a micro-fin tube
    D['finH'] = P['finH']; #fin height in a micro-fin tube
    D['w_b'] = P['w_b'] ; #base width of a single fin
    D['w_e'] = P['w_e']; #top width of a single fin
    D['w_z'] = P['w_z']; #base distance between two neighboring fins
    D['K_T'] = P['K_T'];#400, this is the conductance factor of copper
    D['Ls'] = P['Ls'];#tube unit length
    D['D_b'] = P['D_b'];#tube diameter at the base of the fin
    D['Do'] = P['Do']; #Pipe outside diameter.
    D['D_m'] = P['D_m']; #mean diameter of the micro-fin tube
    D['P_H'] = P['P_H'];# the hydraulical circumference
    D['Acs'] = P['Acs'] ;#cross area of the micro-fin tube, this is the actual cross-section area
    D['Dh_i'] = P['Dh_i'];#inside hydraulical diameter 
    D['Ax'] = P['Ax'];# Inside pipe cross sectional area, based on fin tips
    D['Api'] = P['Api'];# Inside pipe surface area (one tube segment), based on fin tips
     
    if(TXPi['X']>=1.0 or TXPi['X']<=0):
        y, D = FricDPSP_Microfin(TXPi,Gr,D, Ref);
    elif(TXPi['X']<X1):
        TXP1=toTXP(TXPi['T'],0,TXPi['P']);
        TXP2=toTXP(TXPi['T'],X1,TXPi['P']);
        y1, D = FricDPSP_Microfin(TXP1,Gr,D, Ref);
        y2, D = FricDPTP_Microfin(TXP2,Gr,q,D, Ref);
        y=y1+TXPi['X']*(y2-y1)/X1;
    elif(TXPi['X']>X2):
        TXP1=toTXP(TXPi['T'],X2,TXPi['P']);
        TXP2=toTXP(TXPi['T'],1,TXPi['P']);
        y1, D = FricDPTP_Microfin(TXP1,Gr,q,D, Ref);
        y2, D = FricDPSP_Microfin(TXP2,Gr,D, Ref);
        y=y2-(1-TXPi['X'])*(y2-y1)/(1-X2);
    else:
        y, D = FricDPTP_Microfin(TXPi,Gr,q,D, Ref);
     
    return y, P

#The following function is the same as the previous, therefroe it is commented out
# def FricDP(TXP TXPi, #refrigerant state
#          double Gr, #refrigerant mass flux
#          double q,#heat flux
#          CGP * P, Ref):#condenser struct
#     '''
#     #the following function is the interface to call the pressure drop calculation in the condenser
#     '''
#  
#     InsideTube_Dim D;
#     TXP TXP1,TXP2;
#     double y1,y2,y;
#     const double X1=0.1,X2=0.95;
#  
#     D.Microfin = P->Microfin; #microfin type, 0=smooth tube, 1=helical, 2=cross-grooved, 3=herringbone
#     D.finN = P->finN; #fin number in a micro-fin tube
#     D.Di = P->Di;#inside dimater at the fin tip
#     D.gama =P->gama ;#fin apex angle in a micro-fin tube
#     D.beta = P->beta ;    #fin helix angle in a micro-fin tube
#     D.finH = P->finH; #fin height in a micro-fin tube
#     D.w_b = P->w_b ; #base width of a single fin
#     D.w_e = P->w_e; #top width of a single fin
#     D.w_z = P->w_z; #base distance between two neighboring fins
#     D.K_T = P->K_T;#400, this is the conductance factor of copper
#     D.Ls = P->Ls;#tube unit length
#     D.D_b = P->D_b;#tube diameter at the base of the fin
#     D.Do = P->Do; #Pipe outside diameter.
#     D.D_m = P->D_m; #mean diameter of the micro-fin tube
#     D.P_H = P->P_H;# the hydraulical circumference
#     D.Acs = P->Acs ;#cross area of the micro-fin tube, this is the actual cross-section area
#     D.Dh_i = P->Dh_i;#inside hydraulical diameter 
#     D.Ax = P->Ax;# Inside pipe cross sectional area, based on fin tips
#     D.Api =P->Api;# Inside pipe surface area (one tube segment), based on fin tips
#      
#     if(TXPi.X>=1.0 || TXPi.X<=0) {
#         y=FricDPSP_Microfin(TXPi,Gr, &D);
#     } else if(TXPi.X<X1) {
#         TXP1=toTXP(TXPi.T,0,TXPi.P);
#         TXP2=toTXP(TXPi.T,X1,TXPi.P);
#         y1=FricDPSP_Microfin(TXP1,Gr, &D);
#         y2=FricDPTP_Microfin(TXP2,Gr, q, &D);
#         y=y1+TXPi.X*(y2-y1)/X1;
#     } else if(TXPi.X>X2) {
#         TXP1=toTXP(TXPi.T,X2,TXPi.P);
#         TXP2=toTXP(TXPi.T,1,TXPi.P);
#         y1=FricDPTP_Microfin(TXP1,Gr, q, &D);
#         y2=FricDPSP_Microfin(TXP2,Gr, &D);
#         y=y2-(1-TXPi.X)*(y2-y1)/(1-X2);
#     } else {
#         y=FricDPTP_Microfin(TXPi,Gr, q, &D);
#     }
#      
#     return y


def FricDPTP_Microfin(TXPi,#refrigerant state
                      Gr,#mass flux
                      q, #heat flux
                      P, Ref):#microfin tube struct
    '''
    /******************************************
    Friction pressure drop two phase, B.S.
    Kedzierski, M. A., and Choi J. Y., "A generalized pressure drop correlations for 
    evaporation and condensation of alternative refrigerants in smooth and micro-fin 
    tubes" NISTIR 6333, 1999
    ************************************************/
    '''
 
    L_T = P['Ls']; #segment length
    D_H = P['Dh_i'];#hydraulic diameter
    TXP_prop={'T':0,'X':0,'P':0};
 
    #liquid refrigerant properties
    TXP_prop['P']=TXPi['P'];
    TXP_prop['X']=0.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref); #[K]
    DL=PropertyTXPth('D',TXP_prop, Ref); #[kg/m^3] 
    MU_F=PropertyTXPtr('V',TXP_prop, Ref); #[Pa-s]
 
    #vapor refrigerant properties
    TXP_prop['P']=TXPi['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref); #[K]
    DV=PropertyTXPth('D',TXP_prop, Ref); #[kg/m^3] 
    
    G_A=9.81e0;
    KF=1*abs(q)/(L_T*G_A);#this is two-phase flow number, actually, q= (x_in-x_out)*h_fg
     
    RE_D=Gr*D_H/MU_F;#Reynolds number
    F=0.00506e0/pow(RE_D,0.0951e0)*pow(KF,0.1554e0);#two-phase friction factor
    SV=1e0/DL+TXPi['X']*(1e0/DV-1e0/DL);#two-phase specific volume
    DP_FR=SV*pow(Gr,2e0)*2e0*F*L_T/(D_H) #[Pa]
    
    return DP_FR, P


#===============================================================================
# single-phase pressure drop
#===============================================================================
def FricDPSP_Microfin(TXP_loc,#refrigerant state
                      Gr, #mass flux
                      P, Ref):#micro-fin geometry
    '''
    /************************************************************
    single-phase pressure drop in finned tube,
    Haaland, S.S., 1983, "simple and explicit formulas for the friction factor in turbulent
    pipe flow", Journal of fluids engineering, Vol., 105, pp. 89-90
    *************************************************************/
    '''
    
    L_T = P['Ls']; #segment length
    L_IF = P['finH'];#fin height
    R_I = P['D_b']/2.0;#inside radius at the fin bottom
    
    TXP_prop={'T':0,'X':0,'P':0};
 
    if (P['Microfin']<1):#smooth tube
        DP_Smooth = FricDPSP_Smooth(TXP_loc,Gr,2*R_I,L_T, Ref)
        return DP_Smooth, P
     
    #liquid refrigerant properties
    TXP_prop['P']=TXP_loc['P'];
    TXP_prop['X']=0.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    DL=PropertyTXPth('D',TXP_prop, Ref) #[kg/m^3]
 
    #vapor refrigerant properties
    TXP_prop['P']=TXP_loc['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    DV=PropertyTXPth('D',TXP_prop, Ref) #[kg/m^3]
     
    if(TXP_loc['X']>=0.98):
        TXP_loc['X']=1;#to gurantee the single-phase state
    if(TXP_loc['X']<0.1):
        TXP_loc['X']=0;
     
    MU = PropertyTXPtr('V',TXP_loc, Ref) #[Pa-s]

    if(TXP_loc['X']<0.95):
        Density=DL;#to gurantee the single-phase state
    else:
        Density=DV;
     
    RE=Gr*2e0*R_I/MU;#Reunolds number
    F=1e0/pow((-1.8e0*log10(6.9e0/RE+pow((L_IF/(3.7e0*2*R_I)),1.11e0))),2e0)/4e0;#friction factor
 
    DP_SP=pow(Gr,2e0)*F*L_T/R_I/Density #[Pa]
     
    return DP_SP, P


def FricDPSP_Smooth(TXP_loc,#refrigerant state
                    G,#mass flux 
                    D, #inside diameter
                    L, Ref):#segment length
    '''
    #B.S. single-phase pressure drop in smooth tube, this is a common turbulent flow model
    '''
 
    TXP_prop={'T':0,'X':0,'P':0};
    
    #liquid refrigerant properties
    TXP_prop['P']=TXP_loc['P'];
    TXP_prop['X']=0.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    DL=PropertyTXPth('D',TXP_prop, Ref) #[kg/m^3]
 
    #vapor refrigerant properties
    TXP_prop['P']=TXP_loc['P'];
    TXP_prop['X']=1.0;
    TXP_prop['T']=PropertyTXPth('T',TXP_prop, Ref) #[K]
    DV=PropertyTXPth('D',TXP_prop, Ref) #[kg/m^3]
    
    if(TXP_loc['X']>=0.98):
        TXP_loc['X']=1;#to gurantee the single-phase state
    if(TXP_loc['X']<0.1):
        TXP_loc['X']=0;
     
    MU = PropertyTXPtr('V',TXP_loc, Ref) #[Pa-s]

    if(TXP_loc['X']<0.95):
        Density=DL;#to gurantee the single-phase state
    else:
        Density=DV;
     
    RE=G*D/MU;#Reynolds number
    F=0.046e0*pow(RE,(-0.2e0));#turbulent friction factor
    DP=pow(G,2e0)*F*L/(D/2e0)/Density #[Pa]
     
    return DP


def FinEffect_Schmidt(h,#airside heat transfer coefficient
                      k,#fin conductance
                      th,#fin thickness
                      w,#length beween the fin tip and fin base
                      Do):#tube outside diameter
    '''#B.S.-------------------------------------------------------------------------
    #calculating the fin efficiency with Schmidt method
    '''
 
    #Schmidt calculation method
    R_o = Do/2;#tube outside radius
    R_ref=R_o+w;#reference fin radius
 
    beta=R_ref/R_o;
    arg = 2*h/(k*th);
    if(arg<0.0):
        print("CORR::FinEffect_Schmidt","h/(k*yb)<0.0");

    m=pow(arg,0.5);
 
    phi=(beta-1)*(1+0.35*log(beta));
    MRPHI=m*R_o*phi;#this is m*L_F, L_F=R_o*Phi is the equivalent fin length
    eta=tanh(MRPHI)/(MRPHI);
 
    return eta


# def Coef_Hilpert(double G,double Tair,double D):
#     '''
#     /************************************************************
#     convection heat transfer of air flowing across the circular tube
#     *************************************************************/
#     '''
# 
#     double ho=0;
#     const double mu=air.mu(Tair);
#     if(errorLog.IsError()) {
#         errorLog.Add("Coef_Hilpert","mu");
#         return 0;
#     }
# 
#     const double k=air.k(Tair);
#     if(errorLog.IsError()) {
#         errorLog.Add("Coef_Hilpert","k");
#         return 0;
#     }
# 
#     const double Cp=air.Cp(Tair);
#     if(errorLog.IsError()) {
#         errorLog.Add("Coef_Hilpert","Cp");
#         return 0;
#     }
# 
#     const double Pr=Cp*mu/k;
# 
#     const double Re=G*D/mu;
#     double c=0,m=0;
# 
#     if(Re<4)
#     {c=0.989; m=0.33;}
#     else if(Re<40&&Re>=4)
#     {c=0.911; m=0.385;}
#     else if(Re<4000&&Re>=40)
#     {c=0.683; m=0.466;}
#     else if(Re<40000&&Re>=4000)
#     {c=0.193; m=0.618;}
#     else
#     {c=0.027; m=0.805;}
# 
#     const double Nu=c*pow(Re,m)*pow(Pr,1.0/3.0);
# 
#     ho=Nu*k/D;
# 
# 
#     return ho


# def Exercise_Corr():#excise correlation subroutine
# 
#     TXPm = {'T':0.0,'X':0.0,'P':0.0};
#     D=9e-3;
#     G=200;
#     q=15000;
#     T_w=10.49;
#     TXPm['P']=650;
# 
# #    ETdim P;
#     CGP P;
#     P['Microfin']=1;
#     P['beta']=30;
#     P['gama']=45;
#     P['finH']=0.2e-3;
#     P['finN']=60;
#     P['Di'] = 9e-3;
# 
#     std::ofstream outstrm; 
#     outstrm.open("out_corr.xls",std::ios::app); # for outputing the data of all the state points
#     for(int i=1;i<=49;i++)
#     {
#     TXPm.X=0.95+i*0.02;
#     
#     outstrm<<TXPm.X<<"    "<<ConvCoeffEvapTP_Smooth(TXPm,G,D,q)<<std::endl;
# 
# #    outstrm<<TXPm.X<<"    "<<ConvCoeffEvapTP_microfin(TXPm,G,T_w, &P)<<std::endl;
# 
# #    outstrm<<TXPm.X<<"    "<<ConvCoeffCondTP_Microfin(TXPm,G,&P)<<std::endl;
# 
#     }
#     outstrm.close();
# 
# 
#     return 0



#B.S. the following have not been unused
def ConvCoeffEvapInside(TXPi,G,D,T_w,Ref):#ConvCoeffEvapInside(TXP TXPi,double G,double D,double dx,double L)
    '''
    /********************************************************************
    Inside convection coefficient for the refrigerant inside the  #shenbo change
    evaporator.
    ********************************************************************/
    '''
#     double y;
#     const double X1=0.0000001,X2=0.999999;#X2=0.9 0.15
# 
#     if(TXPi.X>=1.0 || TXPi.X<=0) {
#     #    y = ConvCoeffSP(TXPi,G,D);
#     } else if(TXPi.X<X1) {
#         TXP TXP1 = toTXP(TXPi.T,0,TXPi.P);
#         TXP TXP2 = toTXP(TXPi.T,X1,TXPi.P);
# #        double y1 = ConvCoeffSP(TXP1,G,P);
#         if(errorLog.IsError()) {
#             errorLog.Add("ConvCoeffEvapInside");
#             return -1;
#         }
#         double y2 = ConvCoeffEvapTP_Smooth(TXP2,G,D,T_w);#shenbo change  ConvCoeffEvapTP(TXP2,G,D,dx,L);
# #        y = y1+TXPi.X*(y2-y1)/X1;
#     } else if(TXPi.X>X2) {
#         TXP TXP1 = toTXP(TXPi.T,X2,TXPi.P);
#         TXP TXP2 = toTXP(TXPi.T,1,TXPi.P);
#         double y1 = ConvCoeffEvapTP_Smooth(TXP1,G,D,T_w);# shenbo change ConvCoeffEvapTP(TXP1,G,D,dx,L);
#         if(errorLog.IsError()) {
#             errorLog.Add("ConvCoeffEvapInside");
#             return -1;
#         }
# #        double y2 = ConvCoeffSP(TXP2,G,D);
# #        y = y2-(1-TXPi.X)*(y2-y1)/(1-X2);
#     } else {
#         y = ConvCoeffEvapTP_Smooth(TXPi,G,D,T_w);#shenbo change ConvCoeffEvapTP(TXPi,G,D,dx,L);
#     }
# 
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffEvapInside");
#         return -1;
#     }
# 
#     return y;
    return 0;



def IsTwoPhase(X):
    '''
    /********************************************************************
    Tests quality (X) and returns TRUE if it is two phase and
    FALSE otherwise.
    ********************************************************************/
    '''
    
    return 1 if (X-0.5)*(X-0.5)<0.25 else 0


# def ConvCoeffAir(double T,double G,double D,double th,double /*y*/,double z):
#     '''
#     /********************************************************************
#     Air side convection coefficient from "Enhanced fins for air-cooled
#     heat exchangers - heat transfer and friction factor correlations",
#     W. Nakayama, L.P. Xu, 1986 ASME-JSME Thermal Engineering Joint
#     Conference Proceedings, Honolulu, Hawaii, March 20-24, 1983 p.495-
#     p.502.  This correlation is for slotted fins in one or more rows.  
#     T = temperature of entering air
#     G = air mass flux
#     D = hydraulic diameter
#     th = fin thickness
#     y = distance from outside of tube to tip of fin
#     z = fin gap width
#     *status = non-zero indicated error.
#     ********************************************************************/
#     '''
# 
#     # calc Prantl and Reynold number
#     const double mu = air.mu(T);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffAir","mu");
#         return -1;
#     }
#     const double cp = air.Cp(T);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffAir","cp");
#         return -1;
#     }
#     const double k = air.k(T);
#     if(errorLog.IsError()) {
#         errorLog.Add("ConvCoeffAir","k");
#         return -1;
#     }
#     const double Pr = mu*cp/k;
#     if(Pr<0.0) {
#         errorLog.Add("ConvCoeffAir","Pr<0.0");
#         return -1;
#     }
#     const double Re = G*D/mu;
#     if(Re<0.0) {
#         errorLog.Add("ConvCoeffAir","Re<0.0");
#         return -1;
#     }
# 
#     # strip enhancement factor 0.2<fs<0.35
#     # const double fs = 0.32;
#     const double fs = 0.35;
#     if(fs<0.0) {
#         errorLog.Add("ConvCoeffAir","fs<0.0");
#         return -1;
#     }
# 
#     # fin thickness / gap width, geometery factor
#     const double fa = th/z;
#     if(fa<0.0) {
#         errorLog.Add("ConvCoeffAir","fa<0.0");
#         return -1;
#     }
# 
#     # enhancement factor for slotted fins
#     double Fj = 1 + 1.093e3*pow(fa,1.24)*pow(fs,0.944)*pow(Re,-0.58) + 1.097*pow(fa,2.09)*pow(fs,2.26)*pow(Re,0.88);
# 
#     # the 'j-factor' correlation
#     double j = 0.479*pow(Re,-0.644)*Fj;
# 
#     # calc conv coeff (h) from 'j-factor'(j)
#     double h = j*G*cp/pow(Pr,0.666667);
# 
# 
#     return h


# def FinEffect(h,k,th,w,Do):
#     '''
#     /********************************************************************
#     Computes the fin efficiency of a dry circular fin. shenbo change
#     Reference:
#         Gardner, K. A.
#         "Efficiency of Extended Surface"
#         Transactions of the ASME
#         November, 1945
#         pp. 621-631
#     h = air side convection coefficient (W/K/m^2)
#     k = thermal conductivity of the fin material (W/K/m)
#     th = thickness of the fin (m)
#     w = distance from base to tip of fin (m)
#     Do = outside diameter of the pipe (m)
#     ********************************************************************/
#     '''
# 
#     double phi,yb,xe,xb,xr,ub,ue,B,p1,p2,p3;
# 
#     yb=th/2;
#     xb=Do/2;
#     xe=xb+w;
# 
#     xe=xb+w;
#     xr=xe/xb;
#     arg = h/(k*yb);
# 
#     ub=w*pow(arg,0.5)/(xr-1);
#     ue=ub*xr;
#     B=bessi1(ue)/bessk1(ue);
#     p1=bessi1(ub)-B*bessk1(ub);
#     p2=bessi0(ub)+B*bessk0(ub);
#     p3=2/(ub*(1-xr*xr));
#     phi=p3*p1/p2;
# 
#     return phi

 
if __name__=='__main__':
    print('Hello world')
    print ()
    #Exercise_Corr()
