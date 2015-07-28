from __future__ import division
'''
Created on Apr 29, 2015

@author: AmmarBahman
'''

'''This code is for Direct Expansion in Cooling Mode of ECU 18K'''

from Cycle import ECU_DXCycleClass 
from convert_units import in2m, mm2m, cm2m, cfm2cms, F2K, kPa2Pa, C2K, oz2kg, DeltaF2K, cubin2cubm
from ACHPTools import Write2CSV
from CoolProp.Plots import PropsPlot
from CoolProp.Plots import Ts, drawIsoLines
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt
import numpy

def ECUCycle():
    #########################################################################
    ######################     CYCLE INITIALIZATION    ######################
    #########################################################################
    
    ## Here we load parameters that are not a function of operating conditions
    ## They are primarily geometric parameters

    #Instantiate the cycle class
    Cycle=ECU_DXCycleClass()
    
    
    #--------------------------------------
    #--------------------------------------
    #         Cycle parameters
    #--------------------------------------
    #--------------------------------------
    Cycle.Verbosity = 0 #the idea here is to have different levels of debug output 
    Cycle.ImposedVariable = 'Charge' #'Subcooling' # or 'Charge'
    #Cycle.DT_sc_target = 4.638
    Cycle.Charge_target = oz2kg(36.5) #37-44 ounces #kg #uncomment for use with imposed 'Charge'
    Cycle.Mode='AC'
    Cycle.Ref='R407C'
    Cycle.TestName='ECU-18K'  #this and the two next lines can be used to specify exact test conditions
    Cycle.TestDescription='Test#3'
    Cycle.TestDetails='This is the sample cycle for the ECU18K'
    
    
    #--------------------------------------
    #--------------------------------------
    #       Compressor parameters
    #--------------------------------------
    #--------------------------------------
    #A 1.5 ton (18,000 Btu/h) cooling capacity ECU compressor map
    if Cycle.Ref=='R407C':    
        M=[207.31295549477,4.29717231652541,-2.28622302529118,
           0.0347908258163747,-0.0201167288696277,0.0259153689666968,
           9.13169535150059e-5,-6.23623656573478e-5,1.20986974937733e-4,
           -1.07540716639383e-4]
                      
        P=[-511.9893727,-1.867619312,32.35057515,-0.0573,0.0718,-0.2478335,
           -0.000762,0.00116,-0.000798,0.00129]
    
    params={
            'M':M,
            'P':P,
            'Ref':Cycle.Ref, #Refrigerant
            'fp':0, #Fraction of electrical power lost as heat to ambient 
            'Vdot_ratio': 0.9, #Displacement Scale factor
            'Verbosity': 0, # How verbose should the debugging be [0-10]
            }
    
    Cycle.Compressor.Update(**params)
    
    
    #--------------------------------------
    #--------------------------------------
    #      Condenser parameters
    #--------------------------------------
    #--------------------------------------
    Cycle.Condenser.Fins.Tubes.NTubes=30               #Number of tubes (per bank for now!)
    Cycle.Condenser.Fins.Tubes.Nbank=2                 #Number of banks (set to 1 for now!)
    Cycle.Condenser.Fins.Tubes.Npass=2                 #Number of passes (per bank) #averaged if not even
    Cycle.Condenser.Fins.Tubes.Nports=11               #Number of rectangular ports
    Cycle.Condenser.Fins.Tubes.Ltube=in2m(18.24)       #length of a single tube
    Cycle.Condenser.Fins.Tubes.Td=in2m(1)              #Tube outside width (depth)
    Cycle.Condenser.Fins.Tubes.Ht=in2m(0.072)          #Tube outside height (major diameter)
    Cycle.Condenser.Fins.Tubes.b=in2m(0.488)           #Tube spacing   
    Cycle.Condenser.Fins.Tubes.tw=in2m(0.015)          #Tube wall thickness
    Cycle.Condenser.Fins.Tubes.twp=in2m(0.016)         #Port (channel) wall thickness     
    Cycle.Condenser.Fins.Tubes.beta=1.7675             #Port (channel) aspect ratio (=width/height)
    
    Cycle.Condenser.Fins.Fins.FPI=13                   #Fin per inch
    Cycle.Condenser.Fins.Fins.Lf=in2m(1)               #Fin length = tube outside width in this HX
    Cycle.Condenser.Fins.Fins.t=in2m(0.0045)           ##measured## #Fin thickness
    Cycle.Condenser.Fins.Fins.k_fin=117                #Fin thermal conductivity for pure Aluminum
        
    Cycle.Condenser.Fins.Air.Vdot_ha=cfm2cms(1400)     #Air volume flow rate in m^3/s
    Cycle.Condenser.Fins.Air.Tdb=F2K(105)               #Air inlet temperature, K
    Cycle.Condenser.Fins.Air.p=101325                  #Air pressure in Pa
    Cycle.Condenser.Fins.Air.RH=0.2119                 #Air inlet relative humidity
    Cycle.Condenser.Fins.Air.FanPower=884.1              #Fan power, Watts
        
    Cycle.Condenser.Fins.Louvers.Lalpha=25             ##estimated## #Louver angle, in degree
    Cycle.Condenser.Fins.Louvers.lp=mm2m(1.12)         ##measured## #Louver pitch
    
    Cycle.Condenser.Ref=Cycle.Ref
    Cycle.Condenser.Verbosity=0
    params={
            'Ref': Cycle.Ref,   
            'Verbosity':0,
            }
    Cycle.Condenser.Update(**params)
    
    
    #--------------------------------------
    #--------------------------------------
    #     Evaporator Parameters 
    #--------------------------------------
    #--------------------------------------
    Cycle.Evaporator.Fins.Tubes.NTubes_per_bank=9
    Cycle.Evaporator.Fins.Tubes.Nbank=5
    Cycle.Evaporator.Fins.Tubes.Ltube=in2m(19)
    Cycle.Evaporator.Fins.Tubes.OD=in2m(0.31)
    Cycle.Evaporator.Fins.Tubes.ID=Cycle.Evaporator.Fins.Tubes.OD - in2m(0.016)
    Cycle.Evaporator.Fins.Tubes.Pl=in2m(0.645669)
    Cycle.Evaporator.Fins.Tubes.Pt=in2m(1)
    Cycle.Evaporator.Fins.Tubes.Ncircuits=5
    
    Cycle.Evaporator.Fins.Fins.FPI=14.5
    Cycle.Evaporator.Fins.Fins.Pd=in2m(1.0/16.0)
    Cycle.Evaporator.Fins.Fins.xf=in2m(1.0/4.0/2.0)
    Cycle.Evaporator.Fins.Fins.t=in2m(0.006)
    Cycle.Evaporator.Fins.Fins.k_fin=237
    
    Cycle.Evaporator.Fins.Air.Vdot_ha=cfm2cms(300)          #reducing the flow rate from 600cfm to 300cfm
    Cycle.Evaporator.Fins.Air.Tdb=F2K(85)
    Cycle.Evaporator.Fins.Air.p=101325                                              #Evaporator Air pressure in Pa
    Cycle.Evaporator.Fins.Air.RH=0.2833
    Cycle.Evaporator.Fins.Air.FanPower=394.3
    
    Cycle.Evaporator.FinsType = 'WavyLouveredFins'        #WavyLouveredFins, HerringboneFins, PlainFins
    Cycle.Evaporator.Ref=Cycle.Ref
    Cycle.Evaporator.Verbosity=0
    
    params={
            'Ref': Cycle.Ref,
            'Verbosity':0,
            'DT_sh':9.362, #DeltaF2K()
            }
    
    Cycle.Evaporator.Update(**params)
    
    
    # ----------------------------------
    # ----------------------------------
    #       Line Set Return Parameters
    # ----------------------------------
    # ----------------------------------
    params={
            'L':in2m(55),
            'k_tube':0.19,
            't_insul':0.02,
            'k_insul':0.036,
            'T_air':F2K(105),
            'Ref': Cycle.Ref,
            'h_air':10, #0.0000000001 is changed to 10 assumed for forced air convection
            }
     
    Cycle.LineSetReturn.Update(**params)
    Cycle.LineSetReturn.OD=in2m(5.0/8.0)
    Cycle.LineSetReturn.ID=in2m(5.0/8.0)-mm2m(1) #assumed the thickness is 0.0015

    # ----------------------------------
    # ----------------------------------
    #       Line Set Supply Parameters
    # ----------------------------------
    # ----------------------------------
    params={
            'L':in2m(88),                #tube length in m
            'k_tube':0.19,
            't_insul':0, #no insulation
            'k_insul':0.036,
            'T_air':F2K(105),
            'Ref': Cycle.Ref,
            'h_air':10,#0.0000000001 is changed to 10 assumed for forced air convection
            }
     
    Cycle.LineSetSupply.Update(**params)
    Cycle.LineSetSupply.OD=in2m(3.0/8.0)
    Cycle.LineSetSupply.ID=in2m(3.0/8.0)-mm2m(1) #assumed the thickness is 0.0015
    
    # ----------------------------------
    # ----------------------------------
    #       Sight Glass + Filter Drier + MicroMotion Parameters
    # ----------------------------------
    # ----------------------------------
    params={
            'h':in2m(1.370),        #height of sight glass in m
            'D':in2m(1.110),        #diameter of sight glass in m
            'Ref': Cycle.Ref,
            'V': cubin2cubm(11),    #volume of filter drier
            'E': in2m(9.75),        #micromotion width
            'B': in2m(5.12),        #micormotion height
            'F': in2m(2.81),        #micormotion thickness
            }
     
    Cycle.SightGlassFilterDrierMicroMotion.Update(**params)
    Cycle.SightGlassFilterDrierMicroMotion.ID=in2m(3.0/8.0)-mm2m(1) #tube internal diameter
    
    
    #Now solve
    Cycle.PreconditionedSolve()
    
    #Print Cycle outputs
    for id, unit, value in Cycle.OutputList():
        print str(id) + ' = ' + str(value) + ' ' + str(unit)
    
    return Cycle



if __name__=='__main__':
    cycle=ECUCycle()
    
    #Write the outputs to file
    Write2CSV(cycle,open('results/Cycle_Test#3.csv','w'),append=False)
    #Write2CSV(cycle,open('results/Cycle_Test#3.csv','a'),append=True)
    
    #append a second run with different temperauture
    ###Outdoor side###
#     new_outdoor_temp=F2K(115)
#     new_outdoor_RH=0.1103
#     params={
#         'Tin_a':new_outdoor_temp,
#         'RHin_a': new_outdoor_RH
#         }
#     cycle.Condenser.Update(**params)
#     cycle.Condenser.Fins.Air.Tdb=new_outdoor_temp
#     cycle.Condenser.Fins.Air.RH=new_outdoor_RH
#     ###Indoor side###
#     new_indoor_temp=F2K(85)
#     new_indoor_RH=0.2833
#     params={
#         'Tin_a':new_indoor_temp,
#         'RH': new_indoor_RH
#         }
#     cycle.Evaporator.Update(**params)
#     cycle.Evaporator.Fins.Air.Tdb=new_indoor_temp
#     cycle.Evaporator.Fins.Air.RH=new_indoor_RH
#      
#     #file name
#     cycle.TestName='ECU-18K'  #this and the two next lines can be used to specify exact test conditions
#     cycle.TestDescription='Test#3'
#     cycle.TestDetails='Here we changed the air condition on evaporator and condeser'
#     cycle.PreconditionedSolve()  #there seems to be a problem, somewhere
#     Write2CSV(cycle,open('Cycle.csv','a'),append=True)

    
    #P-h & T-s diagrams for R407C
    ref_fluid = 'R407C'
    
    #Experimental results 
    P_exp = [655.8,3108.0,3108.0,3095.0,3095.0,877.8,655.8,655.8,655.8] #in kPa 
    P_exp = numpy.array(P_exp)
    P_exp *= 1000.0 #convert kPa to Pa
    T_exp = [17.92+273.15, 111.0+273.15, PropsSI('T','P',P_exp[2],'Q',1,ref_fluid), PropsSI('T','P',P_exp[3],'Q',0,ref_fluid), 60.57+273.15, 16.25+273.15, PropsSI('T','P',P_exp[6],'Q',1,ref_fluid), 16.66+273.15, 17.92+273.15] #in Kelvin    
    T_exp = numpy.array(T_exp)
    
    #Solve for h_exp and s_exp
    h_exp = [PropsSI('H','P',P_exp[0],'T',T_exp[0],ref_fluid), PropsSI('H','P',P_exp[1],'T',T_exp[1],ref_fluid), 
             PropsSI('H','P',P_exp[2],'Q',1,ref_fluid), PropsSI('H','P',P_exp[3],'Q',0,ref_fluid),
             PropsSI('H','P',P_exp[4],'T',T_exp[4],ref_fluid), PropsSI('H','P',P_exp[4],'T',T_exp[4],ref_fluid),
             PropsSI('H','P',P_exp[6],'Q',1,ref_fluid), PropsSI('H','P',P_exp[7],'T',T_exp[7],ref_fluid), 
             PropsSI('H','P',P_exp[8],'T',T_exp[8],ref_fluid)]
    
    #P_exp[5] = PropsSI('P','Q',0,'T',T_exp[5],ref_fluid)

    s_exp = [PropsSI('S','P',P_exp[0],'T',T_exp[0],ref_fluid), PropsSI('S','P',P_exp[1],'T',T_exp[1],ref_fluid), 
             PropsSI('S','P',P_exp[2],'Q',1,ref_fluid), PropsSI('S','P',P_exp[3],'Q',0,ref_fluid),
             PropsSI('S','P',P_exp[4],'T',T_exp[4],ref_fluid), PropsSI('S','H',h_exp[5],'P',P_exp[5],ref_fluid),
             PropsSI('S','P',P_exp[6],'Q',1,ref_fluid), PropsSI('S','P',P_exp[7],'T',T_exp[7],ref_fluid), 
             PropsSI('S','P',P_exp[8],'T',T_exp[8],ref_fluid)]
    h_exp = numpy.array(h_exp)
    s_exp = numpy.array(s_exp)
    
            
    #convert back to original units kPa, kJ/kg and kJ/kg-K
    P_exp /= 1000.0 
    T_exp = T_exp #keep T in K
    h_exp /= 1000.0
    s_exp /= 1000.0
    
    #Model Results
    P = [cycle.Compressor.pin_r, cycle.Compressor.pout_r, cycle.Condenser.psat_r, cycle.Condenser.psat_r, cycle.Condenser.psat_r + cycle.Condenser.DP_r , cycle.Evaporator.psat_r, cycle.Evaporator.psat_r + cycle.Evaporator.DP_r, cycle.LineSetReturn.pin, cycle.Compressor.pin_r]
    T = [cycle.Compressor.Tin_r, cycle.Compressor.Tout_r, PropsSI('T','P',cycle.Condenser.psat_r,'Q',1,ref_fluid), PropsSI('T','P',cycle.Condenser.psat_r,'Q',0,ref_fluid), cycle.Condenser.Tout_r, cycle.Evaporator.Tin_r, PropsSI('T','P',cycle.Evaporator.psat_r,'Q',1,ref_fluid), cycle.LineSetReturn.Tout, cycle.Compressor.Tin_r]
    h = [cycle.Compressor.hin_r, cycle.Compressor.hout_r, PropsSI('H','P',cycle.Condenser.psat_r,'Q',1,ref_fluid), PropsSI('H','P',cycle.Condenser.psat_r,'Q',0,ref_fluid), cycle.Condenser.hout_r, cycle.Evaporator.hin_r, PropsSI('H','P',cycle.Evaporator.psat_r,'Q',1,ref_fluid), cycle.LineSetReturn.hout, cycle.Compressor.hin_r]
    s = [cycle.Compressor.sin_r, cycle.Compressor.sout_r, PropsSI('S','P',cycle.Condenser.psat_r,'Q',1,ref_fluid), PropsSI('S','P',cycle.Condenser.psat_r,'Q',0,ref_fluid), cycle.Condenser.sout_r, cycle.Evaporator.sin_r, PropsSI('S','P',cycle.Evaporator.psat_r,'Q',1,ref_fluid), PropsSI('S','T',cycle.LineSetReturn.Tout,'P',cycle.LineSetReturn.pin,ref_fluid), cycle.Compressor.sin_r]
    P = numpy.array(P)
    T = numpy.array(T)
    h = numpy.array(h)
    s = numpy.array(s)
    P /= 1000.0 #convert Pa to kPa
    T = T       #keep T in K
    h /= 1000.0 #convert J/kg to kJ/kg
    s /= 1000.0 #convert J/kg-K to kJ/kg-K

    #Plot P-h diagram 
#     ph_plot_R407C = PropsPlot(ref_fluid, 'Ph')
#     ph_plot_R407C.title('$P-h$ $R407C$')
#     ph_plot_R407C.xlabel(r'$h$ $[{kJ}/{kg}]$')
#     ph_plot_R407C.ylabel(r'$P$ $[kPa]$')
#     ph_plot_R407C.axis.set_yscale('log')
#     ph_plot_R407C.grid()
#     plt.plot(h_exp,P_exp, 'bo-', label='Experimental')
#     plt.plot(h,P,'ro--', label='Model')
#     plt.legend(loc='best',fancybox=False)
#     ph_plot_R407C.savefig('images/R407C_Ph_Test3.pdf')    
#     ph_plot_R407C.show()
#      
#     #Plot T-s diagram  
#     ts_plot_R407C = PropsPlot(ref_fluid, 'Ts')
#     ts_plot_R407C.title('$T-s$ $R407C$')
#     ts_plot_R407C.xlabel(r'$s$ $[{kJ}/{kg-K}]$')
#     ts_plot_R407C.ylabel(r'$T$ $[K]$')
#     ts_plot_R407C.grid()
#     plt.plot(s_exp,T_exp, 'bo-', label='Experimental')
#     plt.plot(s,T,'ro--', label='Model')
#     plt.legend(loc='best',fancybox=False)
#     ts_plot_R407C.savefig('images/R407C_Ts_Test3.pdf')    
#     ts_plot_R407C.show()
    
    #Plot T-s and P-h diagrams in one graph
    fig = plt.figure(1, figsize=(16, 8), dpi=100)
    for i, gtype in enumerate(['Ph', 'Ts']):
        ax = plt.subplot(1, 2, i+1)
        if gtype.startswith('P'):
            ax.set_yscale('log')
            plt.grid()
            plt.plot(h_exp,P_exp, 'bo-', label='Experimental')
            plt.plot(h,P,'ro--', label='Model')
            plt.legend(loc='best',fancybox=False)
        if gtype.startswith('T'):
            plt.grid()
            plt.plot(s_exp,T_exp, 'bo-', label='Experimental')
            plt.plot(s,T,'ro--', label='Model')
            plt.legend(loc='best',fancybox=False)
        props_plot = PropsPlot(ref_fluid, gtype, axis=ax)
        props_plot.title(gtype)
        props_plot._draw_graph()
    fig.set_tight_layout(True)
    fig.savefig('images/comined_R407C_Test3.pdf')
    props_plot.show()