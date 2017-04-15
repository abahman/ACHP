from __future__ import division
'''
Created on Dec 6, 2016

@author: AmmarBahman
'''

'''This code is for Cooling Mode of ECU 60K'''

from Cycle import ECU_DXCycleClass, ECU_VICompCycleClass, ECU_VISemiEmpCompCycleClass, ECU_VICompTelloCycleClass
from convert_units import *
from ACHPTools import Write2CSV
import CoolProp
from CoolProp.Plots import PropertyPlot
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt
import numpy as np
from math import pi

#===============================================================================
# Latex render
#===============================================================================
import matplotlib as mpl
#mpl.use('pgf')

def figsize(scale):
    fig_width_pt = 469.755                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size
 
pgf_with_latex = {                      # setup matplotlib to use latex for output
"pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
"text.usetex": True,                # use LaTeX to write all text
"font.family": "serif",
"font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
"font.sans-serif": [],
"font.monospace": [],
"axes.labelsize": 10,               # LaTeX default is 10pt font.
"font.size": 10,
"legend.fontsize": 8,               # Make the legend/label fonts a little smaller
"legend.labelspacing":0.2,
"xtick.labelsize": 8,
"ytick.labelsize": 8,
"figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
"pgf.preamble": [
r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)
#===============================================================================
# END of Latex render
#===============================================================================

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
    Cycle.ImposedVariable = 'Subcooling' #'Subcooling' # or 'Charge'
    Cycle.DT_sc_target = 4.9
    #Cycle.Charge_target = oz2kg(85) #37-44 ounces #kg #uncomment for use with imposed 'Charge'
    Cycle.Mode='AC'
    Cycle.Ref='R407C'
    Cycle.Backend='HEOS' #Backend for refrigerant properties calculation: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    Cycle.TestName='ECU-60K'  #this and the two next lines can be used to specify exact test conditions
    Cycle.TestDescription='Test1'
    Cycle.TestDetails='This is the sample cycle for the ECU60K'
    
    
    #--------------------------------------
    #--------------------------------------
    #       Compressor parameters
    #--------------------------------------
    #--------------------------------------
    #A 5 ton (60,000 Btu/h) cooling capacity ECU compressor map
    if Cycle.Ref=='R407C':    
        M=[493.9060555,14.07268806,-7.345763459,0.112578167,-0.053492912,0.0734872,0.000135982,-7.90201E-05,0.000423591,-0.000313537]
                      
        P=[-2957.179994,3.992088672,118.3575555,0.354681778,-0.0946,-0.88023583,0.00326,-0.00536,0.00137,0.00372]
    
    params={
            'M':M,
            'P':P,
            'Ref':Cycle.Ref, #Refrigerant
            'fp':0.12, #Fraction of electrical power lost as heat to ambient 
            'Vdot_ratio': 1, #Displacement Scale factor
            'Verbosity': 0, # How verbose should the debugging be [0-10]
            'Backend':Cycle.Backend
            }
    
    Cycle.Compressor.Update(**params)
    
    
    #--------------------------------------
    #--------------------------------------
    #      Condenser parameters
    #--------------------------------------
    #--------------------------------------
    Cycle.Condenser.Fins.Tubes.NTubes=52               #Number of tubes (per bank for now!)
    Cycle.Condenser.Fins.Tubes.Nbank=2                 #Number of banks (set to 1 for now!)
    Cycle.Condenser.Fins.Tubes.Npass=2                 #Number of passes (per bank) #averaged if not even
    Cycle.Condenser.Fins.Tubes.Nports=11               #Number of rectangular ports
    Cycle.Condenser.Fins.Tubes.Ltube=in2m(21.26)       #length of a single tube
    Cycle.Condenser.Fins.Tubes.Td=in2m(1)              #Tube outside width (depth)
    Cycle.Condenser.Fins.Tubes.Ht=in2m(0.072)          #Tube outside height (major diameter)
    Cycle.Condenser.Fins.Tubes.b=in2m(0.488)           #Tube spacing   
    Cycle.Condenser.Fins.Tubes.tw=in2m(0.015)          #Tube wall thickness
    Cycle.Condenser.Fins.Tubes.twp=in2m(0.016)         #Port (channel) wall thickness     
    Cycle.Condenser.Fins.Tubes.beta=1.7675             #Port (channel) aspect ratio (=width/height)
    
    Cycle.Condenser.Fins.Fins.FPI=14                   #Fin per inch
    Cycle.Condenser.Fins.Fins.Lf=in2m(1)               #Fin length = tube outside width in this HX
    Cycle.Condenser.Fins.Fins.t=in2m(0.0045)           ##measured## #Fin thickness
    Cycle.Condenser.Fins.Fins.k_fin=117                #Fin thermal conductivity for pure Aluminum
        
    Cycle.Condenser.Fins.Air.Vdot_ha=cfm2cms(3500)     #Air volume flow rate in m^3/s
    Cycle.Condenser.Fins.Air.Tdb=F2K(125)               #Air inlet temperature, K
    Cycle.Condenser.Fins.Air.p=101325                  #Air pressure in Pa
    Cycle.Condenser.Fins.Air.RH=0.199                 #Air inlet relative humidity
    Cycle.Condenser.Fins.Air.FanPower=984.5              #Fan power, Watts
        
    Cycle.Condenser.Fins.Louvers.Lalpha=25             ##estimated## #Louver angle, in degree
    Cycle.Condenser.Fins.Louvers.lp=mm2m(1.12)         ##measured## #Louver pitch
    
    Cycle.Condenser.Ref=Cycle.Ref
    Cycle.Condenser.Verbosity=0
    params={
            'Ref': Cycle.Ref,
            'Backend': Cycle.Backend,
            'Verbosity':0,
            }
    Cycle.Condenser.Update(**params)
    
    
    #--------------------------------------
    #--------------------------------------
    #     Evaporator Parameters 
    #--------------------------------------
    #--------------------------------------
    Cycle.Evaporator.Fins.Tubes.NTubes_per_bank=18
    Cycle.Evaporator.Fins.Tubes.Nbank=4
    Cycle.Evaporator.Fins.Tubes.Ltube=in2m(24.875)
    Cycle.Evaporator.Fins.Tubes.OD=in2m(0.5)
    Cycle.Evaporator.Fins.Tubes.ID=Cycle.Evaporator.Fins.Tubes.OD - 2*in2m(0.019)
    Cycle.Evaporator.Fins.Tubes.Pl=in2m(1.082)
    Cycle.Evaporator.Fins.Tubes.Pt=in2m(1.25)
    Cycle.Evaporator.Fins.Tubes.Ncircuits=6
    
    Cycle.Evaporator.Fins.Fins.FPI=12
    Cycle.Evaporator.Fins.Fins.Pd=in2m(1.0/16.0/2)
    Cycle.Evaporator.Fins.Fins.xf=in2m(1.0/4.0)
    Cycle.Evaporator.Fins.Fins.t=in2m(0.0075)
    Cycle.Evaporator.Fins.Fins.k_fin=237
    
    Cycle.Evaporator.Fins.Air.Vdot_ha=cfm2cms(1700)          
    Cycle.Evaporator.Fins.Air.Tdb=F2K(90)
    Cycle.Evaporator.Fins.Air.p=101325                       #Evaporator Air pressure in Pa
    Cycle.Evaporator.Fins.Air.RH=0.5023
    Cycle.Evaporator.Fins.Air.FanPower=750.2
    
    Cycle.Evaporator.FinsType = 'WavyLouveredFins'        #WavyLouveredFins, HerringboneFins, PlainFins
    Cycle.Evaporator.Verbosity=0
    
    params={
            'Ref': Cycle.Ref,
            'Backend': Cycle.Backend,
            'Verbosity':0,
            'DT_sh':16.42, #DeltaF2K()
            }
    
    Cycle.Evaporator.Update(**params)
    
    
    # ----------------------------------
    # ----------------------------------
    #       Line Set Return Parameters (vapor line)
    # ----------------------------------
    # ----------------------------------
#     params={
#             'L':in2m(85),
#             'k_tube':0.19,
#             't_insul':0.02,
#             'k_insul':0.036,
#             'T_air':F2K(125),
#             'Ref': Cycle.Ref,
#             'h_air':10, #0.0000000001 is changed to 10 assumed for forced air convection
#             }
#      
#     Cycle.LineSetReturn.Update(**params)
#     Cycle.LineSetReturn.OD=in2m(5.0/8.0)
#     Cycle.LineSetReturn.ID=in2m(5.0/8.0)-mm2m(2) #wall thickness is 1mm

    # ----------------------------------
    # ----------------------------------
    #       Line Set Supply Parameters (liquid line)
    # ----------------------------------
    # ----------------------------------
#     params={
#             'L':in2m(183),                #tube length in m
#             'k_tube':0.19,
#             't_insul':0, #no insulation
#             'k_insul':0.036,
#             'T_air':F2K(125),
#             'Ref': Cycle.Ref,
#             'h_air':10,#0.0000000001 is changed to 10 assumed for forced air convection
#             }
#      
#     Cycle.LineSetSupply.Update(**params)
#     Cycle.LineSetSupply.OD=in2m(3.0/8.0)
#     Cycle.LineSetSupply.ID=in2m(3.0/8.0)-mm2m(2) #wall thickness is 1mm
    
    # ----------------------------------
    # ----------------------------------
    #       Sight Glass + Filter Drier + MicroMotion Parameters
    # ----------------------------------
    # ----------------------------------
#     params={
#             'h':in2m(1.370),        #height of sight glass in m
#             'D':in2m(1.110),        #diameter of sight glass in m
#             'Ref': Cycle.Ref,
#             'V': cubin2cubm(30), #volume of filter drier (website = 13.74in^3 calculated) (manual = 16in^3)
#             'D_Micro': in2m(0.21),  #micromotion tube diameter
#             'L_Micro': in2m(14.6),  #micormotion tube length
#             'n_Micro': 2,           #micormotion number of tubes
#             }
#      
#     Cycle.SightGlassFilterDrierMicroMotion.Update(**params)
#     Cycle.SightGlassFilterDrierMicroMotion.ID=in2m(3.0/8.0)-mm2m(2) #wall thickness is 1mm
#     
    
    #Now solve
    Cycle.PreconditionedSolve()
    
    #Print Cycle outputs
    for id, unit, value in Cycle.OutputList():
        print str(id) + ' = ' + str(value) + ' ' + str(unit)
    
    return Cycle

def ECUCycle_NewComp():
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
    Cycle.ImposedVariable = 'Subcooling' #'Subcooling' # or 'Charge'
    Cycle.DT_sc_target = 4.9
    #Cycle.Charge_target = oz2kg(85) #37-44 ounces #kg #uncomment for use with imposed 'Charge'
    Cycle.Mode='AC'
    Cycle.Ref='R407C'
    Cycle.Backend='HEOS' #Backend for refrigerant properties calculation: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    Cycle.TestName='ECU-60K'  #this and the two next lines can be used to specify exact test conditions
    Cycle.TestDescription='Test1'
    Cycle.TestDetails='This is ECU60K with VIComp'
    
    
    #--------------------------------------
    #--------------------------------------
    #       Compressor parameters
    #--------------------------------------
    #--------------------------------------
    #A 5 ton (60,000 Btu/h) cooling capacity ECU compressor map
    if Cycle.Ref=='R407C':    
        M=[358.983723541832,6.32163221180689,-1.45127942206402,0.0629308018032251,0.0215942663002772,0.0138397417834867,4.34652176406554E-4,1.46761858399946E-4,-1.03838833412122E-4,-5.30570544503368E-5]
                      
        P=[-1678.57211390132,-16.2082411263196,101.455598791218,-0.800561879709571,0.582915518150735,-0.94196656452479,0.0015145484818854,0.00639233243871854,-0.00308440373694123,0.00393430167789957]
    
    params={
            'M':M,
            'P':P,
            'Ref':Cycle.Ref, #Refrigerant
            'fp':0.12, #Fraction of electrical power lost as heat to ambient 
            'Vdot_ratio': 1, #Displacement Scale factor
            'Verbosity': 0, # How verbose should the debugging be [0-10]
            'Backend':Cycle.Backend
            }
    
    Cycle.Compressor.Update(**params)
    
    
    #--------------------------------------
    #--------------------------------------
    #      Condenser parameters
    #--------------------------------------
    #--------------------------------------
    Cycle.Condenser.Fins.Tubes.NTubes=52               #Number of tubes (per bank for now!)
    Cycle.Condenser.Fins.Tubes.Nbank=2                 #Number of banks (set to 1 for now!)
    Cycle.Condenser.Fins.Tubes.Npass=2                 #Number of passes (per bank) #averaged if not even
    Cycle.Condenser.Fins.Tubes.Nports=11               #Number of rectangular ports
    Cycle.Condenser.Fins.Tubes.Ltube=in2m(21.26)       #length of a single tube
    Cycle.Condenser.Fins.Tubes.Td=in2m(1)              #Tube outside width (depth)
    Cycle.Condenser.Fins.Tubes.Ht=in2m(0.072)          #Tube outside height (major diameter)
    Cycle.Condenser.Fins.Tubes.b=in2m(0.488)           #Tube spacing   
    Cycle.Condenser.Fins.Tubes.tw=in2m(0.015)          #Tube wall thickness
    Cycle.Condenser.Fins.Tubes.twp=in2m(0.016)         #Port (channel) wall thickness     
    Cycle.Condenser.Fins.Tubes.beta=1.7675             #Port (channel) aspect ratio (=width/height)
    
    Cycle.Condenser.Fins.Fins.FPI=14                   #Fin per inch
    Cycle.Condenser.Fins.Fins.Lf=in2m(1)               #Fin length = tube outside width in this HX
    Cycle.Condenser.Fins.Fins.t=in2m(0.0045)           ##measured## #Fin thickness
    Cycle.Condenser.Fins.Fins.k_fin=117                #Fin thermal conductivity for pure Aluminum
        
    Cycle.Condenser.Fins.Air.Vdot_ha=cfm2cms(3500)     #Air volume flow rate in m^3/s
    Cycle.Condenser.Fins.Air.Tdb=F2K(125)               #Air inlet temperature, K
    Cycle.Condenser.Fins.Air.p=101325                  #Air pressure in Pa
    Cycle.Condenser.Fins.Air.RH=0.199                 #Air inlet relative humidity
    Cycle.Condenser.Fins.Air.FanPower=984.5              #Fan power, Watts
        
    Cycle.Condenser.Fins.Louvers.Lalpha=25             ##estimated## #Louver angle, in degree
    Cycle.Condenser.Fins.Louvers.lp=mm2m(1.12)         ##measured## #Louver pitch
    
    Cycle.Condenser.Ref=Cycle.Ref
    Cycle.Condenser.Verbosity=0
    params={
            'Ref': Cycle.Ref,
            'Backend': Cycle.Backend,
            'Verbosity':0,
            }
    Cycle.Condenser.Update(**params)
    
    
    #--------------------------------------
    #--------------------------------------
    #     Evaporator Parameters 
    #--------------------------------------
    #--------------------------------------
    Cycle.Evaporator.Fins.Tubes.NTubes_per_bank=18
    Cycle.Evaporator.Fins.Tubes.Nbank=4
    Cycle.Evaporator.Fins.Tubes.Ltube=in2m(24.875)
    Cycle.Evaporator.Fins.Tubes.OD=in2m(0.5)
    Cycle.Evaporator.Fins.Tubes.ID=Cycle.Evaporator.Fins.Tubes.OD - 2*in2m(0.019)
    Cycle.Evaporator.Fins.Tubes.Pl=in2m(1.082)
    Cycle.Evaporator.Fins.Tubes.Pt=in2m(1.25)
    Cycle.Evaporator.Fins.Tubes.Ncircuits=6
    
    Cycle.Evaporator.Fins.Fins.FPI=12
    Cycle.Evaporator.Fins.Fins.Pd=in2m(1.0/16.0/2)
    Cycle.Evaporator.Fins.Fins.xf=in2m(1.0/4.0)
    Cycle.Evaporator.Fins.Fins.t=in2m(0.0075)
    Cycle.Evaporator.Fins.Fins.k_fin=237
    
    Cycle.Evaporator.Fins.Air.Vdot_ha=cfm2cms(1700)          
    Cycle.Evaporator.Fins.Air.Tdb=F2K(90)
    Cycle.Evaporator.Fins.Air.p=101325                       #Evaporator Air pressure in Pa
    Cycle.Evaporator.Fins.Air.RH=0.5023
    Cycle.Evaporator.Fins.Air.FanPower=750.2
    
    Cycle.Evaporator.FinsType = 'WavyLouveredFins'        #WavyLouveredFins, HerringboneFins, PlainFins
    Cycle.Evaporator.Verbosity=0
    
    params={
            'Ref': Cycle.Ref,
            'Backend': Cycle.Backend,
            'Verbosity':0,
            'DT_sh':16.42, #DeltaF2K()
            }
    
    Cycle.Evaporator.Update(**params)
    
    
    # ----------------------------------
    # ----------------------------------
    #       Line Set Return Parameters (vapor line)
    # ----------------------------------
    # ----------------------------------
#     params={
#             'L':in2m(85),
#             'k_tube':0.19,
#             't_insul':0.02,
#             'k_insul':0.036,
#             'T_air':F2K(125),
#             'Ref': Cycle.Ref,
#             'h_air':10, #0.0000000001 is changed to 10 assumed for forced air convection
#             }
#      
#     Cycle.LineSetReturn.Update(**params)
#     Cycle.LineSetReturn.OD=in2m(5.0/8.0)
#     Cycle.LineSetReturn.ID=in2m(5.0/8.0)-mm2m(2) #wall thickness is 1mm

    # ----------------------------------
    # ----------------------------------
    #       Line Set Supply Parameters (liquid line)
    # ----------------------------------
    # ----------------------------------
#     params={
#             'L':in2m(183),                #tube length in m
#             'k_tube':0.19,
#             't_insul':0, #no insulation
#             'k_insul':0.036,
#             'T_air':F2K(125),
#             'Ref': Cycle.Ref,
#             'h_air':10,#0.0000000001 is changed to 10 assumed for forced air convection
#             }
#      
#     Cycle.LineSetSupply.Update(**params)
#     Cycle.LineSetSupply.OD=in2m(3.0/8.0)
#     Cycle.LineSetSupply.ID=in2m(3.0/8.0)-mm2m(2) #wall thickness is 1mm
    
    # ----------------------------------
    # ----------------------------------
    #       Sight Glass + Filter Drier + MicroMotion Parameters
    # ----------------------------------
    # ----------------------------------
#     params={
#             'h':in2m(1.370),        #height of sight glass in m
#             'D':in2m(1.110),        #diameter of sight glass in m
#             'Ref': Cycle.Ref,
#             'V': cubin2cubm(30), #volume of filter drier (website = 13.74in^3 calculated) (manual = 16in^3)
#             'D_Micro': in2m(0.21),  #micromotion tube diameter
#             'L_Micro': in2m(14.6),  #micormotion tube length
#             'n_Micro': 2,           #micormotion number of tubes
#             }
#      
#     Cycle.SightGlassFilterDrierMicroMotion.Update(**params)
#     Cycle.SightGlassFilterDrierMicroMotion.ID=in2m(3.0/8.0)-mm2m(2) #wall thickness is 1mm
#     
    
    #Now solve
    Cycle.PreconditionedSolve()
    
    #Print Cycle outputs
    for id, unit, value in Cycle.OutputList():
        print str(id) + ' = ' + str(value) + ' ' + str(unit)
    
    return Cycle

def ECUCycle_VIComp():
    #########################################################################
    ######################     CYCLE INITIALIZATION    ######################
    #########################################################################
    
    ## Here we load parameters that are not a function of operating conditions
    ## They are primarily geometric parameters

    #Instantiate the cycle class
    Cycle=ECU_VICompCycleClass()
    
    
    #--------------------------------------
    #--------------------------------------
    #         Cycle parameters
    #--------------------------------------
    #--------------------------------------
    Cycle.Verbosity = 0 #the idea here is to have different levels of debug output 
    Cycle.ImposedVariable = 'Subcooling' #'Subcooling' # or 'Charge'
    Cycle.DT_sc_target = 4.9
    #Cycle.Charge_target = oz2kg(85) #37-44 ounces #kg #uncomment for use with imposed 'Charge'
    Cycle.Mode='AC'
    Cycle.Ref='R407C'
    Cycle.Backend='HEOS' #Backend for refrigerant properties calculation: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    Cycle.TestName='ECU-60K-VI'  #this and the two next lines can be used to specify exact test conditions
    Cycle.TestDescription='Test1'
    Cycle.TestDetails='This is ECU60K with VIComp new system'
    
    
    #--------------------------------------
    #--------------------------------------
    #       VI Compressor parameters
    #--------------------------------------
    #--------------------------------------
    #Dimensioless compressor map
    if Cycle.Ref=='R407C':    
        M =[-0.08135903196,0.0458453839,0.129432516,0.1325702995,-0.1485129565,9.144668476]
        P = [2102.503527,0.6448139189,0.6561921643,0.1690042408,1.270243282,41.2967932]
        R = [-22.45861341,22.05630275,0.007613286149,-0.01443163273,0.01145168107,0.057749393]
        Te = [317.359575,0.6939083124,0.4131861045,0.6598603942,1.009091398,16.1910005]
    params={
            'M':M,
            'P':P,
            'R':R,
            'Te':Te,
            'Ref':Cycle.Ref, #Refrigerant
            'fp':0.12, #Fraction of electrical power lost as heat to ambient 
            'Vdot_ratio': 1, #Displacement Scale factor
            'Verbosity': 0, # How verbose should the debugging be [0-10]
            'Backend':Cycle.Backend
            }
    
    Cycle.Compressor.Update(**params)
    
    
    #--------------------------------------
    #--------------------------------------
    #      Condenser parameters
    #--------------------------------------
    #--------------------------------------
    Cycle.Condenser.Fins.Tubes.NTubes=52               #Number of tubes (per bank for now!)
    Cycle.Condenser.Fins.Tubes.Nbank=2                 #Number of banks (set to 1 for now!)
    Cycle.Condenser.Fins.Tubes.Npass=2                 #Number of passes (per bank) #averaged if not even
    Cycle.Condenser.Fins.Tubes.Nports=11               #Number of rectangular ports
    Cycle.Condenser.Fins.Tubes.Ltube=in2m(21.26)       #length of a single tube
    Cycle.Condenser.Fins.Tubes.Td=in2m(1)              #Tube outside width (depth)
    Cycle.Condenser.Fins.Tubes.Ht=in2m(0.072)          #Tube outside height (major diameter)
    Cycle.Condenser.Fins.Tubes.b=in2m(0.488)           #Tube spacing   
    Cycle.Condenser.Fins.Tubes.tw=in2m(0.015)          #Tube wall thickness
    Cycle.Condenser.Fins.Tubes.twp=in2m(0.016)         #Port (channel) wall thickness     
    Cycle.Condenser.Fins.Tubes.beta=1.7675             #Port (channel) aspect ratio (=width/height)
    
    Cycle.Condenser.Fins.Fins.FPI=14                   #Fin per inch
    Cycle.Condenser.Fins.Fins.Lf=in2m(1)               #Fin length = tube outside width in this HX
    Cycle.Condenser.Fins.Fins.t=in2m(0.0045)           ##measured## #Fin thickness
    Cycle.Condenser.Fins.Fins.k_fin=117                #Fin thermal conductivity for pure Aluminum
        
    Cycle.Condenser.Fins.Air.Vdot_ha=cfm2cms(3500)     #Air volume flow rate in m^3/s
    Cycle.Condenser.Fins.Air.Tdb=F2K(125)               #Air inlet temperature, K
    Cycle.Condenser.Fins.Air.p=101325                  #Air pressure in Pa
    Cycle.Condenser.Fins.Air.RH=0.199                 #Air inlet relative humidity
    Cycle.Condenser.Fins.Air.FanPower=984.5              #Fan power, Watts
        
    Cycle.Condenser.Fins.Louvers.Lalpha=25             ##estimated## #Louver angle, in degree
    Cycle.Condenser.Fins.Louvers.lp=mm2m(1.12)         ##measured## #Louver pitch
    
    Cycle.Condenser.Ref=Cycle.Ref
    Cycle.Condenser.Verbosity=0
    params={
            'Ref': Cycle.Ref,
            'Backend': Cycle.Backend,
            'Verbosity':0,
            }
    Cycle.Condenser.Update(**params)
    
    
    #--------------------------------------
    #--------------------------------------
    #     Evaporator Parameters 
    #--------------------------------------
    #--------------------------------------
    Cycle.Evaporator.Fins.Tubes.NTubes_per_bank=18
    Cycle.Evaporator.Fins.Tubes.Nbank=4
    Cycle.Evaporator.Fins.Tubes.Ltube=in2m(24.875)
    Cycle.Evaporator.Fins.Tubes.OD=in2m(0.5)
    Cycle.Evaporator.Fins.Tubes.ID=Cycle.Evaporator.Fins.Tubes.OD - 2*in2m(0.019)
    Cycle.Evaporator.Fins.Tubes.Pl=in2m(1.082)
    Cycle.Evaporator.Fins.Tubes.Pt=in2m(1.25)
    Cycle.Evaporator.Fins.Tubes.Ncircuits=6
    
    Cycle.Evaporator.Fins.Fins.FPI=12
    Cycle.Evaporator.Fins.Fins.Pd=in2m(1.0/16.0/2)
    Cycle.Evaporator.Fins.Fins.xf=in2m(1.0/4.0)
    Cycle.Evaporator.Fins.Fins.t=in2m(0.0075)
    Cycle.Evaporator.Fins.Fins.k_fin=237
    
    Cycle.Evaporator.Fins.Air.Vdot_ha=cfm2cms(1700)          
    Cycle.Evaporator.Fins.Air.Tdb=F2K(90)
    Cycle.Evaporator.Fins.Air.p=101325                       #Evaporator Air pressure in Pa
    Cycle.Evaporator.Fins.Air.RH=0.5023
    Cycle.Evaporator.Fins.Air.FanPower=750.2
    
    Cycle.Evaporator.FinsType = 'WavyLouveredFins'        #WavyLouveredFins, HerringboneFins, PlainFins
    Cycle.Evaporator.Verbosity=0
    
    params={
            'Ref': Cycle.Ref,
            'Backend': Cycle.Backend,
            'Verbosity':0,
            'DT_sh':5,#16.42, #DeltaF2K()
            }
    
    Cycle.Evaporator.Update(**params)
    
    # ----------------------------------
    # ----------------------------------
    #       PHEHX Parameters
    # ----------------------------------
    # ----------------------------------
    params={
            'Ref_c':Cycle.Ref,
            'Backend_c': Cycle.Backend,
            'Ref_h':Cycle.Ref,
            'Backend_h': Cycle.Backend,
            'Verbosity':0,
            'DT_sh_target':5,
            
            #Geometric parameters
            'Bp' : in2m(2.875),
            'Lp' : in2m(18), #Center-to-center distance between ports
            'Nplates' : 10,
            'PlateAmplitude' : 0.001, #[m]
            'PlateThickness' : 0.0003, #[m]
            'PlateWavelength' : 0.00626, #[m]
            'InclinationAngle' : 65/180*pi,#[rad]
            'PlateConductivity' : 15.0, #[W/m-K]
            'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
            }
    Cycle.PHEHX.Update(**params)
        
    # ----------------------------------
    # ----------------------------------
    #       Line Set Return Parameters (vapor line)
    # ----------------------------------
    # ----------------------------------
#     params={
#             'L':in2m(85),
#             'k_tube':0.19,
#             't_insul':0.02,
#             'k_insul':0.036,
#             'T_air':F2K(125),
#             'Ref': Cycle.Ref,
#             'h_air':10, #0.0000000001 is changed to 10 assumed for forced air convection
#             }
#      
#     Cycle.LineSetReturn.Update(**params)
#     Cycle.LineSetReturn.OD=in2m(5.0/8.0)
#     Cycle.LineSetReturn.ID=in2m(5.0/8.0)-mm2m(2) #wall thickness is 1mm

    # ----------------------------------
    # ----------------------------------
    #       Line Set Supply Parameters (liquid line)
    # ----------------------------------
    # ----------------------------------
#     params={
#             'L':in2m(183),                #tube length in m
#             'k_tube':0.19,
#             't_insul':0, #no insulation
#             'k_insul':0.036,
#             'T_air':F2K(125),
#             'Ref': Cycle.Ref,
#             'h_air':10,#0.0000000001 is changed to 10 assumed for forced air convection
#             }
#      
#     Cycle.LineSetSupply.Update(**params)
#     Cycle.LineSetSupply.OD=in2m(3.0/8.0)
#     Cycle.LineSetSupply.ID=in2m(3.0/8.0)-mm2m(2) #wall thickness is 1mm
    
    # ----------------------------------
    # ----------------------------------
    #       Sight Glass + Filter Drier + MicroMotion Parameters
    # ----------------------------------
    # ----------------------------------
#     params={
#             'h':in2m(1.370),        #height of sight glass in m
#             'D':in2m(1.110),        #diameter of sight glass in m
#             'Ref': Cycle.Ref,
#             'V': cubin2cubm(30), #volume of filter drier (website = 13.74in^3 calculated) (manual = 16in^3)
#             'D_Micro': in2m(0.21),  #micromotion tube diameter
#             'L_Micro': in2m(14.6),  #micormotion tube length
#             'n_Micro': 2,           #micormotion number of tubes
#             }
#      
#     Cycle.SightGlassFilterDrierMicroMotion.Update(**params)
#     Cycle.SightGlassFilterDrierMicroMotion.ID=in2m(3.0/8.0)-mm2m(2) #wall thickness is 1mm
#     
    
    #Now solve
    Cycle.PreconditionedSolve()
    #Cycle.PreconditionedSolve_new()
    
    #Print Cycle outputs
    for id, unit, value in Cycle.OutputList():
        print str(id) + ' = ' + str(value) + ' ' + str(unit)
    
    return Cycle

def ECUCycle_VISemiEmpComp():
    #########################################################################
    ######################     CYCLE INITIALIZATION    ######################
    #########################################################################
    
    ## Here we load parameters that are not a function of operating conditions
    ## They are primarily geometric parameters

    #Instantiate the cycle class
    Cycle=ECU_VISemiEmpCompCycleClass()
    
    
    #--------------------------------------
    #--------------------------------------
    #         Cycle parameters
    #--------------------------------------
    #--------------------------------------
    Cycle.Verbosity = 0 #the idea here is to have different levels of debug output 
    Cycle.ImposedVariable = 'Subcooling' #'Subcooling' # or 'Charge'
    Cycle.DT_sc_target = 4.9
    #Cycle.Charge_target = oz2kg(85) #37-44 ounces #kg #uncomment for use with imposed 'Charge'
    Cycle.Mode='AC'
    Cycle.Ref='R407C'
    Cycle.Backend='HEOS' #Backend for refrigerant properties calculation: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    Cycle.TestName='ECU-60K-VI-SemiEmp'  #this and the two next lines can be used to specify exact test conditions
    Cycle.TestDescription='Test1'
    Cycle.TestDetails='This is ECU60K with VISemiEmpComp new system'
    
    
    #--------------------------------------
    #--------------------------------------
    #       VI Semi Empricial Compressor parameters
    #--------------------------------------
    #--------------------------------------
    params={
            'N_comp':3600, #compressor speed [RPM]
            'Ref':Cycle.Ref,
            'Backend':Cycle.Backend,
            'Verbosity': 0, # How verbose should the debugging be [0-10]

            #Optimized Compressor Parameters
            'VR1':1.400305991,  #volume ratio before injection
            'VR2':1.431429505, #volume ratio after injection
            'A_suc':8.99E-05, #main suction port area [m^2]
            'A_inj':6.57E-06, #injection port area [m^2]
            'A_dis':1.30E-05, #discharge port area [m^2]
            'A_leak':4.62E-07, #leakage area [m^2]
            'V_suc_comp':6.68E-05, #suction volume [m^3]
            'T_loss':0.616731207, #frictional torque loss [N-m]
            'eta_m':0.99, # motor efficiency
        
            #Extra parameter not included in the original model
            'fp':0.15, #Fraction of electrical power lost as heat to ambient
            'Vdot_ratio': 1.0, #Displacement Scale factor
            }
    
    Cycle.Compressor.Update(**params)
    
    
    #--------------------------------------
    #--------------------------------------
    #      Condenser parameters
    #--------------------------------------
    #--------------------------------------
    Cycle.Condenser.Fins.Tubes.NTubes=52               #Number of tubes (per bank for now!)
    Cycle.Condenser.Fins.Tubes.Nbank=2                 #Number of banks (set to 1 for now!)
    Cycle.Condenser.Fins.Tubes.Npass=2                 #Number of passes (per bank) #averaged if not even
    Cycle.Condenser.Fins.Tubes.Nports=11               #Number of rectangular ports
    Cycle.Condenser.Fins.Tubes.Ltube=in2m(21.26)       #length of a single tube
    Cycle.Condenser.Fins.Tubes.Td=in2m(1)              #Tube outside width (depth)
    Cycle.Condenser.Fins.Tubes.Ht=in2m(0.072)          #Tube outside height (major diameter)
    Cycle.Condenser.Fins.Tubes.b=in2m(0.488)           #Tube spacing   
    Cycle.Condenser.Fins.Tubes.tw=in2m(0.015)          #Tube wall thickness
    Cycle.Condenser.Fins.Tubes.twp=in2m(0.016)         #Port (channel) wall thickness     
    Cycle.Condenser.Fins.Tubes.beta=1.7675             #Port (channel) aspect ratio (=width/height)
    
    Cycle.Condenser.Fins.Fins.FPI=14                   #Fin per inch
    Cycle.Condenser.Fins.Fins.Lf=in2m(1)               #Fin length = tube outside width in this HX
    Cycle.Condenser.Fins.Fins.t=in2m(0.0045)           ##measured## #Fin thickness
    Cycle.Condenser.Fins.Fins.k_fin=117                #Fin thermal conductivity for pure Aluminum
        
    Cycle.Condenser.Fins.Air.Vdot_ha=cfm2cms(3500)     #Air volume flow rate in m^3/s
    Cycle.Condenser.Fins.Air.Tdb=F2K(125)               #Air inlet temperature, K
    Cycle.Condenser.Fins.Air.p=101325                  #Air pressure in Pa
    Cycle.Condenser.Fins.Air.RH=0.199                 #Air inlet relative humidity
    Cycle.Condenser.Fins.Air.FanPower=984.5              #Fan power, Watts
        
    Cycle.Condenser.Fins.Louvers.Lalpha=25             ##estimated## #Louver angle, in degree
    Cycle.Condenser.Fins.Louvers.lp=mm2m(1.12)         ##measured## #Louver pitch
    
    Cycle.Condenser.Ref=Cycle.Ref
    Cycle.Condenser.Verbosity=0
    params={
            'Ref': Cycle.Ref,
            'Backend': Cycle.Backend,
            'Verbosity':0,
            }
    Cycle.Condenser.Update(**params)
    
    
    #--------------------------------------
    #--------------------------------------
    #     Evaporator Parameters 
    #--------------------------------------
    #--------------------------------------
    Cycle.Evaporator.Fins.Tubes.NTubes_per_bank=18
    Cycle.Evaporator.Fins.Tubes.Nbank=4
    Cycle.Evaporator.Fins.Tubes.Ltube=in2m(24.875)
    Cycle.Evaporator.Fins.Tubes.OD=in2m(0.5)
    Cycle.Evaporator.Fins.Tubes.ID=Cycle.Evaporator.Fins.Tubes.OD - 2*in2m(0.019)
    Cycle.Evaporator.Fins.Tubes.Pl=in2m(1.082)
    Cycle.Evaporator.Fins.Tubes.Pt=in2m(1.25)
    Cycle.Evaporator.Fins.Tubes.Ncircuits=6
    
    Cycle.Evaporator.Fins.Fins.FPI=12
    Cycle.Evaporator.Fins.Fins.Pd=in2m(1.0/16.0/2)
    Cycle.Evaporator.Fins.Fins.xf=in2m(1.0/4.0)
    Cycle.Evaporator.Fins.Fins.t=in2m(0.0075)
    Cycle.Evaporator.Fins.Fins.k_fin=237
    
    Cycle.Evaporator.Fins.Air.Vdot_ha=cfm2cms(1700)          
    Cycle.Evaporator.Fins.Air.Tdb=F2K(90)
    Cycle.Evaporator.Fins.Air.p=101325                       #Evaporator Air pressure in Pa
    Cycle.Evaporator.Fins.Air.RH=0.5023
    Cycle.Evaporator.Fins.Air.FanPower=750.2
    
    Cycle.Evaporator.FinsType = 'WavyLouveredFins'        #WavyLouveredFins, HerringboneFins, PlainFins
    Cycle.Evaporator.Verbosity=0
    
    params={
            'Ref': Cycle.Ref,
            'Backend': Cycle.Backend,
            'Verbosity':0,
            'DT_sh':5,#16.42, #DeltaF2K()
            }
    
    Cycle.Evaporator.Update(**params)
    
    # ----------------------------------
    # ----------------------------------
    #       PHEHX Parameters
    # ----------------------------------
    # ----------------------------------
    params={
            'Ref_c':Cycle.Ref,
            'Backend_c': Cycle.Backend,
            'Ref_h':Cycle.Ref,
            'Backend_h': Cycle.Backend,
            'Verbosity':0,
            'DT_sh_target':5,
            
            #Geometric parameters
            'Bp' : in2m(2.875),
            'Lp' : in2m(18), #Center-to-center distance between ports
            'Nplates' : 10,
            'PlateAmplitude' : 0.001, #[m]
            'PlateThickness' : 0.0003, #[m]
            'PlateWavelength' : 0.00626, #[m]
            'InclinationAngle' : 65/180*pi,#[rad]
            'PlateConductivity' : 15.0, #[W/m-K]
            'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
            }
    Cycle.PHEHX.Update(**params)
        
    # ----------------------------------
    # ----------------------------------
    #       Line Set Return Parameters (vapor line)
    # ----------------------------------
    # ----------------------------------
#     params={
#             'L':in2m(85),
#             'k_tube':0.19,
#             't_insul':0.02,
#             'k_insul':0.036,
#             'T_air':F2K(125),
#             'Ref': Cycle.Ref,
#             'h_air':10, #0.0000000001 is changed to 10 assumed for forced air convection
#             }
#      
#     Cycle.LineSetReturn.Update(**params)
#     Cycle.LineSetReturn.OD=in2m(5.0/8.0)
#     Cycle.LineSetReturn.ID=in2m(5.0/8.0)-mm2m(2) #wall thickness is 1mm

    # ----------------------------------
    # ----------------------------------
    #       Line Set Supply Parameters (liquid line)
    # ----------------------------------
    # ----------------------------------
#     params={
#             'L':in2m(183),                #tube length in m
#             'k_tube':0.19,
#             't_insul':0, #no insulation
#             'k_insul':0.036,
#             'T_air':F2K(125),
#             'Ref': Cycle.Ref,
#             'h_air':10,#0.0000000001 is changed to 10 assumed for forced air convection
#             }
#      
#     Cycle.LineSetSupply.Update(**params)
#     Cycle.LineSetSupply.OD=in2m(3.0/8.0)
#     Cycle.LineSetSupply.ID=in2m(3.0/8.0)-mm2m(2) #wall thickness is 1mm
    
    # ----------------------------------
    # ----------------------------------
    #       Sight Glass + Filter Drier + MicroMotion Parameters
    # ----------------------------------
    # ----------------------------------
#     params={
#             'h':in2m(1.370),        #height of sight glass in m
#             'D':in2m(1.110),        #diameter of sight glass in m
#             'Ref': Cycle.Ref,
#             'V': cubin2cubm(30), #volume of filter drier (website = 13.74in^3 calculated) (manual = 16in^3)
#             'D_Micro': in2m(0.21),  #micromotion tube diameter
#             'L_Micro': in2m(14.6),  #micormotion tube length
#             'n_Micro': 2,           #micormotion number of tubes
#             }
#      
#     Cycle.SightGlassFilterDrierMicroMotion.Update(**params)
#     Cycle.SightGlassFilterDrierMicroMotion.ID=in2m(3.0/8.0)-mm2m(2) #wall thickness is 1mm
#     
    
    #Now solve
    Cycle.PreconditionedSolve()
    #Cycle.PreconditionedSolve_new()
    
    #Print Cycle outputs
    for id, unit, value in Cycle.OutputList():
        print str(id) + ' = ' + str(value) + ' ' + str(unit)
    
    return Cycle

def ECUCycle_VICompTello():
    #########################################################################
    ######################     CYCLE INITIALIZATION    ######################
    #########################################################################
    
    ## Here we load parameters that are not a function of operating conditions
    ## They are primarily geometric parameters

    #Instantiate the cycle class
    Cycle=ECU_VICompTelloCycleClass()
    
    
    #--------------------------------------
    #--------------------------------------
    #         Cycle parameters
    #--------------------------------------
    #--------------------------------------
    Cycle.Verbosity = 0 #the idea here is to have different levels of debug output 
    Cycle.ImposedVariable = 'Subcooling' #'Subcooling' # or 'Charge'
    Cycle.DT_sc_target = 4.9
    #Cycle.Charge_target = oz2kg(85) #37-44 ounces #kg #uncomment for use with imposed 'Charge'
    Cycle.Mode='AC'
    Cycle.Ref='R407C'
    Cycle.Backend='HEOS' #Backend for refrigerant properties calculation: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
    Cycle.TestName='ECU-60K-VI'  #this and the two next lines can be used to specify exact test conditions
    Cycle.TestDescription='Test1'
    Cycle.TestDetails='This is ECU60K with VICompTello new system'
    
    
    #--------------------------------------
    #--------------------------------------
    #       VI Compressor parameters
    #--------------------------------------
    #--------------------------------------
    #Dimensioless compressor map
    if Cycle.Ref=='R407C':
        K=[-0.3845,0.3296]
        M=[133.3171898,0.508718380832,-2.15889885692,0.00847246179835,0.009495493298,0.0170835511659,3.65431994895E-05,6.660136064E-06,-4.719716435E-05,-4.61719969253E-05]
        P=[1.003,0.9998 ,1.134 ,0.9032 ,-0.5003 ,0.5136 ,-0.001366 ,-0.009186 ,0.005756 ,-0.00156 ,1.053] 

    params={
            'M':M,
            'P':P,
            'K':K,
            'Ref':Cycle.Ref, #Refrigerant
            'fp':0.12, #Fraction of electrical power lost as heat to ambient 
            'Vdot_ratio': 1, #Displacement Scale factor
            'Verbosity': 0, # How verbose should the debugging be [0-10]
            'Backend':Cycle.Backend
            }
    
    Cycle.Compressor.Update(**params)
    
    
    #--------------------------------------
    #--------------------------------------
    #      Condenser parameters
    #--------------------------------------
    #--------------------------------------
    Cycle.Condenser.Fins.Tubes.NTubes=52               #Number of tubes (per bank for now!)
    Cycle.Condenser.Fins.Tubes.Nbank=2                 #Number of banks (set to 1 for now!)
    Cycle.Condenser.Fins.Tubes.Npass=2                 #Number of passes (per bank) #averaged if not even
    Cycle.Condenser.Fins.Tubes.Nports=11               #Number of rectangular ports
    Cycle.Condenser.Fins.Tubes.Ltube=in2m(21.26)       #length of a single tube
    Cycle.Condenser.Fins.Tubes.Td=in2m(1)              #Tube outside width (depth)
    Cycle.Condenser.Fins.Tubes.Ht=in2m(0.072)          #Tube outside height (major diameter)
    Cycle.Condenser.Fins.Tubes.b=in2m(0.488)           #Tube spacing   
    Cycle.Condenser.Fins.Tubes.tw=in2m(0.015)          #Tube wall thickness
    Cycle.Condenser.Fins.Tubes.twp=in2m(0.016)         #Port (channel) wall thickness     
    Cycle.Condenser.Fins.Tubes.beta=1.7675             #Port (channel) aspect ratio (=width/height)
    
    Cycle.Condenser.Fins.Fins.FPI=14                   #Fin per inch
    Cycle.Condenser.Fins.Fins.Lf=in2m(1)               #Fin length = tube outside width in this HX
    Cycle.Condenser.Fins.Fins.t=in2m(0.0045)           ##measured## #Fin thickness
    Cycle.Condenser.Fins.Fins.k_fin=117                #Fin thermal conductivity for pure Aluminum
        
    Cycle.Condenser.Fins.Air.Vdot_ha=cfm2cms(3500)     #Air volume flow rate in m^3/s
    Cycle.Condenser.Fins.Air.Tdb=F2K(125)               #Air inlet temperature, K
    Cycle.Condenser.Fins.Air.p=101325                  #Air pressure in Pa
    Cycle.Condenser.Fins.Air.RH=0.199                 #Air inlet relative humidity
    Cycle.Condenser.Fins.Air.FanPower=984.5              #Fan power, Watts
        
    Cycle.Condenser.Fins.Louvers.Lalpha=25             ##estimated## #Louver angle, in degree
    Cycle.Condenser.Fins.Louvers.lp=mm2m(1.12)         ##measured## #Louver pitch
    
    Cycle.Condenser.Ref=Cycle.Ref
    Cycle.Condenser.Verbosity=0
    params={
            'Ref': Cycle.Ref,
            'Backend': Cycle.Backend,
            'Verbosity':0,
            }
    Cycle.Condenser.Update(**params)
    
    
    #--------------------------------------
    #--------------------------------------
    #     Evaporator Parameters 
    #--------------------------------------
    #--------------------------------------
    Cycle.Evaporator.Fins.Tubes.NTubes_per_bank=18
    Cycle.Evaporator.Fins.Tubes.Nbank=4
    Cycle.Evaporator.Fins.Tubes.Ltube=in2m(24.875)
    Cycle.Evaporator.Fins.Tubes.OD=in2m(0.5)
    Cycle.Evaporator.Fins.Tubes.ID=Cycle.Evaporator.Fins.Tubes.OD - 2*in2m(0.019)
    Cycle.Evaporator.Fins.Tubes.Pl=in2m(1.082)
    Cycle.Evaporator.Fins.Tubes.Pt=in2m(1.25)
    Cycle.Evaporator.Fins.Tubes.Ncircuits=6
    
    Cycle.Evaporator.Fins.Fins.FPI=12
    Cycle.Evaporator.Fins.Fins.Pd=in2m(1.0/16.0/2)
    Cycle.Evaporator.Fins.Fins.xf=in2m(1.0/4.0)
    Cycle.Evaporator.Fins.Fins.t=in2m(0.0075)
    Cycle.Evaporator.Fins.Fins.k_fin=237
    
    Cycle.Evaporator.Fins.Air.Vdot_ha=cfm2cms(1700)          
    Cycle.Evaporator.Fins.Air.Tdb=F2K(90)
    Cycle.Evaporator.Fins.Air.p=101325                       #Evaporator Air pressure in Pa
    Cycle.Evaporator.Fins.Air.RH=0.5023
    Cycle.Evaporator.Fins.Air.FanPower=750.2
    
    Cycle.Evaporator.FinsType = 'WavyLouveredFins'        #WavyLouveredFins, HerringboneFins, PlainFins
    Cycle.Evaporator.Verbosity=0
    
    params={
            'Ref': Cycle.Ref,
            'Backend': Cycle.Backend,
            'Verbosity':0,
            'DT_sh':5,#16.42, #DeltaF2K()
            }
    
    Cycle.Evaporator.Update(**params)
    
    # ----------------------------------
    # ----------------------------------
    #       PHEHX Parameters
    # ----------------------------------
    # ----------------------------------
    params={
            'Ref_c':Cycle.Ref,
            'Backend_c': Cycle.Backend,
            'Ref_h':Cycle.Ref,
            'Backend_h': Cycle.Backend,
            'Verbosity':0,
            'DT_sh_target':5,
            
            #Geometric parameters
            'Bp' : in2m(2.875),
            'Lp' : in2m(18), #Center-to-center distance between ports
            'Nplates' : 10,
            'PlateAmplitude' : 0.001, #[m]
            'PlateThickness' : 0.0003, #[m]
            'PlateWavelength' : 0.00626, #[m]
            'InclinationAngle' : 65/180*pi,#[rad]
            'PlateConductivity' : 15.0, #[W/m-K]
            'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
            }
    Cycle.PHEHX.Update(**params)
        
    # ----------------------------------
    # ----------------------------------
    #       Line Set Return Parameters (vapor line)
    # ----------------------------------
    # ----------------------------------
#     params={
#             'L':in2m(85),
#             'k_tube':0.19,
#             't_insul':0.02,
#             'k_insul':0.036,
#             'T_air':F2K(125),
#             'Ref': Cycle.Ref,
#             'h_air':10, #0.0000000001 is changed to 10 assumed for forced air convection
#             }
#      
#     Cycle.LineSetReturn.Update(**params)
#     Cycle.LineSetReturn.OD=in2m(5.0/8.0)
#     Cycle.LineSetReturn.ID=in2m(5.0/8.0)-mm2m(2) #wall thickness is 1mm

    # ----------------------------------
    # ----------------------------------
    #       Line Set Supply Parameters (liquid line)
    # ----------------------------------
    # ----------------------------------
#     params={
#             'L':in2m(183),                #tube length in m
#             'k_tube':0.19,
#             't_insul':0, #no insulation
#             'k_insul':0.036,
#             'T_air':F2K(125),
#             'Ref': Cycle.Ref,
#             'h_air':10,#0.0000000001 is changed to 10 assumed for forced air convection
#             }
#      
#     Cycle.LineSetSupply.Update(**params)
#     Cycle.LineSetSupply.OD=in2m(3.0/8.0)
#     Cycle.LineSetSupply.ID=in2m(3.0/8.0)-mm2m(2) #wall thickness is 1mm
    
    # ----------------------------------
    # ----------------------------------
    #       Sight Glass + Filter Drier + MicroMotion Parameters
    # ----------------------------------
    # ----------------------------------
#     params={
#             'h':in2m(1.370),        #height of sight glass in m
#             'D':in2m(1.110),        #diameter of sight glass in m
#             'Ref': Cycle.Ref,
#             'V': cubin2cubm(30), #volume of filter drier (website = 13.74in^3 calculated) (manual = 16in^3)
#             'D_Micro': in2m(0.21),  #micromotion tube diameter
#             'L_Micro': in2m(14.6),  #micormotion tube length
#             'n_Micro': 2,           #micormotion number of tubes
#             }
#      
#     Cycle.SightGlassFilterDrierMicroMotion.Update(**params)
#     Cycle.SightGlassFilterDrierMicroMotion.ID=in2m(3.0/8.0)-mm2m(2) #wall thickness is 1mm
#     
    
    #Now solve
    Cycle.PreconditionedSolve()
    #Cycle.PreconditionedSolve_new()
    
    #Print Cycle outputs
    for id, unit, value in Cycle.OutputList():
        print str(id) + ' = ' + str(value) + ' ' + str(unit)
    
    return Cycle
    
if __name__=='__main__':
    #cycle=ECUCycle()
    #cycle2=ECUCycle_NewComp()
    #cycle3=ECUCycle_VIComp()
    #cycle4=ECUCycle_VISemiEmpComp()
    cycle5=ECUCycle_VICompTello()
    import time
    start = time.time()
    #Write the outputs to file
    #Write2CSV(cycle4,open('results/Cycle_60K_Test1_VISemiEmpComp.csv','w'),append=False)
    #Write2CSV(cycle2,open('results/Cycle_60K_Test1.csv','a'),append=True)
    #Write2CSV(cycle3,open('results/Cycle_60K_Test1.csv','a'),append=True)
    print 'Took '+str(time.time()-start)+' seconds to run Cycle model'
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
#     cycle.TestDescription='Test#2'
#     cycle.TestDetails='Here we changed the air condition on evaporator and condeser'
#     cycle.PreconditionedSolve()  #there seems to be a problem, somewhere
#     Write2CSV(cycle,open('Cycle.csv','a'),append=True)

    
    #P-h & T-s diagrams for R407C
#     ref_fluid = 'HEOS::R407C'
#     
#     #Experimental results
#     P_exp = [625.3,3213.0,3213.0,3145.0,3145.0,990.4,655.6,625.3,625.3,625.3] #in kPa 
#     P_exp = np.array(P_exp)
#     P_exp *= 1000.0 #convert kPa to Pa
#     T_exp = [25.75+273.15, 111.9+273.15, PropsSI('T','P',P_exp[2],'Q',1,ref_fluid), PropsSI('T','P',P_exp[3],'Q',0,ref_fluid), 61.19+273.15, 20.4+273.15,7.122+273.15,PropsSI('T','P',P_exp[6],'Q',1,ref_fluid), 20.75+273.15, 25.75+273.15] #in Kelvin    
#     T_exp = np.array(T_exp)
#      
#     #Solve for h_exp and s_exp
#     h_exp = [PropsSI('H','P',P_exp[0],'T',T_exp[0],ref_fluid), PropsSI('H','P',P_exp[1],'T',T_exp[1],ref_fluid), 
#              PropsSI('H','P',P_exp[2],'Q',1,ref_fluid), PropsSI('H','P',P_exp[3],'Q',0,ref_fluid),
#              PropsSI('H','P',P_exp[4],'T',T_exp[4],ref_fluid), PropsSI('H','P',P_exp[4],'T',T_exp[4],ref_fluid),PropsSI('H','P',P_exp[4],'T',T_exp[4],ref_fluid),
#              PropsSI('H','P',P_exp[7],'Q',1,ref_fluid), PropsSI('H','P',P_exp[8],'T',T_exp[8],ref_fluid), 
#              PropsSI('H','P',P_exp[9],'T',T_exp[9],ref_fluid)]
#      
#     #Recalculate the temperature at inlet of Evap
#     T_exp[6] = PropsSI('T','P',P_exp[6],'H',h_exp[6],ref_fluid)
#  
#     s_exp = [PropsSI('S','P',P_exp[0],'T',T_exp[0],ref_fluid), PropsSI('S','P',P_exp[1],'T',T_exp[1],ref_fluid), 
#              PropsSI('S','P',P_exp[2],'Q',1,ref_fluid), PropsSI('S','P',P_exp[3],'Q',0,ref_fluid),
#              PropsSI('S','P',P_exp[4],'T',T_exp[4],ref_fluid), PropsSI('S','H',h_exp[5],'P',P_exp[5],ref_fluid),PropsSI('S','H',h_exp[6],'P',P_exp[6],ref_fluid),
#              PropsSI('S','P',P_exp[7],'Q',1,ref_fluid), PropsSI('S','P',P_exp[8],'T',T_exp[8],ref_fluid), 
#              PropsSI('S','P',P_exp[9],'T',T_exp[9],ref_fluid)]
#     h_exp = np.array(h_exp)
#     s_exp = np.array(s_exp)
#      
#              
#     #convert back to original units kPa, kJ/kg and kJ/kg-K
#     P_exp /= 1000.0 
#     T_exp = T_exp #keep T in K
#     h_exp /= 1000.0
#     s_exp /= 1000.0
#     
#     #Model Results
#     P = [cycle.Evaporator.psat_r, cycle.Condenser.psat_r, cycle.Condenser.psat_r, cycle.Condenser.psat_r, cycle.Condenser.psat_r, cycle.Evaporator.psat_r, cycle.Evaporator.psat_r, cycle.Evaporator.psat_r]
#     h = [cycle.Compressor.hin_r, cycle.Compressor.hout_r, PropsSI('H','P',cycle.Condenser.psat_r,'Q',1,ref_fluid), PropsSI('H','P',cycle.Condenser.psat_r,'Q',0,ref_fluid), cycle.Condenser.hout_r, cycle.Evaporator.hin_r, PropsSI('H','P',cycle.Evaporator.psat_r,'Q',1,ref_fluid), cycle.Compressor.hin_r]
#     T = [cycle.Compressor.Tin_r, cycle.Compressor.Tout_r, PropsSI('T','P',cycle.Condenser.psat_r,'Q',1,ref_fluid), PropsSI('T','P',cycle.Condenser.psat_r,'Q',0,ref_fluid), cycle.Condenser.Tout_r, cycle.Evaporator.Tin_r, PropsSI('T','P',cycle.Evaporator.psat_r,'Q',1,ref_fluid), cycle.Compressor.Tin_r]
#     s = [cycle.Compressor.sin_r, cycle.Compressor.sout_r, PropsSI('S','P',cycle.Condenser.psat_r,'Q',1,ref_fluid), PropsSI('S','P',cycle.Condenser.psat_r,'Q',0,ref_fluid), cycle.Condenser.sout_r, cycle.Evaporator.sin_r, PropsSI('S','P',cycle.Evaporator.psat_r,'Q',1,ref_fluid), cycle.Compressor.sin_r]
#     P = np.array(P)
#     T = np.array(T)
#     h = np.array(h)
#     s = np.array(s)
#     P /= 1000.0 #convert Pa to kPa
#     T = T       #keep T in K
#     h /= 1000.0 #convert J/kg to kJ/kg
#     s /= 1000.0 #convert J/kg-K to kJ/kg-K
#     
#     #Model Results (VI comp)
#     cycle = cycle2 # just assign the results of cycle2 (VIComp) to cycle so that I don't have to retype everything again
#     P2 = [cycle.Evaporator.psat_r, cycle.Condenser.psat_r, cycle.Condenser.psat_r, cycle.Condenser.psat_r, cycle.Condenser.psat_r, cycle.Evaporator.psat_r, cycle.Evaporator.psat_r, cycle.Evaporator.psat_r]
#     h2 = [cycle.Compressor.hin_r, cycle.Compressor.hout_r, PropsSI('H','P',cycle.Condenser.psat_r,'Q',1,ref_fluid), PropsSI('H','P',cycle.Condenser.psat_r,'Q',0,ref_fluid), cycle.Condenser.hout_r, cycle.Evaporator.hin_r, PropsSI('H','P',cycle.Evaporator.psat_r,'Q',1,ref_fluid), cycle.Compressor.hin_r]
#     T2 = [cycle.Compressor.Tin_r, cycle.Compressor.Tout_r, PropsSI('T','P',cycle.Condenser.psat_r,'Q',1,ref_fluid), PropsSI('T','P',cycle.Condenser.psat_r,'Q',0,ref_fluid), cycle.Condenser.Tout_r, cycle.Evaporator.Tin_r, PropsSI('T','P',cycle.Evaporator.psat_r,'Q',1,ref_fluid), cycle.Compressor.Tin_r]
#     s2 = [cycle.Compressor.sin_r, cycle.Compressor.sout_r, PropsSI('S','P',cycle.Condenser.psat_r,'Q',1,ref_fluid), PropsSI('S','P',cycle.Condenser.psat_r,'Q',0,ref_fluid), cycle.Condenser.sout_r, cycle.Evaporator.sin_r, PropsSI('S','P',cycle.Evaporator.psat_r,'Q',1,ref_fluid), cycle.Compressor.sin_r]
#     P2 = np.array(P2)
#     T2 = np.array(T2)
#     h2 = np.array(h2)
#     s2 = np.array(s2)
#     P2 /= 1000.0 #convert Pa to kPa
#     T2 = T2       #keep T in K
#     h2 /= 1000.0 #convert J/kg to kJ/kg
#     s2 /= 1000.0 #convert J/kg-K to kJ/kg-K
# 
#     #Plot P-h diagram 
#     ph_plot_R407C = PropertyPlot(ref_fluid, 'Ph')
#     ph_plot_R407C.calc_isolines(CoolProp.iQ, num=2)
#     ph_plot_R407C.title('P-h R407C')
#     ph_plot_R407C.xlabel(r'$h$ [{kJ}/{kg}]')
#     ph_plot_R407C.ylabel(r'$P$ [kPa]')
#     ph_plot_R407C.axis.set_yscale('log')
#     ph_plot_R407C.grid()
#     plt.plot(h_exp,P_exp, 'bo-', markersize=4, linewidth=0.75, label='Experimental')
#     plt.plot(h,P,'ro--', markersize=4, linewidth=0.75, label='Model')
#     plt.plot(h2,P2,'ko--', markersize=4, linewidth=0.75, label='Model VIComp')
#     plt.ylim([100,10**4])
#     plt.xlim([150,550])
#     leg=plt.legend(loc='best',fancybox=False,numpoints=1)
#     frame=leg.get_frame()  
#     frame.set_linewidth(0.5)
#     ph_plot_R407C.savefig('images/60K/R407C_Ph_Test1.pdf')    
#     ph_plot_R407C.show()
# #      
# #     #Plot T-s diagram  
#     ts_plot_R407C = PropertyPlot(ref_fluid, 'Ts')
#     ts_plot_R407C.calc_isolines(CoolProp.iQ, num=2)
#     ts_plot_R407C.title('T-s R407C')
#     ts_plot_R407C.xlabel(r'$s$ [{kJ}/{kg-K}]')
#     ts_plot_R407C.ylabel(r'$T$ [K]')
#     ts_plot_R407C.grid()
#     plt.plot(s_exp,T_exp, 'bo-', markersize=4, linewidth=0.75, label='Experimental')
#     plt.plot(s,T,'ro--', markersize=4, linewidth=0.75, label='Model')
#     plt.plot(s2,T2,'ko--', markersize=4, linewidth=0.75, label='Model VIComp')
#     plt.ylim([200,400])
#     plt.xlim([0.8,2.2])
#     leg=plt.legend(loc='best',fancybox=False, numpoints=1)
#     frame=leg.get_frame()  
#     frame.set_linewidth(0.5)
#     ts_plot_R407C.savefig('images/60K/R407C_Ts_Test1.pdf')    
#     ts_plot_R407C.show()
    
    #Plot T-s and P-h diagrams in one graph
#     fig1 = plt.figure(1, figsize=(16, 8), dpi=100)
#     for i, gtype in enumerate(['Ph', 'Ts']):
#         ax = plt.subplot(1, 2, i+1)
#         if gtype.startswith('P'):
#             ax.set_yscale('log')
#             plt.grid()
#             #plt.plot(h_exp,P_exp, 'bo-', label='Experimental')
#             #plt.errorbar(h_exp,P_exp, yerr=0.08*P_exp)
#             plt.plot(h,P,'ro--', label='Model')
#             plt.legend(loc='best',fancybox=False)
#         if gtype.startswith('T'):
#             plt.grid()
#             #plt.plot(s_exp,T_exp, 'bo-', label='Experimental')
#             #plt.errorbar(s_exp,T_exp, yerr=0.005*T_exp)
#             plt.plot(s,T,'ro--', label='Model')
#             plt.legend(loc='best',fancybox=False)
#         props_plot = PropertyPlot(ref_fluid, gtype, axis=ax, figure=fig1)
#         props_plot.calc_isolines(CoolProp.iQ, num=2)
#         props_plot.title(gtype)
#         props_plot.draw()
#     fig1.set_tight_layout(True)
#     fig1.savefig('images/60K/comined_R407C_60K_Test1.pdf')
#     props_plot.show()