'''
Created on Apr 29, 2015

@author: AmmarBahman
'''

'''This code is for Direct Expansion in Cooling Mode of ECU 18K'''

from Cycle import ECU_DXCycleClass 
from convert_units import in2m, mm2m, cm2m, cfm2cms, F2K, kPa2Pa, C2K, oz2kg, DeltaF2K, cubin2cubm
from ACHPTools import Write2CSV


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
    Cycle.ImposedVariable = 'Subcooling' # or 'Charge'
    Cycle.DT_sc_target = 4.638
    #Cycle.Charge_target = oz2kg(37) #37-44 ounces #kg #uncomment for use with imposed 'Charge'
    Cycle.Mode='AC'
    Cycle.Ref='R407C'
    Cycle.TestName='ECU-18K'  #this and the two next lines can be used to specify exact test conditions
    Cycle.TestDescription='Test#1'
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
            'fp':0.15, #Fraction of electrical power lost as heat to ambient 
            'Vdot_ratio': 0.85, #Displacement Scale factor
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
        
    Cycle.Condenser.Fins.Air.Vdot_ha=cfm2cms(1500)     #Air volume flow rate in m^3/s
    Cycle.Condenser.Fins.Air.Tdb=F2K(125)               #Air inlet temperature, K
    Cycle.Condenser.Fins.Air.p=101325                  #Air pressure in Pa
    Cycle.Condenser.Fins.Air.RH=0.199                 #Air inlet relative humidity
    Cycle.Condenser.Fins.Air.FanPower=854.9              #Fan power, Watts
        
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
    
    Cycle.Evaporator.Fins.Air.Vdot_ha=cfm2cms(600)
    Cycle.Evaporator.Fins.Air.Tdb=F2K(90)
    Cycle.Evaporator.Fins.Air.p=101325                                              #Evaporator Air pressure in Pa
    Cycle.Evaporator.Fins.Air.RH=0.5023
    Cycle.Evaporator.Fins.Air.FanPower=392.1
    
    Cycle.Evaporator.FinsType = 'WavyLouveredFins'        #WavyLouveredFins, HerringboneFins, PlainFins
    Cycle.Evaporator.Ref=Cycle.Ref
    Cycle.Evaporator.Verbosity=0
    
    params={
            'Ref': Cycle.Ref,
            'Verbosity':0,
            'DT_sh':6.051, #DeltaF2K()
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
            'T_air':F2K(125),
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
            'T_air':F2K(125),
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
    Write2CSV(cycle,open('Cycle.csv','w'),append=False)
    #Write2CSV(cycle,open('Cycle.csv','a'),append=True)
    
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
