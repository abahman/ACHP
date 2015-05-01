'''
Created on Apr 29, 2015

@author: AmmarBahman
'''

'''This code is for Direct Expansion in Cooling Mode'''

from Cycle import ECU_DXCycleClass 
from convert_units import in2m, mm2m, cfm2cms, F2K, kPa2Pa, C2K

#Instantiate the cycle class
Cycle=ECU_DXCycleClass()

#--------------------------------------
# Cycle parameters
#--------------------------------------
Cycle.Verbosity = 0 #the idea here is to have different levels of debug output 
Cycle.ImposedVariable = 'Subcooling' # or 'Charge'
Cycle.DT_sc_target = 7.0
#Cycle.Charge_target = 1 #kg #uncomment for use with imposed 'Charge'
Cycle.Mode='AC'
Cycle.Ref='R407C'

#--------------------------------------
#       Compressor parameters
#--------------------------------------
#A 3 ton cooling capacity compressor map
M=[217.3163128,5.094492028,-0.593170311,4.38E-02,-2.14E-02,1.04E-02,7.90E-05,-5.73E-05,1.79E-04,-8.08E-05]
P=[-561.3615705,-15.62601841,46.92506685,-0.217949552,0.435062616,-0.442400826,2.25E-04,2.37E-03,-3.32E-03,2.50E-03]

params={
        'M':M,
        'P':P,
        'Ref':Cycle.Ref, #Refrigerant
        'fp':0.2, #Fraction of electrical power lost as heat to ambient 
        'Vdot_ratio': 1.0, #Displacement Scale factor
        'Verbosity': 0, # How verbose should the debugging be [0-10]
        }

Cycle.Compressor.Update(**params)

#--------------------------------------
#      Condenser parameters
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
Cycle.Condenser.Fins.Air.Tmean=F2K(125) 
Cycle.Condenser.Fins.Air.Tdb=F2K(125)              #Air inlet temperature, K
Cycle.Condenser.Fins.Air.p=101325                  #Air pressure in Pa
Cycle.Condenser.Fins.Air.RHmean=0.199
Cycle.Condenser.Fins.Air.RH=0.199                  #Air inlet relative humidity
Cycle.Condenser.Fins.Air.FanPower=855              #Fan power, Watts
    
Cycle.Condenser.Fins.Louvers.Lalpha=25             ##estimated## #Louver angle, in degree
Cycle.Condenser.Fins.Louvers.lp=mm2m(1.12)         ##measured## #Louver pitch

Cycle.Condenser.Ref=Cycle.Ref
Cycle.Condenser.Verbosity=0

#--------------------------------------
# Evaporator Parameters 
#--------------------------------------
Cycle.Evaporator.Fins.Tubes.NTubes_per_bank=32
Cycle.Evaporator.Fins.Tubes.Nbank=3
Cycle.Evaporator.Fins.Tubes.Ltube=0.452
Cycle.Evaporator.Fins.Tubes.OD=0.00913
Cycle.Evaporator.Fins.Tubes.ID=0.00849
Cycle.Evaporator.Fins.Tubes.Pl=0.0191
Cycle.Evaporator.Fins.Tubes.Pt=0.0254
Cycle.Evaporator.Fins.Tubes.Ncircuits=5

Cycle.Evaporator.Fins.Fins.FPI=14.5
Cycle.Evaporator.Fins.Fins.Pd=0.001
Cycle.Evaporator.Fins.Fins.xf=0.001
Cycle.Evaporator.Fins.Fins.t=0.00011
Cycle.Evaporator.Fins.Fins.k_fin=237

Cycle.Evaporator.Fins.Air.Vdot_ha=0.56319
Cycle.Evaporator.Fins.Air.Tmean=297.039
Cycle.Evaporator.Fins.Air.Tdb=297.039
Cycle.Evaporator.Fins.Air.p=101325                                              #Evaporator Air pressure in Pa
Cycle.Evaporator.Fins.Air.RH=0.5
Cycle.Evaporator.Fins.Air.RHmean=0.5
Cycle.Evaporator.Fins.Air.FanPower=438

Cycle.Evaporator.FinsType = 'WavyLouveredFins'        #WavyLouveredFins, HerringboneFins, PlainFins
Cycle.Evaporator.Ref=Cycle.Ref
Cycle.Evaporator.Verbosity=0
Cycle.Evaporator.DT_sh=5

# ----------------------------------
#       Line Set Parameters
# ----------------------------------
# params={
#         'L':7.6,
#         'k_tube':0.19,
#         't_insul':0.02,
#         'k_insul':0.036,
#         'T_air':297,
#         'Ref': Cycle.Ref,
#         'h_air':0.0000000001,
#         }
# 
# Cycle.LineSetSupply.Update(**params)
# Cycle.LineSetReturn.Update(**params)
# Cycle.LineSetSupply.OD=0.009525
# Cycle.LineSetSupply.ID=0.007986
# Cycle.LineSetReturn.OD=0.01905
# Cycle.LineSetReturn.ID=0.017526


#Now solve
from time import time
t1=time()
Cycle.PreconditionedSolve()
print 'Took '+str(time()-t1)+' seconds to run Cycle model'
print 'Cycle COP is '+str(Cycle.COSP)
print 'Cycle refrigerant charge is '+str(Cycle.Charge)+' kg'