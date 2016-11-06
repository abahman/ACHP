'''
Created on Mar 21, 2015

@author: AmmarBahman
'''

from CoolProp.CoolProp import PropsSI
from FinCorrelations import FinInputs
from Evaporator import EvaporatorClass
from convert_units import in2m, cfm2cms, F2K,kPa2Pa, C2K

FinsTubes=FinInputs()

FinsTubes.Tubes.NTubes_per_bank=9
FinsTubes.Tubes.Ncircuits=5
FinsTubes.Tubes.Nbank=5
FinsTubes.Tubes.Ltube=in2m(19)
FinsTubes.Tubes.OD=in2m(0.31)
FinsTubes.Tubes.ID=FinsTubes.Tubes.OD - 2*in2m(0.016)
FinsTubes.Tubes.Pl=in2m(0.645669)
FinsTubes.Tubes.Pt=in2m(1)

FinsTubes.Fins.FPI=14.5
FinsTubes.Fins.Pd=in2m(1.0/16.0)
FinsTubes.Fins.xf=in2m(1.0/4.0/2.0)
FinsTubes.Fins.t=in2m(0.006)
FinsTubes.Fins.k_fin=237

FinsTubes.Air.Vdot_ha=cfm2cms(400)  
FinsTubes.Air.Tdb=F2K(77)
FinsTubes.Air.p=101325
FinsTubes.Air.RH=0.2712
FinsTubes.Air.FanPower=401.9

kwargs={'Ref': 'R407C',
        'mdot_r': 0.0335,
        'psat_r': kPa2Pa(423.4),
        'Fins': FinsTubes,
        'FinsType': 'WavyLouveredFins',                                         #WavyLouveredFins, HerringboneFins, PlainFins
        'hin_r': PropsSI('H','P',1624000,'T',C2K(32.18),'R407C'), #[J/kg]
        'Verbosity': 0,
        'Backend':'TTSE&HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
        }

Evap=EvaporatorClass(**kwargs)
Evap.Update(**kwargs)
Evap.Calculate()

print 'Evaporator heat transfer rate is',Evap.Q,'W'
print 'Evaporator capacity (less fan power) is',Evap.Capacity,'W'
print 'Evaporator fraction of length in two-phase section',Evap.w_2phase,'W'
print 'Evaporator sensible heat ratio',Evap.SHR
for id, unit, value in Evap.OutputList():
    print str(id) + ' = ' + str(value) + ' ' + str(unit)
