from Compressor import CompressorClass
from CoolProp.CoolProp import PropsSI
from convert_units import C2K

kwds={
      'M':[207.31295549477,4.29717231652541,-2.28622302529118,
           0.0347908258163747,-0.0201167288696277,0.0259153689666968,
           9.13169535150059e-5,-6.23623656573478e-5,1.20986974937733e-4,
           -1.07540716639383e-4],
      'P':[-511.9893727,-1.867619312,32.35057515,-0.0573,0.0718,-0.2478335,
           -0.000762,0.00116,-0.000798,0.00129],
      'Ref':'R407C',
      'Tin_r':C2K(7.556),
      'pin_r':423.4*1000,
      'pout_r':1639*1000,
      'fp':0.1, #Fraction of electrical power lost as heat to ambient
      'Vdot_ratio': 1, #Displacement Scale factor
      'Backend':'TTSE&HEOS' #choose between: 'HEOS','TTSE&HEOS','BICUBIC&HEOS','REFPROP','SRK','PR'
      }
Comp=CompressorClass(**kwds)
Comp.Calculate()

print 'Electrical power is: ' + str(Comp.W) + ' W'
print 'Actual mass flow rate is: ' + str(Comp.mdot_r) + ' kg/s'
print 'Isentropic Efficiency is: ' + str(Comp.eta_oi)
print 'Discharge Refrigerant Temperature is: ' + str(Comp.Tout_r) + ' K'
print ' '

'''to print all the output, uncomment the next 2 lines'''
for id, unit, value in Comp.OutputList():
    print str(id) + ' = ' + str(value) + ' ' + str(unit)
