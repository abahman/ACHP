import matplotlib as mpl
import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt
import scipy.interpolate
mpl.style.use('classic')
mpl.style.use('Elsevier.mplstyle')
mpl.rcParams['mathtext.fontset'] = 'custom'


# import matplotlib.pyplot as plt
# from CoolProp.Plots.ConsistencyPlots import ConsistencyFigure
# ff = ConsistencyFigure('R744')
# plt.show()
# plt.close()
# del ff


ref_fluid = CP.AbstractState("HEOS", "R744")
pc = ref_fluid.keyed_output(CP.iP_critical)
Tc = ref_fluid.keyed_output(CP.iT_critical)
ptriple = ref_fluid.keyed_output(CP.iP_triple)
Tmin = 200
Tmax = 1000
pmax = ref_fluid.keyed_output(CP.iP_max)
pt = ref_fluid.keyed_output(CP.iP_triple)
Tt = ref_fluid.keyed_output(CP.iT_triple)
fillcolor = 'g'

fig = plt.figure(figsize = (5,5))
ax = fig.add_subplot(111)
lw = 3

# --------------
# Melting curve
# --------------
melt_args = dict(lw = lw, solid_capstyle = 'round')
TT = []
PP = list(np.logspace(np.log10(pt), np.log10(pmax),1000))
for p in PP:
    TT.append(ref_fluid.melting_line(CP.iT, CP.iP, p))

#Zone VI
for T in np.linspace(max(TT), 350):
    TT.append(T)
    theta = T/273.31
    pi = 1-1.07476*(1-theta**4.6)
    p = pi*632.4e6
    PP.append(p)
 
PPnew = [p/1000 for p in PP]
TTnew = [t-273.15 for t in TT]
plt.plot(TTnew,PPnew,'darkblue',**melt_args)

# ----------------
# Saturation curve
# ----------------
Ts = np.linspace(Tt, Tc, 1000)
ps = CP.CoolProp.PropsSI('P','T',Ts,'Q',0,'R744')

# ------
# Labels
# ------

plt.plot(Ts-273.15,ps/1000,'orange',lw = lw, solid_capstyle = 'round')

# Critical lines
plt.axvline(Tc-273.15, dashes = [2, 2])
plt.axhline(pc/1000, dashes = [2, 2])

# Labels
plt.text(450-273.15, 1e8/1000, 'supercritical',ha= 'center')
plt.text(375-273.15, 2e6/1000, 'supercritical vapor', rotation = 0)
plt.text(275-273.15, 1.5e8/1000, 'supercritical liquid', rotation = 90, ha = 'center')
plt.text(230-273.15, 4.5e6/1000, 'liquid', rotation = 45)
plt.text(255-273.15, 1.5e6/1000, 'vapor', rotation = 45)

plt.ylim(ptriple/1000,1e6)
plt.gca().set_yscale('log')
plt.gca().set_xlim(-70, 300)
plt.ylabel('Pressure [kPa]')
plt.xlabel('Temperature [$\degree$C]')
plt.tight_layout()
plt.savefig('P-T.pdf')
plt.show()