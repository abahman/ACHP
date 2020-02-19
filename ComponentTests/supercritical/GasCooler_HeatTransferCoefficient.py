import matplotlib,os
#matplotlib.use('GTKAgg')
import sys, os
#from FileIO import prep_csv2rec as prep
from matplotlib.mlab import csv2rec

import pylab
import numpy as np
import shutil
from scipy import polyval, polyfit
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.style.use('classic')
mpl.style.use('Elsevier.mplstyle')
mpl.rcParams['mathtext.fontset'] = 'custom'

#===============================================================================
# Latex render
#===============================================================================
# #mpl.use('pgf')
# 
# def figsize(scale):
#     fig_width_pt = 469.755                          # Get this from LaTeX using \the\textwidth
#     inches_per_pt = 1.0/72.27                       # Convert pt to inch
#     golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
#     fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
#     fig_height = fig_width*golden_mean              # height in inches
#     fig_size = [fig_width,fig_height]
#     return fig_size
# 
# pgf_with_latex = {                      # setup matplotlib to use latex for output
# "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
# "text.usetex": True,                # use LaTeX to write all text
# "font.family": "serif",
# "mathtext.fontset" : "custom",
# "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
# "font.sans-serif": [],
# "font.monospace": [],
# "axes.labelsize": 10,               # LaTeX default is 10pt font.
# "font.size": 10,
# "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
# "legend.labelspacing":0.2,
# "xtick.labelsize": 8,
# "ytick.labelsize": 8,
# "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
# "pgf.preamble": [
# r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
# r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
#         ]
#     }
# mpl.rcParams.update(pgf_with_latex)
#===============================================================================
# END of Latex render
#===============================================================================


def rmse(predictions, targets):
    '''
    Root Mean Square Error
    '''
    n = len(predictions)
    RMSE = np.linalg.norm(predictions - targets) / np.sqrt(n) / np.mean(targets) * 100
    return RMSE

def mape(y_pred, y_true):  #maps==mean_absolute_percentage_error
    '''
    Mean Absolute Percentage Error
    '''
    MAPE = np.mean(np.abs((y_true - y_pred) / y_true)) * 100
    return MAPE


    ######RMSE#######
    #rmse_mass = rmse(yuanpei_m_dot_inj_norm,m_dot_inj_norm_exp)
    #print("rmse_mass error is: " + str(rmse_mass) + " %")
    ######MAPE######
    # mape_mass = mape(yuanpei_m_dot_inj_norm,m_dot_inj_norm_exp)
    # print("mape_mass error is: " + str(mape_mass) + " %")
    
from ACHP.Correlations import Petterson_supercritical_average
import CoolProp as CP
from math import pi

AS = CP.AbstractState('HEOS', 'R744')
Tout_r = np.linspace(25+273, 117+273, num=100)
Tin_r = 120+273
T_w = 35+273
psat_r = 10*1000*1000 #10 MPa

mdot_r = [0.02, 0.05, 0.1, 0.2]
Ncircuits = 1
ID = 7.5/1000
Ltube = 0.61
NTubes_per_bank = 18
Nbank = 3
QoverA = 8000 #W/m^2 

TotalLength=Ltube*NTubes_per_bank*Nbank
Lcircuit=TotalLength/Ncircuits
G_r_1 = mdot_r[0]/(Ncircuits*pi*ID**2/4.0)
G_r_2 = mdot_r[1]/(Ncircuits*pi*ID**2/4.0)
G_r_3 = mdot_r[2]/(Ncircuits*pi*ID**2/4.0)
G_r_4 = mdot_r[3]/(Ncircuits*pi*ID**2/4.0)
DoverL = ID/Lcircuit

SC_1 = np.array([]) #empty score array
SC_2 = np.array([]) #empty score array
SC_3 = np.array([]) #empty score array
SC_4 = np.array([]) #empty score array
for i in range(len(Tout_r)):
    h_r_1 = Petterson_supercritical_average(Tout_r[i], Tin_r, T_w, AS, G_r_1, ID, 0, DoverL, mdot_r[0]/Ncircuits, psat_r, QoverA)[0]
    h_r_2 = Petterson_supercritical_average(Tout_r[i], Tin_r, T_w, AS, G_r_2, ID, 0, DoverL, mdot_r[1]/Ncircuits, psat_r, QoverA)[0]
    h_r_3 = Petterson_supercritical_average(Tout_r[i], Tin_r, T_w, AS, G_r_3, ID, 0, DoverL, mdot_r[2]/Ncircuits, psat_r, QoverA)[0]
    h_r_4 = Petterson_supercritical_average(Tout_r[i], Tin_r, T_w, AS, G_r_4, ID, 0, DoverL, mdot_r[3]/Ncircuits, psat_r, QoverA)[0]
    #print (h_r_1,h_r_2,h_r_3,h_r_4)
    SC_1 = np.append(SC_1,h_r_1)
    SC_2 = np.append(SC_2,h_r_2)
    SC_3 = np.append(SC_3,h_r_3)
    SC_4 = np.append(SC_4,h_r_4)
    
#########################
##### plot #######
######################### 
fig=pylab.figure(figsize=(6,4))
x1 = Tout_r-273
y1 = SC_1
y2 = SC_2
y3 = SC_3
y4 = SC_4
plt.plot(x1,y1,':r',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'$\dot m_{r} = 0.02$ kg s$^{-1}$')
plt.plot(x1,y2,'-.b',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'$\dot m_{r} = 0.05$ kg s$^{-1}$')
plt.plot(x1,y3,'--g',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'$\dot m_{r} = 0.1$ kg s$^{-1}$')
plt.plot(x1,y4,'-k',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'$\dot m_{r} = 0.2$ kg s$^{-1}$')

# plt.ylim(0,400)
# plt.xlim(25,120)
plt.xlabel(r'$T_{r}$ [$\degree$C]')
plt.ylabel(r'$\alpha_{avg}$ [W m$^{-2}$ K$^{-1}$]')
leg = plt.legend(loc='best',fancybox=False,numpoints=1)
frame  = leg.get_frame()  
frame.set_linewidth(1.0)
plt.tight_layout(pad=0.2)
plt.savefig('Heat_transfer_coefficient_avg.pdf')
plt.show()
plt.close()