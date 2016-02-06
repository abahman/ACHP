'''
Created on Sep 10, 2015

@author: ammarbahman

Note: you need to have all the test results CSV files executed before run this file.
It will update all the plots based on the latest results from the CSV file.

'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import csv2rec


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

#Experimental Results
T_env = np.array([125,115,105,95,85,75])
m_dot_exp = np.array([0.0447,0.0366,0.0349,0.0379,0.0315,0.0335])
capacity_exp = np.array([5578,5014,5150,5814,5133,5649])
total_power_exp = np.array([4296,3857,3530,3329,2852,2654])
compressor_power_exp = np.array([3049,2591,2252,2030,1779,1582])
COPS_exp = np.array([1.298,1.300,1.459,1.746,1.800,2.129])
charge_exp = np.array([1.2,1.2,1.2,1.2,1.2,1.2])
#Import data from CSV file
data1 = csv2rec('results/Cycle_Test#1.csv',delimiter=',')
data2 = csv2rec('results/Cycle_Test#2.csv',delimiter=',')
data3 = csv2rec('results/Cycle_Test#3.csv',delimiter=',')
data4 = csv2rec('results/Cycle_Test#4.csv',delimiter=',')
data5 = csv2rec('results/Cycle_Test#5.csv',delimiter=',')
data6 = csv2rec('results/Cycle_Test#6.csv',delimiter=',')
#Arrange data in Numpy array for the 6 different tests
m_dot = np.array([data1[2][14],data2[2][14],data3[2][14],data4[2][14],data5[2][14],data6[2][14]])
capacity = np.array([data1[2][11],data2[2][11],data3[2][11],data4[2][11],data5[2][11],data6[2][11]])
total_power = np.array([data1[2][12],data2[2][12],data3[2][12],data4[2][12],data5[2][12],data6[2][12]])
compressor_power = np.array([data1[2][13],data2[2][13],data3[2][13],data4[2][13],data5[2][13],data6[2][13]])
COPS = np.array([data1[2][10],data2[2][10],data3[2][10],data4[2][10],data5[2][10],data6[2][10]])
charge = np.array([data1[2][3],data2[2][3],data3[2][3],data4[2][3],data5[2][3],data6[2][3]])
#to convert string array to integer array
m_dot = m_dot.astype(np.float)
capacity = capacity.astype(np.float)
total_power = total_power.astype(np.float)
compressor_power = compressor_power.astype(np.float)
COPS = COPS.astype(np.float)
charge = charge.astype(np.float)
#plots
#Plot mass flow rate comparison
plt.plot(T_env,m_dot_exp,'-ob',label='Experimental')
plt.errorbar(T_env,m_dot_exp, yerr=0.002*m_dot_exp)
plt.plot(T_env,m_dot,'--or',label='Model')
plt.ylim(0.025,0.055)
plt.xlim(70,130)
plt.legend(loc='best',fancybox=False)
plt.xlabel(r'$T_{env}$ $[\degree\mathrm{F}]$')
plt.ylabel(r'$\dot m_r$ $[\mathrm{kg/s}]$')
#plt.title('Mass flowrate Comparison')
plt.savefig('images/comparison_massflow.pdf')
plt.show()
#Plot Capacity comparison
plt.plot(T_env,capacity_exp,'-ob',label=r'Experimental')
plt.errorbar(T_env,capacity_exp, yerr=0.00905522*capacity_exp)
plt.plot(T_env,capacity,'--or',label=r'Model')
plt.ylim(4000,7000)
plt.xlim(70,130)
plt.legend(loc='best',fancybox=False)
plt.xlabel(r'$T_{env}$ $[\degree\mathrm{F}]$')
plt.ylabel(r'$\dot Q_{evap}$ $[\mathrm{W}]$')
plt.title(r'Capacity Comparison')
plt.savefig('images/comparison_capacity.pdf')
plt.show()
#Plot total power comparison
plt.plot(T_env,total_power_exp,'-ob',label='Experimental')
plt.errorbar(T_env,total_power_exp, yerr=0.02618715*total_power_exp)
plt.plot(T_env,total_power,'--or',label='Model')
plt.ylim(2000,5000)
plt.xlim(70,130)
plt.legend(loc='best',fancybox=False)
plt.xlabel(r'$T_{env}$ $[\degree\mathrm{F}]$')
plt.ylabel(r'$\dot E_t$ $[\mathrm{W}]$')
#plt.title('Total Power Comparison')
plt.savefig('images/comparison_total_power.pdf')
plt.show()
#Plot compressor power comparison
plt.plot(T_env,compressor_power_exp,'-ob',label='Experimental')
plt.errorbar(T_env,compressor_power_exp, yerr=112.5)
plt.plot(T_env,compressor_power,'--or',label='Model')
plt.ylim(1000,3500)
plt.xlim(70,130)
plt.legend(loc='best',fancybox=False)
plt.xlabel(r'$T_{env}$ $[\degree\mathrm{F}]$')
plt.ylabel(r'$\dot W_{comp}$ $[\mathrm{W}]$')
#plt.title('Compressor Power Comparison')
plt.savefig('images/comparison_compressor_power.pdf')
plt.show()
#Plot COPS comparison
plt.plot(T_env,COPS_exp,'-ob',label='Experimental')
plt.errorbar(T_env,COPS_exp, yerr=0.02772727*COPS_exp)
plt.plot(T_env,COPS,'--or',label='Model')
plt.ylim(1,2.4)
plt.xlim(70,130)
plt.legend(loc='best',fancybox=False)
plt.xlabel(r'$T_{env}$ $[\degree\mathrm{F}]$')
plt.ylabel(r'$\mathrm{COP}_{sys}$')
#plt.title('System COP Comparison')
plt.savefig(r'images/comparison_COPS.pdf')
plt.show()
#Plot charge comparison
plt.plot(T_env,charge_exp,'-ob',label='Experimental')
plt.errorbar(T_env,charge_exp, yerr=0.0)
plt.plot(T_env,charge,'--or',label='Model')
plt.ylim(0,1.6)
plt.xlim(70,130)
plt.legend(loc='best',fancybox=False)
plt.xlabel(r'$T_{env}$ $[\degree\mathrm{F}]$')
plt.ylabel(r'$\mathrm{Charge}$ $[\mathrm{kg}]$')
#plt.title('System charge Comparison')
plt.savefig('images/comparison_charge.pdf')
plt.show()
#Combine
fig = plt.figure(1, figsize=(10, 10), dpi=100)
for i, gtype in enumerate(['Mass', 'Capacity', 'Power', 'Compressor', 'COPS','Charge']):
    ax = plt.subplot(3, 2, i+1)
    if gtype.startswith('Mass'):
        plt.plot(T_env,m_dot_exp,'-ob',label='Experimental')
        plt.errorbar(T_env,m_dot_exp, yerr=0.002*m_dot_exp)
        plt.plot(T_env,m_dot,'--or',label='Model')
        plt.ylim(0.02,0.05)
        plt.xlim(70,130)
        plt.legend(loc='best',fancybox=False)
        plt.xlabel(r'$T_{env}$ $[\degree\mathrm{F}]$')
        plt.ylabel(r'$\dot m_r$ $[\mathrm{kg/s}]$')
        #plt.title('Mass flowrate Comparison')
    if gtype.startswith('Capacity'):
        plt.plot(T_env,capacity_exp,'-ob',label='Experimental')
        plt.errorbar(T_env,capacity_exp, yerr=0.00905522*capacity_exp)
        plt.plot(T_env,capacity,'--or',label='Model')
        plt.ylim(4000,7000)
        plt.xlim(70,130)
        plt.legend(loc='best',fancybox=False)
        plt.xlabel(r'$T_{env}$ $[\degree\mathrm{F}]$')
        plt.ylabel(r'$\dot Q_{evap}$ $[\mathrm{W}]$')
        #plt.title('Capacity Comparison')
    if gtype.startswith('Power'):
        plt.plot(T_env,total_power_exp,'-ob',label='Experimental')
        plt.errorbar(T_env,total_power_exp, yerr=0.02618715*total_power_exp)
        plt.plot(T_env,total_power,'--or',label='Model')
        plt.ylim(2000,5000)
        plt.xlim(70,130)
        plt.legend(loc='best',fancybox=False)
        plt.xlabel(r'$T_{env}$ $[\degree\mathrm{F}]$')
        plt.ylabel(r'$\dot E_t$ $[\mathrm{W}]$')
        #plt.title('Total Power Comparison')
    if gtype.startswith('Compressor'):
        plt.plot(T_env,compressor_power_exp,'-ob',label='Experimental')
        plt.errorbar(T_env,compressor_power_exp, yerr=112.5)
        plt.plot(T_env,compressor_power,'--or',label='Model')
        plt.ylim(1000,3500)
        plt.xlim(70,130)
        plt.legend(loc='best',fancybox=False)
        plt.xlabel(r'$T_{env}$ $[\degree\mathrm{F}]$')
        plt.ylabel(r'$\dot W_{comp}$ $[\mathrm{W}]$')
        #plt.title('Compressor Power Comparison')
    if gtype.startswith('COPS'):
        plt.plot(T_env,COPS_exp,'-ob',label='Experimental')
        plt.errorbar(T_env,COPS_exp, yerr=0.02772727*COPS_exp)
        plt.plot(T_env,COPS,'--or',label='Model')
        plt.ylim(1,2.4)
        plt.xlim(70,130)
        plt.legend(loc='best',fancybox=False)
        plt.xlabel(r'$T_{env}$ $[\degree\mathrm{F}]$')
        plt.ylabel(r'$\mathrm{COP}_{sys}$')
        #plt.title('System COP Comparison')
    if gtype.startswith('Charge'):
        plt.plot(T_env,charge_exp,'-ob',label='Experimental')
        plt.errorbar(T_env,charge_exp, yerr=0.0)
        plt.plot(T_env,charge,'--or',label='Model')
        plt.ylim(0,1.6)
        plt.xlim(70,130)
        plt.legend(loc='best',fancybox=False)
        plt.xlabel(r'$T_{env}$ $[\degree\mathrm{F}]$')
        plt.ylabel(r'$\mathrm{Charge}$ $[\mathrm{kg}]$')
        #plt.title('System charge Comparison')
fig.set_tight_layout(True)
plt.savefig('images/comined_comparison.pdf')
plt.show()