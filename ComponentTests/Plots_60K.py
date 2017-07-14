'''
Created on Jul 11, 2017

@author: ammarbahman


'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import csv2rec
import matplotlib as mpl
mpl.style.use('classic')

#===============================================================================
# Latex render
#===============================================================================
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

#Experimental Results
TestNo = np.arange(1,9,1)
m_dot_exp = np.array([0.1049,0.08967,0.08974,0.09389,0.08051,0.08206,0.09311,0.08494]) #[kg/s]
m_dot_tot_exp = np.array([0.1314,0.1134,0.1104,0.1106,0.09419,0.09385,0.1051,0.09789]) #[kg/s]
m_dot_inj_exp = m_dot_tot_exp - m_dot_exp #[kg/s]
capacity_exp = np.array([17.72,16.06,16.5,17.51,15.69,16.09,17.92,16.48]) #[kW]
total_power_exp = np.array([9.097,8.042,7.374,6.795,6.076,5.391,6.125,6.033]) #[kW]
compressor_power_exp = np.array([7.353,6.286,5.544,4.969,4.205,4.16,4.274,4.129]) #[kW]
COPS_exp = np.array([1.948,1.997,2.237,2.576,2.581,2.985,2.927,2.732]) #[-]
charge_exp = np.array([5.01,5.01,5.01,5.01,5.01,5.01,5.01,5.01]) #[kg]
#Import data from CSV file
data1 = csv2rec('results/Cycle_60K_superheat_Test1.csv',delimiter=',')
data2 = csv2rec('results/Cycle_60K_superheat_Test2.csv',delimiter=',')
data3 = csv2rec('results/Cycle_60K_superheat_Test3.csv',delimiter=',')
data4 = csv2rec('results/Cycle_60K_superheat_Test4.csv',delimiter=',')
data5 = csv2rec('results/Cycle_60K_superheat_Test5.csv',delimiter=',')
data6 = csv2rec('results/Cycle_60K_superheat_Test6.csv',delimiter=',')
dataB = csv2rec('results/Cycle_60K_superheat_TestB.csv',delimiter=',')
dataC = csv2rec('results/Cycle_60K_superheat_TestC.csv',delimiter=',')
#Arrange data in Numpy array for the 8 different tests
m_dot = np.array([data1[5][17],data2[5][17],data3[5][17],data4[5][17],data5[5][17],data6[5][17],dataB[5][17],dataC[5][17]])
m_dot_inj = np.array([data1[5][18],data2[5][18],data3[5][18],data4[5][18],data5[5][18],data6[5][18],dataB[5][18],dataC[5][18]])
capacity = np.array([data1[5][14],data2[5][14],data3[5][14],data4[5][14],data5[5][14],data6[5][14],dataB[5][14],dataC[5][14]])
total_power = np.array([data1[5][15],data2[5][15],data3[5][15],data4[5][15],data5[5][15],data6[5][15],dataB[5][15],dataC[5][15]])
compressor_power = np.array([data1[5][16],data2[5][16],data3[5][16],data4[5][16],data5[5][16],data6[5][16],dataB[5][16],dataC[5][16]])
COPS = np.array([data1[5][13],data2[5][13],data3[5][13],data4[5][13],data5[5][13],data6[5][13],dataB[5][13],dataC[5][13]])
charge = np.array([data1[5][3],data2[5][3],data3[5][3],data4[5][3],data5[5][3],data6[5][3],dataB[5][3],dataC[5][3]])
# m_dot = np.array([data1[2][17],data2[2][17],data3[2][17],data4[2][17],data5[2][17],data6[2][17],dataB[2][17],dataC[2][17]])
# m_dot_inj = np.array([data1[2][18],data2[2][18],data3[2][18],data4[2][18],data5[2][18],data6[2][18],dataB[2][18],dataC[2][18]])
# capacity = np.array([data1[2][14],data2[2][14],data3[2][14],data4[2][14],data5[2][14],data6[2][14],dataB[2][14],dataC[2][14]])
# total_power = np.array([data1[2][15],data2[2][15],data3[2][15],data4[2][15],data5[2][15],data6[2][15],dataB[2][15],dataC[2][15]])
# compressor_power = np.array([data1[2][16],data2[2][16],data3[2][16],data4[2][16],data5[2][16],data6[2][16],dataB[2][16],dataC[2][16]])
# COPS = np.array([data1[2][13],data2[2][13],data3[2][13],data4[2][13],data5[2][13],data6[2][13],dataB[2][13],dataC[2][13]])
# charge = np.array([data1[2][3],data2[2][3],data3[2][3],data4[2][3],data5[2][3],data6[2][3],dataB[2][3],dataC[2][3]])
#to convert string array to integer array
m_dot = m_dot.astype(np.float)
m_dot_inj = m_dot_inj.astype(np.float)
capacity = capacity.astype(np.float)
total_power = total_power.astype(np.float)
compressor_power = compressor_power.astype(np.float)
COPS = COPS.astype(np.float)
charge = charge.astype(np.float)

#plots
#Plot mass flow rate comparison
plt.plot(TestNo,m_dot_exp,'-ob',label='Experimental')
plt.errorbar(TestNo,m_dot_exp, yerr=0.001878)
plt.plot(TestNo,m_dot,'--or',label='Model')
plt.ylim(0.0,0.16)
plt.xlim(0,9)
plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
           [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
plt.xlabel(r'Test condition')
plt.ylabel(r'$\dot m_{suc}$ $[\mathrm{kg/s}]$')
leg = plt.legend(loc='best',fancybox=False,numpoints=1)
frame = leg.get_frame()  
frame.set_linewidth(0.5)
plt.tight_layout()
plt.savefig('results/images/60K_massflow.pdf')
plt.show()
#Plot injection mass flow rate comparison
plt.plot(TestNo,m_dot_inj_exp,'-ob',label='Experimental')
plt.errorbar(TestNo,m_dot_inj_exp, yerr=0.002902)
plt.plot(TestNo,m_dot_inj,'--or',label='Model')
plt.ylim(0.0,0.05)
plt.xlim(0,9)
plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
           [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
plt.xlabel(r'Test condition')
plt.ylabel(r'$\dot m_{inj}$ $[\mathrm{kg/s}]$')
leg = plt.legend(loc='best',fancybox=False,numpoints=1)
frame = leg.get_frame()  
frame.set_linewidth(0.5)
plt.tight_layout()
plt.savefig('results/images/60K_massflow_inj.pdf')
plt.show()
#Plot Capacity comparison
plt.plot(TestNo,capacity_exp,'-ob',label=r'Experimental')
plt.errorbar(TestNo,capacity_exp, yerr=0.1679*capacity_exp)
plt.plot(TestNo,capacity/1000,'--or',label=r'Model')
plt.ylim(0,30)
plt.xlim(0,9)
plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
           [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
plt.xlabel(r'Test condition')
plt.ylabel(r'$\dot Q_{evap}$ $[\mathrm{kW}]$')
leg = plt.legend(loc='best',fancybox=False,numpoints=1)
frame = leg.get_frame()  
frame.set_linewidth(0.5)
plt.tight_layout()
plt.savefig('results/images/60K_capacity.pdf')
plt.show()
#Plot total power comparison
plt.plot(TestNo,total_power_exp,'-ob',label='Experimental')
plt.errorbar(TestNo,total_power_exp, yerr=0.2)
plt.plot(TestNo,total_power/1000,'--or',label='Model')
plt.ylim(0,12)
plt.xlim(0,9)
plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
           [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
plt.xlabel(r'Test condition')
plt.ylabel(r'$\dot E_t$ $[\mathrm{kW}]$')
leg = plt.legend(loc='best',fancybox=False,numpoints=1)
frame = leg.get_frame()  
frame.set_linewidth(0.5)
plt.tight_layout()
plt.savefig('results/images/60K_total_power.pdf')
plt.show()
#Plot compressor power comparison
plt.plot(TestNo,compressor_power_exp,'-ob',label='Experimental')
plt.errorbar(TestNo,compressor_power_exp, yerr=0.1125)
plt.plot(TestNo,compressor_power/1000,'--or',label='Model')
plt.ylim(0,10)
plt.xlim(0,9)
plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
           [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
plt.xlabel(r'Test condition')
plt.ylabel(r'$\dot W_{comp}$ $[\mathrm{kW}]$')
leg = plt.legend(loc='best',fancybox=False,numpoints=1)
frame = leg.get_frame()  
frame.set_linewidth(0.5)
plt.tight_layout()
plt.savefig('results/images/60K_compressor_power.pdf')
plt.show()
#Plot COPS comparison
plt.plot(TestNo,COPS_exp,'-ob',label='Experimental')
plt.errorbar(TestNo,COPS_exp, yerr=0.1704*COPS_exp)
plt.plot(TestNo,COPS,'--or',label='Model')
plt.ylim(0,5)
plt.xlim(0,9)
plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
           [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
plt.xlabel(r'Test condition')
plt.ylabel(r'$\mathrm{COP}_{sys}$')
leg = plt.legend(loc='best',fancybox=False,numpoints=1)
frame = leg.get_frame()  
frame.set_linewidth(0.5)
plt.tight_layout()
plt.savefig(r'results/images/60K_COPS.pdf')
plt.show()
#Plot charge comparison
plt.plot(TestNo,charge_exp,'-ob',label='Experimental')
plt.errorbar(TestNo,charge_exp, yerr=0.0)
plt.plot(TestNo,charge,'--or',label='Model')
plt.ylim(0,8)
plt.xlim(0,9)
plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
           [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
plt.xlabel(r'Test condition')
plt.ylabel(r'$\mathrm{Charge}$ $[\mathrm{kg}]$')
leg = plt.legend(loc='best',fancybox=False,numpoints=1)
frame = leg.get_frame()  
frame.set_linewidth(0.5)
plt.tight_layout()
plt.savefig('results/images/60K_charge.pdf')
plt.show()

#Combine
fig = plt.figure(1, figsize=(15,10), dpi=100)
for i, gtype in enumerate(['Mass', 'Injection_Mass', 'Capacity', 'Power', 'Compressor', 'COPS','Charge']):
    ax = plt.subplot(3, 3, i+1)
    if gtype.startswith('Mass'):
        plt.plot(TestNo,m_dot_exp,'-ob',label='Experimental')
        plt.errorbar(TestNo,m_dot_exp, yerr=0.001878)#0.002*m_dot_exp
        plt.plot(TestNo,m_dot,'--or',label='Model')
        plt.ylim(0.0,0.16)
        plt.xlim(0,9)
        plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                   [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
        plt.xlabel(r'Test condition')
        plt.ylabel(r'$\dot m_{suc}$ $[\mathrm{kg/s}]$')
        leg = plt.legend(loc='best',fancybox=False,numpoints=1)
        frame = leg.get_frame()  
        frame.set_linewidth(0.5)
        #plt.title('Mass flowrate Comparison')
    if gtype.startswith('Injection_Mass'):
        plt.plot(TestNo,m_dot_inj_exp,'-ob',label='Experimental')
        plt.errorbar(TestNo,m_dot_inj_exp, yerr=0.002902)#0.002*m_dot_inj_exp
        plt.plot(TestNo,m_dot_inj,'--or',label='Model')
        plt.ylim(0.0,0.05)
        plt.xlim(0,9)
        plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                   [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
        plt.xlabel(r'Test condition')
        plt.ylabel(r'$\dot m_{inj}$ $[\mathrm{kg/s}]$')
        leg = plt.legend(loc='best',fancybox=False,numpoints=1)
        frame = leg.get_frame()  
        frame.set_linewidth(0.5)
        #plt.title('Mass flowrate Comparison')
    if gtype.startswith('Capacity'):
        plt.plot(TestNo,capacity_exp,'-ob',label='Experimental')
        plt.errorbar(TestNo,capacity_exp, yerr=0.1679*capacity_exp)
        plt.plot(TestNo,capacity/1000,'--or',label='Model')
        plt.ylim(0,30)
        plt.xlim(0,9)
        plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                   [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
        plt.xlabel(r'Test condition')
        plt.ylabel(r'$\dot Q_{evap}$ $[\mathrm{kW}]$')
        leg = plt.legend(loc='best',fancybox=False,numpoints=1)
        frame = leg.get_frame()
        frame.set_linewidth(0.5)
        #plt.title('Capacity Comparison')
    if gtype.startswith('Power'):
        plt.plot(TestNo,total_power_exp,'-ob',label='Experimental')
        plt.errorbar(TestNo,total_power_exp, yerr=0.2)#0.03*total_power_exp
        plt.plot(TestNo,total_power/1000,'--or',label='Model')
        plt.ylim(0,12)
        plt.xlim(0,9)
        plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                   [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
        plt.xlabel(r'Test condition')
        plt.ylabel(r'$\dot E_t$ $[\mathrm{kW}]$')
        leg = plt.legend(loc='best',fancybox=False,numpoints=1)
        frame = leg.get_frame()
        frame.set_linewidth(0.5)
        #plt.title('Total Power Comparison')
    if gtype.startswith('Compressor'):
        plt.plot(TestNo,compressor_power_exp,'-ob',label='Experimental')
        plt.errorbar(TestNo,compressor_power_exp, yerr=0.1125)
        plt.plot(TestNo,compressor_power/1000,'--or',label='Model')
        plt.ylim(0,10)
        plt.xlim(0,9)
        plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                   [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
        plt.xlabel(r'Test condition')
        plt.ylabel(r'$\dot W_{comp}$ $[\mathrm{kW}]$')
        leg = plt.legend(loc='best',fancybox=False,numpoints=1)
        frame = leg.get_frame()
        frame.set_linewidth(0.5)
        #plt.title('Compressor Power Comparison')
    if gtype.startswith('COPS'):
        plt.plot(TestNo,COPS_exp,'-ob',label='Experimental')
        plt.errorbar(TestNo,COPS_exp, yerr=0.1704*COPS_exp)
        plt.plot(TestNo,COPS,'--or',label='Model')
        plt.ylim(0,5)
        plt.xlim(0,9)
        plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                   [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
        plt.xlabel(r'Test condition')
        plt.ylabel(r'$\mathrm{COP}_{sys}$')
        leg = plt.legend(loc='best',fancybox=False,numpoints=1)
        frame = leg.get_frame()
        frame.set_linewidth(0.5)
        #plt.title('System COP Comparison')
    if gtype.startswith('Charge'):
        plt.plot(TestNo,charge_exp,'-ob',label='Experimental')
        plt.errorbar(TestNo,charge_exp, yerr=0.0)
        plt.plot(TestNo,charge,'--or',label='Model')
        plt.ylim(0,8)
        plt.xlim(0,9)
        plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                   [r'', r'1', r'2', r'3',r'4/A', r'5', r'6', r'B', r'C', r''])
        plt.xlabel(r'Test condition')
        plt.ylabel(r'$\mathrm{Charge}$ $[\mathrm{kg}]$')
        leg = plt.legend(loc='best',fancybox=False,numpoints=1)
        frame = leg.get_frame()
        frame.set_linewidth(0.5)
        #plt.title('System charge Comparison')
fig.set_tight_layout(True)
plt.savefig('results/images/60K_comined_new.pdf')
plt.show()