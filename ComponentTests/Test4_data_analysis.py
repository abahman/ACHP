'''
Created on Jan 5, 2016

@author: ammarbahman

Note: you need to have all the test data CSV files executed before run this file.
This will calculate and plot the baseline results COP, cooling capacity .. etc.

'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import csv2rec
from CoolProp.CoolProp import PropsSI
from convert_units import C2K, kPa2Pa

plt.rc('font',family='serif',size =11.0)


#Experimental data input
#Import data from CSV file
data = csv2rec('60K_results/ECU_60K_Test4.csv',delimiter=',')

T1_comp_in = np.array(data[561][2])
T2_comp_out= np.array(data[561][3])
T3_cond_in= np.array(data[561][4])
T4_cond_out= np.array(data[561][5])
T5_exv_in= np.array(data[561][6])
T6_exv_out= np.array(data[561][7])
T7_evap_out= np.array(data[561][8])

T10_ballvalve_in= np.array(data[561][9])
T11_ballvalve_out= np.array(data[561][10])
T12_HGB_out= np.array(data[561][11])
T13_after_bypass= np.array(data[561][12])
T14_after_sightglass= np.array(data[561][13])
T15_ballvalve2_in= np.array(data[561][14])
T16_comp_shell_top= np.array(data[561][15])
T17_comp_shell_bottom= np.array(data[561][16])

T20_evap_in= np.array(data[561][18])
T21_evap_in= np.array(data[561][19])
T22_evap_in= np.array(data[561][20])
T23_evap_in= np.array(data[561][21])
T24_evap_in= np.array(data[561][22])
T25_evap_in= np.array(data[561][23])

T40_evap_air_in= np.array(data[561][25])
T41_evap_air_in= np.array(data[561][26])
T42_evap_air_in= np.array(data[561][27])
T43_evap_air_in= np.array(data[561][28])
T44_evap_air_in= np.array(data[561][29])
T45_evap_air_in= np.array(data[561][30])
T46_evap_air_in= np.array(data[561][31])
T47_evap_air_in= np.array(data[561][33])
T48_evap_air_in= np.array(data[561][34])
T_average_evap_air_in= np.array(data[561][35])
T27_evap_out= np.array(data[561][36])
T28_evap_out= np.array(data[561][37])
T29_evap_out= np.array(data[561][38])
T30_evap_out= np.array(data[561][39])
T31_evap_out= np.array(data[561][40])
T32_evap_out= np.array(data[561][41])

T50_evap_air_out= np.array(data[561][43])
T51_evap_air_out= np.array(data[561][44])
T52_evap_air_out= np.array(data[561][45])
T53_evap_air_out= np.array(data[561][46])
T54_evap_air_out= np.array(data[561][47])
T55_evap_air_out= np.array(data[561][48])
T56_evap_air_out= np.array(data[561][49])
T57_evap_air_out= np.array(data[561][50])
T58_evap_air_out= np.array(data[561][51])
T_average_evap_air_out= np.array(data[561][52])
T60_cond_air_in= np.array(data[561][53])
T61_cond_air_in= np.array(data[561][54])
T62_cond_air_in= np.array(data[561][55])
T63_cond_air_in= np.array(data[561][56])
T64_cond_air_in= np.array(data[561][57])
T65_cond_air_in= np.array(data[561][58])
T66_cond_air_in= np.array(data[561][59])
T67_cond_air_in= np.array(data[561][60])
T68_cond_air_in= np.array(data[561][61])
T_average_cond_air_in= np.array(data[561][62])

T64_cond_air_out= np.array(data[561][70])
T65_cond_air_out= np.array(data[561][71])
T66_cond_air_out= np.array(data[561][72])
T67_cond_air_out= np.array(data[561][73])
T68_cond_air_out= np.array(data[561][74])
T69_cond_air_out= np.array(data[561][75])
T70_cond_air_out= np.array(data[561][76])
T71_cond_air_out= np.array(data[561][77])
T72_cond_air_out= np.array(data[561][78])
T_average_cond_air_out= np.array(data[561][79])

T80_evap_new_grid= np.array(data[561][85])
T81_evap_new_grid= np.array(data[561][86])
T82_evap_new_grid= np.array(data[561][87])
T83_evap_new_grid= np.array(data[561][88])
T84_evap_new_grid= np.array(data[561][89])
T85_evap_new_grid= np.array(data[561][90])
T86_evap_new_grid= np.array(data[561][91])
T87_evap_new_grid= np.array(data[561][92])
T88_evap_new_grid= np.array(data[561][93])
T89_evap_new_grid= np.array(data[561][94])
T90_evap_new_grid= np.array(data[561][95])
T91_evap_new_grid= np.array(data[561][96])
T92_evap_new_grid= np.array(data[561][97])
T93_evap_new_grid= np.array(data[561][98])
T94_evap_new_grid= np.array(data[561][99])
T95_evap_new_grid= np.array(data[561][100])
T96_evap_new_grid= np.array(data[561][101])
T97_evap_new_grid= np.array(data[561][102])
T_average_evap_new_grid= np.array(data[561][103])
Air_pressure_out= np.array(data[561][104])
PT_total_power= np.array(data[561][105])

T_dp_chilled_mirror_out= np.array(data[561][112])

m_R= np.array(data[561][114])
T_dp_chilled_mirror_in= np.array(data[561][115])
PT_fan_evap= np.array(data[561][116])
PT_fan_cond= np.array(data[561][117])

VLN_avg= np.array(data[561][120])
VLL_avg= np.array(data[561][121])
kWH= np.array(data[561][122])
AMP_avg= np.array(data[561][123])
WL2ـN= np.array(data[561][124])
WL3ـN= np.array(data[561][125])
Hzsys= np.array(data[561][126])
PFSYS= np.array(data[561][127])
PT_comp_power= np.array(data[561][128])
AL1= np.array(data[561][129])
AL2= np.array(data[561][130])
AL3= np.array(data[561][131])
WL1_N= np.array(data[561][132])

P1_comp_in= np.array(data[561][138])
P2_comp_out= np.array(data[561][139])
P3_cond_out= np.array(data[561][140])
P4_MM_out= np.array(data[561][141])
P5_exv_in= np.array(data[561][142])
P6_exv_out= np.array(data[561][143])
P7_evap_in= np.array(data[561][144])
P8_evap_out= np.array(data[561][145])
Q_nozzle_apx= np.array(data[561][146])
T_nozzle_in= np.array(data[561][147])
P_atm= np.array(data[561][148])
P_nozzle_in= np.array(data[561][149])
P_nozzle_out= np.array(data[561][150])
dP_nozzle= np.array(data[561][151])
T_room_1155= np.array(data[561][152])
RH_room_1155= np.array(data[561][153])
T_room_1151= np.array(data[561][154])
RH_room_1151= np.array(data[561][155])
            
T_env = np.array([125,115,105,95,85,75])
m_dot_exp = np.array([0.0447,0.0366,0.0349,0.0379,0.0315,0.0335])
capacity_exp = np.array([5578,5014,5150,5814,5133,5649])
total_power_exp = np.array([4296,3857,3530,3329,2852,2654])
compressor_power_exp = np.array([3049,2591,2252,2030,1779,1582])
COPS_exp = np.array([1.298,1.300,1.459,1.746,1.800,2.129])
charge_exp = np.array([1.2,1.2,1.2,1.2,1.2,1.2])


#Basic calculation
ref_fluid = 'R407C'
#Pressure Drop Across Evaporator
DELTAP_evap = P8_evap_out - P7_evap_in

#Pressure Drop Across Condensor
DELTAP_cond = kPa2Pa(P3_cond) - kPa2Pa(P_2)

#Subcooled Temperature
T_sat_cond = PropsSI('T','P',P3_cond_out,'Q',0,ref_fluid)
T_sub = T_sat_cond - T4_cond_out

#Condensation temperature for temperature glide refrigerant(R407C)"
T_sat_cond_glide =  temperature(R$,x=0.5,P=P_cond_out)

"Superheat Temperature"
T_sat_evap = temperature(R$,x=1,P=P_comp_in)
T_super = abs(T_sat_evap - T_evap_out)

W_dot_comp = AvgLookup(Text_Pressure$,'KWsys',n_start,n_end)    "!compressor work (meausred)"
W_fan_cond= AvgLookup(Text_Temperature$,'PT_fan_cond',m_start,m_end)    "!condenser fan power consumption (meausred)"
W_fan_evap = AvgLookup(Text_Temperature$,'PT_fan_evap',m_start,m_end)    "!evaporator fan power consumption (meausred)"
E_dot_t = W_dot_comp+W_fan_cond+W_fan_evap        "total electrical power"

m_dot_R = AvgLookup(Text_Temperature$,'m_massflow',m_start,m_end) *convert(g/s, kg/sec)    "!refrigerant mass flow rate (meausred/micromotion)"

"Performance Calculation"
Q_dot_evap = m_dot_R*(h[10]-h[8])
Q_dot_ton = Q_dot_evap*convert(kW,ton)

COP_c = Q_dot_evap/W_dot_comp        "COP w/compressor"
COP_sys = Q_dot_evap/E_dot_t        "COP system"

eta_c = m_dot_R*(h2s-h[1])/W_dot_comp        "isentropic efficiency of the compressor"
h2s=enthalpy(R$,P=P[2],s=s[1])

EER = Q_dot_evap*convert(kW,BTU/h)/E_dot_t/convert(kW,W)          "EER"

f = 60[Hz]
V_dot = m_dot_R*v[1]
V_dot = V *f
V_dot_EN = V_dot*convert(m^3/s,ft^3/hr)
V_dis = V*convert(m^3,in^3)

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
plt.xlabel('$T_{env}$ $[\degree\mathrm{F}]$')
plt.ylabel('$\dot m_r$ $[\mathrm{kg/s}]$')
#plt.title('Mass flowrate Comparison')
plt.savefig('images/comparison_massflow.pdf')
plt.show()
#Plot Capacity comparison
plt.plot(T_env,capacity_exp,'-ob',label='Experimental')
plt.errorbar(T_env,capacity_exp, yerr=0.00905522*capacity_exp)
plt.plot(T_env,capacity,'--or',label='Model')
plt.ylim(4000,7000)
plt.xlim(70,130)
plt.legend(loc='best',fancybox=False)
plt.xlabel('$T_{env}$ $[\degree\mathrm{F}]$')
plt.ylabel('$\dot Q_{evap}$ $[\mathrm{W}]$')
plt.title('Capacity Comparison')
plt.savefig('images/comparison_capacity.pdf')
plt.show()
#Plot total power comparison
plt.plot(T_env,total_power_exp,'-ob',label='Experimental')
plt.errorbar(T_env,total_power_exp, yerr=0.02618715*total_power_exp)
plt.plot(T_env,total_power,'--or',label='Model')
plt.ylim(2000,5000)
plt.xlim(70,130)
plt.legend(loc='best',fancybox=False)
plt.xlabel('$T_{env}$ $[\degree\mathrm{F}]$')
plt.ylabel('$\dot E_t$ $[\mathrm{W}]$')
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
plt.xlabel('$T_{env}$ $[\degree\mathrm{F}]$')
plt.ylabel('$\dot W_{comp}$ $[\mathrm{W}]$')
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
plt.xlabel('$T_{env}$ $[\degree\mathrm{F}]$')
plt.ylabel('$\mathrm{COP}_{sys}$')
#plt.title('System COP Comparison')
plt.savefig('images/comparison_COPS.pdf')
plt.show()
#Plot charge comparison
plt.plot(T_env,charge_exp,'-ob',label='Experimental')
plt.errorbar(T_env,charge_exp, yerr=0.0)
plt.plot(T_env,charge,'--or',label='Model')
plt.ylim(0,1.6)
plt.xlim(70,130)
plt.legend(loc='best',fancybox=False)
plt.xlabel('$T_{env}$ $[\degree\mathrm{F}]$')
plt.ylabel('$\mathrm{Charge}$ $[\mathrm{kg}]$')
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
        plt.xlabel('$T_{env}$ $[\degree\mathrm{F}]$')
        plt.ylabel('$\dot m_r$ $[\mathrm{kg/s}]$')
        #plt.title('Mass flowrate Comparison')
    if gtype.startswith('Capacity'):
        plt.plot(T_env,capacity_exp,'-ob',label='Experimental')
        plt.errorbar(T_env,capacity_exp, yerr=0.00905522*capacity_exp)
        plt.plot(T_env,capacity,'--or',label='Model')
        plt.ylim(4000,7000)
        plt.xlim(70,130)
        plt.legend(loc='best',fancybox=False)
        plt.xlabel('$T_{env}$ $[\degree\mathrm{F}]$')
        plt.ylabel('$\dot Q_{evap}$ $[\mathrm{W}]$')
        #plt.title('Capacity Comparison')
    if gtype.startswith('Power'):
        plt.plot(T_env,total_power_exp,'-ob',label='Experimental')
        plt.errorbar(T_env,total_power_exp, yerr=0.02618715*total_power_exp)
        plt.plot(T_env,total_power,'--or',label='Model')
        plt.ylim(2000,5000)
        plt.xlim(70,130)
        plt.legend(loc='best',fancybox=False)
        plt.xlabel('$T_{env}$ $[\degree\mathrm{F}]$')
        plt.ylabel('$\dot E_t$ $[\mathrm{W}]$')
        #plt.title('Total Power Comparison')
    if gtype.startswith('Compressor'):
        plt.plot(T_env,compressor_power_exp,'-ob',label='Experimental')
        plt.errorbar(T_env,compressor_power_exp, yerr=112.5)
        plt.plot(T_env,compressor_power,'--or',label='Model')
        plt.ylim(1000,3500)
        plt.xlim(70,130)
        plt.legend(loc='best',fancybox=False)
        plt.xlabel('$T_{env}$ $[\degree\mathrm{F}]$')
        plt.ylabel('$\dot W_{comp}$ $[\mathrm{W}]$')
        #plt.title('Compressor Power Comparison')
    if gtype.startswith('COPS'):
        plt.plot(T_env,COPS_exp,'-ob',label='Experimental')
        plt.errorbar(T_env,COPS_exp, yerr=0.02772727*COPS_exp)
        plt.plot(T_env,COPS,'--or',label='Model')
        plt.ylim(1,2.4)
        plt.xlim(70,130)
        plt.legend(loc='best',fancybox=False)
        plt.xlabel('$T_{env}$ $[\degree\mathrm{F}]$')
        plt.ylabel('$\mathrm{COP}_{sys}$')
        #plt.title('System COP Comparison')
    if gtype.startswith('Charge'):
        plt.plot(T_env,charge_exp,'-ob',label='Experimental')
        plt.errorbar(T_env,charge_exp, yerr=0.0)
        plt.plot(T_env,charge,'--or',label='Model')
        plt.ylim(0,1.6)
        plt.xlim(70,130)
        plt.legend(loc='best',fancybox=False)
        plt.xlabel('$T_{env}$ $[\degree\mathrm{F}]$')
        plt.ylabel('$\mathrm{Charge}$ $[\mathrm{kg}]$')
        #plt.title('System charge Comparison')
fig.set_tight_layout(True)
plt.savefig('images/comined_comparison.pdf')
plt.show()