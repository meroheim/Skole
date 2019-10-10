# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 12:43:47 2019

@author: Administrator
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress   # for linear regression
from scipy.optimize import fsolve # to fit function to data
from math import exp

W_al= 5.906       #Vekt til metallbit
T_0 = 23.9 + 273.15        #Romtemperatur
T_f=77              #Temperatur Nitrogen
W_kopper= 3.232 #Vekt isoporkopper
W_mk = 6.049    #Vekt metallbit + tråd
#usikkerhet vekt 0.001 g
#usikkerhet temp 1 K*

Mm_al= 26.98 #molar masse til Al [g/mol]


n_al=(W_al/Mm_al)  #mol
R = 8.314 #(for joule)

L=2.0*(10**5) #Joule/kg




t_før = np.array([0,1,2,3,4,5])
t_før = t_før*60
t_etter = np.array([5.9, 6.9, 7.9, 8.9, 9.9, 10.9, 11.9, 12.9, 13.9, 14.9])          #Skriver inn måledata her direkte
t_etter = t_etter*60
W_tot = [103.8, 100.58, 97.81, 95.07, 92.41, 89.83, 89.13, 86.41, 83.93, 81.52, 79.16, 76.85, 74.57, 72.34, 70.15, 68.00]      #Skriver inn måledata her direkte
W_før = W_tot[0:6]
W_etter = W_tot[6:]

W_etter = np.array(W_etter) - W_al

#FUNKSJON FØR

LinReg_før = linregress(t_før,W_før)

a_før = LinReg_før[0]
b_før = LinReg_før[1]

plt.figure()
plt.title('Mass as a function of time')
t_før_plot = np.linspace(0,t_etter[0])
plt.plot(t_før, W_før, 'xk', Label='Measurement')

t_etter_plot = np.linspace(5*60, 894)
plt.plot(t_etter, W_etter, '+k', label="Measurements")

plt.plot(t_før_plot, a_før*t_før_plot+b_før, '-r', label='Regression before t_1')


#FUNKSJON ETTER

LinReg_etter = linregress(t_etter,W_etter)

a_etter = LinReg_etter[0]
b_etter = LinReg_etter[1]

print('etter', LinReg_etter)
print(LinReg_før)


plt.plot(t_etter_plot, a_etter*t_etter_plot+b_etter, '-b', label='Regression after t_2')

plt.vlines(5*60, 60, 110, colors = "g", linestyle = "--", label = "t_1")
plt.vlines(354, 60, 110, colors = "y", linestyle = "--", label = "t_2")

plt.xlabel("Time, t [s]")
plt.ylabel("Extensive mass, m [g]")

plt.ylim(60, 110)

plt.xticks(ticks = list(range(0,900,120)), fontsize = 11)
plt.yticks(fontsize = 11)


plt.legend()
plt.savefig('maaledata_reg_eksperiment6.png', bbox_inches='tight')
plt.show()

delta_m_1 = abs((a_etter*5*60+b_etter)-(a_før*5*60+b_før))
delta_m_2 = abs((a_etter*354+b_etter)-(a_før*354+b_før))
delta_m =(delta_m_1+delta_m_2)/2

print("Delta_m maks=",delta_m_1)
print("Delta_m min=",delta_m_2)
print("Delta_m =", delta_m)



def f(theta_E):
    delta_Q=(delta_m*0.001)*L
    y_0=theta_E/T_0
    y_f=theta_E/T_f
    epsilon_0 = y_0/(exp(y_0)-1)
    epsilon_f=y_f/(exp(y_f)-1)
    return delta_Q-(3*n_al*R*(((T_0)*epsilon_0)-((T_f)*epsilon_f)))

print('Eksperimentell Theta E=',fsolve(f,290))

m_usikkerhet=(abs(delta_m-delta_m_1)+abs(delta_m-delta_m_2))/2
print('Usikkerhet i Delta_m=', m_usikkerhet)
print('Usikkerhet i T_0 = 1 K')

delta_theta_E = np.sqrt((294.24-287.52)**2 + ((316.19-265.35)*m_usikkerhet)**2)
print('Delta Theta E=', delta_theta_E)
print('Sensitivitet T_0=', 294.24-287.52)
print('Sensitivitet Delta m=', 316.19-265.35)
