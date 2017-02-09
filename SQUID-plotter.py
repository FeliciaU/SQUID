# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 11:31:29 2016

@author: feliciaullstad

Plots the raw SQUID data
"""


import matplotlib.pyplot as plt
import pandas as pd
import os



plt.close("all")    #Closes plots from previous run

#######################################################
#Settings for the SQUID data
filename='F12 MH 5K 20160427.rso.dat'    #Base name of files that need to be read
file_location='/home/feliciaullstad/Desktop/Google Drive Synced files/PhD/SmN data/SQUID/F12 SmN 20160427'
Sample_name='F12'

Fit_value=0.6
trimstart=4
trimend=-16

Scan_program=filename[trimstart:trimend]


######################################################
if not os.path.exists(file_location+'-python'):     #Checks if the python folder exists
    os.makedirs(file_location+'-python')            #If not, it makes it
#######################################################
#Importing and plotting the SQUID data

"""
Scan through datafolder. Find all .dat files. List them.
"""
datfiles_list=[]
print 'You have the following files:'
for file in os.listdir(file_location):
    if file.endswith(".dat"):
        print(file)
        datfiles_list.append(file)

SQUID_data_raw=pd.read_csv(file_location+'/'+filename, header=0, sep=',',skiprows=30)
SQUID_data=SQUID_data_raw[SQUID_data_raw["Long Reg Fit"] >Fit_value]    # Fitlers away data with fit values under Fit_value


"""
Interesting things to plot:
Temperature (K)
Long Moment (emu)
Long Reg Fit
Field (Oe)
"""


fig1=plt.figure()
plot1 = fig1.add_subplot(111)
#plot1=plt.plot(SQUID_data["Temperature (K)"],SQUID_data["Long Moment (emu)"])
plot1=plt.plot(SQUID_data["Field (Oe)"],SQUID_data["Long Moment (emu)"])
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.subplots_adjust(bottom=.1, left=.2)
#plt.xlabel('Temperature (K)')
plt.xlabel('Field (Oe)')
#plt.xlim(-20000,20000)
#plt.ylim(0.0000008,0.0000016)
plt.ylabel('Magnetic moment (emu)')
plt.title(Scan_program+' '+Sample_name+' Fit value: '+str(Fit_value))
plot1=plt.savefig(file_location+'-python/'+Sample_name+'_plot_'+Scan_program+'.pdf', format='pdf', dpi=1200)


#Plot fit?
fig2=plt.figure()
plot2 = fig2.add_subplot(111)
#plot2=plt.plot(SQUID_data["Temperature (K)"],SQUID_data["Long Reg Fit"])
plot2=plt.plot(SQUID_data["Field (Oe)"],SQUID_data["Long Reg Fit"])
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.subplots_adjust(bottom=.1, left=.2)
#plt.xlabel('Temperature (K)')
plt.xlabel('Field (Oe)')
plt.ylabel('Fit')
#plt.xlim(0,300)
plt.title("Fit values for "+ Scan_program+' '+Sample_name)
plot2=plt.savefig(file_location+'-python/'+Sample_name+'_plot_'+Scan_program+'_fit.pdf', format='pdf', dpi=1200)
