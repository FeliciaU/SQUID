# -*- coding: utf-8 -*-
"""
Created on Thu May 26 14:20:35 2016

@author: feliciaullstad

Takes SQUID data and sample dimensions to calculate data in bohr magnetons per ion
"""


import matplotlib.pyplot as plt
import pandas as pd

plt.close("all")


Sample_name='F1'
file_location='/home/feliciaullstad/Desktop/Google_Drive/PhD/SmN data/SQUID/F1 SmN 20160223'
file_name='F1 MH 25 K 20160223.rso.dat'

#select='Temperature'    #Select Temperature or Field to do different x-axis plots
select='Field'
Sample_mass=29.5*10**-3        #SQUID sample mass in g
Sample_thickness=28.1    #Sample thickness in nm
Fit_value=0.65

Sample_Silicon_thickness=230     #Silicon thickness in micrometer

trimstart=3
trimend=-17

Scan_program=file_name[trimstart:trimend]   #String containing the kind of measurement was taken
minmagmom=0    #Magnetisation at 0 field or high temperature
Silicon_rho=2.329       #Density of Si i g/cm³

Sample_volume=Sample_mass/Silicon_rho   #Calculates sample volume in cm³
Sample_area=Sample_volume/(Sample_Silicon_thickness*10**-4)  #Calculates sample area in cm²

SmN_volume=Sample_area*Sample_thickness*10**-7 #SmN volume in cm³
#SmN_density=7.353   #SmN density g/cm³
#SmN_mass=SmN_volume*SmN_density
SmN_latticec=5.035  # SmN lattice parameter in Å
SmN_unit_cell_volume=(SmN_latticec*10**-8)**3
SmN_total_unit_cells=SmN_volume/SmN_unit_cell_volume


myfile=file_location+'/'+file_name
data=pd.read_csv(myfile, header=0, sep=',',skiprows=30)


uB_emu=1.0783*10**20    #multiply your emu value with this to get it in bohr magnetons
Sm_ions=SmN_total_unit_cells    #Total number of Sm ions
#############################################################################

data_cut=data[data["Long Reg Fit"] >Fit_value]

fig1=plt.figure()
if select=='Temperature':
    plot=plt.plot(data_cut['Temperature (K)'],data_cut['Long Moment (emu)'],'.')
    plt.xlabel('Temperature (K)')
elif select=='Field':
    plot=plt.plot(data_cut['Field (Oe)'],data_cut['Long Moment (emu)'],'.')
    plt.xlabel('Magnetic field (Oe)')
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.subplots_adjust(bottom=.1, left=.2)
plt.title(Scan_program+' of '+Sample_name)
#plt.ylim(4e5,3e6)
#plt.xlim(0,140)
#plt.close()

uB_per_Sm=data['Long Moment (emu)']*uB_emu/Sm_ions


data['uB per Sm']=pd.Series(uB_per_Sm)

uB_data_cut=data[data["Long Reg Fit"] >Fit_value]
headers=list(data.columns.values)
print uB_data_cut
minimum=min(uB_data_cut['uB per Sm'])
maximum=max(uB_data_cut['uB per Sm'])
print(str(minimum)+' is smallest moment')

data.to_csv(file_location+'-python/'+Sample_name+'_'+Scan_program+'_'+str(Sample_thickness)+'nm_Magnetic_moment.txt', header=headers, index=None, sep=' ', mode='w')


Magnetic_moment=maximum-minimum
Manual_moment=maximum-minimum-minmagmom
print('Maximum uB value in data set is '+str(max(uB_data_cut['uB per Sm']))+' uB/Sm')
print('Minimum uB value in data set is '+str(min(uB_data_cut['uB per Sm']))+' uB/Sm')
print ''
print('Max magnetic moment of sample '+Sample_name+' is '+str(Magnetic_moment)+' u_B/Sm')
print ''
print('Manual magnetic moment of sample '+Sample_name+' is '+str(Manual_moment)+' u_B/Sm')

fig2=plt.figure()
if select=='Temperature':
    plot2=plt.plot(uB_data_cut['Temperature (K)'],(uB_data_cut['uB per Sm']-minimum),'.')
    plt.xlabel('Temperature (K)')
elif select=='Field':
    plot2=plt.plot(uB_data_cut['Field (Oe)'],(uB_data_cut['uB per Sm']),'.')
    plt.xlabel('Magnetic field (Oe)')
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.subplots_adjust(bottom=.1, left=.2)
plt.ylabel('Magnetic moment ($\mu _{B}$/Sm)')
plt.title(Scan_program+' of '+Sample_name)
#ax.text(40, 0.02, 'Max moment: '+str(Magnetic_moment)+' u$_B$/Sm', style='italic',
#        bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
#ax.text(40, 0.015, 'Manual moment: '+str(Manual_moment)+' u$_B$/Sm', style='italic',
#        bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
plot=plt.savefig(file_location+'-python/'+Sample_name+'_'+Scan_program+'_'+str(Sample_thickness)+'nm_Magnetic_moment_plot.pdf', format='pdf', dpi=1200)
plt.show()
