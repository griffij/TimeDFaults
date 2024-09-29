"""Make plots of slip rate and paleoearthquake data
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Rectangle

fig = plt.figure(1) 

data = np.genfromtxt('../../OtagoCosmo_publish/data/age_offset_rock_creek_sampled.csv', delimiter=',', skip_header=1) [:,1:]
data_dunstan = np.genfromtxt('../../OtagoCosmo_publish/data/age_offset_neds_creek_sampled.csv', delimiter=',', skip_header=1) [:,1:]
paleo_data = np.genfromtxt('../data/inputs/Hyde_paleodata_simple.csv', delimiter=',', skip_header=1)
paleo_data_dunstan = np.genfromtxt('../data/inputs/Dunstan_paleodata_simple.csv', delimiter=',', skip_header=1) 
#Convert ages to ka
data[:,0] = data[:,0]/1000
data[:,1] = data[:,1]/1000
data_dunstan[:,0] = data_dunstan[:,0]/1000
data_dunstan[:,1] = data_dunstan[:,1]/1000   
# Convert age uncertainties to 2 sigma
data[:,1] = data[:,1]*2.
data_dunstan[:,1] = data_dunstan[:,1]*2.
print(data)
print(paleo_data)
paleo_offsets = [1.8,4,6,8]
paleo_offsets_dunstan = [1.5,3.0,4.5,6.,7.5]


plt.scatter(data[:,0], data[:,2], marker='s', c='mediumblue', label='Hyde displacement')
plt.errorbar(data[:,0], data[:,2],xerr=data[:,1], yerr=data[:,3], fmt='o', c='mediumblue',)
plt.scatter(paleo_data[:,0], paleo_offsets, marker='s', c='cornflowerblue', label='Hyde paleoearthquake')
plt.errorbar(paleo_data[:,0], paleo_offsets, xerr=((paleo_data[:,0]-paleo_data[:,1]), (paleo_data[:,2]-paleo_data[:,0])), yerr=[0.5,0.5,0.5,0.5],fmt='o', c='cornflowerblue')
plt.scatter(data_dunstan[:,0], data_dunstan[:,2], marker='o', c='saddlebrown', label='Dunstan displacement')
plt.errorbar(data_dunstan[:,0], data_dunstan[:,2],xerr=data_dunstan[:,1], yerr=data_dunstan[:,3], fmt='o', c='saddlebrown')
plt.scatter(paleo_data_dunstan[:,0], paleo_offsets_dunstan, marker='o', c='sandybrown', label='Dunstan paleoearthquake')
plt.errorbar(paleo_data_dunstan[:,0], paleo_offsets_dunstan,
             xerr=((paleo_data_dunstan[:,0]-paleo_data_dunstan[:,1]),(paleo_data_dunstan[:,2]-paleo_data_dunstan[:,0])), yerr=[0.5,0.5,0.5,0.5,0.5],
             fmt='o', c='sandybrown')
plt.legend()
ax = plt.gca()
# Add rectangle for inset 
rect = Rectangle((0, 0), 50, 9,                                                                           
                 linewidth=1., edgecolor='0.3', linestyle='dashed',                                                            
                 facecolor='none')
ax.add_patch(rect)
#ax = plt.gca()  
ax.set_xlabel('Age (ka)', fontsize=14)
ax.set_ylabel('Vertical displacement (m)', fontsize=14)
ax.set_xlim(0,350)
ax.set_ylim(0,35)
# Now do inset
fig.add_axes([0.51, 0.17, 0.37, 0.30])
ax = plt.gca()

plt.scatter(data[:,0], data[:,2],marker='s', c='mediumblue')
plt.errorbar(data[:,0], data[:,2],xerr=data[:,1], yerr=data[:,3], fmt='o', c='mediumblue')
plt.scatter(paleo_data[:,0], paleo_offsets, marker='s', c='cornflowerblue')
plt.errorbar(paleo_data[:,0], paleo_offsets, xerr=((paleo_data[:,0]-paleo_data[:,1]), (paleo_data[:,2]-paleo_data[:,0])), yerr=[0.5,0.5,0.5,0.5],fmt='o', c='cornflowerblue')
plt.scatter(data_dunstan[:,0], data_dunstan[:,2], marker='o', c='saddlebrown')
plt.errorbar(data_dunstan[:,0], data_dunstan[:,2],xerr=data_dunstan[:,1], yerr=data_dunstan[:,3], fmt='o', c='saddlebrown')
plt.scatter(paleo_data_dunstan[:,0], paleo_offsets_dunstan, marker='o', c='sandybrown')
plt.errorbar(paleo_data_dunstan[:,0], paleo_offsets_dunstan,
             xerr=((paleo_data_dunstan[:,0]-paleo_data_dunstan[:,1]),(paleo_data_dunstan[:,2]-paleo_data_dunstan[:,0])), yerr=[0.5,0.5,0.5,0.5,0.5],
             fmt='o', c='sandybrown') 
plt.xlim(0,50)
plt.ylim(0,9)


plt.savefig('Hyde_Dunstan_data.png')
