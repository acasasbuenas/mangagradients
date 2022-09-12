"""
Created on Wed 31 Aug 2022
@author: Alba Casasbuenas
INFO SIG -> AGN O NO 
"""

from marvin.tools import Maps
from marvin.tools import Image
from marvin.utils.general.general import get_drpall_table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits, ascii
from astropy.table import Table
import time

t=Table.read('tabla_mangaagn_updated2.ecsv',format='ascii.ecsv')

for i in range(len(t)):
	try:
		t['NNEIGH'][i]=t['NNEIGH'][i][0]
		#print(i)
	except:
		pass

print('total de galaxias', len(t['PLATEIFU']))

#LISTA DE GALAXIAS EDGE-ON
print('Galaxias edge-on:', np.sum(t['P_EDGEON']>0.75))
plateifu_edgeon=t['PLATEIFU'][t['P_EDGEON']>0.75]

#LISTA DE GALAXIAS CON EL NÚCLEO ÓPTICO LEJOS DEL PÍXEL CENTRAL
print('Galaxias a revisar:', np.sum(t['FLAG']=='REVISAR'))
plateifu_revisar=t['PLATEIFU'][t['FLAG']=='REVISAR']

#LISTA DE GALAXIAS PEQUEÑAS
mask_px=((t['PX_V']>160) & (t['PX_VHA']>160) & (t['PX_VOIII']>160))
print('Galaxias grandes:', np.sum(mask_px))

#LISTA DE GALAXIAS SIN FLAG:
print('Galaxias sin flag: ', np.sum((t['FLAG']=='OK') ))
mask_noflag=((t['FLAG']=='OK') & mask_px)
plateifu_noflag=t['PLATEIFU'][mask_noflag]

#LISTA DE GALAXIAS CON BARRA:
mask_barA=(t['P_BAR']>2/3)
mask_barAB=((t['P_BAR']>1/3) & (t['P_BAR']<2/3))
mask_barB=(t['P_BAR']<1/3)

#LISTA DE GALAXIAS AGN:
mask_agn=(t['AGN']>0)
mask_noagn=~mask_agn
print('AGN: ', np.sum((t['AGN'])>0 ))

#LISTA DE GALAXIAS AISLADAS/NO AISLADAS:
mask_nosig=(t['NNEIGH']>0)
mask_sig=~mask_nosig
print('ISOLATED GALAXIES: ', np.sum(mask_sig))

#MÁSCARAS DE TIPO MORFOLÓGICO

mask_s2=(t['MORPH']=='S2')
mask_s1=(t['MORPH']=='S1')
mask_s0=(t['MORPH']=='S0')
mask_e=(t['MORPH']=='E')
mask_irr=(t['MORPH']=='-')

#MÁSCARAS DE MASA

mask_mass1=(t['LOGMASS']<=10)
mask_mass2=((t['LOGMASS']>10) & (t['LOGMASS']<=10.5))
mask_mass3=((t['LOGMASS']>10.5) & (t['LOGMASS']<=11))
mask_mass4=((t['LOGMASS']>11) & (t['LOGMASS']<=11.5))
mask_mass5=(t['LOGMASS']>11.5)

#MÁSCARAS DE REDSHIFT

mask_z1=(t['NSA_Z']<=0.03)
mask_z2=((t['NSA_Z']>0.03) & (t['NSA_Z']<=0.05))
mask_z3=(t['NSA_Z']>0.05)

print('\n Total de galaxias', np.sum(mask_sig))
print('Galaxias con flag: ', np.sum(mask_sig & (t['FLAG']=='REVISAR')))
print('Galaxias peques: ', np.sum(mask_sig & ~mask_px))
print('Galaxias utilizadas: ', np.sum(mask_sig & mask_noflag))
print('Galaxias utilizadas con AGN: ', np.sum(mask_sig & mask_noflag & mask_agn))

#CÍRCULO CON EL TAMAÑO DE LA FIBRA DE MANGA
fiber_size=2 #arcsec
psf_size=2.5 #arcsec

b=[1, 2, 3]

#WE'LL ONLY STUDY SIG GALAXIES -> ISOLATED ONES

####################################################################################
###########################  ON-KC ESTELAR  ########################################
####################################################################################

circle1_st=plt.Circle((0,0),fiber_size,fill=False,color='darkviolet', linewidth=2)
circle2_st=plt.Circle((0,0),fiber_size*2,fill=False,color='darkviolet', linestyle='dashed', linewidth=2)

d_stellar=np.sqrt((t['KCRA_ST'])**2+(t['KCDEC_ST'])**2)

#BAR A + NOT AGN
barA_noagn_b1=np.sum(d_stellar[mask_noflag & mask_sig & mask_barA & mask_noagn]<fiber_size)
barA_noagn_b2=np.sum((d_stellar[mask_noflag & mask_sig & mask_barA & mask_noagn]>fiber_size) & (d_stellar[mask_noflag & mask_sig & mask_barA & mask_noagn]<2*fiber_size)) 
barA_noagn_b3=np.sum(d_stellar[mask_noflag & mask_sig & mask_barA & mask_noagn]>2*fiber_size)
#BAR A + ISOLATED
barA_agn_b1=np.sum(d_stellar[mask_noflag & mask_sig & mask_barA & mask_agn]<fiber_size)
barA_agn_b2=np.sum((d_stellar[mask_noflag & mask_sig & mask_barA & mask_agn]>fiber_size) & (d_stellar[mask_noflag & mask_sig & mask_barA & mask_agn]<2*fiber_size)) 
barA_agn_b3=np.sum(d_stellar[mask_noflag & mask_sig & mask_barA & mask_agn]>2*fiber_size)
#BAR AB + INTERACTING
barAB_noagn_b1=np.sum(d_stellar[mask_noflag & mask_sig & mask_barAB & mask_noagn]<fiber_size)
barAB_noagn_b2=np.sum((d_stellar[mask_noflag & mask_sig & mask_barAB & mask_noagn]>fiber_size) & (d_stellar[mask_noflag & mask_sig & mask_barAB & mask_noagn]<2*fiber_size)) 
barAB_noagn_b3=np.sum(d_stellar[mask_noflag & mask_sig & mask_barAB & mask_noagn]>2*fiber_size)
#BAR AB + ISOLATED
barAB_agn_b1=np.sum(d_stellar[mask_noflag & mask_sig & mask_barAB & mask_agn]<fiber_size)
barAB_agn_b2=np.sum((d_stellar[mask_noflag & mask_sig & mask_barAB & mask_agn]>fiber_size) & (d_stellar[mask_noflag & mask_sig & mask_barAB & mask_agn]<2*fiber_size)) 
barAB_agn_b3=np.sum(d_stellar[mask_noflag & mask_sig & mask_barAB & mask_agn]>2*fiber_size)
#BAR B + INTERACTING
barB_noagn_b1=np.sum(d_stellar[mask_noflag & mask_sig & mask_barB & mask_noagn]<fiber_size)
barB_noagn_b2=np.sum((d_stellar[mask_noflag & mask_sig & mask_barB & mask_noagn]>fiber_size) & (d_stellar[mask_noflag & mask_sig & mask_barB & mask_noagn]<2*fiber_size)) 
barB_noagn_b3=np.sum(d_stellar[mask_noflag & mask_sig & mask_barB & mask_noagn]>2*fiber_size)
#BAR B + ISOLATED
barB_agn_b1=np.sum(d_stellar[mask_noflag & mask_sig & mask_barB & mask_agn]<fiber_size)
barB_agn_b2=np.sum((d_stellar[mask_noflag & mask_sig & mask_barB & mask_agn]>fiber_size) & (d_stellar[mask_noflag & mask_sig & mask_barB & mask_agn]<2*fiber_size)) 
barB_agn_b3=np.sum(d_stellar[mask_noflag & mask_sig & mask_barB & mask_agn]>2*fiber_size)

barA_noagn_ST=[barA_noagn_b1,barA_noagn_b2,barA_noagn_b3]
barA_agn_ST=[barA_agn_b1,barA_agn_b2,barA_agn_b3]
barAB_noagn_ST=[barAB_noagn_b1,barAB_noagn_b2,barAB_noagn_b3]
barAB_agn_ST=[barAB_agn_b1,barAB_agn_b2,barAB_agn_b3]
barB_noagn_ST=[barB_noagn_b1,barB_noagn_b2,barB_noagn_b3]
barB_agn_ST=[barB_agn_b1,barB_agn_b2,barB_agn_b3]

bottom1_ST=np.zeros(3)
bottom2_ST=bottom1_ST+barB_agn_ST
bottom3_ST=bottom2_ST+barB_noagn_ST
bottom4_ST=bottom3_ST+barAB_agn_ST
bottom5_ST=bottom4_ST+barAB_noagn_ST
bottom6_ST=bottom5_ST+barA_agn_ST

top_ST=bottom6_ST+barA_noagn_ST
agns_ST=np.array(barA_agn_ST)+np.array(barAB_agn_ST)+np.array(barB_agn_ST)

porcentajes_ST=np.round(np.nan_to_num((top_ST/np.sum(top_ST))*100, nan=0))
porcentajes_agn_ST=np.round(np.nan_to_num((agns_ST/top_ST)*100,nan=0))

barlabels_ST=[]
for i in range(len(b)):
	a_ST=str(int(porcentajes_ST[i]))+'% (AGN:' + str(int(porcentajes_agn_ST[i]))+'% )'
	barlabels_ST.append(a_ST)

####################################################################################
###########################  ON-KC HALPHA  #########################################
####################################################################################

circle1_ha=plt.Circle((0,0),fiber_size,fill=False,color='darkviolet', linewidth=2)
circle2_ha=plt.Circle((0,0),fiber_size*2,fill=False,color='darkviolet', linewidth=2,linestyle='dashed')

d_ha=np.sqrt((t['KCRA_HA'])**2+(t['KCDEC_HA'])**2)

#BAR A + INTERACTING
barA_noagn_b1=np.sum(d_ha[mask_noflag & mask_sig & mask_barA & mask_noagn]<1.3)
barA_noagn_b2=np.sum((d_ha[mask_noflag & mask_sig & mask_barA & mask_noagn]>1.3) & (d_ha[mask_noflag & mask_sig & mask_barA & mask_noagn]<2.7)) 
barA_noagn_b3=np.sum(d_ha[mask_noflag & mask_sig & mask_barA & mask_noagn]>2.7)
#BAR A + ISOLATED
barA_agn_b1=np.sum(d_ha[mask_noflag & mask_sig & mask_barA & mask_agn]<1.3)
barA_agn_b2=np.sum((d_ha[mask_noflag & mask_sig & mask_barA & mask_agn]>1.3) & (d_ha[mask_noflag & mask_sig & mask_barA & mask_agn]<2.7)) 
barA_agn_b3=np.sum(d_ha[mask_noflag & mask_sig & mask_barA & mask_agn]>2.7)
#BAR AB + INTERACTING
barAB_noagn_b1=np.sum(d_ha[mask_noflag & mask_sig & mask_barAB & mask_noagn]<1.3)
barAB_noagn_b2=np.sum((d_ha[mask_noflag & mask_sig & mask_barAB & mask_noagn]>1.3) & (d_ha[mask_noflag & mask_sig & mask_barAB & mask_noagn]<2.7)) 
barAB_noagn_b3=np.sum(d_ha[mask_noflag & mask_sig & mask_barAB & mask_noagn]>2.7)
#BAR AB + ISOLATED
barAB_agn_b1=np.sum(d_ha[mask_noflag & mask_sig & mask_barAB & mask_agn]<1.3)
barAB_agn_b2=np.sum((d_ha[mask_noflag & mask_sig & mask_barAB & mask_agn]>1.3) & (d_ha[mask_noflag & mask_sig & mask_barAB & mask_agn]<2.7)) 
barAB_agn_b3=np.sum(d_ha[mask_noflag & mask_sig & mask_barAB & mask_agn]>2.7)
#BAR B + INTERACTING
barB_noagn_b1=np.sum(d_ha[mask_noflag & mask_sig & mask_barB & mask_noagn]<1.3)
barB_noagn_b2=np.sum((d_ha[mask_noflag & mask_sig & mask_barB & mask_noagn]>1.3) & (d_ha[mask_noflag & mask_sig & mask_barB & mask_noagn]<2.7)) 
barB_noagn_b3=np.sum(d_ha[mask_noflag & mask_sig & mask_barB & mask_noagn]>2.7)
#BAR B + ISOLATED
barB_agn_b1=np.sum(d_ha[mask_noflag & mask_sig & mask_barB & mask_agn]<1.3)
barB_agn_b2=np.sum((d_ha[mask_noflag & mask_sig & mask_barB & mask_agn]>1.3) & (d_ha[mask_noflag & mask_sig & mask_barB & mask_agn]<2.7)) 
barB_agn_b3=np.sum(d_ha[mask_noflag & mask_sig & mask_barB & mask_agn]>2.7)

barA_noagn_HA=[barA_noagn_b1,barA_noagn_b2,barA_noagn_b3]
barA_agn_HA=[barA_agn_b1,barA_agn_b2,barA_agn_b3]
barAB_noagn_HA=[barAB_noagn_b1,barAB_noagn_b2,barAB_noagn_b3]
barAB_agn_HA=[barAB_agn_b1,barAB_agn_b2,barAB_agn_b3]
barB_noagn_HA=[barB_noagn_b1,barB_noagn_b2,barB_noagn_b3]
barB_agn_HA=[barB_agn_b1,barB_agn_b2,barB_agn_b3]

bottom1_HA=np.zeros(3)
bottom2_HA=bottom1_HA+barB_agn_HA
bottom3_HA=bottom2_HA+barB_noagn_HA
bottom4_HA=bottom3_HA+barAB_agn_HA
bottom5_HA=bottom4_HA+barAB_noagn_HA
bottom6_HA=bottom5_HA+barA_agn_HA

top_HA=bottom4_HA+barA_noagn_HA
agns_HA=np.array(barA_agn_HA)+np.array(barAB_agn_HA)+np.array(barB_agn_HA)

porcentajes_HA=np.round(np.nan_to_num((top_HA/np.sum(top_HA))*100, nan=0))
porcentajes_agn_HA=np.round( np.nan_to_num((agns_HA/top_HA)*100,nan=0))

barlabels_HA=[]
for i in range(len(b)):
	a_HA=str(int(porcentajes_HA[i]))+'% (AGN:' + str(int(porcentajes_agn_HA[i]))+'% )'
	barlabels_HA.append(a_HA)

####################################################################################
###########################  ON-KC OIII  ###########################################
####################################################################################

circle1_oiii=plt.Circle((0,0),fiber_size,fill=False,color='darkviolet', linewidth=2)
circle2_oiii=plt.Circle((0,0),fiber_size*2,fill=False,color='darkviolet', linewidth=2,linestyle='dashed')

d_oiii=np.sqrt((t['KCRA_OIII'])**2+(t['KCDEC_OIII'])**2)

#BAR A + INTERACTING
barA_noagn_b1=np.sum(d_oiii[mask_noflag & mask_sig & mask_barA & mask_noagn]<1.3)
barA_noagn_b2=np.sum((d_oiii[mask_noflag & mask_sig & mask_barA & mask_noagn]>1.3) & (d_oiii[mask_noflag & mask_sig & mask_barA & mask_noagn]<2.7)) 
barA_noagn_b3=np.sum(d_oiii[mask_noflag & mask_sig & mask_barA & mask_noagn]>2.7)
#BAR A + ISOLATED
barA_agn_b1=np.sum(d_oiii[mask_noflag & mask_sig & mask_barA & mask_agn]<1.3)
barA_agn_b2=np.sum((d_oiii[mask_noflag & mask_sig & mask_barA & mask_agn]>1.3) & (d_oiii[mask_noflag & mask_sig & mask_barA & mask_agn]<2.7)) 
barA_agn_b3=np.sum(d_oiii[mask_noflag & mask_sig & mask_barA & mask_agn]>2.7)
#BAR AB + INTERACTING
barAB_noagn_b1=np.sum(d_oiii[mask_noflag & mask_sig & mask_barAB & mask_noagn]<1.3)
barAB_noagn_b2=np.sum((d_oiii[mask_noflag & mask_sig & mask_barAB & mask_noagn]>1.3) & (d_oiii[mask_noflag & mask_sig & mask_barAB & mask_noagn]<2.7)) 
barAB_noagn_b3=np.sum(d_oiii[mask_noflag & mask_sig & mask_barAB & mask_noagn]>2.7)
#BAR AB + ISOLATED
barAB_agn_b1=np.sum(d_oiii[mask_noflag & mask_sig & mask_barAB & mask_agn]<1.3)
barAB_agn_b2=np.sum((d_oiii[mask_noflag & mask_sig & mask_barAB & mask_agn]>1.3) & (d_oiii[mask_noflag & mask_sig & mask_barAB & mask_agn]<2.7)) 
barAB_agn_b3=np.sum(d_oiii[mask_noflag & mask_sig & mask_barAB & mask_agn]>2.7)
#BAR B + INTERACTING
barB_noagn_b1=np.sum(d_oiii[mask_noflag & mask_sig & mask_barB & mask_noagn]<1.3)
barB_noagn_b2=np.sum((d_oiii[mask_noflag & mask_sig & mask_barB & mask_noagn]>1.3) & (d_oiii[mask_noflag & mask_sig & mask_barB & mask_noagn]<2.7)) 
barB_noagn_b3=np.sum(d_oiii[mask_noflag & mask_sig & mask_barB & mask_noagn]>2.7)
#BAR B + ISOLATED
barB_agn_b1=np.sum(d_oiii[mask_noflag & mask_sig & mask_barB & mask_agn]<1.3)
barB_agn_b2=np.sum((d_oiii[mask_noflag & mask_sig & mask_barB & mask_agn]>1.3) & (d_oiii[mask_noflag & mask_sig & mask_barB & mask_agn]<2.7)) 
barB_agn_b3=np.sum(d_oiii[mask_noflag & mask_sig & mask_barB & mask_agn]>2.7)

barA_noagn_oiii=[barA_noagn_b1,barA_noagn_b2,barA_noagn_b3]
barA_agn_oiii=[barA_agn_b1,barA_agn_b2,barA_agn_b3]
barAB_noagn_oiii=[barAB_noagn_b1,barAB_noagn_b2,barAB_noagn_b3]
barAB_agn_oiii=[barAB_agn_b1,barAB_agn_b2,barAB_agn_b3]
barB_noagn_oiii=[barB_noagn_b1,barB_noagn_b2,barB_noagn_b3]
barB_agn_oiii=[barB_agn_b1,barB_agn_b2,barB_agn_b3]

bottom1_oiii=np.zeros(3)
bottom2_oiii=bottom1_oiii+barB_agn_oiii
bottom3_oiii=bottom2_oiii+barB_noagn_oiii
bottom4_oiii=bottom3_oiii+barAB_agn_oiii
bottom5_oiii=bottom4_oiii+barAB_noagn_oiii
bottom6_oiii=bottom5_oiii+barA_agn_oiii
top_oiii=bottom6_oiii+barA_noagn_oiii
agns_oiii=np.array(barA_agn_oiii)+np.array(barAB_agn_oiii)+np.array(barB_agn_oiii)


porcentajes_oiii=np.round(np.nan_to_num((top_oiii/np.sum(top_oiii))*100, nan=0))
porcentajes_agn_oiii=np.round( np.nan_to_num((agns_oiii/top_oiii)*100,nan=0))

barlabels_oiii=[]
for i in range(len(b)):
	a_oiii=str(int(porcentajes_oiii[i]))+'% (AGN:' + str(int(porcentajes_agn_oiii[i]))+'%)'
	barlabels_oiii.append(a_oiii)

####################################################################################
###########################  Hacemos los gráficos ##################################
####################################################################################

fig,axs = plt.subplots(2,3,figsize=(16,10))

axs[0,0].set_title('Stellar Kinematics')
axs[0,0].set_xlim(-10.5, 10.5)
axs[0,0].set_ylim(-10.5, 10.5)
axs[0,0].set_xlabel('x-offset [ON-KC]')
axs[0,0].set_ylabel('y-offset [ON-KC]')
axs[0,0].scatter(0,0,marker='+',color='darkviolet',zorder=10)
axs[0,0].add_artist(circle1_st)
axs[0,0].add_artist(circle2_st)
axs[0,0].errorbar(t['KCRA_ST'][mask_noagn],  t['KCDEC_ST'][mask_noagn], yerr=t['KCDEC_ST_ERR'][mask_noagn], xerr=t['KCRA_ST_ERR'][mask_noagn], fmt='ks',fillstyle='none',zorder=-1,label='Not Confirmed AGN')
axs[0,0].errorbar(t['KCRA_ST'][mask_agn],  t['KCDEC_ST'][mask_agn], yerr=t['KCDEC_ST_ERR'][mask_agn], xerr=t['KCRA_ST_ERR'][mask_agn], fmt='bs',fillstyle='none',zorder=0,label='Confirmed AGN')
#axs[0,0].errorbar(t['KCRA_ST'][mask_agn],  t['KCDEC_ST'][mask_agn], yerr=t['KCDEC_ST_ERR'][mask_agn], xerr=t['KCRA_ST_ERR'][mask_agn], fmt='gs',fillstyle='none',zorder=0)
axs[0,0].legend(loc='upper right')

p6_ST=axs[1,0].bar(b,barA_noagn_ST, bottom=bottom6_ST,color='steelblue',label='Not AGN + Bar A')
p5_ST=axs[1,0].bar(b,barA_agn_ST, bottom=bottom5_ST,color='rebeccapurple',label='AGN + Bar A')
p4_ST=axs[1,0].bar(b,barAB_noagn_ST, bottom=bottom4_ST,color='lightskyblue',label='Not AGN + Bar AB')
p3_ST=axs[1,0].bar(b,barAB_agn_ST, bottom=bottom3_ST,color='blueviolet',label='AGN + Bar AB')
p2_ST=axs[1,0].bar(b,barB_noagn_ST, bottom=bottom2_ST,color='paleturquoise',label='Not AGN + Bar B')
p1_ST=axs[1,0].bar(b,barB_agn_ST,bottom=bottom1_ST,color='mediumpurple',label='AGN + Bar B')

axs[1,0].bar_label(p6_ST,labels=barlabels_ST)
axs[1,0].set_xticks(b, labels=['d<2 \" ','2 \" d < 4 \"','d > 4 \"'])
axs[1,0].legend(loc='upper right')
axs[1,0].set_xlabel('Number of galaxies')
axs[1,0].set_ylabel(r'ON-KC distance')

axs[0,1].set_title(r'H$\alpha$ Kinematics')
axs[0,1].set_xlim(-10.5, 10.5)
axs[0,1].set_ylim(-10.5, 10.5)
axs[0,1].set_xlabel('x-offset [ON-KC]')
axs[0,1].set_ylabel('y-offset [ON-KC]')
axs[0,1].scatter(0,0,marker='+',color='darkviolet',zorder=10)
axs[0,1].add_artist(circle1_ha)
axs[0,1].add_artist(circle2_ha)
axs[0,1].errorbar(t['KCRA_HA'][mask_noagn],  t['KCDEC_HA'][mask_noagn], yerr=t['KCDEC_HA_ERR'][mask_noagn], xerr=t['KCRA_HA_ERR'][mask_noagn], fmt='ks',fillstyle='none',zorder=-1,label='Not Confirmed AGN')
axs[0,1].errorbar(t['KCRA_HA'][mask_agn],  t['KCDEC_HA'][mask_agn], yerr=t['KCDEC_HA_ERR'][mask_agn], xerr=t['KCRA_HA_ERR'][mask_agn], fmt='bs',fillstyle='none',zorder=0,label='Confirmed AGN')
#axs[0,1].errorbar(t['KCRA_HA'][mask_agn],  t['KCDEC_HA'][mask_agn], yerr=t['KCDEC_HA_ERR'][mask_agn], xerr=t['KCRA_HA_ERR'][mask_agn], fmt='gs',fillstyle='none',zorder=0)
axs[0,1].legend(loc='upper right')

p6_HA=axs[1,1].bar(b,barA_noagn_HA, bottom=bottom6_HA,color='steelblue',label='Not AGN + Bar A')
p5_HA=axs[1,1].bar(b,barA_agn_HA, bottom=bottom5_HA,color='rebeccapurple',label='AGN + Bar A')
p4_HA=axs[1,1].bar(b,barAB_noagn_HA,  bottom=bottom4_HA, color='lightskyblue', label='Not AGN + Bar AB')
p3_HA=axs[1,1].bar(b,barAB_agn_HA,bottom=bottom3_HA,color='blueviolet',label='AGN + Bar AB')
p2_HA=axs[1,1].bar(b,barB_noagn_HA,  bottom=bottom2_HA, color='paleturquoise', label='Not AGN + Bar B')
p1_HA=axs[1,1].bar(b,barB_agn_HA,bottom=bottom1_HA,color='mediumpurple',label='AGN + Bar B')

axs[1,1].bar_label(p6_HA,labels=barlabels_HA)
axs[1,1].set_xticks(b, labels=['d<2 \" ','2 \" d < 4 \"','d > 4 \"'])
axs[1,1].legend(loc='upper right')
axs[1,1].set_xlabel('Number of galaxies')
axs[1,1].set_ylabel('ON-KC distance')

axs[0,2].set_title('[OIII] Kinematics')
axs[0,2].set_xlim(-10.5, 10.5)
axs[0,2].set_ylim(-10.5, 10.5)
axs[0,2].set_xlabel('x-offset [ON-KC]')
axs[0,2].set_ylabel('y-offset [ON-KC]')
axs[0,2].scatter(0,0,marker='+',color='darkviolet',zorder=10)
axs[0,2].add_artist(circle1_oiii)
axs[0,2].add_artist(circle2_oiii)
axs[0,2].errorbar(t['KCRA_OIII'][mask_noagn],  t['KCDEC_OIII'][mask_noagn], yerr=t['KCDEC_OIII_ERR'][mask_noagn], xerr=t['KCRA_OIII_ERR'][mask_noagn], fmt='ks',fillstyle='none',zorder=-1,label='Not Confirmed AGN')
axs[0,2].errorbar(t['KCRA_OIII'][mask_agn],  t['KCDEC_OIII'][mask_agn], yerr=t['KCDEC_OIII_ERR'][mask_agn], xerr=t['KCRA_OIII_ERR'][mask_agn], fmt='bs',fillstyle='none',zorder=0,label='Confirmed AGN')
#axs[0,2].errorbar(t['KCRA_OIII'][mask_agn],  t['KCDEC_OIII'][mask_agn], yerr=t['KCDEC_OIII_ERR'][mask_agn], xerr=t['KCRA_OIII_ERR'][mask_agn], fmt='gs',fillstyle='none',zorder=0)
axs[0,2].legend(loc='upper right')

p6_oiii=axs[1,2].bar(b,barA_noagn_oiii,  bottom=bottom6_oiii, color='steelblue', label='Not AGN + Bar A')
p5_oiii=axs[1,2].bar(b,barA_agn_oiii, bottom=bottom5_oiii,color='rebeccapurple',label='AGN + Bar A')
p4_oiii=axs[1,2].bar(b,barAB_noagn_oiii,  bottom=bottom4_oiii, color='lightskyblue', label='Not AGN + Bar AB')
p3_oiii=axs[1,2].bar(b,barAB_agn_oiii,bottom=bottom3_oiii,color='blueviolet',label='AGN + Bar AB')
p2_oiii=axs[1,2].bar(b,barB_noagn_oiii,  bottom=bottom2_oiii, color='paleturquoise', label='Not AGN + Bar B')
p1_oiii=axs[1,2].bar(b,barB_agn_oiii,bottom=bottom1_oiii,color='mediumpurple',label='AGN + Bar B')

axs[1,2].bar_label(p6_oiii,labels=barlabels_oiii)
axs[1,2].set_xticks(b, labels=['d<2 \" ','2 \" d < 4 \"','d > 4 \"'])
axs[1,2].legend(loc='upper right')
axs[1,2].set_xlabel('Number of galaxies')
axs[1,2].set_ylabel('ON-KC distance')
#fig.suptitle('SDSS-BASED ISOLATED GALAXIES')
fig.savefig('sig_offset_01.png')
plt.show()

####################################################################################
################## LO REPETIMOS PARA LOS KC ########################################
####################################################################################

####################################################################################
###########################  KC ESTELAR - KC HALPHA ################################
####################################################################################

circle1_st_ha=plt.Circle((0,0),fiber_size,fill=False,color='darkviolet', linewidth=2)
circle2_st_ha=plt.Circle((0,0),fiber_size*2,fill=False,color='darkviolet', linewidth=2,linestyle='dashed')

d_st_ha=np.sqrt((t['KCRA_ST'] - t['KCRA_HA'])**2+(t['KCDEC_ST']- t['KCDEC_HA'])**2)
#BAR A + INTERACTING
barA_noagn_b1=np.sum(d_st_ha[mask_noflag & mask_sig & mask_barA & mask_noagn]<fiber_size)
barA_noagn_b2=np.sum((d_st_ha[mask_noflag & mask_sig & mask_barA & mask_noagn]>fiber_size) & (d_st_ha[mask_noflag & mask_sig & mask_barA & mask_noagn]<2*fiber_size)) 
barA_noagn_b3=np.sum(d_st_ha[mask_noflag & mask_sig & mask_barA & mask_noagn]>2*fiber_size)
#BAR A + ISOLATED
barA_agn_b1=np.sum(d_st_ha[mask_noflag & mask_sig & mask_barA & mask_agn]<fiber_size)
barA_agn_b2=np.sum((d_st_ha[mask_noflag & mask_sig & mask_barA & mask_agn]>fiber_size) & (d_st_ha[mask_noflag & mask_sig & mask_barA & mask_agn]<2*fiber_size)) 
barA_agn_b3=np.sum(d_st_ha[mask_noflag & mask_sig & mask_barA & mask_agn]>2*fiber_size)
#BAR AB + INTERACTING
barAB_noagn_b1=np.sum(d_st_ha[mask_noflag & mask_sig & mask_barAB & mask_noagn]<fiber_size)
barAB_noagn_b2=np.sum((d_st_ha[mask_noflag & mask_sig & mask_barAB & mask_noagn]>fiber_size) & (d_st_ha[mask_noflag & mask_sig & mask_barAB & mask_noagn]<2*fiber_size)) 
barAB_noagn_b3=np.sum(d_st_ha[mask_noflag & mask_sig & mask_barAB & mask_noagn]>2*fiber_size)
#BAR AB + ISOLATED
barAB_agn_b1=np.sum(d_st_ha[mask_noflag & mask_sig & mask_barAB & mask_agn]<fiber_size)
barAB_agn_b2=np.sum((d_st_ha[mask_noflag & mask_sig & mask_barAB & mask_agn]>fiber_size) & (d_st_ha[mask_noflag & mask_sig & mask_barAB & mask_agn]<2*fiber_size)) 
barAB_agn_b3=np.sum(d_st_ha[mask_noflag & mask_sig & mask_barAB & mask_agn]>2*fiber_size)
#BAR B + INTERACTING
barB_noagn_b1=np.sum(d_st_ha[mask_noflag & mask_sig & mask_barB & mask_noagn]<fiber_size)
barB_noagn_b2=np.sum((d_st_ha[mask_noflag & mask_sig & mask_barB & mask_noagn]>fiber_size) & (d_st_ha[mask_noflag & mask_sig & mask_barB & mask_noagn]<2*fiber_size)) 
barB_noagn_b3=np.sum(d_st_ha[mask_noflag & mask_sig & mask_barB & mask_noagn]>2*fiber_size)
#BAR B + ISOLATED
barB_agn_b1=np.sum(d_st_ha[mask_noflag & mask_sig & mask_barB & mask_agn]<fiber_size)
barB_agn_b2=np.sum((d_st_ha[mask_noflag & mask_sig & mask_barB & mask_agn]>fiber_size) & (d_st_ha[mask_noflag & mask_sig & mask_barB & mask_agn]<2*fiber_size)) 
barB_agn_b3=np.sum(d_st_ha[mask_noflag & mask_sig & mask_barB & mask_agn]>2*fiber_size)

barA_noagn_st_ha=[barA_noagn_b1,barA_noagn_b2,barA_noagn_b3]
barA_agn_st_ha=[barA_agn_b1,barA_agn_b2,barA_agn_b3]
barAB_noagn_st_ha=[barAB_noagn_b1,barAB_noagn_b2,barAB_noagn_b3]
barAB_agn_st_ha=[barAB_agn_b1,barAB_agn_b2,barAB_agn_b3]
barB_noagn_st_ha=[barB_noagn_b1,barB_noagn_b2,barB_noagn_b3]
barB_agn_st_ha=[barB_agn_b1,barB_agn_b2,barB_agn_b3]

bottom1_st_ha=np.zeros(3)
bottom2_st_ha=bottom1_st_ha+barB_agn_st_ha
bottom3_st_ha=bottom2_st_ha+barB_noagn_st_ha
bottom4_st_ha=bottom3_st_ha+barAB_agn_st_ha
bottom5_st_ha=bottom4_st_ha+barAB_noagn_st_ha
bottom6_st_ha=bottom5_st_ha+barA_agn_st_ha

top_st_ha=bottom4_st_ha+barA_noagn_st_ha
agns_st_ha=np.array(barA_agn_st_ha)+np.array(barAB_agn_st_ha)+np.array(barB_agn_st_ha)

porcentajes_st_ha=np.round(np.nan_to_num((top_st_ha/np.sum(top_st_ha))*100, nan=0))
porcentajes_agn_st_ha=np.round( np.nan_to_num((agns_st_ha/top_st_ha)*100,nan=0))

barlabels_st_ha=[]
for i in range(len(b)):
	a_st_ha=str(int(porcentajes_st_ha[i])) +'% (AGN:' + str(int(porcentajes_agn_st_ha[i])) +'% )'
	barlabels_st_ha.append(a_st_ha)

####################################################################################
###########################   KC ESTELAR - KC OIII #################################
####################################################################################

circle1_st_oiii=plt.Circle((0,0),fiber_size,fill=False,color='darkviolet', linewidth=2)
circle2_st_oiii=plt.Circle((0,0),fiber_size*2,fill=False,color='darkviolet', linewidth=2,linestyle='dashed')

d_st_oiii=np.sqrt((t['KCRA_ST'] - t['KCRA_OIII'])**2+(t['KCDEC_ST']- t['KCDEC_OIII'])**2)
#BAR A + INTERACTING
barA_noagn_b1=np.sum(d_st_oiii[mask_noflag & mask_sig & mask_barA & mask_noagn]<fiber_size)
barA_noagn_b2=np.sum((d_st_oiii[mask_noflag & mask_sig & mask_barA & mask_noagn]>fiber_size) & (d_st_oiii[mask_noflag & mask_sig & mask_barA & mask_noagn]<2*fiber_size)) 
barA_noagn_b3=np.sum(d_st_oiii[mask_noflag & mask_sig & mask_barA & mask_noagn]>2*fiber_size)
#BAR A + ISOLATED
barA_agn_b1=np.sum(d_st_oiii[mask_noflag & mask_sig & mask_barA & mask_agn]<fiber_size)
barA_agn_b2=np.sum((d_st_oiii[mask_noflag & mask_sig & mask_barA & mask_agn]>fiber_size) & (d_st_oiii[mask_noflag & mask_sig & mask_barA & mask_agn]<2*fiber_size)) 
barA_agn_b3=np.sum(d_st_oiii[mask_noflag & mask_sig & mask_barA & mask_agn]>2*fiber_size)
#BAR AB + INTERACTING
barAB_noagn_b1=np.sum(d_st_oiii[mask_noflag & mask_sig & mask_barAB & mask_noagn]<fiber_size)
barAB_noagn_b2=np.sum((d_st_oiii[mask_noflag & mask_sig & mask_barAB & mask_noagn]>fiber_size) & (d_st_oiii[mask_noflag & mask_sig & mask_barAB & mask_noagn]<2*fiber_size)) 
barAB_noagn_b3=np.sum(d_st_oiii[mask_noflag & mask_sig & mask_barAB & mask_noagn]>2*fiber_size)
#BAR AB + ISOLATED
barAB_agn_b1=np.sum(d_st_oiii[mask_noflag & mask_sig & mask_barAB & mask_agn]<fiber_size)
barAB_agn_b2=np.sum((d_st_oiii[mask_noflag & mask_sig & mask_barAB & mask_agn]>fiber_size) & (d_st_oiii[mask_noflag & mask_sig & mask_barAB & mask_agn]<2*fiber_size)) 
barAB_agn_b3=np.sum(d_st_oiii[mask_noflag & mask_sig & mask_barAB & mask_agn]>2*fiber_size)
#BAR B + INTERACTING
barB_noagn_b1=np.sum(d_st_oiii[mask_noflag & mask_sig & mask_barB & mask_noagn]<fiber_size)
barB_noagn_b2=np.sum((d_st_oiii[mask_noflag & mask_sig & mask_barB & mask_noagn]>fiber_size) & (d_st_oiii[mask_noflag & mask_sig & mask_barB & mask_noagn]<2*fiber_size)) 
barB_noagn_b3=np.sum(d_st_oiii[mask_noflag & mask_sig & mask_barB & mask_noagn]>2*fiber_size)
#BAR B + ISOLATED
barB_agn_b1=np.sum(d_st_oiii[mask_noflag & mask_sig & mask_barB & mask_agn]<fiber_size)
barB_agn_b2=np.sum((d_st_oiii[mask_noflag & mask_sig & mask_barB & mask_agn]>fiber_size) & (d_st_oiii[mask_noflag & mask_sig & mask_barB & mask_agn]<2*fiber_size)) 
barB_agn_b3=np.sum(d_st_oiii[mask_noflag & mask_sig & mask_barB & mask_agn]>2*fiber_size)

barA_noagn_st_oiii=[barA_noagn_b1,barA_noagn_b2,barA_noagn_b3]
barA_agn_st_oiii=[barA_agn_b1,barA_agn_b2,barA_agn_b3]
barAB_noagn_st_oiii=[barAB_noagn_b1,barAB_noagn_b2,barAB_noagn_b3]
barAB_agn_st_oiii=[barAB_agn_b1,barAB_agn_b2,barAB_agn_b3]
barB_noagn_st_oiii=[barB_noagn_b1,barB_noagn_b2,barB_noagn_b3]
barB_agn_st_oiii=[barB_agn_b1,barB_agn_b2,barB_agn_b3]

bottom1_st_oiii=np.zeros(3)
bottom2_st_oiii=bottom1_st_oiii+barB_agn_st_oiii
bottom3_st_oiii=bottom2_st_oiii+barB_noagn_st_oiii
bottom4_st_oiii=bottom3_st_oiii+barAB_agn_st_oiii
bottom5_st_oiii=bottom4_st_oiii+barAB_noagn_st_oiii
bottom6_st_oiii=bottom5_st_oiii+barA_agn_st_oiii

top_st_oiii=bottom4_st_oiii+barA_noagn_st_oiii
agns_st_oiii=np.array(barA_agn_st_oiii)+np.array(barAB_agn_st_oiii)+np.array(barB_agn_st_oiii)

porcentajes_st_oiii=np.round(np.nan_to_num((top_st_oiii/np.sum(top_st_oiii))*100, nan=0))
porcentajes_agn_st_oiii=np.round( np.nan_to_num((agns_st_oiii/top_st_oiii)*100, nan=0))

barlabels_st_oiii=[]
for i in range(len(b)):
	a_st_oiii=str(int(porcentajes_st_oiii[i])) +'% (AGN:' + str(int(porcentajes_agn_st_oiii[i])) +'% )'
	barlabels_st_oiii.append(a_st_oiii)
	
	
####################################################################################
###########################   KC HA - KC OIII #################################
####################################################################################

circle1_ha_oiii=plt.Circle((0,0),fiber_size,fill=False,color='darkviolet', linewidth=2)
circle2_ha_oiii=plt.Circle((0,0),fiber_size*2,fill=False,color='darkviolet', linewidth=2,linestyle='dashed')

d_ha_oiii=np.sqrt((t['KCRA_HA'] - t['KCRA_OIII'])**2+(t['KCDEC_HA']- t['KCDEC_OIII'])**2)
#BAR A + INTERACTING
barA_noagn_b1=np.sum(d_ha_oiii[mask_noflag & mask_sig & mask_barA & mask_noagn]<fiber_size)
barA_noagn_b2=np.sum((d_ha_oiii[mask_noflag & mask_sig & mask_barA & mask_noagn]>fiber_size) & (d_ha_oiii[mask_noflag & mask_sig & mask_barA & mask_noagn]<2*fiber_size)) 
barA_noagn_b3=np.sum(d_ha_oiii[mask_noflag & mask_sig & mask_barA & mask_noagn]>2*fiber_size)
#BAR A + ISOLATED
barA_agn_b1=np.sum(d_ha_oiii[mask_noflag & mask_sig & mask_barA & mask_agn]<fiber_size)
barA_agn_b2=np.sum((d_ha_oiii[mask_noflag & mask_sig & mask_barA & mask_agn]>fiber_size) & (d_ha_oiii[mask_noflag & mask_sig & mask_barA & mask_agn]<2*fiber_size)) 
barA_agn_b3=np.sum(d_ha_oiii[mask_noflag & mask_sig & mask_barA & mask_agn]>2*fiber_size)
#BAR AB + INTERACTING
barAB_noagn_b1=np.sum(d_ha_oiii[mask_noflag & mask_sig & mask_barAB & mask_noagn]<fiber_size)
barAB_noagn_b2=np.sum((d_ha_oiii[mask_noflag & mask_sig & mask_barAB & mask_noagn]>fiber_size) & (d_ha_oiii[mask_noflag & mask_sig & mask_barAB & mask_noagn]<2*fiber_size)) 
barAB_noagn_b3=np.sum(d_ha_oiii[mask_noflag & mask_sig & mask_barAB & mask_noagn]>2*fiber_size)
#BAR AB + ISOLATED
barAB_agn_b1=np.sum(d_ha_oiii[mask_noflag & mask_sig & mask_barAB & mask_agn]<fiber_size)
barAB_agn_b2=np.sum((d_ha_oiii[mask_noflag & mask_sig & mask_barAB & mask_agn]>fiber_size) & (d_ha_oiii[mask_noflag & mask_sig & mask_barAB & mask_agn]<2*fiber_size)) 
barAB_agn_b3=np.sum(d_ha_oiii[mask_noflag & mask_sig & mask_barAB & mask_agn]>2*fiber_size)
#BAR B + INTERACTING
barB_noagn_b1=np.sum(d_ha_oiii[mask_noflag & mask_sig & mask_barB & mask_noagn]<fiber_size)
barB_noagn_b2=np.sum((d_ha_oiii[mask_noflag & mask_sig & mask_barB & mask_noagn]>fiber_size) & (d_ha_oiii[mask_noflag & mask_sig & mask_barB & mask_noagn]<2*fiber_size)) 
barB_noagn_b3=np.sum(d_ha_oiii[mask_noflag & mask_sig & mask_barB & mask_noagn]>2*fiber_size)
#BAR B + ISOLATED
barB_agn_b1=np.sum(d_ha_oiii[mask_noflag & mask_sig & mask_barB & mask_agn]<fiber_size)
barB_agn_b2=np.sum((d_ha_oiii[mask_noflag & mask_sig & mask_barB & mask_agn]>fiber_size) & (d_ha_oiii[mask_noflag & mask_sig & mask_barB & mask_agn]<2*fiber_size)) 
barB_agn_b3=np.sum(d_ha_oiii[mask_noflag & mask_sig & mask_barB & mask_agn]>2*fiber_size)

barA_noagn_ha_oiii=[barA_noagn_b1,barA_noagn_b2,barA_noagn_b3]
barA_agn_ha_oiii=[barA_agn_b1,barA_agn_b2,barA_agn_b3]
barAB_noagn_ha_oiii=[barAB_noagn_b1,barAB_noagn_b2,barAB_noagn_b3]
barAB_agn_ha_oiii=[barAB_agn_b1,barAB_agn_b2,barAB_agn_b3]
barB_noagn_ha_oiii=[barB_noagn_b1,barB_noagn_b2,barB_noagn_b3]
barB_agn_ha_oiii=[barB_agn_b1,barB_agn_b2,barB_agn_b3]

bottom1_ha_oiii=np.zeros(3)
bottom2_ha_oiii=bottom1_ha_oiii+barB_agn_ha_oiii
bottom3_ha_oiii=bottom2_ha_oiii+barB_noagn_ha_oiii
bottom4_ha_oiii=bottom3_ha_oiii+barAB_agn_ha_oiii
bottom5_ha_oiii=bottom4_ha_oiii+barAB_noagn_ha_oiii
bottom6_ha_oiii=bottom5_ha_oiii+barA_agn_ha_oiii

top_ha_oiii=bottom4_ha_oiii+barA_noagn_ha_oiii
agns_ha_oiii=np.array(barA_agn_ha_oiii)+np.array(barAB_agn_ha_oiii)+np.array(barB_agn_ha_oiii)

porcentajes_ha_oiii=np.round(np.nan_to_num((top_ha_oiii/np.sum(top_ha_oiii))*100, nan=0))
porcentajes_agn_ha_oiii=np.round( np.nan_to_num((agns_ha_oiii/top_ha_oiii)*100, nan=0))

barlabels_ha_oiii=[]
for i in range(len(b)):
	a_ha_oiii=str(int(porcentajes_ha_oiii[i])) +'% (AGN:' + str(int(porcentajes_agn_ha_oiii[i])) +'% )'
	barlabels_ha_oiii.append(a_ha_oiii)

####################################################################################
###########################  Hacemos los gráficos ##################################
####################################################################################

fig,axs = plt.subplots(2,3,figsize=(16,10))

axs[0,0].set_title(r'Stellar and H$\alpha$ Kinematics')
axs[0,0].set_xlim(-11.5, 11.5)
axs[0,0].set_ylim(-11.5, 11.5)
axs[0,0].set_xlabel(r'x-offset [KC(st)-KC(H$\alpha$)]')
axs[0,0].set_ylabel(r'y-offset [KC(st)-KC(H$\alpha$)]')
axs[0,0].scatter(0,0,marker='+',color='darkviolet',zorder=10)
axs[0,0].add_artist(circle1_st_ha)
axs[0,0].add_artist(circle2_st_ha)
axs[0,0].errorbar((t['KCRA_ST']-t['KCRA_HA'])[mask_noagn],  (t['KCDEC_ST']-t['KCDEC_HA'])[mask_noagn], yerr=np.sqrt((t['KCDEC_ST_ERR'])**2 + (t['KCDEC_HA_ERR'])**2)[mask_noagn],  xerr=np.sqrt((t['KCRA_ST_ERR'])**2 + (t['KCRA_HA_ERR'])**2)[mask_noagn],  fmt='ks',fillstyle='none',zorder=-1,label='Not confirmed AGN')
axs[0,0].errorbar((t['KCRA_ST']-t['KCRA_HA'])[mask_agn],  (t['KCDEC_ST']-t['KCDEC_HA'])[mask_agn], yerr=np.sqrt((t['KCDEC_ST_ERR'])**2 + (t['KCDEC_HA_ERR'])**2)[mask_agn],  xerr=np.sqrt((t['KCRA_ST_ERR'])**2 + (t['KCRA_HA_ERR'])**2)[mask_agn],  fmt='bs',fillstyle='none',zorder=0,label='Confirmed AGN')
axs[0,0].legend(loc='upper right')

p6_st_ha=axs[1,0].bar(b,barA_noagn_st_ha, bottom=bottom6_st_ha, color='steelblue',label='Not AGN + Bar A')
p5_st_ha=axs[1,0].bar(b,barA_agn_st_ha, bottom=bottom5_st_ha, color='rebeccapurple', label='AGN + Bar A')
p4_st_ha=axs[1,0].bar(b,barAB_noagn_st_ha, bottom=bottom4_st_ha, color='lightskyblue', label='Not AGN + Bar AB')
p3_st_ha=axs[1,0].bar(b,barAB_agn_st_ha, bottom=bottom3_st_ha, color='blueviolet',label='AGN + Bar AB')
p2_st_ha=axs[1,0].bar(b,barB_noagn_st_ha, bottom=bottom2_st_ha, color='paleturquoise',label='Not AGN + Bar B')
p1_st_ha=axs[1,0].bar(b,barB_agn_st_ha,bottom=bottom1_st_ha,color='mediumpurple',label='AGN + Bar B')

axs[1,0].bar_label(p6_st_ha,labels=barlabels_st_ha)
axs[1,0].set_xticks(b, labels=['d<2 \" ','2 \" d < 4 \"','d > 4 \"'])
axs[1,0].legend(loc='upper right')
axs[1,0].set_xlabel('Number of galaxies')
axs[1,0].set_ylabel(r'KC(st)-KC(H$\alpha$) distance')

axs[0,1].set_title('Stellar and [OIII] Kinematics')
axs[0,1].set_xlim(-11.5, 11.5)
axs[0,1].set_ylim(-11.5, 11.5)
axs[0,1].set_xlabel('x-offset [KC(st)-KC([OIII])]')
axs[0,1].set_ylabel('y-offset [KC(st)-KC([OIII])]')
axs[0,1].scatter(0,0,marker='+',color='darkviolet',zorder=10)
axs[0,1].add_artist(circle1_st_oiii)
axs[0,1].add_artist(circle2_st_oiii)
axs[0,1].errorbar((t['KCRA_ST']-t['KCRA_OIII'])[mask_noagn],  (t['KCDEC_ST']-t['KCDEC_OIII'])[mask_noagn], yerr=np.sqrt((t['KCDEC_ST_ERR'])**2 + (t['KCDEC_OIII_ERR'])**2)[mask_noagn],  xerr=np.sqrt((t['KCRA_ST_ERR'])**2 + (t['KCRA_OIII_ERR'])**2)[mask_noagn],  fmt='ks',fillstyle='none',zorder=-1,label='Not Confirmed AGN')
axs[0,1].errorbar((t['KCRA_ST']-t['KCRA_OIII'])[mask_agn],  (t['KCDEC_ST']-t['KCDEC_OIII'])[mask_agn], yerr=np.sqrt((t['KCDEC_ST_ERR'])**2 + (t['KCDEC_OIII_ERR'])**2)[mask_agn],  xerr=np.sqrt((t['KCRA_ST_ERR'])**2 + (t['KCRA_OIII_ERR'])**2)[mask_agn],  fmt='bs',fillstyle='none',zorder=0,label='Confirmed AGN')
#axs[0,1].errorbar((t['KCRA_ST']-t['KCRA_OIII'])[mask_agn],  (t['KCDEC_ST']-t['KCDEC_OIII'])[mask_agn], yerr=np.sqrt((t['KCDEC_ST_ERR'])**2 + (t['KCDEC_OIII_ERR'])**2)[mask_agn],  xerr=np.sqrt((t['KCRA_ST_ERR'])**2 + (t['KCRA_OIII_ERR'])**2)[mask_agn],  fmt='gs',fillstyle='none',zorder=0)
axs[0,1].legend(loc='upper right')

p6_st_oiii=axs[1,1].bar(b,barA_noagn_st_oiii, bottom=bottom6_st_oiii, color='steelblue', label='Not AGN + Bar A')
p5_st_oiii=axs[1,1].bar(b,barA_agn_st_oiii, bottom=bottom5_st_oiii, color='rebeccapurple', label='AGN + Bar A')
p4_st_oiii=axs[1,1].bar(b,barAB_noagn_st_oiii, bottom=bottom4_st_oiii, color='lightskyblue', label='Not AGN + Bar AB')
p3_st_oiii=axs[1,1].bar(b,barAB_agn_st_oiii, bottom=bottom3_st_oiii, color='blueviolet', label='AGN + Bar AB')
p2_st_oiii=axs[1,1].bar(b,barB_noagn_st_oiii,  bottom=bottom2_st_oiii, color='paleturquoise', label='Not AGN + Bar B')
p1_st_oiii=axs[1,1].bar(b,barB_agn_st_oiii,bottom=bottom1_st_oiii,color='mediumpurple',label='AGN + Bar B')

axs[1,1].bar_label(p6_st_oiii,labels=barlabels_st_oiii)
axs[1,1].set_xticks(b, labels=['d<2 \" ','2 \" d < 4 \"','d > 4 \"'])
axs[1,1].legend(loc='upper right')
axs[1,1].set_xlabel('Number of galaxies')
axs[1,1].set_ylabel('KC(st)-KC([OIII]) distance')

axs[0,2].set_title(r'H$\alpha$ and [OIII] Kinematics')
axs[0,2].set_xlim(-9.5, 9.5)
axs[0,2].set_ylim(-9.5, 9.5)
axs[0,2].set_xlabel(r'x-offset [KC(H$\alpha$)-KC(H$\alpha$)]')
axs[0,2].set_ylabel(r'y-offset [KC(H$\alpha$)-KC(H$\alpha$)]')
axs[0,2].scatter(0,0,marker='+',color='darkviolet',zorder=10)
axs[0,2].add_artist(circle1_ha_oiii)
axs[0,2].add_artist(circle2_ha_oiii)
axs[0,2].errorbar((t['KCRA_HA']-t['KCRA_OIII'])[mask_noagn],  (t['KCDEC_HA']-t['KCDEC_OIII'])[mask_noagn], yerr=np.sqrt((t['KCDEC_HA_ERR'])**2 + (t['KCDEC_OIII_ERR'])**2)[mask_noagn],  xerr=np.sqrt((t['KCRA_HA_ERR'])**2 + (t['KCRA_OIII_ERR'])**2)[mask_noagn],  fmt='ks',fillstyle='none',zorder=-1,label='Not Confirmed AGN')
axs[0,2].errorbar((t['KCRA_HA']-t['KCRA_OIII'])[mask_agn],  (t['KCDEC_HA']-t['KCDEC_OIII'])[mask_agn], yerr=np.sqrt((t['KCDEC_HA_ERR'])**2 + (t['KCDEC_OIII_ERR'])**2)[mask_agn],  xerr=np.sqrt((t['KCRA_HA_ERR'])**2 + (t['KCRA_OIII_ERR'])**2)[mask_agn],  fmt='bs',fillstyle='none',zorder=0,label='Confirmed AGN')
#axs[0,2].errorbar((t['KCRA_HA']-t['KCRA_OIII'])[mask_agn],  (t['KCDEC_HA']-t['KCDEC_OIII'])[mask_agn], yerr=np.sqrt((t['KCDEC_HA_ERR'])**2 + (t['KCDEC_OIII_ERR'])**2)[mask_agn],  xerr=np.sqrt((t['KCRA_HA_ERR'])**2 + (t['KCRA_OIII_ERR'])**2)[mask_agn],  fmt='gs',fillstyle='none',zorder=0)
axs[0,2].legend(loc='upper right')

p6_ha_oiii=axs[1,2].bar(b,barA_noagn_ha_oiii,  bottom=bottom6_ha_oiii, color='steelblue', label='Not AGN + Bar A')
p5_ha_oiii=axs[1,2].bar(b,barA_agn_ha_oiii, bottom=bottom5_ha_oiii,color='rebeccapurple',label='AGN + Bar A')
p4_ha_oiii=axs[1,2].bar(b,barAB_noagn_ha_oiii,  bottom=bottom4_ha_oiii, color='lightskyblue', label='Not AGN + Bar AB')
p3_ha_oiii=axs[1,2].bar(b,barAB_agn_ha_oiii, bottom=bottom3_ha_oiii,color='blueviolet',label='AGN + Bar AB')
p2_ha_oiii=axs[1,2].bar(b,barB_noagn_ha_oiii, bottom=bottom2_ha_oiii, color='paleturquoise', label='Not AGN + Bar B')
p1_ha_oiii=axs[1,2].bar(b,barB_agn_ha_oiii, bottom=bottom1_ha_oiii, color='mediumpurple',label='AGN + Bar B')

axs[1,2].bar_label(p6_ha_oiii,labels=barlabels_ha_oiii)
axs[1,2].set_xticks(b, labels=['d<2 \" ','2 \" d < 4 \"','d > 4 \"'])
axs[1,2].legend(loc='upper right')
axs[1,2].set_xlabel('Number of galaxies')
axs[1,2].set_ylabel(r'KC(H$\alpha$)-KC([OIII]) distance')

#fig.suptitle('SDSS-BASED ISOLATED GALAXIES')
fig.savefig('sig_offset_02.png')

plt.show()

