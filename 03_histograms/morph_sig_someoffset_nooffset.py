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

t=Table.read('tabla_mangaagn_updated.ecsv',format='ascii.ecsv')

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

b3=[1, 2, 3]
b5=[1,2,3,4,5]

d_stellar=np.sqrt((t['KCRA_ST'])**2+(t['KCDEC_ST'])**2)
d_ha=np.sqrt((t['KCRA_HA'])**2+(t['KCDEC_HA'])**2)
d_oiii=np.sqrt((t['KCRA_OIII'])**2+(t['KCDEC_OIII'])**2)
d_st_ha=np.sqrt((t['KCRA_ST'] - t['KCRA_HA'])**2+(t['KCDEC_ST']- t['KCDEC_HA'])**2)
d_st_oiii=np.sqrt((t['KCRA_ST'] - t['KCRA_OIII'])**2+(t['KCDEC_ST']- t['KCDEC_OIII'])**2)
d_ha_oiii=np.sqrt((t['KCRA_HA'] - t['KCRA_OIII'])**2+(t['KCDEC_HA']- t['KCDEC_OIII'])**2)

mask_nooffset=((d_stellar<fiber_size) & (d_ha<fiber_size) & (d_oiii<fiber_size) & (d_st_ha<fiber_size) & (d_st_oiii<fiber_size) & (d_ha_oiii<fiber_size) )
mask_someoffset=((d_stellar>2*fiber_size) | (d_ha>2*fiber_size) | (d_oiii>2*fiber_size) | (d_st_ha>2*fiber_size) | (d_st_oiii>2*fiber_size) | (d_ha_oiii>2*fiber_size) ) 

print('Galaxias sin offset: ', np.sum(mask_nooffset & mask_noflag & mask_sig))
print('Galaxias con algún offset: ', np.sum(mask_someoffset & mask_noflag & mask_sig))

####################################################################################
####################### SAMPLE: SIG GALAXIES WITH NO OFFSET ########################
####################################################################################

####################################################################################
####################### MORPHOLOGICAL TYPE  ########################################
####################################################################################

#BAR A + NOT AGN
barA_noagn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_e)
barA_noagn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_s0)
barA_noagn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_s1)
barA_noagn_b4=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_s2)
barA_noagn_b5=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_irr)

#BAR A + AGN
barA_agn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_e)
barA_agn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_s0) 
barA_agn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_s1)
barA_agn_b4=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_s2) 
barA_agn_b5=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_irr)

#BAR AB + NOT AGN
barAB_noagn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_e)
barAB_noagn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_s0) 
barAB_noagn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_s1)
barAB_noagn_b4=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_s2) 
barAB_noagn_b5=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_irr)

#BAR AB + AGN
barAB_agn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_e)
barAB_agn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_s0) 
barAB_agn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_s1)
barAB_agn_b4=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_s2) 
barAB_agn_b5=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_irr)

#BAR B + NOT AGN
barB_noagn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_e)
barB_noagn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_s0) 
barB_noagn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_s1)
barB_noagn_b4=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_s2) 
barB_noagn_b5=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_irr)

#BAR B + AGN
barB_agn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_e)
barB_agn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_s0) 
barB_agn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_s1)
barB_agn_b4=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_s2) 
barB_agn_b5=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_irr)

barA_noagn_MORPH=[barA_noagn_b1,barA_noagn_b2,barA_noagn_b3,barA_noagn_b4,barA_noagn_b5]
barA_agn_MORPH=[barA_agn_b1,barA_agn_b2,barA_agn_b3,barA_agn_b4,barA_agn_b5]
barAB_noagn_MORPH=[barAB_noagn_b1,barAB_noagn_b2,barAB_noagn_b3,barAB_noagn_b4,barAB_noagn_b5]
barAB_agn_MORPH=[barAB_agn_b1,barAB_agn_b2,barAB_agn_b3,barAB_agn_b4,barAB_agn_b5]
barB_noagn_MORPH=[barB_noagn_b1,barB_noagn_b2,barB_noagn_b3,barAB_noagn_b4,barAB_noagn_b5]
barB_agn_MORPH=[barB_agn_b1,barB_agn_b2,barB_agn_b3,barAB_agn_b4,barAB_agn_b5]

bottom1_MORPH=np.zeros(5)
bottom2_MORPH=bottom1_MORPH+barB_agn_MORPH
bottom3_MORPH=bottom2_MORPH+barB_noagn_MORPH
bottom4_MORPH=bottom3_MORPH+barAB_agn_MORPH
bottom5_MORPH=bottom4_MORPH+barAB_noagn_MORPH
bottom6_MORPH=bottom5_MORPH+barA_agn_MORPH

top_MORPH=bottom6_MORPH+barA_noagn_MORPH
agns_MORPH=np.array(barA_agn_MORPH)+np.array(barAB_agn_MORPH)+np.array(barB_agn_MORPH)

porcentajes_MORPH=np.round(np.nan_to_num((top_MORPH/np.sum(top_MORPH))*100, nan=0))
porcentajes_agn_MORPH=np.round(np.nan_to_num((agns_MORPH/top_MORPH)*100,nan=0))

barlabels_MORPH=[]
for i in range(len(b5)):
	a_MORPH=str(int(porcentajes_MORPH[i]))+'% (AGN:' + str(int(porcentajes_agn_MORPH[i]))+'%)'
	barlabels_MORPH.append(a_MORPH)
	
#NOW WE'LL ONLY STUDY SIG GALAXIES -> ISOLATED ONES

####################################################################################
####################### STELLAR MASS  ##############################################
####################################################################################

#BAR A + NOT AGN
barA_noagn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_mass1)
barA_noagn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_mass2)
barA_noagn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_mass3)
barA_noagn_b4=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_mass4)
barA_noagn_b5=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_mass5)

#BAR A + ISOLATED
barA_agn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_mass1)
barA_agn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_mass2) 
barA_agn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_mass3)
barA_agn_b4=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_mass4) 
barA_agn_b5=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_mass5)

#BAR AB + INTERACTING
barAB_noagn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_mass1)
barAB_noagn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_mass2) 
barAB_noagn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_mass3)
barAB_noagn_b4=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_mass4) 
barAB_noagn_b5=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_mass5)

#BAR AB + ISOLATED
barAB_agn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_mass1)
barAB_agn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_mass2) 
barAB_agn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_mass3)
barAB_agn_b4=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_mass4) 
barAB_agn_b5=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_mass5)

#BAR B + INTERACTING
barB_noagn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_mass1)
barB_noagn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_mass2) 
barB_noagn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_mass3)
barB_noagn_b4=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_mass4) 
barB_noagn_b5=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_mass5)

#BAR B + ISOLATED
barB_agn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_mass1)
barB_agn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_mass2) 
barB_agn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_mass3)
barB_agn_b4=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_mass4) 
barB_agn_b5=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_mass5)

barA_noagn_MASS=[barA_noagn_b1,barA_noagn_b2,barA_noagn_b3,barA_noagn_b4,barA_noagn_b5]
barA_agn_MASS=[barA_agn_b1,barA_agn_b2,barA_agn_b3,barA_agn_b4,barA_agn_b5]
barAB_noagn_MASS=[barAB_noagn_b1,barAB_noagn_b2,barAB_noagn_b3,barAB_noagn_b4,barAB_noagn_b5]
barAB_agn_MASS=[barAB_agn_b1,barAB_agn_b2,barAB_agn_b3,barAB_agn_b4,barAB_agn_b5]
barB_noagn_MASS=[barB_noagn_b1,barB_noagn_b2,barB_noagn_b3,barAB_noagn_b4,barAB_noagn_b5]
barB_agn_MASS=[barB_agn_b1,barB_agn_b2,barB_agn_b3,barAB_agn_b4,barAB_agn_b5]

bottom1_MASS=np.zeros(5)
bottom2_MASS=bottom1_MASS+barB_agn_MASS
bottom3_MASS=bottom2_MASS+barB_noagn_MASS
bottom4_MASS=bottom3_MASS+barAB_agn_MASS
bottom5_MASS=bottom4_MASS+barAB_noagn_MASS
bottom6_MASS=bottom5_MASS+barA_agn_MASS

top_MASS=bottom6_MASS+barA_noagn_MASS
agns_MASS=np.array(barA_agn_MASS)+np.array(barAB_agn_MASS)+np.array(barB_agn_MASS)

porcentajes_MASS=np.round(np.nan_to_num((top_MASS/np.sum(top_MASS))*100, nan=0))
porcentajes_agn_MASS=np.round(np.nan_to_num((agns_MASS/top_MASS)*100,nan=0))

barlabels_MASS=[]
for i in range(len(b5)):
	a_MASS=str(int(porcentajes_MASS[i]))+'% (AGN:' + str(int(porcentajes_agn_MASS[i]))+'%)'
	barlabels_MASS.append(a_MASS)
	
####################################################################################
####################### REDSHIFT Z  ################################################
####################################################################################

#BAR A + NOT AGN
barA_noagn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_z1)
barA_noagn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_z2)
barA_noagn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_z3)

#BAR A + ISOLATED
barA_agn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_z1)
barA_agn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_z2) 
barA_agn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_z3)

#BAR AB + INTERACTING
barAB_noagn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_z1)
barAB_noagn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_z2) 
barAB_noagn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_z3)

#BAR AB + ISOLATED
barAB_agn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_z1)
barAB_agn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_z2) 
barAB_agn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_z3)

#BAR B + INTERACTING
barB_noagn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_z1)
barB_noagn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_z2) 
barB_noagn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_z3)

#BAR B + ISOLATED
barB_agn_b1=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_z1)
barB_agn_b2=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_z2) 
barB_agn_b3=np.sum(mask_nooffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_z3)

barA_noagn_Z=[barA_noagn_b1,barA_noagn_b2,barA_noagn_b3]
barA_agn_Z=[barA_agn_b1,barA_agn_b2,barA_agn_b3]
barAB_noagn_Z=[barAB_noagn_b1,barAB_noagn_b2,barAB_noagn_b3]
barAB_agn_Z=[barAB_agn_b1,barAB_agn_b2,barAB_agn_b3]
barB_noagn_Z=[barB_noagn_b1,barB_noagn_b2,barB_noagn_b3]
barB_agn_Z=[barB_agn_b1,barB_agn_b2,barB_agn_b3]

bottom1_Z=np.zeros(3)
bottom2_Z=bottom1_Z+barB_agn_Z
bottom3_Z=bottom2_Z+barB_noagn_Z
bottom4_Z=bottom3_Z+barAB_agn_Z
bottom5_Z=bottom4_Z+barAB_noagn_Z
bottom6_Z=bottom5_Z+barA_agn_Z

top_Z=bottom6_Z+barA_noagn_Z
agns_Z=np.array(barA_agn_Z)+np.array(barAB_agn_Z)+np.array(barB_agn_Z)

porcentajes_Z=np.round(np.nan_to_num((top_Z/np.sum(top_Z))*100, nan=0))
porcentajes_agn_Z=np.round(np.nan_to_num((agns_Z/top_Z)*100,nan=0))

barlabels_Z=[]
for i in range(len(b3)):
	a_Z=str(int(porcentajes_Z[i]))+'% (AGN:' + str(int(porcentajes_agn_Z[i]))+'%)'
	barlabels_Z.append(a_Z)
		
####################################################################################
###########################  Hacemos los gráficos ##################################
####################################################################################

fig,axs = plt.subplots(1,3,figsize=(30,8))

p6_MORPH=axs[0].bar(b5,barA_noagn_MORPH, bottom=bottom6_MORPH,color='steelblue',label='Not AGN + Bar A')
p5_MORPH=axs[0].bar(b5,barA_agn_MORPH, bottom=bottom5_MORPH,color='rebeccapurple',label='AGN + Bar A')
p4_MORPH=axs[0].bar(b5,barAB_noagn_MORPH, bottom=bottom4_MORPH,color='lightskyblue',label='Not AGN + Bar AB')
p3_MORPH=axs[0].bar(b5,barAB_agn_MORPH, bottom=bottom3_MORPH,color='blueviolet',label='AGN + Bar AB')
p2_MORPH=axs[0].bar(b5,barB_noagn_MORPH, bottom=bottom2_MORPH,color='paleturquoise',label='Not AGN + Bar B')
p1_MORPH=axs[0].bar(b5,barB_agn_MORPH,bottom=bottom1_MORPH,color='mediumpurple',label='AGN + Bar B')

axs[0].bar_label(p6_MORPH,labels=barlabels_MORPH)
axs[0].set_xticks(b5, labels=['E','S0','S1','S2','Irr'])
axs[0].legend(loc='upper left')
axs[0].set_ylabel('Number of galaxies')
axs[0].set_xlabel(r'Morphological type')
axs[0].set_title('Classification by Morphological type')

p6_MASS=axs[1].bar(b5,barA_noagn_MASS, bottom=bottom6_MASS,color='steelblue',label='Not AGN + Bar A')
p5_MASS=axs[1].bar(b5,barA_agn_MASS, bottom=bottom5_MASS,color='rebeccapurple',label='AGN + Bar A')
p4_MASS=axs[1].bar(b5,barAB_noagn_MASS, bottom=bottom4_MASS,color='lightskyblue',label='Not AGN + Bar AB')
p3_MASS=axs[1].bar(b5,barAB_agn_MASS, bottom=bottom3_MASS,color='blueviolet',label='AGN + Bar AB')
p2_MASS=axs[1].bar(b5,barB_noagn_MASS, bottom=bottom2_MASS,color='paleturquoise',label='Not AGN + Bar B')
p1_MASS=axs[1].bar(b5,barB_agn_MASS,bottom=bottom1_MASS,color='mediumpurple',label='AGN + Bar B')

axs[1].bar_label(p6_MASS,labels=barlabels_MASS)
axs[1].set_xticks(b5, labels=['M<10','10<M<10.5','10.5<M<11','11<M<11.5','M>11.5'])
axs[1].legend(loc='upper right')
axs[1].set_ylabel('Number of galaxies')
axs[1].set_xlabel(r'log(Mass) [M_$\odot$]')
axs[1].set_title('Classification by Stellar Mass')

p6_Z=axs[2].bar(b3,barA_noagn_Z, bottom=bottom6_Z,color='steelblue',label='Not AGN + Bar A')
p5_Z=axs[2].bar(b3,barA_agn_Z, bottom=bottom5_Z,color='rebeccapurple',label='AGN + Bar A')
p4_Z=axs[2].bar(b3,barAB_noagn_Z, bottom=bottom4_Z,color='lightskyblue',label='Not AGN + Bar AB')
p3_Z=axs[2].bar(b3,barAB_agn_Z, bottom=bottom3_Z,color='blueviolet',label='AGN + Bar AB')
p2_Z=axs[2].bar(b3,barB_noagn_Z, bottom=bottom2_Z,color='paleturquoise',label='Not AGN + Bar B')
p1_Z=axs[2].bar(b3,barB_agn_Z,bottom=bottom1_Z,color='mediumpurple',label='AGN + Bar B')

axs[2].bar_label(p6_Z,labels=barlabels_Z)
axs[2].set_xticks(b3, labels=['z<0.03','0.03<z<0.05','z>0.05'])
axs[2].legend(loc='upper left')
axs[2].set_ylabel('Number of galaxies')
axs[2].set_xlabel(r'Redshift (z)')
axs[2].set_title('Classification by Redshift')

#fig.suptitle('SDSS-Based Isolated Galaxies with no kinematic offset')
fig.savefig('nooffset_morph.png')

plt.show()



####################################################################################
####################### SAMPLE: SIG GALAXIES WITH SOME OFFSET ######################
####################################################################################

####################################################################################
####################### MORPHOLOGICAL TYPE  ########################################
####################################################################################

#BAR A + NOT AGN
barA_noagn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_e)
barA_noagn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_s0)
barA_noagn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_s1)
barA_noagn_b4=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_s2)
barA_noagn_b5=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_irr)

#BAR A + AGN
barA_agn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_e)
barA_agn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_s0) 
barA_agn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_s1)
barA_agn_b4=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_s2) 
barA_agn_b5=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_irr)

#BAR AB + NOT AGN
barAB_noagn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_e)
barAB_noagn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_s0) 
barAB_noagn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_s1)
barAB_noagn_b4=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_s2) 
barAB_noagn_b5=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_irr)

#BAR AB + AGN
barAB_agn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_e)
barAB_agn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_s0) 
barAB_agn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_s1)
barAB_agn_b4=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_s2) 
barAB_agn_b5=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_irr)

#BAR B + NOT AGN
barB_noagn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_e)
barB_noagn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_s0) 
barB_noagn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_s1)
barB_noagn_b4=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_s2) 
barB_noagn_b5=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_irr)

#BAR B + AGN
barB_agn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_e)
barB_agn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_s0) 
barB_agn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_s1)
barB_agn_b4=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_s2) 
barB_agn_b5=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_irr)

barA_noagn_MORPH=[barA_noagn_b1,barA_noagn_b2,barA_noagn_b3,barA_noagn_b4,barA_noagn_b5]
barA_agn_MORPH=[barA_agn_b1,barA_agn_b2,barA_agn_b3,barA_agn_b4,barA_agn_b5]
barAB_noagn_MORPH=[barAB_noagn_b1,barAB_noagn_b2,barAB_noagn_b3,barAB_noagn_b4,barAB_noagn_b5]
barAB_agn_MORPH=[barAB_agn_b1,barAB_agn_b2,barAB_agn_b3,barAB_agn_b4,barAB_agn_b5]
barB_noagn_MORPH=[barB_noagn_b1,barB_noagn_b2,barB_noagn_b3,barAB_noagn_b4,barAB_noagn_b5]
barB_agn_MORPH=[barB_agn_b1,barB_agn_b2,barB_agn_b3,barAB_agn_b4,barAB_agn_b5]

bottom1_MORPH=np.zeros(5)
bottom2_MORPH=bottom1_MORPH+barB_agn_MORPH
bottom3_MORPH=bottom2_MORPH+barB_noagn_MORPH
bottom4_MORPH=bottom3_MORPH+barAB_agn_MORPH
bottom5_MORPH=bottom4_MORPH+barAB_noagn_MORPH
bottom6_MORPH=bottom5_MORPH+barA_agn_MORPH

top_MORPH=bottom6_MORPH+barA_noagn_MORPH
agns_MORPH=np.array(barA_agn_MORPH)+np.array(barAB_agn_MORPH)+np.array(barB_agn_MORPH)

porcentajes_MORPH=np.round(np.nan_to_num((top_MORPH/np.sum(top_MORPH))*100, nan=0))
porcentajes_agn_MORPH=np.round(np.nan_to_num((agns_MORPH/top_MORPH)*100,nan=0))

barlabels_MORPH=[]
for i in range(len(b5)):
	a_MORPH=str(int(porcentajes_MORPH[i]))+'% (AGN:' + str(int(porcentajes_agn_MORPH[i]))+'%)'
	barlabels_MORPH.append(a_MORPH)
	
#NOW WE'LL ONLY STUDY SIG GALAXIES -> ISOLATED ONES

####################################################################################
####################### STELLAR MASS  ##############################################
####################################################################################

#BAR A + NOT AGN
barA_noagn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_mass1)
barA_noagn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_mass2)
barA_noagn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_mass3)
barA_noagn_b4=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_mass4)
barA_noagn_b5=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_mass5)

#BAR A + AGN
barA_agn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_mass1)
barA_agn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_mass2) 
barA_agn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_mass3)
barA_agn_b4=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_mass4) 
barA_agn_b5=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_mass5)

#BAR AB + NOT AGN
barAB_noagn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_mass1)
barAB_noagn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_mass2) 
barAB_noagn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_mass3)
barAB_noagn_b4=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_mass4) 
barAB_noagn_b5=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_mass5)

#BAR AB + AGN
barAB_agn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_mass1)
barAB_agn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_mass2) 
barAB_agn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_mass3)
barAB_agn_b4=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_mass4) 
barAB_agn_b5=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_mass5)

#BAR B + NOT AGN
barB_noagn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_mass1)
barB_noagn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_mass2) 
barB_noagn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_mass3)
barB_noagn_b4=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_mass4) 
barB_noagn_b5=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_mass5)

#BAR B + AGN
barB_agn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_mass1)
barB_agn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_mass2) 
barB_agn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_mass3)
barB_agn_b4=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_mass4) 
barB_agn_b5=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_mass5)

barA_noagn_MASS=[barA_noagn_b1,barA_noagn_b2,barA_noagn_b3,barA_noagn_b4,barA_noagn_b5]
barA_agn_MASS=[barA_agn_b1,barA_agn_b2,barA_agn_b3,barA_agn_b4,barA_agn_b5]
barAB_noagn_MASS=[barAB_noagn_b1,barAB_noagn_b2,barAB_noagn_b3,barAB_noagn_b4,barAB_noagn_b5]
barAB_agn_MASS=[barAB_agn_b1,barAB_agn_b2,barAB_agn_b3,barAB_agn_b4,barAB_agn_b5]
barB_noagn_MASS=[barB_noagn_b1,barB_noagn_b2,barB_noagn_b3,barAB_noagn_b4,barAB_noagn_b5]
barB_agn_MASS=[barB_agn_b1,barB_agn_b2,barB_agn_b3,barAB_agn_b4,barAB_agn_b5]

bottom1_MASS=np.zeros(5)
bottom2_MASS=bottom1_MASS+barB_agn_MASS
bottom3_MASS=bottom2_MASS+barB_noagn_MASS
bottom4_MASS=bottom3_MASS+barAB_agn_MASS
bottom5_MASS=bottom4_MASS+barAB_noagn_MASS
bottom6_MASS=bottom5_MASS+barA_agn_MASS

top_MASS=bottom6_MASS+barA_noagn_MASS
agns_MASS=np.array(barA_agn_MASS)+np.array(barAB_agn_MASS)+np.array(barB_agn_MASS)

porcentajes_MASS=np.round(np.nan_to_num((top_MASS/np.sum(top_MASS))*100, nan=0))
porcentajes_agn_MASS=np.round(np.nan_to_num((agns_MASS/top_MASS)*100,nan=0))

barlabels_MASS=[]
for i in range(len(b5)):
	a_MASS=str(int(porcentajes_MASS[i]))+'% (AGN:' + str(int(porcentajes_agn_MASS[i]))+'%)'
	barlabels_MASS.append(a_MASS)
	
####################################################################################
####################### REDSHIFT Z  ################################################
####################################################################################

#BAR A + NOT AGN
barA_noagn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_z1)
barA_noagn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_z2)
barA_noagn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_noagn & mask_z3)

#BAR A + AGN
barA_agn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_z1)
barA_agn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_z2) 
barA_agn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barA & mask_agn & mask_z3)

#BAR AB + NOT AGN
barAB_noagn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_z1)
barAB_noagn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_z2) 
barAB_noagn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_noagn & mask_z3)

#BAR AB + AGN
barAB_agn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_z1)
barAB_agn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_z2) 
barAB_agn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barAB & mask_agn & mask_z3)

#BAR B + NOT AGN
barB_noagn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_z1)
barB_noagn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_z2) 
barB_noagn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_noagn & mask_z3)

#BAR B + AGN
barB_agn_b1=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_z1)
barB_agn_b2=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_z2) 
barB_agn_b3=np.sum(mask_someoffset & mask_noflag & mask_sig & mask_barB & mask_agn & mask_z3)

barA_noagn_Z=[barA_noagn_b1,barA_noagn_b2,barA_noagn_b3]
barA_agn_Z=[barA_agn_b1,barA_agn_b2,barA_agn_b3]
barAB_noagn_Z=[barAB_noagn_b1,barAB_noagn_b2,barAB_noagn_b3]
barAB_agn_Z=[barAB_agn_b1,barAB_agn_b2,barAB_agn_b3]
barB_noagn_Z=[barB_noagn_b1,barB_noagn_b2,barB_noagn_b3]
barB_agn_Z=[barB_agn_b1,barB_agn_b2,barB_agn_b3]

bottom1_Z=np.zeros(3)
bottom2_Z=bottom1_Z+barB_agn_Z
bottom3_Z=bottom2_Z+barB_noagn_Z
bottom4_Z=bottom3_Z+barAB_agn_Z
bottom5_Z=bottom4_Z+barAB_noagn_Z
bottom6_Z=bottom5_Z+barA_agn_Z

top_Z=bottom6_Z+barA_noagn_Z
agns_Z=np.array(barA_agn_Z)+np.array(barAB_agn_Z)+np.array(barB_agn_Z)

porcentajes_Z=np.round(np.nan_to_num((top_Z/np.sum(top_Z))*100, nan=0))
porcentajes_agn_Z=np.round(np.nan_to_num((agns_Z/top_Z)*100,nan=0))

barlabels_Z=[]
for i in range(len(b3)):
	a_Z=str(int(porcentajes_Z[i]))+'% (AGN:' + str(int(porcentajes_agn_Z[i]))+'%)'
	barlabels_Z.append(a_Z)
		
####################################################################################
###########################  Hacemos los gráficos ##################################
####################################################################################

fig,axs = plt.subplots(1,3,figsize=(30,8))

p6_MORPH=axs[0].bar(b5,barA_noagn_MORPH, bottom=bottom6_MORPH,color='steelblue',label='Not AGN + Bar A')
p5_MORPH=axs[0].bar(b5,barA_agn_MORPH, bottom=bottom5_MORPH,color='rebeccapurple',label='AGN + Bar A')
p4_MORPH=axs[0].bar(b5,barAB_noagn_MORPH, bottom=bottom4_MORPH,color='lightskyblue',label='Not AGN + Bar AB')
p3_MORPH=axs[0].bar(b5,barAB_agn_MORPH, bottom=bottom3_MORPH,color='blueviolet',label='AGN + Bar AB')
p2_MORPH=axs[0].bar(b5,barB_noagn_MORPH, bottom=bottom2_MORPH,color='paleturquoise',label='Not AGN + Bar B')
p1_MORPH=axs[0].bar(b5,barB_agn_MORPH,bottom=bottom1_MORPH,color='mediumpurple',label='AGN + Bar B')

axs[0].bar_label(p6_MORPH,labels=barlabels_MORPH)
axs[0].set_xticks(b5, labels=['E','S0','S1','S2','Irr'])
axs[0].legend(loc='upper left')
axs[0].set_ylabel('Number of galaxies')
axs[0].set_xlabel(r'Morphological type')
axs[0].set_title('Classification by Morphological type')

p6_MASS=axs[1].bar(b5,barA_noagn_MASS, bottom=bottom6_MASS,color='steelblue',label='Not AGN + Bar A')
p5_MASS=axs[1].bar(b5,barA_agn_MASS, bottom=bottom5_MASS,color='rebeccapurple',label='AGN + Bar A')
p4_MASS=axs[1].bar(b5,barAB_noagn_MASS, bottom=bottom4_MASS,color='lightskyblue',label='Not AGN + Bar AB')
p3_MASS=axs[1].bar(b5,barAB_agn_MASS, bottom=bottom3_MASS,color='blueviolet',label='AGN + Bar AB')
p2_MASS=axs[1].bar(b5,barB_noagn_MASS, bottom=bottom2_MASS,color='paleturquoise',label='Not AGN + Bar B')
p1_MASS=axs[1].bar(b5,barB_agn_MASS,bottom=bottom1_MASS,color='mediumpurple',label='AGN + Bar B')

axs[1].bar_label(p6_MASS,labels=barlabels_MASS)
axs[1].set_xticks(b5, labels=['M<10','10<M<10.5','10.5<M<11','11<M<11.5','M>11.5'])
axs[1].legend(loc='upper right')
axs[1].set_ylabel('Number of galaxies')
axs[1].set_xlabel(r'log(Mass) [M_$\odot$]')
axs[1].set_title('Classification by Stellar Mass')

p6_Z=axs[2].bar(b3,barA_noagn_Z, bottom=bottom6_Z,color='steelblue',label='Not AGN + Bar A')
p5_Z=axs[2].bar(b3,barA_agn_Z, bottom=bottom5_Z,color='rebeccapurple',label='AGN + Bar A')
p4_Z=axs[2].bar(b3,barAB_noagn_Z, bottom=bottom4_Z,color='lightskyblue',label='Not AGN + Bar AB')
p3_Z=axs[2].bar(b3,barAB_agn_Z, bottom=bottom3_Z,color='blueviolet',label='AGN + Bar AB')
p2_Z=axs[2].bar(b3,barB_noagn_Z, bottom=bottom2_Z,color='paleturquoise',label='Not AGN + Bar B')
p1_Z=axs[2].bar(b3,barB_agn_Z,bottom=bottom1_Z,color='mediumpurple',label='AGN + Bar B')

axs[2].bar_label(p6_Z,labels=barlabels_Z)
axs[2].set_ylim(0,32)
axs[2].set_xticks(b3, labels=['z<0.03','0.03<z<0.05','z>0.05'])
axs[2].legend(loc='upper left')
axs[2].set_ylabel('Number of galaxies')
axs[2].set_xlabel(r'Redshift (z)')
axs[2].set_title('Classification by Redshift')

#fig.suptitle('SDSS-Based Isolated Galaxies with some kinematic offset')
fig.savefig('someoffset_morph.png')

plt.show()



