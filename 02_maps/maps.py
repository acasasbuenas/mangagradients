"""
Created on Tue 30 Aug 2022
@author: Alba Casasbuenas
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
import glob

#######################################################################
############ Definimos las funciones que harán falta ##################
#######################################################################

pix_size=0.5 #arcsec
PSF=2.5 #arcsec
fiber_size=2 #arcsec


#Interquartile Range Method: finds the outlier values in an array
#We use this later to mask outlier values in the velocity maps.
def IQR(array):
	Q1,Q3=np.nanpercentile(array,[25,75])
	IQR=Q3-Q1
	ul=Q3+1.5*IQR
	ll=Q1-1.5*IQR
	return([ul,ll])

#######################################################################
############# Muestra =  todas las galaxias de MaNGA ##################
#######################################################################

#FULL MANGA SAMPLE 
data = get_drpall_table()
plateifus_data=data['plateifu']

#Finds the basic parameters of a MaNGA galaxy given its plate-ifu combination.
def basic_parameters(plateifu_ex):
	mangaid_ex = (data['mangaid'][plateifus_data==plateifu_ex])[0]
	redshift_ex = ((data['nsa_z']) [plateifus_data==plateifu_ex])[0]
	dist_ex	= ((data['nsa_zdist']) [plateifus_data==plateifu_ex])[0] #*c/H0 para Mpc
	objra_ex = ((data['objra']) [plateifus_data==plateifu_ex])[0]
	objdec_ex = ((data['objdec']) [plateifus_data==plateifu_ex])[0]
	return(plateifu_ex,mangaid_ex,objra_ex,objdec_ex,redshift_ex,dist_ex)

##############################################################
########### MORPHOLOGICAL CLASSIFICATION #####################
##############################################################

hdul_morph=fits.open('manga-morphology-dl-DR17.fits')
data_morph=hdul_morph[1].data
hdul_morph.close()

hdul_gema=fits.open('manga_GEMA_2.0.2.fits')
data_gema=hdul_gema[2].data
data_gema1=hdul_gema[12].data
hdul_gema.close()

hdul_agn=fits.open('manga_agn-v1_0_1.fits')
data_agn=hdul_agn[1].data
hdul_agn.close()

hdul_pipe3d=fits.open('SDSS17Pipe3D_v3_1_1.fits')
data_pipe3d=hdul_pipe3d[1].data
hdul_pipe3d.close()

#Finds the morphological and agn classification of a MaNGA galaxy given its plate-ifu combination
def morph_class(plateifu_ex):

	ttype=data_morph['T-Type'][data_morph['PLATEIFU']==plateifu_ex]
	pltg=data_morph['P_LTG'][data_morph['PLATEIFU']==plateifu_ex]
	ps0=data_morph['P_S0'][data_morph['PLATEIFU']==plateifu_ex]
	pbar=data_morph['P_bar'][data_morph['PLATEIFU']==plateifu_ex]
	pedgeon=data_morph['P_edge'][data_morph['PLATEIFU']==plateifu_ex]
	vc=data_morph['Visual_Class'][data_morph['PLATEIFU']==plateifu_ex]
	pmerger=data_gema1['p_merger'][data_gema1['mangaid']==str(mangaid_ex_i)]
	mangaid_ex_agn='manga-'+plateifu_ex
	wise=data_agn['WISE_AGN'][data_agn['MANGAID']==mangaid_ex_agn]
	bat=data_agn['BAT_AGN'][data_agn['MANGAID']==mangaid_ex_agn]
	radio=data_agn['RADIO_AGN'][data_agn['MANGAID']==mangaid_ex_agn]
	radioclass=data_agn['RADIO_CLASS'][data_agn['MANGAID']==mangaid_ex_agn] #high excitation radio galaxy (quasar mode) or low excitation radio galaxy (radio mode)
	broad=data_agn['BROAD_AGN'][data_agn['MANGAID']==mangaid_ex_agn]
	
	# MORPHOLOGY:
	if ttype>0 and ttype<3 and pltg >= 0.5 and vc==3:
		morph='S1'
	elif ttype>3 and pltg >= 0.5 and vc==3:
		morph='S2'
	elif ttype<0 and ps0>0.5 and pltg<0.5 and vc==2:
		morph='S0'
	elif ttype<0 and ps0<0.5 and pltg<0.5 and vc==2:
		morph='E'
	else:
		morph='-'	
	#BAR:
	if pbar.size<1:
		pbar=np.nan
	else:
		pbar=pbar[0]
	#EDGE ON:
	if pedgeon.size<1:
		pedgeon=np.nan
	else:
		pedgeon=pedgeon[0]
	#MERGER	
	if pmerger.size>0:
		merger=pmerger[0]
	else:
		merger= np.nan
	#AGN
	agn=0
	if wise.size<1 and bat.size<1 and radio.size<1 and broad.size<1:
		agn=np.nan
	#AGN CLASS			
	agn_class=''
	if np.sum(wise)>0: 
		agn=agn+1
		agn_class=agn_class+' WISE_AGN '
	if np.sum(bat)>0:
		agn=agn+1
		agn_class=agn_class+' BAT_AGN '
	if np.sum(radio)>0:
		agn=agn+1
		agn_class=agn_class+' RADIO_AGN '
	if np.sum(broad)>0:
		agn=agn+1
		agn_class=agn_class+' BROAD_AGN '
	
	return(morph,pbar,pedgeon,merger,agn,agn_class)

def isolation(plateifu_ex):
	mangaid_ex = (data['mangaid'][plateifus_data==plateifu_ex])[0]
	nneigh=data_gema['Nneigh'][data_gema['mangaid']==str(mangaid_ex)]
	qnn=data_gema['Q_nn'][data_gema['mangaid']==str(mangaid_ex)]
	qlss=data_gema['Q_lss'][data_gema['mangaid']==str(mangaid_ex)]
	if nneigh.size<1:
		nneigh=0
	if qnn.size<1:
		qnn=0
	if qlss.size<1:
		qlss=0
	return nneigh,qnn,qlss #sdss-based isolated galaxy -> nneigh=0
	
def mass(plateifu_ex):
	logmass=data_pipe3d['log_Mass'][data_pipe3d['plateifu']==plateifu_ex]
	return logmass

##############################################################
############### KINEMATIC ANALYSIS ###########################
##############################################################

#definimos el array dist_arr para hacer el analisis con varios valores de dist
#luego, el centro cinematico sera la media y el error la desviacion tipica

dist_arr=np.arange((PSF/2)/pix_size, (PSF*2)/pix_size, 0.5)
#vamos a poner el de la PSF el último para hacer los gráficos con ese valor
dist_arr=np.delete(dist_arr,np.where(dist_arr==PSF/pix_size))
dist_arr=np.append(dist_arr,PSF/pix_size)

#performs the kinematic analysis of the MaNGA galaxy given its plate-ifu combination.
def kin_analysis(plateifu_ex):

	maps = Maps(plateifu=plateifu_ex)
	image = Image(plateifu=plateifu_ex) #para ver la imagen de sloan
	
	##############################################################
	######## Cargamos las variables a estudiar ###################
	##############################################################
	
	#Mean g-band weighted signal-to-noise ratio per pixel
	#we only choose pixels with snr>6
	snr=maps.spx_snr
	mask_snr=(snr.value < 6)	
	
	#stellar continuum (para el plot)
	flux_cont=maps.emline_sew_cnt_ha_6564
	mask_fluxcont=(flux_cont.value<1e-20)
	
	#halpha flux (para el plot)
	flux_ha=maps.emline_sflux_Ha_6564
	mask_fluxha=flux_ha.pixmask.get_mask('DONOTUSE',dtype=bool) 

	#oiii flux 
	flux_oiii=maps.emline_sflux_oiii_5008
	mask_fluxoiii=flux_oiii.pixmask.get_mask('DONOTUSE',dtype=bool) 

	#stellar velocity field:
	v_map=maps.stellar_vel
	mask_v=v_map.pixmask.get_mask('DONOTUSE',dtype=bool) #| mask_snr
	#mask_v=mask_v | mask_fluxcont 
	v=v_map.value
	v[mask_v]=np.nan
	v[mask_snr]=np.nan
	v[mask_fluxcont]=np.nan
	sat1,sat2=IQR(np.abs(v))
	v[np.abs(v)>sat1]=np.nan
	px_v=np.sum(~np.isnan(v))
	
	#halpha velocity field:
	vha_map=maps.emline_gvel_ha_6564
	mask_vha=vha_map.pixmask.get_mask('DONOTUSE',dtype=bool) # | mask_snr
	#mask_vha=mask_vha | mask_fluxha
	v_ha=vha_map.value
	v_ha[mask_vha]=np.nan
	v_ha[mask_snr]=np.nan
	v_ha[mask_fluxha]=np.nan
	sat1_ha,sat2_ha=IQR(np.abs(v_ha))
	v_ha[np.abs(v_ha)>sat1_ha]=np.nan
	px_vha=np.sum(~np.isnan(v_ha))
	
	#oiii velocity field:
	voiii_map=maps.emline_gvel_oiii_5008
	mask_voiii=voiii_map.pixmask.get_mask('DONOTUSE',dtype=bool) # | mask_snr 
	#mask_voiii=mask_voiii | mask_fluxoiii
	v_oiii=voiii_map.value
	v_oiii[mask_voiii]=np.nan
	v_oiii[mask_snr]=np.nan
	v_oiii[mask_fluxoiii]=np.nan
	sat1_oiii,sat2_oiii=IQR(np.abs(v_oiii))
	v_oiii[np.abs(v_oiii)>sat1_oiii]=np.nan
	px_voiii=np.sum(~np.isnan(v_oiii))
	
	#stellar continuum contours
	fcont=flux_cont.value
	fcont_ha=flux_cont.value
	fcont_oiii=flux_cont.value
	fcont[mask_v]=np.nan
	fcont_ha[mask_vha]=np.nan
	fcont_oiii[mask_voiii]=np.nan
	
	#funciones para cambiar de posición en píxeles a posición en arcsec
	#tomamos como origen el centro óptico (se calcula más tarde)
	def xpix_to_arcsec(A):
		return (A-on[1][0])*pix_size
	def ypix_to_arcsec(A):
		return (A-on[0][0])*pix_size
	
	##############################################################
	##### Empezamos a calcular los gradientes ####################
	##############################################################

	n=v.shape[0]
	nn=v.shape
	X=np.zeros(nn)
	Y=np.zeros(nn)

	for i in np.arange(n):
		X[:,i]=i
		Y[i,:]=i
		
	px_cen=np.array(nn)/2
	on=np.where(flux_cont.value == np.nanmax(flux_cont.value)) #optical nucleus
	dist_flag=np.sqrt((on[0]-px_cen[0])**2+(on[1]-px_cen[1])**2)[0]
	if dist_flag>5:
		flag='REVISAR'
		print(flag, dist_flag)
	else:
		flag='OK'

	#array de las posiciones del centro cinemático para cada valor de dist
	KC_x=[]
	KC_y=[]
	KCha_x=[]
	KCha_y=[]
	KCoiii_x=[]
	KCoiii_y=[]

	for dist_ind in range(len(dist_arr)):
		dist=dist_arr[dist_ind]
		diff=np.zeros(nn)
		diff_ha=np.zeros(nn)
		diff_oiii=np.zeros(nn)
		for j in np.arange(n):
			for k in np.arange(n):
				#distancia de cada píxel a los demás:
				cpix=[j,k]
				d=np.sqrt((X-X[j,k])**2+(Y-Y[j,k])**2)
				#seleccionamos solo los que están a una distancia dist del pixel j
				t=((d>0.1) & (d<dist))
				count0=np.count_nonzero(t)
				if count0>0:
	          			resta=np.abs(v[j,k]-v[t])
	          			diff[j,k]=np.nanmean(resta)
	          			resta_ha=np.abs(v_ha[j,k]-v_ha[t])
	          			diff_ha[j,k]=np.nanmean(resta_ha)
	          			resta_oiii=np.abs(v_oiii[j,k]-v_oiii[t])
	          			diff_oiii[j,k]=np.nanmean(resta_oiii)  			           		
	
		vmaps=np.array([v,diff,v_ha,diff_ha,v_oiii,diff_oiii])
	
		##############################################################
		############ FINDING THE KINEMATIC CENTER ####################
		##############################################################

		#KC with the stellar velocity
		
		rad=20 #20 pix -> 10 arcsec, como en GL2015 (igual hay que adaptarlo)
		diff_peak=np.where(diff == np.nanmax(diff))
		diff_peak=on
		d=np.sqrt((X-X[diff_peak])**2+(Y-Y[diff_peak])**2)
		diff_kc=diff
		diff_kc[d>rad] = np.nan
		pix_largediff=np.where(diff_kc > (np.nanmean(diff_kc)+np.nanstd(diff_kc)))
		if np.count_nonzero(pix_largediff)>0:
			pix_weights=diff[pix_largediff]
			xpix_average=round(np.average(pix_largediff[1],weights=pix_weights))
			ypix_average=round(np.average(pix_largediff[0],weights=pix_weights))
			pix_average=([xpix_average,ypix_average])
			KC_x.append(xpix_average)
			KC_y.append(ypix_average)
		else:
			print('Cannot compute the kinematic center with stellar velocities') 

		#KC with the Halpha gas

		diff_peak_ha=np.where(diff_ha == np.nanmax(diff_ha))
		diff_peak_ha=on
		d_ha=np.sqrt((X-X[diff_peak_ha])**2+(Y-Y[diff_peak_ha])**2)
		diff_ha_kc=diff_ha
		diff_ha_kc[d_ha>rad]=np.nan
		pix_largediff_ha=np.where(diff_ha_kc > (np.nanmean(diff_ha_kc)+np.nanstd(diff_ha_kc)))	
		if np.count_nonzero(pix_largediff_ha)>0:
			pix_weights_ha=diff_ha[pix_largediff_ha]
			xpix_average_ha=round(np.average(pix_largediff_ha[1],weights=pix_weights_ha))
			ypix_average_ha=round(np.average(pix_largediff_ha[0],weights=pix_weights_ha))
			pix_average_ha=([xpix_average_ha,ypix_average_ha])
			KCha_x.append(xpix_average_ha)
			KCha_y.append(ypix_average_ha)
		else:
			print('Cannot compute the kinematic center with gas velocities (Ha)')	

		#KC with the OIII gas
		
		diff_peak_oiii=np.where(diff_oiii == np.nanmax(diff_oiii))
		diff_peak_oiii=on
		d_oiii=np.sqrt((X-X[diff_peak_oiii])**2+(Y-Y[diff_peak_oiii])**2)
		diff_oiii_kc=diff_oiii
		diff_oiii_kc[d_oiii>rad]=np.nan
		pix_largediff_oiii=np.where(diff_oiii_kc > (np.nanmean(diff_oiii_kc)+np.nanstd(diff_oiii_kc)))	
		if np.count_nonzero(pix_largediff_oiii)>0:
			pix_weights_oiii=diff_oiii[pix_largediff_oiii]
			xpix_average_oiii=round(np.average(pix_largediff_oiii[1],weights=pix_weights_oiii))
			ypix_average_oiii=round(np.average(pix_largediff_oiii[0],weights=pix_weights_oiii))
			pix_average_oiii=([xpix_average_oiii,ypix_average_oiii])
			KCoiii_x.append(xpix_average_oiii)
			KCoiii_y.append(ypix_average_oiii)
		else:
			print('Cannot compute the kinematic center with gas velocities (OIII) ')	
	
		#Centro cinemático: media + desviación estándar de los valores calculados con diferentes aperturas. 

	xpix_average=np.nanmean(np.array(KC_x))
	ypix_average=np.nanmean(np.array(KC_y))
	xpix_std=np.nanstd(np.array(KC_x))
	ypix_std=np.nanstd(np.array(KC_y))

	xpix_average_ha=np.nanmean(np.array(KCha_x))
	ypix_average_ha=np.nanmean(np.array(KCha_y))
	xpix_std_ha=np.nanstd(np.array(KCha_x))
	ypix_std_ha=np.nanstd(np.array(KCha_y))

	xpix_average_oiii=np.nanmean(np.array(KCoiii_x))
	ypix_average_oiii=np.nanmean(np.array(KCoiii_y))
	xpix_std_oiii=np.nanstd(np.array(KCoiii_x))
	ypix_std_oiii=np.nanstd(np.array(KCoiii_y))
	
	KCRA_ST=xpix_to_arcsec(xpix_average)
	KCDEC_ST=ypix_to_arcsec(ypix_average)
	KCRA_ST_ERR=xpix_std*pix_size
	KCDEC_ST_ERR=ypix_std*pix_size
	
	KCRA_HA=xpix_to_arcsec(xpix_average_ha)
	KCDEC_HA=ypix_to_arcsec(ypix_average_ha)
	KCRA_HA_ERR=xpix_std_ha*pix_size
	KCDEC_HA_ERR=ypix_std_ha*pix_size
	
	KCRA_OIII=xpix_to_arcsec(xpix_average_oiii)
	KCDEC_OIII=ypix_to_arcsec(ypix_average_oiii)
	KCRA_OIII_ERR=xpix_std_oiii*pix_size
	KCDEC_OIII_ERR=ypix_std_oiii*pix_size
	
	##############################################################	
	######### Hacemos los gráficos ###############################
	##############################################################
	
	ftitle='PLATE-IFU '+ str(plateifu_ex) + ' - MaNGA ID ' + str(mangaid_ex_i) + '\n  Morph=' + morph_class_ex[0] + ' P_bar=' + str(round(morph_class_ex[1],2)) + ' P_edgeon=' + str(round(morph_class_ex[2],2)) + ' Nneigh=' + str(int(nneigh)) + ' \nFLAG=' + flag + ' AGN=' + str(morph_class_ex[4]) + str(morph_class_ex[5]) 
	vlabels=np.array([r'Stellar velocity (km/s)',r'Stellar velocity gradient(km/s)',r'H$\alpha$ velocity (km/s)',r'H$\alpha$ velocity gradient (km/s)',r'[OIII] velocity (km/s)',r'[OIII] velocity gradient (km/s)'])
	cmaps=['RdBu_r','plasma','RdBu_r','plasma','RdBu_r','plasma']
	cmap_pix = colors.ListedColormap(['None', 'black'])

	extent_arcsec=[xpix_to_arcsec(0),xpix_to_arcsec(X[-1,-1]), ypix_to_arcsec(0),ypix_to_arcsec(Y[-1,-1])]

	Z=np.zeros(nn)
	Z_ha=np.zeros(nn)
	Z_oiii=np.zeros(nn)
	Z[pix_largediff]=1
	Z_ha[pix_largediff_ha]=1
	Z_oiii[pix_largediff_oiii]=1

	#imagen de SDSS, continuum and halpha flux
	image.plot()
	imagetitle='PLATEIFU_'+ str(plateifu_ex) + '_img.png'
	plt.savefig(imagetitle)
	
	fig0,axs0=plt.subplots(1,3,figsize=(18,4))
	flux_cont.plot(fig=fig0,ax=axs0[0])
	flux_ha.plot(fig=fig0,ax=axs0[1])
	flux_oiii.plot(fig=fig0,ax=axs0[2])
	fig0title='PLATEIFU_'+ str(plateifu_ex) + '_fmap.png'
	fig0.savefig(fig0title)
		
	#velocity and gradient maps
	fig,axs = plt.subplots(3,2,figsize=(9,9),dpi=100)
	c=0
	for row in range(3):
		for col in range(2):
			ax=axs[row,col]
			vmap=ax.imshow(vmaps[c,:],cmap=cmaps[c], origin='lower', interpolation='nearest',extent=extent_arcsec,zorder=1)
			fig.colorbar(vmap,ax=ax,label=vlabels[c])
			#ax.contour(flux_cont.value,extent=extent_arcsec,colors='black',linewidths=0.5)
			ax.set_xlabel(r'$\Delta \delta$ (arcsec)')
			ax.set_ylabel(r'$\Delta \alpha$ (arcsec)')
			ax.set_xlim(extent_arcsec[0],extent_arcsec[1])
			ax.set_ylim(extent_arcsec[2],extent_arcsec[3])
			c=c+1
	axs[0,0].scatter(xpix_to_arcsec(on[1]),ypix_to_arcsec(on[0]),marker='+',color='white',zorder=4)
	axs[0,0].contour(fcont,extent=extent_arcsec,linewidths=1,zorder=2,colors='forestgreen',alpha=0.5)
	axs[0,1].contour(fcont,extent=extent_arcsec,linewidths=1,zorder=2,colors='forestgreen',alpha=0.5)	
	axs[0,1].contour(Z,extent=extent_arcsec,levels=0,colors='black',zorder=3)		
	axs[0,1].scatter(xpix_to_arcsec(on[1]),ypix_to_arcsec(on[0]),marker='+',color='white',zorder=4)
	axs[0,1].scatter(xpix_to_arcsec(xpix_average),ypix_to_arcsec(ypix_average), color='black',marker=',',zorder=6)
	axs[0,1].errorbar(xpix_to_arcsec(xpix_average),ypix_to_arcsec(ypix_average), xerr=xpix_std*pix_size,yerr=ypix_std*pix_size,fmt='o',color='black',zorder=5)
	axs[1,0].scatter(xpix_to_arcsec(on[1]),ypix_to_arcsec(on[0]),marker='+',color='white',zorder=4)
	axs[1,0].contour(fcont_ha,extent=extent_arcsec,linewidths=1,zorder=2,colors='forestgreen',alpha=0.5)
	axs[1,1].contour(fcont_ha,extent=extent_arcsec,linewidths=1,zorder=2,colors='forestgreen',alpha=0.5)
	axs[1,1].contour(Z_ha,extent=extent_arcsec,levels=0,colors='black',zorder=3)
	axs[1,1].scatter(xpix_to_arcsec(on[1]),ypix_to_arcsec(on[0]),marker='+',color='white',zorder=4)		
	axs[1,1].scatter(xpix_to_arcsec(xpix_average_ha),ypix_to_arcsec(ypix_average_ha), color='black',marker=',',zorder=6)
	axs[1,1].errorbar(xpix_to_arcsec(xpix_average_ha),ypix_to_arcsec(ypix_average_ha), xerr=xpix_std_ha*pix_size,yerr=ypix_std_ha*pix_size,fmt='o',color='black',zorder=5)
	axs[2,0].scatter(xpix_to_arcsec(on[1]),ypix_to_arcsec(on[0]),marker='+',color='white',zorder=4)
	axs[2,0].contour(fcont_oiii,extent=extent_arcsec,linewidths=1,zorder=2,colors='forestgreen',alpha=0.5)
	axs[2,1].contour(fcont_oiii,extent=extent_arcsec,linewidths=1,zorder=2,colors='forestgreen',alpha=0.5)
	axs[2,1].contour(Z_oiii,extent=extent_arcsec,levels=0,colors='black',zorder=3)
	axs[2,1].scatter(xpix_to_arcsec(on[1]),ypix_to_arcsec(on[0]),marker='+',color='white',zorder=4)		
	axs[2,1].scatter(xpix_to_arcsec(xpix_average_oiii),ypix_to_arcsec(ypix_average_oiii), color='black',marker=',',zorder=6)
	axs[2,1].errorbar(xpix_to_arcsec(xpix_average_oiii),ypix_to_arcsec(ypix_average_oiii), xerr=xpix_std_oiii*pix_size,yerr=ypix_std_oiii*pix_size,fmt='o',color='black',zorder=5)
	fig.suptitle(ftitle)
	figtitle='PLATEIFU_'+ str(plateifu_ex) + '_vmap.png'
	fig.savefig(figtitle)
	
	#plt.show()
	
	return(KCRA_ST, KCDEC_ST, KCRA_ST_ERR, KCDEC_ST_ERR, KCRA_HA, KCDEC_HA, KCRA_HA_ERR, KCDEC_HA_ERR, KCRA_OIII, KCDEC_OIII, KCRA_OIII_ERR, KCDEC_OIII_ERR,px_v,px_vha,px_voiii,flag)


##############################################################	
########### Definimos las galaxias a estudiar ################
##############################################################

'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
#primero elegimos qué tabla usamos: la de todo MaNGA o la de AGNs: 
t=Table.read('tabla_manga_full_updated.ecsv',format='ascii.ecsv')
#t=Table.read('tabla_mangaagn_updated.ecsv',format='ascii.ecsv')
'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

#definimos los offsets del KC al ON, y los offsets entre KC de las componentes
d_stellar=np.sqrt((t['KCRA_ST'])**2+(t['KCDEC_ST'])**2)
d_ha=np.sqrt((t['KCRA_HA'])**2+(t['KCDEC_HA'])**2)
d_oiii=np.sqrt((t['KCRA_OIII'])**2+(t['KCDEC_OIII'])**2)
d_st_ha=np.sqrt((t['KCRA_ST'] - t['KCRA_HA'])**2+(t['KCDEC_ST']- t['KCDEC_HA'])**2)
d_st_oiii=np.sqrt((t['KCRA_ST'] - t['KCRA_OIII'])**2+(t['KCDEC_ST']- t['KCDEC_OIII'])**2)
d_ha_oiii=np.sqrt((t['KCRA_HA'] - t['KCRA_OIII'])**2+(t['KCDEC_HA']- t['KCDEC_OIII'])**2)

##############################################################	
########### Definimos algunas máscaras útiles ################
##############################################################

#GALAXIAS GRANDES
mask_px=((t['PX_V']>160) & (t['PX_VHA']>160) & (t['PX_VOIII']>160))
#GALAXIAS CON FLAG REVISAR (ON LEJOS DEL CENTRO DE LA IMAGEN)
mask_revisar=(t['FLAG']=='REVISAR')
#GALAXIAS QUE USAREMOS PARA LOS ANÁLISIS (NO SE DESCARTAN) 
mask_noflag=((t['FLAG']=='OK') & mask_px)
#GALAXIAS CON DISTINTAS PROBABILIDADES DE BARRA
mask_barA=(t['P_BAR']>2/3)
mask_barAB=((t['P_BAR']>1/3) & (t['P_BAR']<2/3))
mask_barB=(t['P_BAR']<1/3)
#GALAXIAS CONFIRMADAS COMO AGN
mask_agn=(t['AGN']>0)
mask_noagn=~mask_agn
#LISTA DE GALAXIAS AISLADAS (SIG) / NO AISLADAS:
mask_nosig=(t['NNEIGH']>0)
mask_sig=~mask_nosig
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
#MÁSCARAS DE SI EXISTE O NO ALGÚN OFFSET 
mask_nooffset=((d_stellar<fiber_size) & (d_ha<fiber_size) & (d_oiii<fiber_size) & (d_st_ha<fiber_size) & (d_st_oiii<fiber_size) & (d_ha_oiii<fiber_size) )
mask_someoffset=((d_stellar>2*fiber_size) | (d_ha>2*fiber_size) | (d_oiii>2*fiber_size) | (d_st_ha>2*fiber_size) | (d_st_oiii>2*fiber_size) | (d_ha_oiii>2*fiber_size) ) 

print('\nTAMAÑO DE LAS MUESTRAS') #añadir los parámetros relevantes
print('Número de galaxias en la muestra total:', len(t['PLATEIFU']))
print('Número de galaxias a revisar:', np.sum(mask_revisar))
print('Número de galaxias pequeñas:', np.sum(~mask_px))
print('Número de galaxias no descartadas para el análisis:', np.sum(mask_noflag))
print('Número de galaxias que son AGN',np.sum(mask_agn))
print('Número de galaxias sin ningún offset',np.sum(mask_noflag & mask_nooffset))
print('Número de galaxias con algún offset',np.sum(mask_noflag & mask_someoffset))

'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

##############################################################	
##### Definimos las galaxias que queremos representar ########
##############################################################

#Ejemplo: galaxias grandes, sin flag, aisladas, confirmadas como AGN que también presentan algún offset:
plateifu_list=t['PLATEIFU'][mask_noflag & mask_sig & mask_agn & mask_someoffset]
#Ejemplo: galaxias barradas marcadas con el flag "revisar" que indica que el ON está lejos del centro de la imagen:
plateifu_list=t['PLATEIFU'][mask_revisar]

'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

#(podemos definir la muestra que queramos según lo que nos interese)

################################################################	
## Si queremos representar todas las galaxias de la muestra: ###
################################################################

#primero elegimos qué tabla usamos: la de todo MaNGA o la de AGNs: 
t=Table.read('tabla_manga_full_updated.ecsv',format='ascii.ecsv') #todo MaNGA
#t=Table.read('tabla_mangaagn_updated.ecsv',format='ascii.ecsv') #catálogo de AGNs

plateifus=t['PLATEIFU']

#Hacemos el código de 200 en 200 y vigilamos que no dé error (si da error se vuelve a correr)
print('GALAXIAS 0000-0200')
plateifu_list=plateifus[0:200] #lista de los nombres de las galaxias que vamos a representar

for i in range(len(plateifu_list)):
	start=time.time()
	plateifu_i=plateifu_list[i]
	#plateifu_i=np.random.choice(plateifu_list)
	print('PLATEIFU' , plateifu_i)
	mangaid_ex_i=basic_parameters(plateifu_i)[1]
	print('mangaid', mangaid_ex_i)
	nneigh=isolation(plateifu_i)[0]
	qnn=isolation(plateifu_i)[1]
	qlss=isolation(plateifu_i)[2]
	print('nnegih', nneigh, 'qnn', qnn, 'qlss', qlss)
	print('P(merge):', morph_class(plateifu_i)[3])
	basic_parameters_ex=basic_parameters(plateifu_i)
	mangaid_ex_i=basic_parameters_ex[1]
	morph_class_ex=morph_class(plateifu_i)
	kin_an_ex=kin_analysis(plateifu_i)
	
	print('Galaxy ', i+1, '/', len(plateifu_list), 'finished \n ') 
	
	
