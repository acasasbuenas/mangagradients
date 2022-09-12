"""
Created on Mon 8 Aug 2022
@author: Alba Casasbuenas
:)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits, ascii
from astropy.table import Table, vstack
import time
import glob

listatablas=glob.glob('/path/*') #path de la carpeta donde están todas las tablas separadas
t=[]

for i in range(len(listatablas)):
	t.append(Table.read(listatablas[i],format='ascii.ecsv'))

tablafinal=t[0]
for i in range(len(t)-1):
	tablafinal=vstack([tablafinal,t[i+1]])

nombretabla='tabla_xxxxx.ecsv'
tablafinal.write(nombretabla, format='ascii.ecsv', overwrite=True)  

