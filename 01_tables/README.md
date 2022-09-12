# Generate tables with the relevant parameters

The program **table.py** does the relevant kinematic analysis for all the galaxies in a sample.

The part we need to change is a section framed with exclamation symbols, and it goes from line 455 to 465.

- **plateifus**: list of the plateifu combinations of the galaxies to study. 
  
  For example, in order to choose the full MaNGA sample, we would just get all the plateifu combinations of all the galaxies in the drpall MaNGA galaxies: 
 ```
data = get_drpall_table()
plateifus_data=data['plateifu']
plateifus=plateifus_data
 ```
  To choose the first six thousand galaxies observed with MaNGA in order to study the AGN subsample, we have generated a table with all the observation dates of the galaxies (tabla_fechas.ecsv). 
 ```  
t=Table.read('tabla_fechas.ecsv',format='ascii.ecsv')
t.sort('DATEOBS')
plateifus=t['PLATEIFU'][0:6522]
 ```
 
 - **plateifu_list**: list of galaxies to run at one time. Because this code is slow and Marvin sometimes gives errors when trying to retrieve the maps, we chose to run the code 200 galaxies at once, and then put all the tables together in the same folder, and put them together with the program **stacktables.py** 

This program has the code to make the graphs for each galaxy commented, because it takes too long and we're using a big sample of galaxies. In order to see the stellar maps, check the folder **02_maps**. 
