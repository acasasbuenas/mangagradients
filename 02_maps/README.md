# Create the velocity and gradient maps 

The program **maps.py** does the relevant kinematic analysis for a list of galaxies, and shows three plots for each galaxy: one with the SDSS optical image of the galaxy, another with the continuum, halpha and [OIII] flux, and a final one with the velocity and gradient maps for each of the components. 

The sections we need to change are framed with exclamation symbols, and they go from line 465-470, and from line 430 to 550 approximately.

- **t** is the table we read from which to get the relevant galactic parameters, important if we want to obtain the maps of a subsample of galaxies with certain properties. We can choose the full MaNGA sample or the AGN subsample (first six thousand galaxies observed through mid-2018)

- We define afterwards some useful masks to choose the galaxy subsamples we might be interested in. We can add any masks we'd consider useful.

- **plateifu_list** is the list of the plate-ifu combination of the subsample of galaxies we want to represent. We can define it through the combination of the relevant masks. Afterwards, inside the for loop, we can choose through the **plateifu_i** definition if we want to represent all the galaxies in order `plateifu_i=plateifu_list[i]` or at random `plateifu_i=np.random.choice(plateifu_list)`

- NOTE: this programs shows and saves all the figures represented. If we don't want this, we can comment the `plt.show()` or the `savefig` lines in the kin_analysis program 
