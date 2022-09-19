# Inferring dynamical distinct components in nearby galaxies using velocity gradients
## Alba Casasbuenas Corral (Beca de Verano IAC, 2022)
### Supervisors: Begoña García-Lorenzo & Natacha Zanon Dametto

The stellar and ionized gas kinematics trace the dynamical structure of galaxies and constrain the processes driving their evolution. [García-Lorenzo et al. (2015)](https://www.aanda.org/articles/aa/full_html/2015/01/aa23485-14/aa23485-14.html) developed a simple approach for estimating the kinematic center location of galaxies using the velocity maps of the ionized gas of a small sample of galaxies using the CALIFA galaxy survey. In this work, we have adapted this method to the MaNGA galaxy survey, using a much larger sample of galaxies in the Local Universe. Moreover, we have completed the analysis by using both stellar and ionized gas velocity fields. We have found that there are several factors that cause departures from pure rotational motions in a galaxy, that manifest as an offset between their kinematic center and optical nucleus. Two of these factors are galaxy mergers and interactions and stellar mass. 

This is the repository of the python programs used in this work for the kinematic analysis of the MaNGA sample, as well as the tables generated. To use them, it's necessary to have the Marvin tools package installed (available at the Denso machine at IAC). 

The python programs included in this repository also use the following catalogs: the MaNGA Pipe3D value-added catalog (Sánchez et al., 2016), the MaNGA Morphology Deep Learning DR17 catalog (Domínguez Sánchez et al., 2022), the GEMA-VAC: Galaxy Environment for MaNGA Value-Added catalog (Argudo-Fernández et al., 2015) and the MaNGA AGN catalog (Comerford et al., 2020). The links to download these value-added catalogs as .fits tables are: 

- https://data.sdss.org/sas/dr17/manga/spectro/pipe3d/v3_1_1/3.1.1/
- https://data.sdss.org/sas/dr17/env/MANGA_MORPHOLOGY/deep_learning/1.1.1/
- https://data.sdss.org/sas/dr17/env/MANGA_GEMA/2.0.2/
- https://data.sdss.org/sas/dr17/env/MANGA_AGN/v1_0_1/
