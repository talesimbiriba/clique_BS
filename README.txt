************************************
************************************

# clique BS
Matlab code for Band Selection in RKHS of hyperspectral images. 

This Package has all programs and third part software to run the simulations presented in [1]. 
All third part software used here has free distribution and noncomercial use licenses. 
The routine setup.m sets up the matlab path for our and third party software and other 
configuration needed.

All third party software is listed below. 

The use and redistribution of this code is allowed for noncomercial purposes as long as the 
copyrights presented here are replicated. 

Author: Tales Imbiriba.

Date: March 2016. 

Ref.:

[1]  T. Imbiriba, J.C.M. Bermudez, C. Richard, "Band selection for nonlinear unmixing of hyperspectral images as maximal click problem", submitted to the IEEE Transaction on Image Processing Magazine.


************************************
************************************

SYSTEM REQUIREMENTS

This matlab code use third party software that were compiled in a linux plataform. Thus, for the proper operation we strongly recomend a linux environment.  


************************************
************************************


Main Files:

setup.m 
	This program adds all the needed code (ours and third part) to matlabs path. Furthermore,
	it also creates goToCliqueDir.m script needed by the clique BS function 'clique_coherence_bandselection()' in 'clique_BS/Source/BandSelection/'. This script is needed to write and read files in the directory where the MaxCQL executable is. 


compareBSMethods.m
	Performs the simulations with synthetic data used in [1]. It compares BS algorithms using clique, greedy, KKM and SK-Hype. 


Cuprite_BS_And_SKHype.m
	Performs the simulations presented in [1] using the Cuprite scene aquired by the AVIRIS instrument.


PaviaU_BS_And_SKHype.m
	Performs the simulations presented in [1] using the Pavia scene aquired by the ROSIS specctrometer.


************************************
************************************


Third part software:


	SK-Hype (Matlab code avaiable at http://cedric-richard.fr/pub.html)

	MaxCQL	(linux executable available at http://home.mis.u-picardie.fr/~cli/EnglishPage.html)		 
	


