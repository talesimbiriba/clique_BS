************************************
************************************

# clique BS
Matlab code for Band Selection in RKHS of hyperspectral images. 

This Package has all programs and third part software to run the simulations presented in [1] and [2]. 
All third part software used here has free distribution and noncomercial use licenses. 
The routine setup.m sets up the matlab path for our and third party software and other 
configuration needed.

All third party software is listed below. 

The use and redistribution of this code is allowed for noncomercial purposes as long as the 
copyrights presented here are replicated. 

Author: Tales Imbiriba.

Date: March 2016. 

Ref.:

[1]  T. Imbiriba, J.C.M. Bermudez, C. Richard, "Band selection for nonlinear unmixing of hyperspectral images as maximal click problem", submitted to the IEEE Transaction on Image Processing Magazine. (http://ieeexplore.ieee.org/document/7867872/)

[2] T. Imbiriba, J.C.M. Bermudez, C. Richard, "Technical Report: Band selection for nonlinear unmixing of hyperspectral images as a maximal clique problem", [Online] Available: https://arxiv.org/abs/1603.00437


************************************
************************************

SYSTEM REQUIREMENTS

This matlab code use third party software that were compiled in a linux plataform. Thus, for the proper operation we strongly recomend a linux-like environment.  


************************************
************************************


Main Files:

setup.m 
	This program adds all the needed code (ours and third part) to matlabs path. Furthermore,
	it also creates goToCliqueDir.m script needed by the clique BS function 'clique_coherence_bandselection()' in 'clique_BS/Source/BandSelection/'. This script is needed to write and read files in the directory where the MaxCQL executable is. 


compareBSMethods_And_ProgBS.m
	Performs the simulations with synthetic data used in [1] and [2]. It compares BS algorithms using clique, greedy, KKM, SK-Hype, and Progressive Band selection. 


Cuprite_BS_And_SKHype_recError.m
	Performs the unmix using the proposed methods and SK-Hype and display the Reconstruction error (not abndance RMSE) of each method. This simulation corresponds to Table III in [1].


realLab_JFM89_BS_Unmix.m 
	Performs simulations using the RealLab data set provided by prof. John Mustard. The program outputs the RMSE for the abundances (Tables presented in tables IV and V in [1] and VIII to XI in [2]). 


Cuprite_BS_And_SKHype.m
	Performs the simulations presented in [2] using the Cuprite scene aquired by the AVIRIS instrument.


PaviaU_BS_And_SKHype.m
	Performs the simulations presented in [2] using the Pavia scene aquired by the ROSIS specctrometer.


compareBSMethods.m
	Performs the simulations with synthetic data used in [2]. It compares BS algorithms using clique, greedy, KKM, and SK-Hype (does NOT include PBS). 


CompareBSMethods_LargeM.m
	Performs simulations presented in [2] Table VII. Reconstruction error for the Cuprite image considering large values for M (preselected number of bands).


************************************
************************************


Third part software (also included in this package):


	SK-Hype (Matlab code avaiable at http://cedric-richard.fr/pub.html)

	MaxCQL	(linux executable available at http://home.mis.u-picardie.fr/~cli/EnglishPage.html)		 
	


