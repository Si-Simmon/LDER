# LDER
A classification for RSVP-BCI 
This repository contains code to improve the classification performance for RSVP-BCI by 
enhancing ERP features. These codes are for the paper “LDER: A classification framework 
based on ERP enhancement in RSVP task”. Users employing the code in a publication 
should cite this paper. 

All codes were written in MATLAB 2018b 
“ERP_enhance_2” contains the code used to enhance ERP components 
“ERP_classification_Main” contains the code used to classify targets and non-targets by 
using STHCP method. This code implements the STHCP algorithm outlined in the following 
article. Users employing the code in a publication should cite this paper. 

Cui Y, Xie S, Xie X, Zhang X and Liu X 2022 Dynamic probability integration for 
electroencephalography-based rapid serial visual presentation performance enhancement: 
Application in nighttime vehicle detection Front. Comput. Neurosc. 16 
(https://doi.org/10.3389/fncom.2022.1006361) 

The other codes are functions contained in these two main methods. 


Notes: 
1. Our codes need to run based on the EEGLAB, which can be found at 
https://sccn.ucsd.edu/eeglab/download.php. Data need to be imported through eeglab and 
the format is N*T, before running “ERP_enhance_2”. 
2. After running “ERP_enhance_2”, the data should be segmented into the format of N*T*S, 
and loaded in “ERP_classification_Main”. 
