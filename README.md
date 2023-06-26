# seak_seo_dsr
# Southeast Alaska Outside Demersal Shelf Rockfish Assessment

This repository will hold data, files and r script for the SEO DSR assessment that ADF&G performs as part of the NPFMC stock assessments.  Resultant assessments associated with this code can be found here:

2022 assessment: https://www.fisheries.noaa.gov/resource/data/2022-assessment-demersal-shelf-rockfish-stock-complex-southeast-outside-subdistrict

This repository will contain the following sub directories that contain different analysis that are part of the assessment.  They are,

1. **Data processing:** This folder contains script and raw data files for putting together data for use in the analysis contained in the folders below.  This includes catch data and biological data.
2. **Density estimation:**  Script and data for estimating yelloweye rockfish density using distance sampling methods and transect data from ROV surveys.
3. **Unreported discard estimation:**  Script and data for estimating unreported discards of yelloweye rockfish in the halibut fishery
4. **REMA analysis:** This folder contains the scripts for running the random effects model that is the current model used in the federal assessment.
5. **Production models:** This folder contains scripts for developing and running the bayesian state-space surplus production model (SS-SPM).

## Code path to follow for assessment:

1. [`r_helper/`] Folder with pre-canned functions referenced throughout the assessment

2. **Update catch and portside sampleing data**: blah blah

3. **Estimate yelloweye density for new ROV surveys** If new areas have been sampled since the last assessment, estimate density and then biomass.

  a. [`YE_density_estimation/Code/ROV_Processing.R`]: Script for processing raw ROV data, prepping it for GIS work, and then complete the processing
  
  b. [`YE_density_estimation/Code/ROV_base_distance_modelling.R`]: Examine covariates in the data and run the basic distance models
  
  c. [`YE_density_estimation/Code/ROV_dist_model_avg.R`]: Code for averaging multiple distance models for estimate of yelloweye density.
  
  d. [`Data_processing/Code/Biomass_Subdistrict_Calcs.R`]: Code for converting yelloweye density to biomass estimates using the latest available weight data from portside samples. 

4. **Calculate yelloweye cpue in the IPHC longline survey (FISS)**: This will be a moving target as the IPHC website is undergoing continues revisions. In addition, these methods are currently underoing some revisions with the help of IPHC staff to better align with best practices.  Specifically, cpue estimates need to be adjusted for hook saturation.  To update the data, go to

  a. https://www.iphc.int/data/fiss-data-query : This page will allow you to download the NPUE data for non-halibut and save it to [`Data_processing/Data/IPHC_raw/`] folder. Follow the filters to select IPHC regulatory areas 3A and 2C and yelloweye rockfish in the "Select Non-Pacific halibut species" tab. Just above the "tableau" bar, there is the download button.  Press it, select "Crosstab" and then "CSV" and Non-Pacific halibut data.  Make sure you label this file with the year when saved to the proper folder.
  b. https://www.iphc.int/datatest/fiss-pacific-halibut-data : This page will allow you to download the WPUE data for halibut and save it to the same folder. Follow the same instructions (except you don't have to select species) to download and save in the same location. 
  c. [`IPHC_Survey_CPUE_index`]: Code for calculating IPHC CPUE index of abundance for use in the REMA model and the SS-SPM.  
  
5.   
