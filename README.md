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

2. **Update catch and portside sampling data**: Download the most recent data from these sites and merge with older data using [`Data_processing/Code/Harvest reconstruction.R`]

   a. *Portside sampling for biological data* is available here: https://oceanak.dfg.alaska.local/analytics/saw.dll?Answers&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2Fport%20sampling%20bio%20data
   
   b. *Yelloweye harvest data* is available here:

      ---*Groundfish fish tickets:* https://oceanak.dfg.alaska.local/analytics/saw.dll?Answers&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2FYelloweye%20harvest%20SEO%20in%20groundfish%20fish%20ticket%20data%202015-current

      ---*Salmon fish tickets:* https://oceanak.dfg.alaska.local/analytics/saw.dll?Answers&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2FYelloweye%20harvest%20SEO%20in%20salmon%20fish%20ticket%20data%202015-current
   
   c. *Halibut harvest data* is available here.  (This is necessary for estimating yelloweye bycatch.) https://oceanak.dfg.alaska.local/analytics/saw.dll?Answers&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2FHalibut%20harvest%20SEO%20in%20fish%20ticket%20data%202007-2022 

3. **Estimate yelloweye density for new ROV surveys** If new areas have been sampled since the last assessment, estimate density and then biomass.
   
   a. [`YE_density_estimation/Code/ROV_Processing.R`]: Script for processing raw **ROV** data, prepping it for GIS work, and then complete the processing.
   
   b. [`YE_density_estimation/Code/ROV_base_distance_modelling.R`]: Examine covariates in the data and run the basic distance models.

   c. [`YE_density_estimation/Code/ROV_dist_model_avg.R`]: Code for averaging multiple distance models for estimate of yelloweye density.

   d. [`Data_processing/Code/Biomass_Subdistrict_Calcs.R`]: Code for converting yelloweye density to biomass estimates using the latest available weight data from portside samples.

4. **Calculate yelloweye cpue in the IPHC longline survey (FISS)**: This will be a moving target as the IPHC website is undergoing continues revisions. In addition, these methods are currently underoing some revisions with the help of IPHC staff to better align with best practices.  Specifically, cpue estimates need to be adjusted for hook saturation.  To update the data, go to

   a. https://www.iphc.int/data/fiss-data-query : This page will allow you to download the NPUE data for non-halibut and save it to [`Data_processing/Data/IPHC_raw/`] folder. Follow the filters to select IPHC regulatory areas 3A and 2C and yelloweye rockfish in the "Select Non-Pacific halibut species" tab. Just above the "tableau" bar, there is the download button.  Press it, select "Crosstab" and then "CSV" and *Non-Pacific halibut data* and then download *Set and Pacific halibut data*.  Make sure you label this file with the year when saved to the proper folder.  In 2023 I had problem downloading and importing the csv files and had to download the excel file and then save as a csv file.  

   b. https://www.iphc.int/datatest/fiss-pacific-halibut-data : As of 2023 all the data you need is described above, but halibut data is also available here so I'm including it in case things change again.  In this arrangement of the data base you have to download the set data and halibut data separately and then merge them in the code.  

   c. [`Data_processing/Code/IPHC_Survey_CPUE_index.R`]: Code for calculating IPHC CPUE index of abundance for use in the REMA model and the SS-SPM.  Going into the 2023 CIE review I am debating what stations to include in the CPUE calculations.  Right now I have only been using stations that encounter yelloweye 40% of the time, but will discuss at review...
   
5.  **Calculate expected yelloweye bycatch in the halibut fishery** [`Data_processing/Code/IPHC_Survey_Expected_Bycatch.R`]: Estimate the expected bycatch from the halibut fishery using the IPHC survey data and historical halibut harvests.  CAS estimates can be retrieved from INSERT HERE  

6. **REMA model** Run the random effects model using the biomass estimates and IPHC cpue estimates.  This is the model used in the 2022 SAFE report.
