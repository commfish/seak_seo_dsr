# seak_seo_dsr

<img align="left" src="https://github.com/commfish/cr_brf/blob/main/figures/SRIlogo.jfif" width="100">

# Southeast Alaska Outside Demersal Shelf Rockfish Assessment

<br clear="left"/>

The Southeast Outside subdistrict (SEO) of the Gulf of Alaska (GOA) demersal shelf rockfish (DSR) complex (yelloweye, quillback, copper, rosethorn, China, canary, and tiger rockfish) is assessed on a biennial cycle, with a full stock assessment conducted every second year. The SEO area is comprised of four management sections: East Yakutat (EYKT), Northern Southeast Outside (NSEO), Central Southeast Outside (CSEO), and Southern Southeast Outside (SSEO) (Figure 1). Yelloweye rockfish (*Sebastes ruberrimus*) have been managed as a Tier 4 stock; however, given that the current assessment models do not produce reliable estimates of $F_{35}$ or $F_{40}$, the stock has been determined to be more appropriately in Tier 5, where $F_{OFL}$ = *M* and $maxF_{ABC}$ = *0.75* * *M*. The other species in the complex are managed as Tier 6.   Prior to 2022, the stock assessment was based on biomass estimates of yelloweye rockfish derived from Alaska Department of Fish and Game (ADF&G) submersible and ROV surveys.  Harvest recommendations were established by applying an assumed mortality rate to the lower 90% confidence interval of the estimate as a hedge against uncertainty.  The recommended acceptable biological catch (ABC) and overfishing level (OFL) for the current assessment are based on a random effects models (REMA) applied to the ROV biomass estimates and CPUE estimates of yelloweye rockfish in the International Pacific Halibut Commission (IPHC) longline survey.  A new Bayesian state-space surplus production model (SS-SPM) is in development and schedule for CIE review in September 2023.

<img src="https://github.com/commfish/seak_seo_dsr/blob/main/readme/seomap.png" width="600" />

*Figure 1. The Southeast Outside Subdistrict of the Gulf of Alaska (SEO) with the Alaska Department 
of Fish and Game groundfish management areas used for managing the demersal shelf rockfish fishery: 
East Yakutat (EYKT), Northern Southeast Outside (NSEO), Central Southeast Outside (CSEO), and 
Southern Southeast Outside (SSEO) Sections.*

The assessments associated with this code can be found here:

* [2022 SAFE report](https://www.fisheries.noaa.gov/resource/data/2022-assessment-demersal-shelf-rockfish-stock-complex-southeast-outside-subdistrict)
* [2024 SAFE report](https://www.npfmc.org/wp-content/PDFdocuments/SAFE/2024/GOAdsr.pdf)

### Repo directory: 

This repository will hold data, files and r script for the SEO DSR assessment that ADF&G performs as part of the NPFMC stock assessments.  This repository contains the following sub directories that contain different analysis that are part of the assessment.  They are:

1. **Data processing:** This folder contains script and raw data files for putting together data for use in the analysis contained in the folders below.  This includes catch data and biological data.
2. **Density estimation:**  Script and data for estimating yelloweye rockfish density using distance sampling methods and transect data from ROV surveys. The ROV survey program was suspended in 2023 due to lack of funding.
4. **REMA analysis:** This folder contains the scripts for running the random effects model that is the current model used in the federal assessment.
5. **Production models_JAGS:** This folder contains scripts for developing and running the bayesian state-space surplus production model (SS-SPM) in JAGS.  This model was reviewed by a CIE panel in September 2023 and was transitioned to STAN for further development.
6. **Production models_STAN:** This folder contains scripts for developing and running the bayesian state-space surplus production model (SS-SPM) in STAN.  This modeling is the working directory and follows the recomendations of the CIE panel.
7. **Output:** This folder will contain output files from various analysis.
8. **Figures:** This folder will contain figures associated with the assessment.
9. **r_helper:** This folder contains functions built specific to the assessment. 

## Code path and data sources for the assessment:

1. [`r_helper/`]: Folder with pre-canned functions referenced throughout the assessment

2. **Update catch and portside sampling data**: Download the most recent data from these sites and merge with older data using [`Data_processing/Code/Harvest reconstruction.R`]

    1. [*OceanAK Portside sampling for biological data*](https://oceanak.adfg.alaska.gov/analytics/saw.dll?Answers&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2Fport%20sampling%20bio%20data)  
   
    2. *Yelloweye harvest data* is available here:

        1. [*OceanAK Groundfish fish tickets*](https://oceanak.dfg.alaska.local/analytics/saw.dll?Answers&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2FYelloweye%20harvest%20SEO%20in%20groundfish%20fish%20ticket%20data%202015-current)

        2. [*OceanAK Salmon fish tickets*](https://oceanak.dfg.alaska.local/analytics/saw.dll?Answers&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2FYelloweye%20harvest%20SEO%20in%20salmon%20fish%20ticket%20data%202015-current)
      
        3. [*CAS bycatch in halibut fishery (NOAA observer estimates)*](https://akfinbi.psmfc.org/analytics/saw.dll?Dashboard): Make sure you are signed into the dahsboard [here](https://akfinbi.psmfc.org/bi-security-login/login.jsp;jsessionid=8WlrKptywRXX_FwL4FJEorojht32nDFuLjRKqWNOyIM1xZDc9y4M!2006146099?msi=false&redirect=L2FuYWx5dGljcy9zYXcuZGxsP0Rhc2hib2FyZCZoYXNoPXRTMUVwWDR1TUNmN3dTbjV6UGlXUWh0Vk9xcHZDZmgyR0FROU5UdDRZaHBRTkpWUWZDTE1ZR2FkTXV0U3JhdVo=) to make sure you can use the link. There are two prepared reports currently titled *Joy - Yelloweye Bycatch* and *Joy - Yelloweye Bycatch with ADFG stat area* from AKFIN database under the tab *Analyst Specific Reports*.  Save these files in [`Data_processing/Data/Harvests/CONF/`] and rename the files as *CAS Yelloweye Rockfish Bycatch in Halibut Fishery 2013-present.csv* and *CAS Yelloweye Bycatch with Stat 2013-present.csv*, respectively.  This report was completed in September 2023 by Jean Lee (jean.lee@noaa.gov) and is currently being reviewed.  This data contains confidential information and will not be stored on the publicly available repository.  
   
    3. [*Halibut harvest data*](https://oceanak.dfg.alaska.local/analytics/saw.dll?Answers&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2FHalibut%20harvest%20SEO%20in%20fish%20ticket%20data%202007-2022). This is necessary for estimating yelloweye bycatch. 

3. **Compile harvest data**: [`Data_processing/Code/YE Harvest reconstruction.R`] Code for compiling fish ticket data, sport fish data and subsistence data to create a catch history for use in assessment.
   
5. **Estimate catch in the foreign trawl fleet**: A Pacific ocean perch trawl fishery in the Gulf of Alaska developed in the early 1960’s with large effort by the U.S.S.R and Japanese fleets. At the height of the fishery in 1965, the catches of all rockfish, including Pacific ocean perch, exceeded 370,000 metric tons (mt). Generally, catches declined following this peak until foreign fishing was abolished from the Gulf of Alaska in 1987. During the early period of this foreign fishery (1961-1974), foreign catches of rockfish were often reported in crude management groups, including “Pacific ocean perch” or “other rockfish”, with no differentiating between species. With implementation of a fishery observer program in 1975 and 1977 onward, species composition, including Pacific ocean perch and Yelloweye Rockfish, of foreign catches became available. Donnie Arther with ADF&G Sport Fish Division compiled records of this harvest to understand and estimate yelloweye removals in this fishery.  The data is extremely limited but provides "ball park" estimates for use in a bayesian stock assessment.  The raw data is in [`Data_processing/Data/YERF_Foreign_reconstruction.csv`] and the code for estimating catch is in [`Data_processing/Code/SEAK_YERF_Foreign_fleet_reconstruction.R`].  The estimates for use in the production model are saved as [`Data_processing/Data/Harvrests/Foreign_YERF_SEAK.csv`].  

6. **Biological analysis**: [`Data_processing/Code/biological.R`] Examine and make graphics of biological data (age, sex, length, weight and maturity) for yelloweye rockfish using portside samples. The data is fomratted for the use in an age structured assessment if one is developed.  Biological trends can be incorporated into the assessment regardless of the assessment model being used.  

7. **Estimate yelloweye density for new ROV surveys**: If new areas have been sampled since the last assessment, estimate density and then biomass.
   
    1. [`YE_density_estimation/Code/ROV_Processing.R`]: Script for processing raw **ROV** data, prepping it for GIS work, and then complete the processing.
   
    2. [`YE_density_estimation/Code/ROV_base_distance_modelling.R`]: Examine covariates in the data and run the basic distance models.

    3. [`YE_density_estimation/Code/ROV_dist_model_avg.R`]: Code for averaging multiple distance models for estimate of yelloweye density.

    4. [`Data_processing/Code/Biomass_Subdistrict_Calcs.R`]: Code for converting yelloweye density to biomass estimates using the latest available weight data from portside samples.

8. **Calculate yelloweye cpue in the IPHC longline survey (FISS)**: This will be a moving target as the IPHC website is undergoing continues revisions. In addition, these methods are currently underoing some revisions with the help of IPHC staff to better align with best practices.  Specifically, cpue estimates need to be adjusted for hook saturation.  To update the data, go to

    1. [FISS data query](https://www.iphc.int/data/fiss-data-query): This page will allow you to download the NPUE data for non-halibut and save it to the [`Data_processing/Data/IPHC_raw/`] folder. Follow the filters to select IPHC regulatory areas 3A and 2C and yelloweye rockfish in the "Select Non-Pacific halibut species" tab. Just above the "tableau" bar, there is the download button.  Press it, select "Crosstab" and then "CSV" and *Non-Pacific halibut data* and then download *Set and Pacific halibut data*.  Make sure you label this file with the year when saved to the proper folder.  In 2023 I had problem downloading and importing the csv files and had to download the excel file and then save as a csv file.  
    
    2. Get the hook adjustment factor data to correct CPUE estimators for hook saturation. This data is not posted regularly and may require a [data request](https://www.iphc.int/form/data-request).  If someone has already requested the most recent data it may be posted at [IPHC hook adjustments](https://www.iphc.int/data/fiss-survey-raw-survey-data/), Document IPHC-YEAR-FISS-HADJ.  The online data is a large, heavily formated excel sheet that will need some work to convert to csv file.      

    3. [FISS Halibut data](https://www.iphc.int/datatest/fiss-pacific-halibut-data): As of 2023 all the data you need is described above, but halibut data is also available here so I'm including it in case things change again.  In this arrangement of the data base you have to download the set data and halibut data separately and then merge them in the code.  

    4. [`Data_processing/Code/IPHC_Survey_CPUE_index.R`]: Code for calculating IPHC CPUE index of abundance for use in the REMA model and the SS-SPM.  Use all stations that have encountered a yelloweye at least once in the time series. 
   
9.  **Calculate expected yelloweye bycatch in the halibut fishery**: This section involves estimating historical discards of yelloweye in the halibut fishery by applying bycatch rates estimated during contemporary times to historical catches and propogating uncertainty.  

    1. [`Data_processing/Code/Halibut_harvest_reconstruction.R`]: Code for compiling halibut harvest data and apportion catches to SEO management areas.
   
    2. [`Data_processing/Code/IPHC_Survey_Expected_Bycatch.R`]: Estimate the expected bycatch from the halibut fishery using the IPHC survey data and historical halibut harvests.  CAS estimates of bycatch from the NOAA observer program are stored in [`Data_processing/Data/Harvests/CONF/CAS Yelloweye Rockfish Bycatch in Halibut Fishery 2013-present.csv`] and [`CAS Yelloweye Bycatch with Stat 2013-present.csv`].
   
10. **REMA model**: Run the random effects model using the biomass estimates and IPHC cpue estimates [`REMA/Code/REMA_run.R`].  This is the model used in the 2022 SAFE report.  **Model 22.2** is the current assessment model for the NPFMC as of 2022.  This code contains multiple rema models, methods for comparisons and figures. 

11. **Production_models_STAN**: This folder contains code for running the state-space surplus production model (SS-SPM) in STAN, plotting results, saving output and comparing different models. Details of this model are included in the appendices of the 2022 SAFE report linked above but have been modified according to the CIE review in 2023. As of this update (January 2024) several parameterizations of the model have been built in STAN and coded to accommodate a number of different scenarios and allows for changing of key priors without having to modify the model files (.stan).  

    1. [`Production_models_STAN/Code/stan_helper.R`]: Functions for working with the SPMs. 

    2. [`Production_models_STAN/Code/spm_modeling.R`] Streamlined code for simulating data, running models and examining results.  There is preliminary code to set up a full simulation study and output results. Functions used in this script are contained in the stan_helper.R code.
    
    3. [`Production_models_STAN/Code/divergent_exam.R`] This code provides a loop to examine different models, priors and parameterizations to examine solutions for dealing with divergent transitions.
    
    4. [`Production_models_STAN/Code/simple_schaefer_eg.R`] This code provides a simple schaefer model that was used to examine divergent transitions. I was able to reduce divergences by marginalizing catchability and changing the prior on sigma but have not had success in getting rid of divergent transitions in the more complex models.
   
    5. [`Production_models_STAN/Models/`]:  This folder contains the models (.stan files) in a number of different parameterizations and with different levels of flexibility and functionality.  Please see the code in [`Production_models_STAN/Code/spm_modeling.R`] for details on the particulars of the model constructions. 

12. **Production_models_JAGS**: This folder contains code for running the state-space surplus production model (SS-SPM) in JAGS, plotting results, saving output and comparing different models. Details of this model are included in the appendices of the 2022 SAFE report linked above. There are 3 stages to running the model, all of which involve a Pella-Tomlinson surplus production model currently configured to maximize the production curve at $B_{40}$ (biomass at 40% of virgin biomass).  The first two stages involves running the models under minimal assumptions (i.e., priors) and using those results to construct informative priors for *K* (carrying capacity for the SEO as a whole) and $\phi_{80}$ (the proportion of K present in 1980).  Part of the reason for this is the difficulty in the long term data series (back to 1888) and the spatial resolution available in the data.  Since 1980 data is available for the 4 management areas in the SEO.  Prior to that data is only available at the scale of the SEO.  The ROV derived biomass estimates that form the foundation of this assessment are at the scale of the 4 management areas during different years making it necessary to come up with a statistical way to estimate SEO-wide biomass.  I have made efforts to combine this into one model, but have so far failed to come up with a better solution that takes advantage of all of the data sources. This analysis also included a risk analysis to determine how under- and overestimating bycatch in the halibut fishery affects the assessment.

    1. [`Production_models_JAGS/Code/SPM_helper.R`]: Functions for working with the SPMs. 

    2. *Stage 1*: This stage involves running the spatially stratified SPM using data going back to 1980.  This model has minimal priors placed on key parameters and the goal of this stage is to produce posterior estimates of SEO biomass to use in Stage 2.  The model can be constructed and run using [`Production_models_JAGS/Code/Stage1_run_multiple_models.R`]. The code  allows you to load and run multiple models but I use it to run single models. The original model construction code [`Production_models_JAGS/Code/Stage1_run_model.R`] is no longer used, but included in the repo for model building purposes, if needed (I use the multiple_model scripts and modify the models in the [`Production_models_JAGS/Models/`] folder in a text editor) The model for stage one is currently [`Production_models_JAGS/Models/v22.3_Stage1`].
   
    3. *Stage 2*:  This stage involves  simple, non-stratified, production model that uses the posterior distributions of SEO biomass from 1996 through now and catch estimates going back to 1881.  Catch data consist of broad estimates of removals by the foreign fleet during the 60's and estimates of discards in the halibut fishery derived from 7.ii [`Data_processing_JAGS/Code/IPHC_Survey_Expected_Bycatch.R`].  Extra variance is included in estimating historical bycatch and discarding in the halibut fishery.  As with stage 1 there are two available scripts,[`Production_models_JAGS/Code/Stage2_run_multiple_models.R`] and [`Production_models_JAGS/Code/Stage2_run_model.R`], with the *multiple_models* script being actively used and the single model script degraded, but available if needed. The stage 2 model used for v22.3 is [`Production_models_JAGS/Models/v22.3_Stage2`].  
   
    4. *K and* $\phi_{80}$ *priors*:  With the results from Stage 1 and 2, use [`Production_models/Code/K and phi prior development.R`] to extract prior distributions for *K* and $\phi_{80}$.
   
    5. *R prior*:  [`Production_models_JAGS/Code/r prior development.R`] provides methods for developing a prior for *r*, the rate of intrinsic growth at low population density.  This script follows the methods of McAllister et al. 2001 and produces beta and gamma distributions based on the results from Stage 2.  The final model used for v22.3 used an uninformative beta prior, $Beta$(1,1), which resulted in beta distributions similar to the priors.  
   
    6. *Stage 3*:  This model is the same as Stage 1, a spatially stratified SPM, but now includes priors on *K* and $\phi_{80}$.  The *K* for the 4 management areas are estimated using a dirichlet distribution such that the sum of the management area *K's* add up to *K*.  Results from this model are used to estimate *MSY*, $F_{msy}$, stock status ($\phi_{now}$), and other biological reference points as well as project the population into the future and perform risk analysis.  This model may be run with [`Production_models_JAGS/Code/Stage3_run_multiple_models.R`] and [`Production_models_JAGS/Code/Stage3_run_model.R`].  As with stage 1 and 2, the *multiple_model* script is the one actively used while the single model script is somewhat degraded but available if needed.  The stage 3 model used for v22.3 is [`Production_models_JAGS/Models/v22.3_Stage3`].
   
    7. *Projections*: [`Production_models_JAGS/Code/Projections.R`]  This is code to project the results of Stage 3 models into the future under various harvest policies.  
   
    8. *Simulations*: This folder [`Production_models_JAGS/Code/Simulations/.R`] contains code for running simulations to determine the models ability to accruately estimate model parameters and associated biological reference points.  

        1. [`Production_models_JAGS/Code/Simulations/PHASE3_simulator.R`]: Modify the desired model to simulate data from known parameter values.
      
        2. [`Production_models_JAGS/Code/Simulations/Data_simulator.R`]: Create simulated data sets from a given model output and compare simulated values to posteriors from the original model.This is mostly a visual tool to compare simulated data to posteriors from the model.
      
        3. [`Production_models_JAGS/Code/Simulations/Model_simulated_data.R`]: Run the simulations.  Given that these models were taking over 12 hours to converge, this is a long term commitment.  
      
        4. [`Production_models_JAGS/Code/Simulations/Sim_exam.R`]: Examine the results of the simulations to see how well the model predicted parameter values.  Results stored in [`Production_models_JAGS/Code/Output/Sims/`] and [`Production_models_JAGS/Code/Figures/Sims/`].

    9. ***Outstanding issues with the SS-SPMs***: *for CIE reviewers* This model was reviewed by the CIE committee in September 2023.  Reviews are available upon request.  This assessment and code has been modified and updated to accomodate CIE recomendations.  


   
    
    

   

  
