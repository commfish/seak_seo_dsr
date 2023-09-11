# seak_seo_dsr
# Southeast Alaska Outside Demersal Shelf Rockfish Assessment

The Southeast Outside subdistrict (SEO) of the Gulf of Alaska (GOA) demersal shelf rockfish (DSR) complex (yelloweye, quillback, copper, rosethorn, China, canary, and tiger rockfish) is assessed on a biennial cycle, with a full stock assessment conducted every second year. The SEO area is comprised of four management sections: East Yakutat (EYKT), Northern Southeast Outside (NSEO), Central Southeast Outside (CSEO), and Southern Southeast Outside (SSEO) (Figure 1). Yelloweye rockfish (*Sebastes ruberrimus*) have been managed as a Tier 4 stock; however, given that the current assessment models do not produce reliable estimates of $F_{35}$ or $F_{40}$, the stock has been determined to be more appropriately in Tier 5, where $F_{OFL}$ = *M* and $maxF_{ABC}$ = *0.75* * *M*. The other species in the complex are managed as Tier 6.   Prior to 2022, the stock assessment was based on biomass estimates of yelloweye rockfish derived from Alaska Department of Fish and Game (ADF&G) submersible and ROV surveys.  Harvest recommendations were established by applying an assumed mortality rate to the lower 90% confidence interval of the estimate as a hedge against uncertainty.  The recommended acceptable biological catch (ABC) and overfishing level (OFL) for the current assessment are based on a random effects models (REMA) applied to the ROV biomass estimates and CPUE estimates of yelloweye rockfish in the International Pacific Halibut Commission (IPHC) longline survey.  A new Bayesian state-space surplus production model (SS-SPM) is in development and schedule for CIE review in September 2023.

<img src="https://github.com/commfish/seak_seo_dsr/blob/main/readme/seomap.png" width="600" />

*Figure 1. The Southeast Outside Subdistrict of the Gulf of Alaska (SEO) with the Alaska Department 
of Fish and Game groundfish management areas used for managing the demersal shelf rockfish fishery: 
East Yakutat (EYKT), Northern Southeast Outside (NSEO), Central Southeast Outside (CSEO), and 
Southern Southeast Outside (SSEO) Sections.*

The assessments associated with this code can be found here:

* [2022 SAFE report](https://www.fisheries.noaa.gov/resource/data/2022-assessment-demersal-shelf-rockfish-stock-complex-southeast-outside-subdistrict)

### Repo directory: 

This repository will hold data, files and r script for the SEO DSR assessment that ADF&G performs as part of the NPFMC stock assessments.  This repository contains the following sub directories that contain different analysis that are part of the assessment.  They are:

1. **Data processing:** This folder contains script and raw data files for putting together data for use in the analysis contained in the folders below.  This includes catch data and biological data.
2. **Density estimation:**  Script and data for estimating yelloweye rockfish density using distance sampling methods and transect data from ROV surveys.
4. **REMA analysis:** This folder contains the scripts for running the random effects model that is the current model used in the federal assessment.
5. **Production models:** This folder contains scripts for developing and running the bayesian state-space surplus production model (SS-SPM).
6. **Output:** This folder will contain output files from various analysis.
7. **Figures:** This folder will contain figures associated with the assessment.
8. **r_helper:** This folder contains functions built specific to the assessment. 

## Code path and data sources for the assessment:

1. [`r_helper/`]: Folder with pre-canned functions referenced throughout the assessment

2. **Update catch and portside sampling data**: Download the most recent data from these sites and merge with older data using [`Data_processing/Code/Harvest reconstruction.R`]

    1. [*OceanAK Portside sampling for biological data*](https://oceanak.dfg.alaska.local/analytics/saw.dll?Answers&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2Fport%20sampling%20bio%20data)  
   
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
    
    2. Get the hook adjustment factor data to correct CPUE estimators for hook saturation. This data is not posted regularly and may require a [data request](https://www.iphc.int/form/data-request).  If someone has already requested the most recent data it may be posted at [IPHC hook adjustments](https://www.iphc.int/data/fiss-data-query), Document IPHC-YEAR-FISS-HADJ.  The online data is a large, heavily formated excel sheet that will need some work to convert to csv file.      

    3. [FISS Halibut data](https://www.iphc.int/datatest/fiss-pacific-halibut-data): As of 2023 all the data you need is described above, but halibut data is also available here so I'm including it in case things change again.  In this arrangement of the data base you have to download the set data and halibut data separately and then merge them in the code.  

    4. [`Data_processing/Code/IPHC_Survey_CPUE_index.R`]: Code for calculating IPHC CPUE index of abundance for use in the REMA model and the SS-SPM.  Going into the 2023 CIE review I am debating what stations to include in the CPUE calculations.  Right now I have only been using stations that encounter yelloweye 40% of the time, but will discuss at review...
   
9.  **Calculate expected yelloweye bycatch in the halibut fishery**: This section involves estimating historical discards of yelloweye in the halibut fishery by applying bycatch rates estimated during contemporary times to historical catches and propogating uncertainty.  

    1. [`Data_processing/Code/Halibut_harvest_reconstruction.R`]: Code for compiling halibut harvest data and apportion catches to SEO management areas.
   
    2. [`Data_processing/Code/IPHC_Survey_Expected_Bycatch.R`]: Estimate the expected bycatch from the halibut fishery using the IPHC survey data and historical halibut harvests.  CAS estimates of bycatch from the NOAA observer program are stored in [`Data_processing/Data/Harvests/CONF/CAS Yelloweye Rockfish Bycatch in Halibut Fishery 2013-present.csv`] and [`CAS Yelloweye Bycatch with Stat 2013-present.csv`].
   
10. **REMA model**: Run the random effects model using the biomass estimates and IPHC cpue estimates [`REMA/Code/REMA_run.R`].  This is the model used in the 2022 SAFE report.  **Model 22.2** is the current assessment model for the NPFMC as of 2022.  This code contains multiple rema models, methods for comparisons and figures. 

11. **Surplus Production Models**: This folder contains code for running the state-space surplus production model (SS-SPM) in JAGS, plotting results, saving output and comparing different models. Details of this model are included in the appendices of the 2022 SAFE report linked above. There are 3 stages to running the model, all of which ivole a Pella-Tomlinson surplus production model currently configured to maximize the production curve at $B_{40}$ (biomass at 40% of virgin biomass).  The first two stages involves running the models under minimal assumptions (i.e., priors) and using those results to construct informative priors for *K* (carrying capacity for the SEO as a whole) and $\phi_{80}$ (the proportion of K present in 1980).  Part of the reason for this is the difficulty in the long term data series (back to 1888) and the spatial resolution available in the data.  Since 1980 data is available for the 4 management areas in the SEO.  Prior to that data is only available at the scale of the SEO.  The ROV derived biomass estimates that form the foundation of this assessment are at the scale of the 4 management areas during different years making it necessary to come up with a statistical way to estimate SEO-wide biomass.  I have made efforts to combine this into one model, but have so far failed to come up with a better solution that takes advantage of all of the data sources. This analysis also included a risk analysis to determine how under- and overestimating bycatch in the halibut fishery affects the assessment.

    1. [`Production_models/Code/SPM_helper.R`]: Functions for working with the SPMs. 

    2. *Stage 1*: This stage involves running the spatially stratified SPM using data going back to 1980.  This model has minimal priors placed on key parameters and the goal of this stage is to produce posterior estimates of SEO biomass to use in Stage 2.  The model can be constructed and run using [`Production_models/Code/Stage1_run_multiple_models.R`]. The code  allows you to load and run multiple models but I use it to run single models. The original model construction code [`Production_models/Code/Stage1_run_model.R`] is no longer used, but included in the repo for model building purposes, if needed (I use the multiple_model scripts and modify the models in the [`Production_models/Models/`] folder in a text editor) The model for stage one is currently [`Production_models/Models/v22.3_Stage1`].
   
    3. *Stage 2*:  This stage involves  simple, non-stratified, production model that uses the posterior distributions of SEO biomass from 1996 through now and catch estimates going back to 1881.  Catch data consist of broad estimates of removals by the foreign fleet during the 60's and estimates of discards in the halibut fishery derived from 7.ii [`Data_processing/Code/IPHC_Survey_Expected_Bycatch.R`].  Extra variance is included in estimating historical bycatch and discarding in the halibut fishery.  As with stage 1 there are two available scripts,[`Production_models/Code/Stage2_run_multiple_models.R`] and [`Production_models/Code/Stage2_run_model.R`], with the *multiple_models* script being actively used and the single model script degraded, but available if needed. The stage 2 model used for v22.3 is [`Production_models/Models/v22.3_Stage2`].  
   
    4. *K and* $\phi_{80}$ *priors*:  With the results from Stage 1 and 2, use [`Production_models/Code/K and phi prior development.R`] to extract prior distributions for *K* and $\phi_{80}$.
   
    5. *R prior*:  [`Production_models/Code/r prior development.R`] provides methods for developing a prior for *r*, the rate of intrinsic growth at low population density.  This script follows the methods of McAllister et al. 2001 and produces beta and gamma distributions based on the results from Stage 2.  The final model used for v22.3 used an uninformative beta prior, $Beta$(1,1), which resulted in beta distributions similar to the priors.  
   
    6. *Stage 3*:  This model is the same as Stage 1, a spatially stratified SPM, but now includes priors on *K* and $\phi_{80}$.  The *K* for the 4 management areas are estimated using a dirichlet distribution such that the sum of the management area *K's* add up to *K*.  Results from this model are used to estimate *MSY*, $F_{msy}$, stock status ($\phi_{now}$), and other biological reference points as well as project the population into the future and perform risk analysis.  This model may be run with [`Production_models/Code/Stage3_run_multiple_models.R`] and [`Production_models/Code/Stage3_run_model.R`].  As with stage 1 and 2, the *multiple_model* script is the one actively used while the single model script is somewhat degraded but available if needed.  The stage 3 model used for v22.3 is [`Production_models/Models/v22.3_Stage3`].
   
    7. *Projections*: [`Production_models/Code/Projections.R`]  This is code to project the results of Stage 3 models into the future under various harvest policies.  
   
    8. *Simulations*: This folder [`Production_models/Code/Simulations/.R`] contains code for running simulations to determine the models ability to accruately estimate model parameters and associated biological reference points.  

        1. [`Production_models/Code/Simulations/PHASE3_simulator.R`]: Modify the desired model to simulate data from known parameter values.
      
        2. [`Production_models/Code/Simulations/Data_simulator.R`]: Create simulated data sets from a given model output and compare simulated values to posteriors from the original model.This is mostly a visual tool to compare simulated data to posteriors from the model.
      
        3. [`Production_models/Code/Simulations/Model_simulated_data.R`]: Run the simulations.  Given that these models were taking over 12 hours to converge, this is a long term commitment.  
      
        4. [`Production_models/Code/Simulations/Sim_exam.R`]: Examine the results of the simulations to see how well the model predicted parameter values.  Results stored in [`Production_models/Code/Output/Sims/`] and [`Production_models/Code/Figures/Sims/`].

    9. ***Outstanding issues with the SS-SPMs***: *for CIE reviewers*

        1. *Model underestimates r:* Assuming that the simulations are set up correctly, the simulations indicate that the model is underestimating the r parameter in the production model and hence underestimating overall productivity of the stock.  My hunch is that this may be due to misspecification in the r hyper-prior, or issues with constraints on process error. ***I am very eager to work with CIE reviewers on this issue!*** My thoughts going forward are:
           
            1. Try a single r parameter shared by the four management areas.  But the stage 2 model is a simple, unstratifed production model and it appears to underestimate r as well?
            2. Loosening bounds on the process error. Generally, the looser the process error bounds, the higher the estimate of r.  Left unbounded the model will try to "connect the dots" of the biomass point estimates, so there is a balance to this approach and process error *should* be small for such a long-lived species.
            3. Apply the process error to only the biological portion of the model and exclude the catch part of the equation.  The model is currently set up according to Ono et al. (2012) with process error applied to the entire equation (biological production minus catch).  But I also model a lot of error into the catch portion of the model because of uncertainty coming from the bycatch estimates in the haibut fishery.  Given how much uncertainty is applied to the catch part of the equation, perhaps this is conflating the model's ability to estimate r?  

        3. *Long run times:* One Stage of the model can take over 12 hours to converge.  I will try reformulating the model such that $B_t$ = $\phi_t$ * *K*.  Supposedly the algorythm has an easier time finding $\phi$ than it does *B*.  Switching from JAGS to STAN may also allow faster processing as well.  
      
        4. *3 stage approach:* I know this part seems like it could be controversial and raises eyebrows.  My thoughts are that stage 1 and 2 are done to produce some vaguelly informative priors so it's OK.  At the very least, I think labelling this as a 3-stage approach is misleading and there should just be a section on how the *K* and $\phi_{80}$ priors are developed.  That being said, I am very eager to work with the CIE reviewers on this.
      
        5. *IPHC data:* I have met with IPHC staff to discuss how I am using the survey data in this assessment and they were, by and large, very supportive of what is presented above.  They had several recomendations for fine tuning the approach:
      
            1. Filtering out "ineffective" stations. I wasn't aware of this before, but this is a very simple task.
          
            2. Correcting CPUE/NPUE/WCPUE for hook saturation effects.  I have gotten the requested data on correcting for hook saturation and will be incorporating this into the analysis.
          
            3. Filtering stations for only those that are in yelloweye habitat.  The station locations can vary by up to 3 nm per year, so a given station may drift in and out of yelloweye habitat.  This may be best accomplished with a habitat map... likely using the one defined for the ROV survey.  The ADF&G State Rockfish Initiative (SRI) has assigned various people to rework the habitat maps but a timeline for completing that task is not available. I think the best approach for now is to filter for the stations that are within the habitat polygons used in the ROV survey and will look to adopt that approach going forward.  I have also considered an alternative approach whereby I use stations that encounted yelloweye at least X% of the time for the cpue index and have played with values from 10-40% for selecting stations.  *I am very eager to work with CIE reviewers on this issue!*
 
    

   

  
