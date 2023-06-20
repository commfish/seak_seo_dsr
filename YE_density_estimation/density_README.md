# seak_seo_dsr
## DSR ROV distance analysis and density and biomass estimation

TThis project entails the current stock assessment methodology for Demersal Shelf Rockfish (DSR) in Southeast Alaska utilizing Remote Operated Vehicle (ROV) surveys to produce biomass estimates of Yelloweye Rockfish (which comprise the majority of the DSR complex) and assess species composition of the DSR.  The methods uses data from the DSR ROV surveys in Southeast Alaska to estimate Yelloweye density in the SEAK subdistrict that has been surveyed.  Randomly determined ROVs transects of approximately 1 km in length are performed in the subdistrict of interest.  DSR and other groundfish are counted along the transect and the distance from the ROV are measured using stereoscopic methods.  The analysis entails basic distance analysis whereby the perpendicular distances of fish from transect lines are used to model detection probabilities and those detection probabilities are then used to estimate density of fish in the overall study area.  Producing biomass estimates entails combining estimates of density, average weight of Yelloweye, and the total area of Yelloweye habitat in the subdistrict.  Mean weights are derived from harvested fish sampled in ports and also by estimating weights of Yelloweye whose lengths were estimated in the ROV surveys.  Furthermore, we consider the proportion of Yelloweye in the ROV survey whose lengths are greater than the smallest Yelloweyes in the harvests to modify the biomass estimate.  Variances from these components are combined and the lower 90% confidence interval of the biomass estimate are used to set harvest policy.  The major unknown in this analysis is the total amount of Yelloweye habitat in each subdistrict, which is a major component of the biomass calculations.  The estimate of Yelloweye habitat has been developed from commercial longline data combined with NOAA charts and currently has no variance associated with the estimate.  Yelloweye habitat is currently being assessed and future biomass estimates may involve a stratified sampling design based around a more complete understanding of Yelloweye habitat and its associated variance.  
These scripts are recorded in a logical order such that they often produce csv files for use in subsequent scripts:

1.	[`DSR_ROV_SEAK_ROV_Processing.r`]:  This script takes the raw ROV data and processes it for analysis with the distance package in R which produces density estimates.  Files produced:
      <ol type="a">
         <li>Output/new_transect_qc.csv</li>
         <li>Output/SSEO_smooth_transect_output.csv</li>
         <li>Output/SSEO_2020_smooth_predict_lengths.csv</li>
         <li>Output/2020_SSEO_smooth_predict.csv</li>
         <li>Data/SSEO_distance_data_GIStran_for_analysis.csv</li>
      </ol>
      
2.	[`DSR_ROV_SEAK_base_distance_modelling.R`]:  This script explores the processed ROV data and runs a suite of candidate models to estimate yelloweye density.  There is code for comparing model output, ranking them, plotting output, and producing a model list for model averaging and bootstrapping.
      <ol type="a">
         <li>Data/Bootstrap_Model_List.csv</li>
      </ol>
	
3.	[`DSR_ROV_SEAK_LWexam.R`]:  This script compares lengths of ROV and Port sampled fish, examines LW relationships in the port samples, removes possible outliers and organizes the LW port data for computing biomass and for estimating weight from ROV lengths.
	   <ol type="a">
         <li>Data/Data/YE_LW_rel.csv</li>
      </ol>

4.	[`DSR_ROV_SEAK_basic_biomass_calc.R`]:  This script generates biomass estimates from one distance model (presumably the best one) using mean weights from the port samples and estimated mean weights from the ROV measured fish.  Bootstrap used to estimate weights of ROV fish and the proportion of ROV fish above a certain cutoff length.
   
5.	[`DSR_ROV_SEAK_dist_model_avg.R`]:  This script provides the bootstrap for generating a density estimate by averaging models from DSR_ROV_SEAK_base_distance_modelling.R
      <ol type="a">
         <li>Output/Model_Avg_boot_GIStran_1000.csv</li>
      </ol>

6.	[`DSR_ROV_SEAK_Full_bootstrap_biomass_estimation.R`]:  This script combines the ROV bootstrap, the model averaging bootstrap and resampling of the selected port weight samples to generate biomass estimates.
      <ol type="a">
         <li>Output/Model_Avg_biomass_relP_p4_5K.csv</li>
      </ol>



   
