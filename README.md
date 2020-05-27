# 2015-Coral-bleaching-recovery
This repository contains the data and code for the analyses for the manuscript:

Matsuda, S.B., Huffmyer, A.S., Lenz, E., Davidson, J., Hancock, J., Przybylowski, A., Innis, T., Gates, R., Barott, K. (2020). Coral bleaching susceptibility is predictive of subsequent mortality within but not between coral species, Front. Ecol. Evol., doi: 10.3389/fevo.2020.00178

Abstract:
Marine heat waves instigated by anthropogenic climate change are causing increasingly frequent and severe coral bleaching events that often lead to widespread coral mortality. While community-wide increases in coral mortality following bleaching events have been documented on reefs around the world, the ecological consequences for conspecific individual colonies exhibiting contrasting phenotypes during heat stress (e.g. bleached vs. not bleached) are not well understood. Here we describe the ecological outcomes of the two dominant reef-building coral species in K?ne?ohe Bay, Hawai?i, Montipora capitata and Porites compressa, by monitoring the fates of individuals that exhibited either a bleaching susceptible phenotype (bleached) or resistant phenotype (non-bleached) following the second of two consecutive coral bleaching events in Hawai?i in 2015. Conspecific pairs of adjacent bleaching susceptible vs. resistant corals were tagged on patch reefs in two regions of Kaneohe Bay with different seawater residence times and terrestrial influence. The ecological consequences (symbiont recovery and mortality) were monitored for two years following the peak of the bleaching event. Bleaching susceptible corals suffered higher partial mortality than bleaching resistant corals of the same species in the first 6 months following heat stress. Surprisingly, P. compressa had greater resilience following bleaching (faster pigment recovery and lower post-bleaching mortality) than M. capitata, despite having less resistance to bleaching (higher bleaching prevalence and severity). These differences indicate that bleaching susceptibility of a species is not always a good predictor of mortality following a bleaching event. By tracking the fate of individual colonies of resistant and susceptible phenotypes, contrasting ecological consequences of heat stress were revealed that were undetectable at the population level. Furthermore, this approach revealed individuals that underwent particularly rapid recovery from mortality, including some colonies over a meter in diameter that recovered all live tissue cover from >60% partial mortality within just one year. These coral pairs (44 pairs of each species) continue to be maintained and monitored in the field, serving as a 'living library' for future investigations on the ecology and physiology of coral bleaching.

In this repository, please find data sheets and .Rmd for all data analyses in the manuscript.

This manuscript was first released as a preprint (https://doi.org/10.1101/2019.12.17.880161). For copies of the .Rmd used for this release, please contact SMatsuda.   


**Color pigmentation analysis:**    
Script: ColorScore.Rmd.     
Data: Bleaching_FINAL.csv.     
columns    
*TagID*: Unique identification number for each coral colony  
*Pair*: Unique identification number for each set of bleached and non-bleached adjacent pairs of colonies  
*Species*: Coral species
*Site*: Lagoon: 4=PR4, Inner Lagoon; 13=PR13, Outer Lagoon    
*B2015*: Bleaching phenotype (if a colony bleached during the 2015 bleaching event) Y=Yes, N=No    
*Month0*: Color scores at the time of peack bleaching, October 2015    
*Month1.5-Month24*: Color scores the number of months post peak bleaching (columns= months post peak bleaching)    

**Partial mortality analysis:**  
Script: Partial_Mortality.Rmd  
Data: Mortality_Scores_Percent.csv. 
columns. 
*TagID*: Unique identification number for each coral colony    
*Pair*: Unique identification number for each set of bleached and non-bleached adjacent pairs of colonies    
*Species*: Coral species. 
*Site*: Lagoon: 4=PR4, Inner Lagoon; 13=PR13, Outer Lagoon      
*B2015*: Bleaching phenotype (if a colony bleached during the 2015 bleaching event) Y=Yes, N=No      
*Month0*: Partial mortality scores (20% intervals) at peack bleaching, October 2015      
*Month1.5-Month24*: Partial mortality scores the number of months post peak bleaching (columns= months post peak bleaching)    

**Benthic transect surveys analysis:**. 
Script: Transect_Data.Rmd. 

Data:   
Benthic_Transect_Coral_Cover.csv   
*Reef*: Inner = Inner Lagoon, Outer = Outer Lagoon  
*Depth*: Depth of transect (meters) 
*Month*: 0 (peak bleaching) and number of months post peak bleaching    
*Survey*: Survey replicate    
*CoralCover*: % coral cover 


Benthic_Transect_df.csv    
columns    
*Species*: Coral species (Mcap = *Montipora capitata*, Pcomp = *Porites compressa*)    
*Bleaching_status*: Bleaching phenotype (Bleached, Pale, Healthy)    
*BNB*: Bleaching color score (1=Bleached, 2=Pale, 3=Healthy)    
*Date*: Date of transect    
*Depth*: Depth of transect (meters)    
*Reef*: Inner = Inner Lagoon, Outer = Outer Lagoon    
*Percent*: % of each bleaching phenotype out of total    
*Rep*: Transect replicate by depth     

**Temperature analysis:**        
Script: Temp_Analysis.Rmd        
Data:         
Temp_Data_PR12.csv    
columns    
*Date*: Date    
*Time*: Time 12hr clock    
*Time24*: Date 24hr clock    
*Temp F*: Temperature in degrees F    
*Temp C*: Temperature in degrees C    
*Corrected Temp C*: Temp C logger calibrated data    
    
Temp_Data_PR1.csv (PACIOOS sensor on Coconut Island)    
columns    
*Date*: Date    
*Time*: Time 12hr clock    
*Time24*: Date 24hr clock    
*Temp C*: Temperature in degrees C    

Temp_Data_PR4_PR13.csv      
*Date*: Date    
*Time*: Time 12hr clock    
*Time24*: Date 24hr clock    
*Temp C*: Temperature in degrees C    
*Corrected Temp C*: Temp C logger calibrated data    

