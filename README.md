# Pandemic-System-Identification

## Description

This repository contains the codes and data for the system identification part related to the Model Predictive Pandemic Control (MPPC) research project by TU/e MOVEMENT Research Group.

## Julia Files

`DataCooking.jl` and `DataCollection.jl`: Julia script to obtain data in order to fit an age-stratified model. The data are not directly available, thus raw datasets have been manipulated in order to be employed. 

`TestingSIDTHE.jl`: Code to test the behaviour of the system to parameters obtained from moving horizon estimation (MHE).

## Data Folders 

<b><u>Age Compartments Distributions</u></b>: 
Data from [Istituto Superiore di Sanità](https://covid19.infn.it/iss/) (ISS) in collaboration with INFN.

`issITAdeceased.csv`: Trend of percentages by age groups related to DECEASED INDIVIDUALS (Deaths, reported on the day of death) in Italy.

`issITAhospitalised.csv`: Trend of percentages by age groups related to HOSPITALISED INDIVIDUALS (New hospitalizations, reported on the date of admission) in Italy.

`issITAicu.csv`: Trend of percentages by age groups related to INTENSIVE CARE INDIVIDUALS (New ICU admissions, referred to the date of admission) in Italy.

`issITApositivi.csv`: Trend of percentages by age groups related to POSITIVE INDIVIDUALS (New Positive Tests, referred to the date of Testing) in Italy.

<b><u>Age New Cases Data</u></b>: 
Data from [Istituto Superiore di Sanità](https://covid19.infn.it/iss/) (ISS) in collaboration with INFN.

`issITAnew_ICU.csv`: Trend of daily numbers by age group of new ICU cases in Italy. Data reported refers to the 7-day interval moving average.

`issITAnew_ricoveri.csv`: Trend of daily numbers by age group of new Hospitalisations cases in Italy. Data reported refers to the 7-day interval moving average.

`issITAnewdeceased.csv`: Trend of daily numbers by age group of new Deceased (Deaths) cases in Italy. Data reported refers to the 7-day interval moving average.

`issITAnewpos.csv`: Trend of daily numbers by age group of new Positive cases in Italy. Data reported refers to the 7-day interval moving average.

<b><u>National Trends</u></b>: 
Data from [Protezione Civile Github Dataset](https://github.com/pcm-dpc/COVID-19) in collaboration with INFN.

`WeekTrend_ITA.csv`: Trend of daily numbers of CURRENTLY POSITIVE, HEALED, DECEASED and TOTAL CASES. Data in 7-day moving average.

`DailyTrend_ITA.csv`: Trend of daily numbers of:
 - Currently Hospitalised
 - Currently in ICU
 - Currently Positive
 - Δ variation Positive
 - New Positive Cases 
 - Currently Healed
 - Currently Deceased
 - Cumulated Number of Positive Cases
 - New Δ(+) ICUs 
 - Δ variation ICUs
 - Δ variation Hospitals

    & extra categories like, Δ variation in total positive cases, Δ variation in new positive cases etc...

<b><u>Reconstructed Datasets</u></b>: 
Reconstructed datasets of the Italian population divided in four age groups:
## Age Groups

1. **_u40_**:
   Individuals aged **0 to 39 years**. It includes children, young adults, and those in early adulthood.
  
2. **_mid_**:
   Individuals aged **40 to 59 years**, considered the largest amount of working-age population.

3. **_old_**:
   Individuals aged **60 to 79 years**. This group includes senior individuals, Heightened risks due to age-related vulnerabilitie, more prone to experience severe outcomes from the virus.

4. **_ger_**:
   Individuals aged **80 years and above**. Commonly referred to as geriatric population. Heavily affected by the COVID-19 pandemic, facing high mortality rates.

   All the data are extrapolated from the code `DataCooking.jl` and `DataCollection.jl`.
   The dataset are used to inform and identify the model.

## MHE Matlab Folders 

<b><u>MHE Age Stratified - Single alpha</u></b>: 
Folder containing all the codes to run the simulation of the SIDTHE age-stratified model, with single $\alpha$ (virus infectivity).
Main files in the folder: 

`testScriptMHE.m`: Main test script to run the MHE on the whole horizon.
`bayesMHEObj.m`: Function performing Hyperparameters Autotuning of the MHE objective funztion by implementig bayesian optimization. 
`runMHE.m`: Function containing MHE optimization problem in CasADi framework. 
`PlottingMHE.m`: Code to plot the results of the MHE. 

<b><u>MHE Age Stratified</u></b>: 
Folder containing all the codes to run the simulation of the SIDTHE age-stratified model, with multiple $\alpha$ (virus infectivity).
Main files in the folder: 

`testScriptMHE_a.m`: Main test script to run the MHE on the whole horizon.
`bayesMHEObj_a.m`: Function performing Hyperparameters Autotuning of the MHE objective funztion by implementig bayesian optimization. 
`runMHE_a.m`: Function containing MHE optimization problem in CasADi framework. 
`PlottingMHE_a.m`: Code to plot the results of the MHE. 

