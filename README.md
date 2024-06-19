# Pandemic-System-Identification

## Description

This repository contains the code and data for the system identification part related to the Model Predictive Pandemic Control (MPPC) research project by TU/e MOVEMENT Research Group.

## Files

`MHEsidthe.jl`: Initial Julia script containing the main code for system identification. The method used for the identification is Moving Horizon Estimation (MHE)

`SIDTTHE_data_DEF.mat`: Dataset used in the system identification process. These data have been collected from the [Italian COVID-19 repository](https://github.com/pcm-dpc/COVID-19) and pre-processed. Ensure this file is in the same directory as the Julia script for the code to run correctly

`DataCooking.jl` and `DataCollection.jl`: Julia script to obtain data in order to fit an age-stratified model. The data are not directly available, thus raw datasets have been manipulated n order to be employed. 

`PlotCode.jl`: Code to plot results of interest in order to keep other codes more clean

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