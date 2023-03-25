# RSVSystematicReview

This repository contains the data and code used to produce the results
presented in “Age-dependent risk of respiratory syncytial virus
infection: A systematic review and hazard modeling from serological
data”.

The raw data and script for data processing for all the different
studies is located in the folder “Data”:

-   `Data_seroprevalence.csv` (data from reports on seropositivity)
-   `data_under5.txt` (pre-F data of individuals &lt; 5 years)
-   `data_over5.txt` (pre-F data of groups &gt; 5 years)
-   `cleanedData.R` (scripts used to clean the data)
-   `cleanedData.RDS` (cleaned data)

The main scripts are located in the folder “MainModel”:

-   `FoIEstimation.R` (estimation of FoIs and decay rate assuming SIS
    catalytic model)
-   `MixtureModel.R` (estimation of the means and variances for the
    susceptible group and the ever-infected/recovered group assuming SIS
    model)

Two sensitivity analyses are located in the folder
“SensitivityAnalysis”:

-   `SensitivityAnalysis_1.R` (estimation of FoIs using non-informative
    prior)
-   `SensitivityAnalysis_2.R` (mixture model using SIW catalytic model)
