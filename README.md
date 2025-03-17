# Spatial Loadings
## Spatially Varying Loadings in Modeling the Opioid Syndemic

In this project our goal is to investigate whether spatially modeling the loadings provides a more nuanced approach to understanding the opioid syndemic. There are six different outcomes considered in this project: 
- Death counts: county-level counts; included are the counts of unintentional overdose deaths involving illicit opioids.
- ED visits: county-level counts; included are all ED visits related to medication and drug overdoses.
- Treatment counts: county-level counts; included are the number of uninsured individuals and Medicaid beneficiaries with opioid use disorder served by treatment programs.
- Buprenorphine prescriptions: county-level counts; counts include residents who receive buprenorphine prescriptions.
- HCV infections: county-level counts; included acute and chronic cases of HCV cases.
- HIV infections: county-level counts; included are newly diagnosed HIV infections.

Description of the files:
- **NCdata.xlsx** includes the data;
- **Adj_Matrix_NC.RData** is the adjacency matrix used in our ICAR model;
- **SpatialLoadings.R** is the R-code we used to estimate the factor with spatially varying loadings.
