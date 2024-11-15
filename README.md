# Idealized_leaftrait_responses

### Repository description
This repository contains the R functions necessary for simulating idealized leaf-trait responses to environmental change scenarios in R as first described in the manuscript of Odé et al. (2024) "Temporal constraints on leaf-level trait plasticity for next-generation land surface models". 

### Summary
This model can be used to simulate idealized leaf-trait respones in ci:ca as well as stomatal conductance to environmental changes, according to our new conceptual framework of temporal leaf trait dynamics. Using the EEO-based P-model (see Prentice et al., 2014) and principles from stomatal anatomy (Franks et al., 2012; Mcelwain et al., 2016), our conceptual framework considers (i) temporal separation of dynamics in the leaf interior to atmospheric CO2 concentration (ci:ca) from the optimal ci:ca ratio (χ(optimal)), which is represented by the y-axis of the framework, and (ii) dynamics in the operational stomatal conductance (gs(operational)) within the constraint of anatomical maximum stomatal conductance (gsmax), which is represented by the x-axis of the framework. To simulate leaf trait responses this model uses the R ‘plantecophys’ package (version 1.4-6), with the Photosyn function (Duursma, 2015; DOI 10.1371/journal.pone.0143346), combined with the EEO-based ‘P-model’ for acclimation in leaf biochemistry, using the  “optimal_vcmax_R” repository as reported in (Smith et al., 2019; DOI 10.5281/zenodo.14026447). 

Model variables and abbreviations are provided in the model code.

#### References
- Duursma, R. A. (2015). Plantecophys - An R package for analysing and modelling leaf gas exchange data. PLoS ONE, 10(11). https://doi.org/10.1371/journal.pone.0143346

- Franks, P. J., Leitch, I. J., Ruszala, E. M., Hetherington, A. M., & Beerling, D. J. (2012). Physiological framework for adaptation of stomata to CO2 from glacial to future concentrations. Philosophical Transactions of the Royal Society B: Biological Sciences, 367(1588), 537–546. https://doi.org/10.1098/rstb.2011.0270

- Mcelwain, J. C., Yiotis, C., & Lawson, T. (2016). Using modern plant trait relationships between observed and theoretical maximum stomatal conductance and vein density to examine patterns of plant macroevolution. New Phytologist, 209(1), 94–103. https://doi.org/10.1111/nph.13579

- Prentice, I. C., Dong, N., Gleason, S. M., Maire, V., & Wright, I. J. (2014). Balancing the costs of carbon gain and water transport: Testing a new theoretical framework for plant functional ecology. Ecology Letters, 17(1), 82–91. https://doi.org/10.1111/ele.12211

- Smith, N. G., Keenan, T. F., Colin Prentice, I., Wang, H., Wright, I. J., Niinemets, Ü., Crous, K. Y., Domingues, T. F., Guerrieri, R., Yoko Ishida, F., Kattge, J., Kruger, E. L., Maire, V., Rogers, A., Serbin, S. P., Tarvainen, L., Togashi, H. F., Townsend, P. A., Wang, M., … Zhou, S. X. (2019). Global photosynthetic capacity is optimized to the environment. In Ecology Letters (Vol. 22, Issue 3, pp. 506–517). Blackwell Publishing Ltd. https://doi.org/10.1111/ele.13210


### Contact
Any questions or issues can be submitted via GitHub or directed to Astrid Odé (<a.ode@uu.nl>).
