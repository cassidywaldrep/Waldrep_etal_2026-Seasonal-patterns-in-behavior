# Time–activity budgets of eastern mallards: Seasonal patterns and differences among migratory strategies

This repository contains the data and code used in the following publication:

All analyses were completed using R. See manuscript for details about analyses.

Please contact the corresponding author with questions about this data package or to seek potential collaborations using these data.
___
### Description



___
### Authors
Cassidy L. Waldrep    
Department of Biology, University of Saskatchewan, Saskatoon, SK, Canada

John M. Coluccy     
Ducks Unlimited, Inc., Great Lakes/Atlantic Region, Ann Arbor, MI, USA

Nathaniel R. Huck     
Minnesota Department of Natural Resources, Brainerd, MN, USA

Josh C. Stiller     
New York State Department of Environmental Conservation, Albany, NY, USA

Jacob N. Straub    
SUNY Brockport, Brockport, NY, USA

Mathieu Tétreault         
Canadian Wildlife Service, Environment and Climate Change Canada, Québec, QC, Canada

Mitch D. Weegman     
Department of Biology, University of Saskatchewan, Saskatoon, SK, Canada


___
### Files
`github_modelcode.R`: NIMBLE model for multinomial logistic regression for modeling feeding, preening, resting, and flying, indexing on substage and migratory strategyy. Note that the data used for this model is not on the GitHub. To protect the integrity of the project's data, instead, a summarized version per substage is attached.

### Data file column names
`bandnum_index`: 1021 unique individuals
`substage_index`: 8 substages used for analysis. 1 = early breeding, 2 = fall migration, 3 = late spring migration, 4 = prefall migration, 5 = prespring migration, 6 = spring migration, 7 = summer, 8 = winter 
`study_year_index`: years 2022-2025
`pop_index`: migratory status. 1 = migratory, 2 = partial migrant, 3 = resident
`total_num_fixes`: total number of behavioural fixes per individual per substage
`num_fixes_fly`: total number of flying fixes per individual per substage
`num_fixes_preen`: total number of preening fixes per individual per substage
`num_fixes_feed`: total number of feeding fixes per individual per substage
`num_fixes_stationary`: total number of stationary fixes per individual per substage

