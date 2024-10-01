## What `HZZ` framework do?
This framework can be used to process NanoAOD samples, including processing various variables, correction, and selection. 
1. 5 regions including Dilepton, EGammaFromMisid, SingleElectron, SinglePhoton, and ZGamma have been designed for signal selection and estimation of DY Data Driven in the VBS ZZ analysis.
- config/samples_*.txt provides samples.
- config/REGION(eg.dilepton)/*.yaml provides triggers,corrections and selections.
- src/*Trees.cc defines the variables to be stored.
- You can design new regions based on your own analysis by adding the above files.
2. 