# RNAseq_pipeline
Pipeline to process and analyse RNAseq data.

To download repository including submodules, please use the command:
```
git clone --recursive https://github.com/lmb-cell-biology/RNAseq_pipeline.git
```

Please don't forget to update submodules everytime you update the repository:
```
git pull
git submodule update
```

### Create environmental variables in your bashrc file
Add the following lines in your bashrc file
```
export RNAseq_analysis="[path_to_repository]/RNAseq_analysis.R"
export cummeRbund="[path_to_repository]/exploratory_analysis_cummeRbund.R"
```
