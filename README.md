# Rapid-SL
Rapid-SL performs embarrassingly parallel computations to find all synthetic lethal reactions or genes with a predefined maximum cardinality. 

### Requirements
The following tools are needed to use RapidSL for finding synthetic lethal sets:
1. [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/)
2. CPLEX v12.9.0 or higher.

### Rapid-SL overview:

![alt text](https://github.com/CSBLaboratory/RapidSL/blob/main/RapidSL_abstract_flowchart.png)


### To get the synthetic lethal sets follow the steps below:
##### Initialize COBRA toolbox:
```
initCobraToolbox
``` 

##### Load the COBRA model:
``` 
load ecoli_core_model.mat
```

##### Specify the medium condition (if needed):
```
model = changeRxnBounds(model, 'EX_glc(e)', -10, 'l');
model = changeRxnBounds(model, 'EX_o2(e)', -20, 'l');
```


##### Specify the list of reactions to be ignored for lethality analysis (if exists). For example, use "findExcRxns" function from COBRA toolbox to exclude exchange reactions from the analysis:
``` 
eliList = model.rxns(findExcRxns(model));
```

##### Specify the maximum desired cardinality:
```
MaxCardinality = 3;
```

##### Specify the cutoff-ratio (set it to an empty matrix to use the default value of 1% of the maximum growth-rate of the wild-type strain):
```
cutOff = [];
```

##### Pass the inputs to Rapid-SL:
```
LethalSets = RapidSLOuterLoop(model, MaxCardinality, cutOff, eliList)
```

### Citing Rapid-SL
If you use Rapid-SL in your work, please cite:
> Manshadi et al.,"Rapid-SL: An efficient implementation of Fast-SL to identify higher order synthetic lethals", in prepration.
