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


##### Based on the type of lethality analysis, specify the list of reactions or genes to be ignored (if exists). For example, use `findExcRxns` function from COBRA toolbox to exclude exchange reactions from the analysis:
``` 
eliList = model.rxns(findExcRxns(model));
```

##### Specify the maximum desired cardinality:
```
MaxCardinality = 3;
```

##### Specify the cutoff-ratio. Set it to an empty matrix to use the default value of 1% of the maximum growth-rate of the wild-type strain:
```
cutOff = [];
```

##### Choose the type of analysis. Use `'Rxn'` for obtaining synthetic lethal reactions and `'Gene'` for synthetic lethal genes. Passing the analysis mode to the RapidSL function is optional. If no variable is specified, the analysis will be run for finding synthetic lethal reactions:
```
Mode = 'Rxn';
```

##### Pass the inputs to Rapid-SL:
```
LethalSets = RapidSL(model, MaxCardinality, cutOff, eliList, Mode)
```

### Citing Rapid-SL
If you use Rapid-SL in your work, please cite:
> Manshadi et al.,"[Rapid-SL](https://www.nature.com/articles/s41598-022-18177-w) identifies synthetic lethal sets with an arbitrary cardinality".
