clc
close
clear
load COAD_SGD.mat
rules = model_cancer.rules;
n = 1;
for i = 1 : length(rules)
   str = rules{i};
   if sum((str == '|'))>3
      RxnId(n, 1) = i;
      n = n + 1;
   end
end
mat = model_cancer.rxnGeneMat;
offGenesIDs = find(sum(mat) == 0)';
mat(RxnId, :) = zeros(length(RxnId), length(model_cancer.genes));
orGenesIDs = setdiff(find(sum(mat) == 0)', offGenesIDs);

