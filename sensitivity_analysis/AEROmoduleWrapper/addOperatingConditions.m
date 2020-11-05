function Inputs = addOperatingConditions(Inputs,OperatingCondition)

n_dim = length(Inputs.Marginals);
n_op  = length(OperatingCondition);

for i=1:n_op
    j = i + n_dim;
   
    Inputs.Marginals(j).Type = 'Constant';
    Inputs.Marginals(j).Name = OperatingCondition(i).Name;
    Inputs.Marginals(j).Parameters = OperatingCondition(i).Parameters; 
end