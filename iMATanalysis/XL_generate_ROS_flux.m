%readme
%1- iMAT
%2- FLUX minimization to generate reference state
function [Flux_min_iMAT,solution,offRxn,HRxnCandidates] = XL_generate_ROS_flux(model_sinked, GeneNameExpressed,epsilon_f,epsilon_r)
%% iMAT
ecoli = model_sinked;
expressionMimicforROSproject = ones(size(GeneNameExpressed));%convert the hit set to a "expression value" for iMAT algorithm
expression_gene=struct;
expression_gene.value = expressionMimicforROSproject;%all set to 1 for Yans Gene in model 
expression_gene.gene = GeneNameExpressed;%all the Yans Gene in model
[expressionRxns , ~] = mapExpressionToReactions_xl(ecoli, expression_gene);
HRxnCandidates = ecoli.rxns(expressionRxns == 1);
%do iMAT
[solution, MILPproblem] = iMAT_xl(ecoli, expressionRxns,epsilon_f,epsilon_r, 0.25, 0.75);
%% MINIMIZATION of flux
MILPproblem = solution2constraint(MILPproblem,solution);
A = sparse(size(MILPproblem.A,1)+2*(length(ecoli.rxns)),size(MILPproblem.A,2)+length(ecoli.rxns));
[m,n,s] = find(MILPproblem.A);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end
% add flux measuremnet 
for i = 1:length(ecoli.rxns)
    A(size(MILPproblem.A,1)+i,i) = -1;
    A(size(MILPproblem.A,1)+i,size(MILPproblem.A,2)+i) = 1;
end
for i = 1:length(ecoli.rxns)
    A(size(MILPproblem.A,1)+length(ecoli.rxns)+i,i) = 1;
    A(size(MILPproblem.A,1)+length(ecoli.rxns)+i,size(MILPproblem.A,2)+i) = 1;
end
MILPproblem.A = A;
%create other inputs
MILPproblem.lb = [MILPproblem.lb;zeros(length(ecoli.rxns),1)];
MILPproblem.ub = [MILPproblem.ub;1000*ones(length(ecoli.rxns),1)];
MILPproblem.b = [MILPproblem.b;zeros(2*length(ecoli.rxns),1)];
if size(MILPproblem.csense,1) == 1 %for some model, the csense is a string instead of vector
    csense1(1:2*length(ecoli.rxns)) = 'G';
    csense = [MILPproblem.csense,csense1];
else %is a vector
    csense1(1:(2*length(ecoli.rxns)),1) = 'G';
    csense = [MILPproblem.csense; csense1];
end
MILPproblem.csense = csense;
%modify vartype if it is an MILP
if isfield(MILPproblem,'vartype')
    vartype1(1:length(ecoli.rxns),1) = 'C';
    vartype = [MILPproblem.vartype;vartype1];
    MILPproblem.vartype = vartype;
end
% Creating c (objective function)
c_old = zeros(length(MILPproblem.c),1);
c_minFlux = ones(length(ecoli.rxns),1);
c = [c_old;c_minFlux];
MILPproblem.c = c;
% set the objective sense as minimize 
MILPproblem.osense = 1;
% sovle it
solution = solveCobraMILP(MILPproblem, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', 3);
if solution.stat == 1
    Flux_min_iMAT = solution.full(1:length(ecoli.rxns));
    % off rxns 
    RHindex = find(expressionRxns >= 0.75);
    RHrxns = ecoli.rxns(RHindex);
    offRxn = RHrxns(abs(Flux_min_iMAT(RHindex)) <= 1e-8);
else
    Flux_min_iMAT = 0;
    offRxn = 0;
end


end
