addpath /opt/cobra/cobratoolbox/
initCobraToolbox;
load('iJO1366.mat');
ecoli = iJO1366;
%add sink
[num,txt,all] = xlsread('sinks_to_fix_rxn.xlsx','met_list');
sinksRxns = all(2:end,1);
ecoli = addSinkReactions(ecoli, sinksRxns);
%% LB constraints
%add constrains; biomass unconstrained; sinks and non_lb_exchange are 0.1
ecoli = xlsConst2model('./open_exchange.xlsx', ecoli);
ecoli = xlsConst2model('DeletedGenes.xlsx',ecoli);
%%
%generate GPR
ecoli = changeObjective (ecoli, 'BIOMASS_Ec_iJO1366_core_53p95M'); 
ecoli = generateRules(ecoli);
ecoli = buildRxnGeneMat(ecoli);   

%so that generate the three set
[num,txt,all] = xlsread('Yans_gene_in_model.xlsx','geneYans_all');
wormSize_all = num;
GeneNamesFull = all(2:end,1);
[num,txt,all] = xlsread('Yans_gene_in_model.xlsx','geneYans_hit');
wormSize_hit = num;
GeneNameExpressed = all(2:end,1);

%in this trial, all the genes are opened regardness of OR and AND gate
%also, use dynamic epsilon to open all hit genes
[num,txt,all] = xlsread('Yans_gene_in_model.xlsx','ROS_rxns');
ROSrxns = all(1:end,1);
ROSrxns(end+1,1) = {'wormSize' };
GeneNameExpressed0 = GeneNameExpressed;
%% generate the iMAT flux
[epsilon_f,epsilon_r] = makeEpsilonSeq(ecoli,ecoli.rxns,0.05,0.5);
[Flux_min_iMAT,solution,offRxn, HRxnCandidates] = XL_generate_ROS_flux(ecoli,GeneNameExpressed,epsilon_f,epsilon_r);
%% 1. exam ROS reactions
onRxns = ecoli.rxns(abs(Flux_min_iMAT) > 1e-7);
ROSon = intersect(ROSrxns,onRxns);
ROSon(:,2) = ecoli.subSystems(ismember(ecoli.rxns,ROSon));
ROSflux = sum(abs(Flux_min_iMAT(ismember(ecoli.rxns,ROSrxns))))
ROSflux_percentatage0 = ROSflux / sum(abs(Flux_min_iMAT))
ROSrxnsPercentage0 = length(intersect(onRxns,ROSrxns)) / length(ROSrxns)
%% 2. exam ETC rxns
ETCrxns = ecoli.rxns(strcmp(ecoli.subSystems,'Oxidative Phosphorylation'));
ETCflux = sum(abs(Flux_min_iMAT(strcmp(ecoli.subSystems,'Oxidative Phosphorylation'))))
ETCflux_percentatage0 = ETCflux / sum(abs(Flux_min_iMAT))
ETCrxnsPercentage0 = length(intersect(onRxns,ETCrxns)) / length(ETCrxns)
%% 3. exam folate 
folateRxns = ecoli.rxns(strcmp(ecoli.subSystems,'Folate Metabolism'));
folateFlux = sum(abs(Flux_min_iMAT(strcmp(ecoli.subSystems,'Folate Metabolism'))))
Folateflux_percentatage0 = folateFlux / sum(abs(Flux_min_iMAT))
FolateRxnPercentage0 = length(intersect(onRxns,folateRxns)) / length(folateRxns)
