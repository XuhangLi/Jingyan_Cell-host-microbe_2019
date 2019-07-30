addpath /opt/cobra/cobratoolbox/
initCobraToolbox;
load('iJO1366.mat');
ecoli = iJO1366;
%add sink
[num,txt,all] = xlsread('sinks_to_fix_rxn.xlsx','met_list');
sinksRxns = all(2:end,1);
ecoli = addSinkReactions(ecoli, sinksRxns);
%% LB constraints
%add constrains; biomass unconstrained(but core is forbid); sinks and non_lb_exchange are 0.1
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
[epsilon_f,epsilon_r] = makeEpsilonSeq(ecoli,ecoli.rxns,0.05,0.5);
%%
parpool(20);
pctRunOnAll initCobraToolbox;
Flux_min_iMAT_rand = zeros(length(ecoli.rxns),10100);
offRxnNum = zeros(10100,1);
HRxnCandidatesNum = zeros(10100,1);
ROSflux_percentatage = zeros(10100,1);
ROSrxnsPercentage = zeros(10100,1);
ETCflux_percentatage = zeros(10100,1);
ETCrxnsPercentage = zeros(10100,1);
Folateflux_percentatage = zeros(10100,1);
FolateRxnPercentage = zeros(10100,1);
NaNindx = zeros(10100,1);
%randomization
randseq = 1:10100;  %perform more than 10000 in case of infeasible
parfor j = 1:10100
   rng(randseq(j));
   ind1 = randperm(length(ecoli.genes))';
   indx = ind1(1:length(GeneNameExpressed0),1);
   GeneNameExpressed_permuted = ecoli.genes(indx);
   GeneNameExpressed = GeneNameExpressed_permuted;
   %% generate the iMAT flux
   [Flux_min_iMAT,solution,offRxn, HRxnCandidates] = XL_generate_ROS_flux(ecoli,GeneNameExpressed,epsilon_f,epsilon_r);
   if solution.stat == 1 %skip numeric problem
       Flux_min_iMAT_rand(:,j) = Flux_min_iMAT;
       offRxnNum(j) = length(offRxn);
       HRxnCandidatesNum(j) = length(HRxnCandidates);
        %% 1. exam ROS reactions
        onRxns = ecoli.rxns(abs(Flux_min_iMAT) > 1e-7);
        ROSon = intersect(ROSrxns,onRxns);
        ROSflux = sum(abs(Flux_min_iMAT(ismember(ecoli.rxns,ROSrxns))))
        ROSflux_percentatage(j) = ROSflux / sum(abs(Flux_min_iMAT))
        ROSrxnsPercentage(j) = length(ROSon) / length(ROSrxns)
        %% 2. exam ETC rxns
        ETCrxns = ecoli.rxns(strcmp(ecoli.subSystems,'Oxidative Phosphorylation'));
        ETCflux = sum(abs(Flux_min_iMAT(strcmp(ecoli.subSystems,'Oxidative Phosphorylation'))))
        ETCflux_percentatage(j) = ETCflux / sum(abs(Flux_min_iMAT))
        ETCrxnsPercentage(j) = length(intersect(onRxns,ETCrxns)) / length(ETCrxns)
        %% 3. exam folate 
        folateRxns = ecoli.rxns(strcmp(ecoli.subSystems,'Folate Metabolism'));
        folateFlux = sum(abs(Flux_min_iMAT(strcmp(ecoli.subSystems,'Folate Metabolism'))))
        Folateflux_percentatage(j) = folateFlux / sum(abs(Flux_min_iMAT))
        FolateRxnPercentage(j) = length(intersect(onRxns,folateRxns)) / length(folateRxns) 
   else
        NaNindx(j) = 1;
   end
end
    save('randGeneFlux.mat','Flux_min_iMAT_rand');
    save('offRxnNum.mat','offRxnNum');
    save('HRxnCandidatesNum.mat','HRxnCandidatesNum');
    save('ROSflux_percentatage.mat','ROSflux_percentatage');
    save('ROSrxnsPercentage.mat','ROSrxnsPercentage');
    save('ETCflux_percentatage.mat','ETCflux_percentatage');
    save('ETCrxnsPercentage.mat','ETCrxnsPercentage');
    save('Folateflux_percentatage.mat','Folateflux_percentatage');
    save('FolateRxnPercentage.mat','FolateRxnPercentage');
    save('NaNindx.mat','NaNindx');
