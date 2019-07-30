load NaNindx.mat
load offRxnNum.mat
offRxnNum(NaNindx==1)=[];
offRxnNum = offRxnNum(1:10000);
load HRxnCandidatesNum.mat
HRxnCandidatesNum(NaNindx==1)=[];
HRxnCandidatesNum = HRxnCandidatesNum(1:10000);
load ROSflux_percentatage.mat
ROSflux_percentatage(NaNindx==1)=[];
ROSflux_percentatage = ROSflux_percentatage(1:10000);
load ROSrxnsPercentage.mat
ROSrxnsPercentage(NaNindx==1)=[];
ROSrxnsPercentage = ROSrxnsPercentage(1:10000);
load ETCflux_percentatage.mat
ETCflux_percentatage(NaNindx==1)=[];
ETCflux_percentatage = ETCflux_percentatage(1:10000);
load ETCrxnsPercentage.mat
ETCrxnsPercentage(NaNindx==1)=[];
ETCrxnsPercentage = ETCrxnsPercentage(1:10000);
load Folateflux_percentatage.mat
Folateflux_percentatage(NaNindx==1)=[];
Folateflux_percentatage = Folateflux_percentatage(1:10000);
load FolateRxnPercentage.mat
FolateRxnPercentage(NaNindx==1)=[];
FolateRxnPercentage = FolateRxnPercentage(1:10000);
load randGeneFlux.mat
Flux_min_iMAT_rand(:,NaNindx==1)=[];
Flux_min_iMAT_rand = Flux_min_iMAT_rand(:,1:10000);
%% QC of the randomization
figure(1)
hist(offRxnNum);
figure(2)
hist(HRxnCandidatesNum);
figure(3)
hist(offRxnNum ./ HRxnCandidatesNum);
%% p value of ROS flux; note the ROSflux_percentatage0 is calculated in step#2
figure(1)
hold on
histogram(ROSflux_percentatage);
line([ROSflux_percentatage0, ROSflux_percentatage0], ylim, 'LineWidth', 2, 'Color', 'r');
%figure(2)
%hist(ROSrxnsPercentage);
p1 = (1+ sum(ROSflux_percentatage >= ROSflux_percentatage0)) / 10001
p2 = (1+ sum(ROSrxnsPercentage >= ROSrxnsPercentage0)) / 10001

%% same for oxphos
figure(2)
hold on
histogram(ETCflux_percentatage);
line([ETCflux_percentatage0, ETCflux_percentatage0], ylim, 'LineWidth', 2, 'Color', 'r');
%figure(2)
%hist(ETCrxnsPercentage);
p1 = (1+ sum(ETCflux_percentatage >= ETCflux_percentatage0)) / 10001
p2 = (1+ sum(ETCrxnsPercentage >= ETCrxnsPercentage0)) / 10001
%% same for ETC
randGeneFlux = Flux_min_iMAT_rand;
[num txt all] = xlsread('ETC_set.xlsx');
redoxRxnList = ismember(ecoli.rxns,all(:,1));
for i = 1:size(randGeneFlux,2)
    myFlux_min_iMAT = randGeneFlux(:,i);
    myonRxns = ecoli.rxns(abs(myFlux_min_iMAT) > 1e-7);
    myredoxFlux = sum(abs(myFlux_min_iMAT(redoxRxnList)));
    myredoxflux_percentatage(i) = myredoxFlux / sum(abs(myFlux_min_iMAT));
    myredoxrxnsPercentage(i) = length(intersect(myonRxns,ecoli.rxns(redoxRxnList))) / length(ecoli.rxns(redoxRxnList));
end
myperc = sum(abs(Flux_min_iMAT(redoxRxnList)))/sum(abs(Flux_min_iMAT));
figure(3)
histogram(myredoxflux_percentatage);
line([myperc, myperc], ylim, 'LineWidth', 2, 'Color', 'r');
%figure(2)
%hist(myredoxrxnsPercentage);
p1 = (1+ sum(myredoxflux_percentatage >= sum(abs(Flux_min_iMAT(redoxRxnList)))/sum(abs(Flux_min_iMAT)))) / 10001
p2 = (1+ sum(myredoxrxnsPercentage >= sum(abs(Flux_min_iMAT(redoxRxnList))>1e-7)/sum(redoxRxnList))) / 10001

%% same for folate metabolism
hold on
figure(4)
histogram(Folateflux_percentatage);
line([Folateflux_percentatage0, Folateflux_percentatage0], ylim, 'LineWidth', 2, 'Color', 'r');

%figure(2)
%hist(FolateRxnPercentage);
p1 = (1+ sum(Folateflux_percentatage >= Folateflux_percentatage0)) / 10001
p2 = (1+ sum(FolateRxnPercentage >= FolateRxnPercentage0)) / 10001
%% same for all the subsystem
[num txt all] = xlsread('subsystems.xlsx');
subsystems = all(2:end,1);
randGeneFlux = Flux_min_iMAT_rand;
for zz = 1:length(subsystems)
    mysubsys = subsystems(zz);
    redoxRxnList = strcmp(ecoli.subSystems,mysubsys);
    for i = 1:size(randGeneFlux,2)
        myFlux_min_iMAT = randGeneFlux(:,i);
        myonRxns = ecoli.rxns(abs(myFlux_min_iMAT) > 1e-7);%numeric cutoff for valid flux
        myredoxFlux = sum(abs(myFlux_min_iMAT(redoxRxnList)));
        myredoxflux_percentatage(i) = myredoxFlux / sum(abs(myFlux_min_iMAT));
        myredoxrxnsPercentage(i) = length(intersect(myonRxns,ecoli.rxns(redoxRxnList))) / length(ecoli.rxns(redoxRxnList));
    end
    %hist(myredoxflux_percentatage);
    %hist(myredoxrxnsPercentage);
    myonRxns0 = ecoli.rxns(abs(Flux_min_iMAT) > 1e-7);
    p1(zz) = (1+ sum(myredoxflux_percentatage >= (sum(abs(Flux_min_iMAT(redoxRxnList)))/sum(abs(Flux_min_iMAT))))) / 10001;
    p2(zz) = (1+ sum(myredoxrxnsPercentage >= (length(intersect(myonRxns0,ecoli.rxns(redoxRxnList))) / length(ecoli.rxns(redoxRxnList))))) / 10001;
    zz
end

%% same for all redox 
randGeneFlux = Flux_min_iMAT_rand;
redoxList = {'nad_c','nadh_c','fad_c','fadh2_c','nadp_c','nadph_c','gthox_c','gthrd_c','mqn8_c','mql8_c','q8_c','q8h2_c','2dmq8','2dmql8'};
redoxRxnList = any(ecoli.S(ismember(ecoli.mets,redoxList),:),1);
for i = 1:size(randGeneFlux,2)
    myFlux_min_iMAT = randGeneFlux(:,i);
    myonRxns = ecoli.rxns(abs(myFlux_min_iMAT) > 1e-7);
    myredoxFlux = sum(abs(myFlux_min_iMAT(redoxRxnList)));
    myredoxflux_percentatage(i) = myredoxFlux / sum(abs(myFlux_min_iMAT));
    myredoxrxnsPercentage(i) = length(intersect(myonRxns,ecoli.rxns(redoxRxnList))) / length(ecoli.rxns(redoxRxnList));
end
%figure(1)
%hist(myredoxflux_percentatage);
%figure(2)
%hist(myredoxrxnsPercentage);
p1 = (1+ sum(myredoxflux_percentatage >= redoxflux_percentatage0)) / 10001
p2 = (1+ sum(myredoxrxnsPercentage >= redoxrxnsPercentage0)) / 10001

%% negative control - enrichment in random rxn set 
NrandRxn = 9;
indx = randperm(length(ecoli.rxns));
randRxnList = zeros(length(ecoli.rxns),1);
randRxnList(indx(1:NrandRxn)) = 1;
randRxnList = any(randRxnList,2);

myonRxns0 = ecoli.rxns(abs(Flux_min_iMAT) > 1e-7);
myrandFlux0 = sum(abs(Flux_min_iMAT(randRxnList)));
myrandflux_percentatage0 = myrandFlux0 / sum(abs(Flux_min_iMAT));
myrandrxnsPercentage0 = length(intersect(myonRxns0,ecoli.rxns(randRxnList))) / length(ecoli.rxns(randRxnList));
for i = 1:size(randGeneFlux,2)
    myFlux_min_iMAT = randGeneFlux(:,i);
    myonRxns = ecoli.rxns(abs(myFlux_min_iMAT) > 1e-7);
    myrandFlux = sum(abs(myFlux_min_iMAT(randRxnList)));
    myrandflux_percentatage(i) = myrandFlux / sum(abs(myFlux_min_iMAT));
    myrandrxnsPercentage(i) = length(intersect(myonRxns,ecoli.rxns(randRxnList))) / length(ecoli.rxns(randRxnList));
end
%figure(1)
%hist(myrandflux_percentatage);
%figure(2)
%hist(myrandrxnsPercentage);
p1 = (1+ sum(myrandflux_percentatage >= myrandflux_percentatage0)) / 10001
p2 = (1+ sum(myrandrxnsPercentage >= myrandrxnsPercentage0)) / 10001

%% calculate p-value for each reaction
% analysis the nulll distribution of "on" reaction flux (percentatge of total flux)
totalflux = sum(abs(randGeneFlux),1);
relativerandGeneFlux = abs(randGeneFlux) ./ totalflux .* 100; %as percentage

fluxAve = mean(relativerandGeneFlux,2);
fluxStd = std(relativerandGeneFlux,0,2);
[fluxAve_sorted I] = sort(fluxAve);
fluxStd_sorted = fluxStd(I);
figure
hold on
bar(log10(abs(fluxAve_sorted)))
errorbar(log10(abs(fluxAve_sorted)),log10(abs(fluxStd_sorted)) ,'.')

totalFlux0 = sum(abs(Flux_min_iMAT));
for i = 1: length(Flux_min_iMAT)
    pValue(i) = ttest(relativerandGeneFlux(i,:),abs(Flux_min_iMAT(i))/totalFlux0);
end
bar(sort(pValue)); 
ecoli.rxns(pValue<0.05)
%significant rxns:DHORDfum,4PCPpp,MALt3pp,sink_dgslnt_c
%% report
for i = 1: length(Flux_min_iMAT)
    pValue(i) = (1+sum(relativerandGeneFlux(i,:) >= 100*abs(Flux_min_iMAT(i))/totalFlux0)) / 10001;
end
SigRxn = ecoli.rxns(pValue<0.05)  
pValue1 = pValue(ismember(ecoli.rxns,SigRxn))