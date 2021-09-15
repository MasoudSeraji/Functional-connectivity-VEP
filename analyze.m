% Script to check the OSF shared data
% - extract latency from preprocessed data
% - summarize latency differences and statistical analysis
% - summarize connectivity differences and statistical analysis
% - Use an SVN to assess diagnostic power of Latency,L, CC
%
% Needs
%  ../data/preprocessedData.mat   (preprocessed table with preprocessed EEG
%                                   signal matrices)
% ../data/results.mat             (results table, with latency, L, CC per subject)
%
% Dependencies:
%   Uses linear model utilities of the BayesFactor toolbox avaiable here:
%   https://github.com/klabhub/bayesFactor
%
% BK/MS -2021
%% Load pre-processed data
clear all;
load ../data/preprocessedData.mat

%%  Estimate latency from the preprocessed data.
N = height(preprocessed);
pathways ={'Magno','Parvo','Konio'};
time= 1:512;
nrSamples= numel(time);
Latency = nan(N*3,1);
Subject = cell(N*3,1);
Status = cell(N*3,1);
Pathway = cell(N*3,1);
mVep = cell(N*3,1);
SR=512;
% T=T(101:end);
% Loop over subjects and pathways
for i=1:N
    if strcmpi(preprocessed.Status{i},'HC')
        window = [(75*SR)/1000 (120*SR)/1000];
    else
        window = [(120*512)/1000 (160*512)/1000];
    end
    for p = 1:3
        thisRow  = (i-1)*3+p;
        thisVEP = permute(preprocessed.(pathways{p}){i},[2 1 3]); % Put Time dimension first
        thisVEP = nanmean(nanmean(thisVEP,3),2); % Average over electrodes and trials
        thisVEP=thisVEP(101:end);
        inWindow = time>=window(1) & time <=window(2);
        [~,Latency(thisRow)] =  max(thisVEP(inWindow));
        Latency(thisRow) = ((Latency(thisRow)+window(1)-1)*1000)/512;
        Subject{thisRow}  = preprocessed.Subject{i};
        Status{thisRow} = preprocessed.Status{i};
        Pathway{thisRow} = lower(pathways{p});
    end
end
% This is the Latency part of the results table; for the summary below it
% gives the same results as the full results table on OSF.
results = table(Subject,Status,Latency,Pathway);

%% Statistical analysis of the results table shared on OSF.
load ../data/results.mat  % Load OSF stored results table
results.Subject = categorical(results.Subject);
% Latency Results
% Does Latency change with disease status in a pathway specific manner?
summary = groupsummary(results,{'Status','Pathway'},{'mean','std'},'Latency')
isHC = strcmpi(summary.Status,'HC');
g =hedgesG(summary.mean_Latency(~isHC),summary.mean_Latency(isHC),summary.std_Latency(~isHC),summary.std_Latency(isHC),summary.GroupCount(~isHC),summary.GroupCount(isHC));
pathway = summary.Pathway(isHC);
latencyG  = table(pathway,g)
% Statistics
lmm = fitlme(results,'Latency~Status*Pathway + (1|Subject)','DummyVarCoding','Effects');

anova(lmm)
R = residuals(lmm);
normal = makedist('norm','mu',0,'sigma',nanstd(R));
[notNormal,p] = kstest(R,'cdf',normal)

if isempty(which('lm.contrast'))
    warning('Skipping contrasts calculations that depend on the bayesFactor toolbox');
else
    % pairwise linear contrasts: is magno effect larger than parvo or konio for MS?
    msMagno = lm.contrast(lmm,{'Status','MS','Pathway','magno'},{'Status','HC','Pathway','magno'});
    msParvo = lm.contrast(lmm,{'Status','MS','Pathway','parvo'},{'Status','HC','Pathway','parvo'});
    msKonio = lm.contrast(lmm,{'Status','MS','Pathway','konio'},{'Status','HC','Pathway','konio'});
    
    [HLat.mIsP,~,~,~,~,strLatMP] = lm.posthoc(lmm,msMagno,msParvo);
    [HLat.mIsK,~,~,~,~,strLatMK] = lm.posthoc(lmm,msMagno,msKonio);
    [HLat.pIsK,~,~,~,~,strLatKP] = lm.posthoc(lmm,msParvo,msKonio);
    [HLat.pM,~,~,~,~,strLatM] = lm.posthoc(lmm,msMagno);
    [HLat.pP,~,~,~,~,strLatP] = lm.posthoc(lmm,msParvo);
    [HLat.pK,~,~,~,~,strLatK] = lm.posthoc(lmm,msKonio);
    HLat
end
%%  Synchrony stats
% Does L or CC change with disease status in a pathway specific manner?
summary = groupsummary(results,{'Status','Pathway'},{'mean','std'},{'L','CC'})
isHC = strcmpi(summary.Status,'HC');
g =hedgesG(summary.mean_L(~isHC),summary.mean_L(isHC),summary.std_L(~isHC),summary.std_L(isHC),summary.GroupCount(~isHC),summary.GroupCount(isHC));
pathway = summary.Pathway(isHC);
LG  = table(pathway,g)

lmmL = fitlme(results,'L~Status*Pathway + (1|Subject)','DummyVarCoding','Effects');
anova(lmmL)
if isempty(which('lm.contrast'))
    warning('Skipping contrasts calculations that depend on the bayesFactor toolbox');
else
    % pairwise linear contrasts: is magno effect larger than parvo or konio for MS?
    
    msMagno = lm.contrast(lmmL,{'Status','MS','Pathway','magno'},{'Status','HC','Pathway','magno'});
    msParvo = lm.contrast(lmmL,{'Status','MS','Pathway','parvo'},{'Status','HC','Pathway','parvo'});
    msKonio = lm.contrast(lmmL,{'Status','MS','Pathway','konio'},{'Status','HC','Pathway','konio'});
    
    [HL.mIsP,~,~,~,~,strLMP] = lm.posthoc(lmmL,msMagno,msParvo);
    [HL.mIsK,~,~,~,~,strLMK] = lm.posthoc(lmmL,msMagno,msKonio);
    [HL.pIsK,~,~,~,~,strLKP] = lm.posthoc(lmmL,msParvo,msKonio);
    [HL.pM,~,~,~,~,strLM] = lm.posthoc(lmmL,msMagno);
    [HL.pP,~,~,~,~,strLP] = lm.posthoc(lmmL,msParvo);
    [HL.pK,~,~,~,~,strLK] = lm.posthoc(lmmL,msKonio);
    HL
end

g =hedgesG(summary.mean_CC(~isHC),summary.mean_CC(isHC),summary.std_CC(~isHC),summary.std_CC(isHC),summary.GroupCount(~isHC),summary.GroupCount(isHC));
pathway = summary.Pathway(isHC);
CCG  = table(pathway,g)

lmmC = fitlme(results,'CC~Status*Pathway + (1|Subject)','DummyVarCoding','Effects');
if isempty(which('lm.contrast'))
    warning('Skipping contrasts calculations that depend on the bayesFactor toolbox');
else
    % pairwise linear contrasts: is magno effect larger than parvo or konio for MS?
    msMagno = lm.contrast(lmmC,{'Status','MS','Pathway','magno'},{'Status','HC','Pathway','magno'});
    msParvo = lm.contrast(lmmC,{'Status','MS','Pathway','parvo'},{'Status','HC','Pathway','parvo'});
    msKonio = lm.contrast(lmmC,{'Status','MS','Pathway','konio'},{'Status','HC','Pathway','konio'});
    
    [HC.mIsP,~,~,~,~,strCMP,cMP] = lm.posthoc(lmmC,msMagno,msParvo);
    [HC.mIsK,~,~,~,~,strCMK,cMK] = lm.posthoc(lmmC,msMagno,msKonio);
    [HC.pIsK,~,~,~,~,strCKP,cKP] = lm.posthoc(lmmC,msParvo,msKonio);
    [HC.pM,~,~,~,~,strCM] = lm.posthoc(lmmC,msMagno);
    [HC.pP,~,~,~,~,strCP] = lm.posthoc(lmmC,msParvo);
    [HC.pK,~,~,~,~,strCK] = lm.posthoc(lmmC,msKonio);
    HC
end

anova(lmmC)

%% Classification
% Can we predict disease status from Latency, L, and CC (and combinations
% thereof?)
T = unstack(results,{'Latency','L','CC'},'Pathway');
models = {'Status~Latency_magno',...
    'Status~L_magno'....
    'Status~CC_magno',...
    'Status~L_magno+CC_magno',...
    'Status~L_magno+CC_magno+L_parvo+CC_parvo+L_konio+CC_konio',...
    'Status~Latency_magno+Latency_parvo+Latency_konio'};
summary = table;  % Construct summary
for m =1:numel(models)
    % Fit a SVM with linear kernel
    model = models{m};
    svm = fitcsvm(T,model,'Standardize',true,'KernelFunction','linear','KernelScale','auto');
    cv = crossval(svm,'Leaveout','on');
    
    predicted = kfoldPredict(cv);
    pctMiss = mean(~strcmpi(predicted,svm.Y));
    sensitivity = sum(strcmpi(svm.Y,'MS')&strcmpi(predicted,'MS'))/sum(strcmpi(svm.Y,'MS')); % True positives
    specificity = sum(strcmpi(svm.Y,'HC')&strcmpi(predicted,'HC'))/sum(strcmpi(svm.Y,'HC')); % True negatives
    
    summary = [summary; table({model},pctMiss,sensitivity,specificity,'VariableNames',{'model','pctMiss','sensitivity','specifity'})]; %#ok<AGROW>
end
summary =sortrows(summary,'pctMiss')
