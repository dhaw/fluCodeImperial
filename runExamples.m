function f=runExamples
%This code will generate all figures from the paper

%All data files requird can be downloaded here:
%https://github.com/dhaw/fluCodeImperial 

load('vaxparams.mat','vaxparams')
load('cellNN.mat','cellNN','NNsums','NNbarTotal')
%NNsums: state totals; NNbarTotal: age-groups summed by states
%load('cellOtherParams.mat','cellOtherParams')
load('cellx.mat','cellx')
load('cellxsto.mat','cellxsto')
load('celly.mat','celly');
load('cellyScaled.mat','cellyScaled');
load('cellZ2.mat','cellZ2','realZ2');

%% Notes
%Some codes inclde flags for various model variants that were used in
%devising the study, or for other purposes, but not included in the manuscript

%% Load/select objects
%Select state for single-state plots:
state=1;%1=CA, 2=CO, 3=CT, 4=GA, 5=MD, 6=MN, 7=NM, 8=NY, 9=OR, 10=TN
names={'California','Colorado','Connecticut','Georgia','Maryland','Minnesota','New Mexico','New York','Oregon','Tennessee'};
stateName=names{state};
%
NNbar=cellNN{state}';
xdata=cellx{state};
ydata=celly{state};
ydataScaled=cellyScaled{state};
xsto=cellxsto{state};
params=xsto(end,:);
%nstates=length(cellx);
%nages=size(celly{1},2);

%% List of input parameters ("params"/rows of mcmc output "xsto"): %xxToDo - finish list/naming
%1 - 
%2 - 
%3-27 - age-mixing matrix coefficients (listed vertically)
%28 - proportion immune at start (of those born pre-1957)
%29 - proportion of cases asymptomatic
%30 - relative infectiousness of asymptomatic cases
%31 - t_0 (seed time
%32 - R0
%33 - gamma=1/TR (recovery rate)

tswitch=243;
plotComp=1;%Plots data and simulation aggregated by week (comparison)
plotEpis=0;%Plots ODE output
%Single simulation (final parameter set in MCMC output):
[f,g,z2]=subPandemicSimulationVax(NNbar,params,xdata,plotComp,plotEpis,ydata,tswitch,vaxparams);

%% FIGURE 1 - plot data:
plotDataAllStates(cellNN,cellx,cellyScaled);

%% FIGURE 2 - trajectories with uncertainty from MCMC (selected state):
subWave2Z1(NNbar,xdata,ydata,ydataScaled,0,vaxparams,0,xsto);

%% FIGURE 3 - fall-wave cumulative cases (Z2) for all states:
plotz2(cellZ2,realZ2,NNsums)%xxToDo - check

%% FIGURE 4 - effect of delay in school openings (selected state):
vaxTestPlot(NNbar,xdata,xsto,vaxparams);

%% FIGURE 5 - decision framework (selected state):
subSeverity(NNbar,xsto,xdata,ydata,vaxparams,stateName);

%% FIGURE S1 - first and second wave cumulative cases for all states:
plotDataBars(cellx,cellyScaled,cellNN);

%% FIGURE S2 - plots as figure 2 for all states in a single figure:
subUncertaintyPlotAll(cellNN,cellxsto,cellx,celly,vaxparams);

%% FIGURE S3 - MCMC traces for parameters a and b (see above for parameter indices):
a=32; b=33;
plotMCMCproj(xsto,a,b);

%% Figure S4
plotVax09(vaxparams,NNbarTotal);

%% Run MCMC:
[xsto_new, outsto, history, accept_rate,covmat]=subMCMC(NNbar,xdata,ydataScaled,xsto(1,:),vaxparams);
%Note a subtle change in threshold for each in the likelihood function
%subLhood.m was used. It is currently set of CA and is commented should
%the user wish to experiment

%% Comment in to make cellZ2/realZ2 (saved objects):
%cellZ2 is a cell array of simulated secondary attack rates (fall wave)
%realZ2 is an array of secondary attack rates in data
%{
nstates=length(cellx);
nages=size(celly{1},2);
cellZ2=cell(nstates,1);
realZ2=zeros(nages,nstates);
NNsums=zeros(nstates,1);
for i=1:length(cellx)
    [Z2i,realZ2i]=subWave2Z1(cellNN{i}',cellx{i},celly{i},cellyScaled{i},0,vaxparams,0,cellxsto{i});%xxToDo - add flag for no plot (keep on above)
    cellZ2{i}=Z2i;
    realZ2(:,i)=realZ2i;
    NNsums(i)=sum(cellNN{i});
end
%}
%% Comment in to cable of mean/95% CB values for marginals
%{
colPerState=4;
lc=length(cellxsto);
[a,b]=size(cellxsto{1});
C=[1.9200    0.4268    0.5260    0.2554    0.1665;
    1.7600    8.7522    2.2855    1.0876    1.2190;
    4.0700    4.5939    6.6160    4.5939    2.9494;
    0.9000    0.8885    1.6180    2.3847    1.6919;
    0.2300    0.2975    0.5712    0.8756    1.8930]';
C=reshape(C,25,1);
clim=[.8*C,1.2*C];
%plim=[.5,0;90,-30;.5,0;.7,.4;.6,.4;200,0;3,1;2,.001];
%plim=[plim(:,2),plim(:,1)];
plim=[0,.3;-31,31;0,.4;.495,.605;.4,.6;-100,100;1.305,1.595;1/4,1/1.5];
%plim=[plim(1:2,:);clim;plim(3:end,:)];
table=zeros(b,colPerState*lc);
for i=1:lc
    xstoi=cellxsto{i};
    ind=(i-1)*colPerState+1;
    table(:,ind)=mean(xsto,1);
    %table(:,ind+1:ind+2)=prctile(xstoi,[2.5,97.5],1)';
    table(:,ind+1:ind+3)=prctile(xstoi,[50,2.5,97.5],1)';
end
plim=[plim;clim];
table=[plim,[table([1,2,28:33],:);table(3:27,:)]];%Just numbers for latex table
writematrix(round(table,2),'mcmcSummary.csv')
%}