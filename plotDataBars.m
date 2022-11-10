function f=plotDataBars(cellx,cellyScaled,cellNN)
nstates=length(cellx);
A=zeros(nstates,2,5);%state/wave/age
thresh=35;
for i=1:nstates
    x=cellx{i};
    y=cellyScaled{i};
    NN=sum(cellNN{i});
    A(i,1,:)=nansum(y(x<thresh,:),1)/NN*1e3;
    A(i,2,:)=nansum(y(x>=thresh,:),1)/NN*1e3;
end

names={'CA','CO','CT','GA','MD','MN','NM','NY','OR','TN'};
plotBarStackGroups(A,names);
%Axis labes hard coded in "plotBarStackGroups"