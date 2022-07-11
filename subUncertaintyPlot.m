function [f,Z2]=subUncertaintyPlot(NNbar,xsto,xdata,ydata,theta,vaxparams)%,sd)%,simTheta)
thisage=0;%=0 for total
burn=20000;
int=200;
plotStuff=1;
mumult=1;
cmap=lines(7);
%
col=[0.2422    0.1504    0.6603;
    0.2504    0.1650    0.7076;
    0.2578    0.1818    0.7511;
    0.2647    0.1978    0.7952;
    0.2706    0.2147    0.8364];%lines(7);
if thisage>0
    col1=col(thisage,:);
else
    col1=.5*[1,1,1];
    col2=cmap(6,:);%.5*[1,0,0];
end
vlinecol=0*[1,1,1];%col(1,:);

if thisage==0
    plotData=nansum(ydata,2);
else
    plotData=ydata(:,thisage);
end

lx=size(xsto,1);
%xdata=1:36;
%xbounds=15;%[5,9,21,28];
thresh=15;%Match with cdcLhoodsW5
xbounds=[1,35];%[5,9,21,28];

mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
%
mu=mu/mumult;
%}
nbar=size(mu,2);

toRun=floor((lx-burn)/int);
pands=zeros(length(xdata),toRun);
Z2=zeros(toRun,nbar);
j=1;
for i=burn+1:int:lx
    %[pandsimi,~,zi]=subPandemicSimulation(NNbar,xsto(i,:),xdata,0,0,ydata,243);
    %With vaccination:
    try
        [pandsimi,~,zi]=subPandemicSimulationVax(NNbar,xsto(i,:),xdata,0,0,ydata,243,vaxparams);
    %catch ME
    %end
    %Add distribution on to "fixed" parameters:
    %{
    p12=normrnd(theta([1,2]),sd*theta([1,2]));
    [pandsimi,~,zi]=subPandemicSimulationVax(NNbar,[p12,xsto(i,:)],xdata,0,0,ydata,243,vaxparams);
    %}    
%Sum or specific age group
    if thisage==0
        pands(:,j)=sum(pandsimi,2);
    else
        pands(:,j)=pandsimi(:,thisage);
    end
    Z2(j,:)=zi;
    j=j+1;
    catch ME
    end
end
%h2mean=sum(mean(Z2,1)./mu);
%}
pandsq=prctile(pands,[5,50,95],2);
y1=pandsq(:,1);
y2=pandsq(:,2);
y3=pandsq(:,3);

maxy=max(max(pands));%max([max(y3),max(plotData)]);
%Plot a simulation:
%{
%[mle,~,~]=subPandemicSimulationVax(NNbar,theta(3:end),xdata,0,0,0,243,vaxparams);
[mle,~,~]=subPandemicSimulationVax(NNbar,theta,xdata,0,0,0,243,vaxparams);
if thisage==0
    mleSum=sum(mle,2);
else
    mleSum=mle(:,thisage);
end
maxy=max([max(y3),max(plotData),max(mleSum)]);
%}
if plotStuff==1
    fs=10; lw=2;
    figure
    hold on
    plot([xbounds(2),xbounds(2)],[0,maxy],'--','linewidth',2,'color',vlinecol)
    
    t1=xdata(xdata<36);
    plot_distribution_prctile(t1,pands(xdata<36,:)','color',col1,'prctile',(0:5:100));%Need hold off commented out in this function
    %h1=plot(-1,-1,'-','color',col1,'linewidth',lw);
    %{
    h1=plot(t1,y2(xdata<36,:),'-','linewidth',2,'color',col1);
    h2=plot(t1,plotData(xdata<36,:),'k-','linewidth',2,'color','k');%col(2,:));
    %}
    %
    t2=xdata(xdata>34);
    plot_distribution_prctile(t2,pands(xdata>34,:)','color',col2,'prctile',(0:5:100));%Need hold off commented out in this function
    %h2=plot(-1,-1,'-','color',col2,'linewidth',lw);
    
    h1=plot(xdata,y2,'-','linewidth',2,'color',col1);
    h2=plot(xdata,plotData,'k-','linewidth',2,'color','k');%col(2,:));
    %}
    set(gca,'fontsize',fs)
    xlabel('MMWR week')
    ylabel('Weekly symptomatic incidence')
    axis([xdata(1),xdata(end),0,maxy])
    %legend([h2,h1,h3],'Data','MCMC','MLE','location','NW')
    legend([h2,h1],'Data (adjusted with case-hospitalisation multipliers)','Model projections (aggregated over all states)','location','NW')
    grid on
    grid minor
    box on
end
f=pands;

%{
fs=10; lw=.5;
    %maxy=max(max(y3),max(plotData));%max(max(pands)),max(sum(ydata,2)));%,ydata(1:thresh,:)]));
    figure
    hold on
    %plot([thresh+xcut,thresh+xcut],[0,maxy],'--','linewidth',2,'color',col(1,:))
    %plot([xbounds(1),xbounds(1)],[0,maxy],'--','linewidth',2,'color',col(1,:))
    plot([xbounds(2),xbounds(2)],[0,maxy],'--','linewidth',2,'color',vlinecol)
    %{
    plot(xdata(1:thresh)+xcut,pands(1:thresh,:),'linewidth',lw,'color',[.5,.5,.5]);
    plot(xdata(thresh:end)+xcut,pands(thresh:end,:),'linewidth',lw,'color','k');
    %}
    %
    %{
    %All black:
    tvec=xdata;
    tvec2=[tvec;flipud(tvec)];
    inBetween=[y1;flipud(y3)];
    fill(tvec2,inBetween,col1,'facealpha',.2);
    plot(tvec,y1,'-','linewidth',1,'color',col1);
    plot(tvec,y3,'-','linewidth',1,'color',col1);
    h1=plot(tvec,y2,'-','linewidth',2,'color',col1);
    h2=plot(xdata,plotData,'k-','linewidth',2,'color','k');%col(2,:));
    %}
    %
    %Black and red:
    tvec=xdata;
    t1=xdata(xdata<36);
    t1b=[t1;flipud(t1)];
    ib1=[y1(xdata<36);flipud(y3(xdata<36))];
    fill(t1b,ib1,col1,'facealpha',.2);
    plot(t1,y1(xdata<36),'-','linewidth',1,'color',col1);
    plot(t1,y3(xdata<36),'-','linewidth',1,'color',col1);
    %
    t2=xdata(xdata>34);
    t2b=[t2;flipud(t2)];
    ib2=[y1(xdata>34);flipud(y3(xdata>34))];
    fill(t2b,ib2,col2,'facealpha',.2);
    plot(t2,y1(xdata>34),'-','linewidth',1,'color',col2);
    plot(t2,y3(xdata>34),'-','linewidth',1,'color',col2);
    h1=plot(tvec,y2,'-','linewidth',2,'color',col1);
    h2=plot(xdata,plotData,'k-','linewidth',2,'color','k');%col(2,:));
    %}
    %
    %plot((thresh:xdata(end))+xcut,plotData(thresh:end),'k--','linewidth',2,'color',[0,0,0])%col(2,:))
    %}
    %h3=plot(xdata,mleSum,'--','linewidth',2,'color',.5*[1,1,1]);%col(3,:));
%}