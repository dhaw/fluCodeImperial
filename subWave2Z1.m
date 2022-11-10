function [Z2,realZ2]=subWave2Z1(NNbar,xdata,ydata,ydataNX,thetac,vaxparams,i,xsto)
%Plot data and simulation outputs for a given state
%sampling from MCMC output (xsto) from spring wave calibration

NNtot=sum(NNbar);
realZ2=nansum(ydataNX(xdata>34,:),1);
realZ2=realZ2'/NNtot*1000;
[~,z2]=subUncertaintyPlotOne(NNbar,xsto,xdata,ydata,thetac,vaxparams);
Z2=z2;
end

function [f,Z2]=subUncertaintyPlotOne(NNbar,xsto,xdata,ydata,theta,vaxparams)
thisage=0;%=0 for total, 1/2/3/4/5 for relevant age group only
burn=20000;
int=200;
plotStuff=1;
mumult=1;
cmap=lines(7);
col=[0.2422    0.1504    0.6603;
    0.2504    0.1650    0.7076;
    0.2578    0.1818    0.7511;
    0.2647    0.1978    0.7952;
    0.2706    0.2147    0.8364];
if thisage>0
    col1=col(thisage,:);
else
    col1=.5*[1,1,1];
    col2=cmap(6,:);
end
vlinecol=0*[1,1,1];

if thisage==0
    plotData=nansum(ydata,2);
else
    plotData=ydata(:,thisage);
end
lx=size(xsto,1);
xbounds=[1,35];
mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
mu=mu/mumult;
nbar=size(mu,2);
toRun=floor((lx-burn)/int);
pands=zeros(length(xdata),toRun);
Z2=zeros(toRun,nbar);
j=1;
for i=burn+1:int:lx
    try
        [pandsimi,~,zi]=subPandemicSimulationVax(NNbar,xsto(i,:),xdata,0,0,ydata,243,vaxparams);
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
pandsq=prctile(pands,[5,50,95],2);
y1=pandsq(:,1);
y2=pandsq(:,2);
y3=pandsq(:,3);
maxy=max(max(pands));
if plotStuff==1
    fs=10; lw=2;
    figure
    hold on
    plot([xbounds(2),xbounds(2)],[0,maxy],'--','linewidth',2,'color',vlinecol)
    t1=xdata(xdata<36);
    plot_distribution_prctile(t1,pands(xdata<36,:)','color',col1,'prctile',(0:5:100));
    t2=xdata(xdata>34);
    plot_distribution_prctile(t2,pands(xdata>34,:)','color',col2,'prctile',(0:5:100));
    %
    h1=plot(xdata,y2,'-','linewidth',2,'color',col1);
    h2=plot(xdata,plotData,'k-','linewidth',2,'color','k');
    %
    set(gca,'fontsize',fs)
    xlabel('MMWR week')
    ylabel('Weekly symptomatic incidence')
    axis([xdata(1),xdata(end),0,maxy])
    legend([h2,h1],'Data (adjusted with case-hospitalisation multipliers)','Model projections (aggregated over all states)','location','NW')
    grid on
    grid minor
    box on
end
f=pands;
end