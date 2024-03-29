function f=subUncertaintyPlotAll(cellNN,cellxsto,cellx,celly,vaxparams)
%Plot data and simulation outputs for all states
%sampling from MCMC output (xsto) from spring wave calibration

thisage=0;%=0 for total, 1/2/3/4/5 for relevant age group only
burn=20000;
int=200;
mumult=1;
names={'California','Colorado','Connecticut','Georgia','Maryland','Minnesota','New Mexico','New York','Oregon','Tennessee'};
cmap=lines(7);
col1=.7*[1,1,1];
col2=cmap(6,:);
vlinecol=0*[1,1,1];
lc=length(cellNN);
fs=10; lw=2;
figure;
tiledlayout(lc/2,2,'TileSpacing','compact')
for i=1:lc
    NNbar=cellNN{i};
    NNbar=NNbar';
    xsto=cellxsto{i};
    xdata=cellx{i};
    ydata=celly{i};
    plotData=nansum(ydata,2);
    lx=size(xsto,1);
    xbounds=[1,35];
    mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
    mu=mu/mumult;
    nbar=size(mu,2);
    toRun=floor((lx-burn)/int);
    pands=zeros(length(xdata),toRun);
    Z2=zeros(toRun,nbar);
    j=1;
    for ii=burn+1:int:lx
        try
            [pandsimi,~,zi]=subPandemicSimulationVax(NNbar,xsto(ii,:),xdata,0,0,ydata,243,vaxparams);
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
    maxy=max(max(max(pands)),max(plotData));
    nexttile
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
        title(names{i})
        axis([xdata(1),xdata(end),0,maxy])
        grid on
        grid minor
        box on
        if i>lc-2
            xlabel('MMWR week')
        end
        if i==5
            ylabel('Weekly symptomatic incidence');
        end
end
lg=legend(nexttile(9),[h2,h1],'Data (adjusted with case-hospitalisation multipliers)','Model projections (aggregated over all ages)');
lg.Location='southoutside';
end