function f=vaxTestPlot(NNbar,xdata,xsto,vaxparams)
%Plot cumulative cases in projected fall wave with varied delay in school
%openings
%Samples from MCMC output from selected state

mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
sig=[129.6984 329.7708 231.8866 134.0024 58.4746];
na=length(mu);
tswitchvec=(0:12);%Weeks
lt=length(tswitchvec);
%
burn=20000;
int=200;
NNtot=sum(NNbar);
prior=xsto(burn+1:int:end,:);
lp=size(prior,1);

xfrom=prctile(xsto(burn+1:end,:),[2.5],1);
xto=prctile(xsto(burn+1:end,:),[97.5],1);

rands=normrnd(repmat(mu,400,1),repmat(sig,400,1));

Z=zeros(na,lt,lp);
for i=1:lt
    tswitch=7*tswitchvec(i)+243;%Days
    Zi=nan(na,lp);
    parfor j=1:lp
        pj=prior(j,:);
        [~,~,z2]=subPandemicSimulationVax(NNbar,pj,xdata,0,0,0,tswitch,vaxparams);
        ri=rands(j,:);
        Zi(:,j)=z2'./ri'/NNtot*1000;
    end
    Z(:,i,:)=Zi;
end
ZZ=nansum(Z,1);
ZZ=(ZZ(1,1,:)-ZZ(1,11,:))./ZZ(1,1,:);
ZZ=reshape(ZZ,length(ZZ),1,1);
f=prctile(ZZ,[2.5,50,97.5]);
%
pandsq=prctile(Z,[.25,50,75],3);
y1=pandsq(:,:,1);
y2=pandsq(:,:,2);
y3=pandsq(:,:,3);
figure
fs=10; lw=2;
cmap=lines(na);
tvec=tswitchvec;
tvec2=[tvec,fliplr(tvec)];
inBetween=[y1,fliplr(y3)];
h=zeros(1,na);
hold on
for i=1:na
    Zi=reshape(Z(i,:,:),lt,lp);
    plot_distribution_prctile(tvec,Zi','color',cmap(i,:),'prctile',(0:5:50));
    h(i)=plot(-1,-1,'-','color',cmap(i,:),'linewidth',lw);
end
maxY=max(max(y3));
plot([10,10],[0,maxY],'k--','linewidth',2)
xlabel('Delay in opening school terms (weeks)','FontSize',fs);
ylabel('Cumulative hospitalisations in fall wave (per 1,000)');
set(gca,'FontSize',fs);
axis([tswitchvec(1),tswitchvec(end),0,maxY])
legend(h,{'0-4','5-17','18-49','50-64','65+'},'location','NE')
grid on
grid minor
box on
hold off