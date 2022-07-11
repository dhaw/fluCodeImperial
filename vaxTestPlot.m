function f=vaxTestPlot(NNbar,xdata,xsto,vaxparams)
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

xfrom=prctile(xsto(burn+1:end,:),[2.5],1);%(burn+1:int:end,:)
xto=prctile(xsto(burn+1:end,:),[97.5],1);%(burn+1:int:end,:)

%rands=normrnd(1,0,lp,na);
%rands=repmat(rands,1,na).*repmat(sig,lp,1)+repmat(mu,lp,1);
rands=normrnd(repmat(mu,400,1),repmat(sig,400,1));

Z=zeros(na,lt,lp);
for i=1:lt
    tswitch=7*tswitchvec(i)+243;%Days
    Zi=nan(na,lp);
    parfor j=1:lp
        pj=prior(j,:);
        %if min(pj-xfrom)>0 && min(xto-pj)>0
        [~,~,z2]=subPandemicSimulationVax(NNbar,pj,xdata,0,0,0,tswitch,vaxparams);
        %Zi(:,j)=z2'./mu'/NNtot*1000;
        ri=rands(j,:);%(normrnd(mu,sig);
        Zi(:,j)=z2'./ri'/NNtot*1000;
        %end
    end
    Z(:,i,:)=Zi;
end
%}
%
%Z(Z<0)=0;
ZZ=nansum(Z,1);
ZZ=(ZZ(1,1,:)-ZZ(1,11,:))./ZZ(1,1,:);
ZZ=reshape(ZZ,length(ZZ),1,1);
f=prctile(ZZ,[2.5,50,97.5]);%[min(ZZ),max(ZZ)];
%}

pandsq=prctile(Z,[.25,50,75],3);%25,50,75
y1=pandsq(:,:,1);
y2=pandsq(:,:,2);
y3=pandsq(:,:,3);
figure
fs=10; lw=2;
cmap=lines(na);
tvec=tswitchvec;
tvec2=[tvec,fliplr(tvec)];
inBetween=[y1,fliplr(y3)];
%plot(tswitchvec,Z,'-','linewidth',lw)%,'color',col1);
h=zeros(1,na);
hold on
for i=1:na
    %{
    fill(tvec2,inBetween(i,:),cmap(i,:),'facealpha',.2);
    plot(tvec,y1(i,:),'-','linewidth',1,'color',cmap(i,:));
    plot(tvec,y3(i,:),'-','linewidth',1,'color',cmap(i,:));
    h(i)=plot(tvec,y2(i,:),'-','linewidth',2,'color',cmap(i,:));
    %}
    Zi=reshape(Z(i,:,:),lt,lp);
    plot_distribution_prctile(tvec,Zi','color',cmap(i,:),'prctile',(0:5:50));%Need hold off commented out in this function
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
%}