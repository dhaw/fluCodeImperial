function f=plotVax09(vaxparams,NNbar)
xString={'Oct 09','Nov 09','Dec 09','Jan 10','Feb 10','Mar 10','Apr 10','May10'};

y=vaxparams(:,2:end,1);
ly=size(y,2);
y=y.*repmat(NNbar,1,ly)/sum(NNbar);
f=sum(sum(y));

x=(10:10+ly-1);
fs=10; lw=2;
figure
bar(x,y,'stacked')%'-','linewidth',lw)
set(gca,'fontsize',fs)
xticks(x)
xticklabels(xString)
xtickangle(45)
xlabel('Month')
ylabel('Monthly vaccine allocation per capita')
legend('0-4','5-19','20-49','50-64','65+','location','NE');
grid on
grid minor
box on
