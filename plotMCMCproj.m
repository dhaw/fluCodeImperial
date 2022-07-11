function f=plotMCMCproj(xsto,a,b)
fs=10; lw=2; ms=5;
burn=5e4;
int=10;
figure
scatter(xsto(burn+1:int:end,a),xsto(burn+1:int:end,b),'o','markeredgecolor','k','markerfacecolor','k','markerfacealpha',.2)
%hist2d(xsto(:,a),xsto(:,b),50)
xlabel(strcat('param ',num2str(a)))
ylabel(strcat('param ',num2str(b)))
set(gca,'fontsize',fs)
axis tight
grid on
grid minor
box on
