function f=subSeverity(NNbar,xsto,xdata,ydata,vaxparams,stateName)
%{
Plot decision framework for a given state:
Compare 2009 with a hypothetical flu pandemic with inlated clinical
severity, with and wihout school closures triggered. 
%}

burn=20000;
int=200;
na=length(NNbar);
pstar=0:.05:1;
Hvec=0:5e1:3e3;
numIter=100;

lp=length(pstar);
lh=length(Hvec);

meanHosp=[393.0256  999.3054  702.6867  406.0680  177.1958];
sdHosp=[129.6984 329.7708 231.8866 134.0024 58.4746];
%
%Counterfactual scenario with clinical severity inflated by "factor":
factor=2;
uncertaintyThresh=.1;%Pre-determined uncertainty threshold
ploc=1500;%Maximum allowed admissions
M=[meanHosp;meanHosp/factor;meanHosp/factor];
S=[sdHosp;sdHosp/factor;sdHosp/factor];
tswitch=243+[0,0,70];
lm=size(M,1);
X=zeros(lm,lh);
%
xsto=xsto(burn+1:int:end,:);
lx=size(xsto,1);
for m=1:lm
    meanHosp=M(m,:);
    sdHosp=S(m,:);
    hospVals=zeros(lx,numIter);
    for i=1:lx
        [~,~,z2]=subPandemicSimulationVax(NNbar,xsto(i,:),xdata,0,0,ydata,tswitch(m),vaxparams);
        for j=1:numIter
            mult=normrnd(meanHosp,sdHosp);
            hospVals(i,j)=sum(z2./mult);
        end
    end
    hospVals=reshape(hospVals,lx*numIter,1);
    for h=1:lh
        X(m,h)=length(hospVals(hospVals>Hvec(h)))/lx/numIter;
    end
end
f=X;
%
fs=10; lw=2; ms=7;
cmap=lines(7);
figure;
hold on
h1=plot(X(1,:),Hvec,'linewidth',lw,'color',cmap(1,:));
h2=plot(X(2,:),Hvec,'linewidth',lw,'color',cmap(2,:));
h3=plot(X(3,:),Hvec,'--','linewidth',lw,'color',cmap(2,:));
plot([.1,.1],[Hvec(1),Hvec(end)],'k--','linewidth',1.5)
plot(uncertaintyThresh,ploc,'ko','markersize',ms,'markerfacecolor','k')
set(gca,'fontsize',fs)
title(strcat('Probabalistic risk score: ',stateName))
xlabel('Proportion of simulations having fall-wave hospitalisations in excess of h')
ylabel('h')
axis([0,1,Hvec(1),Hvec(end)])
legend([h1,h2,h3],'2009 pandemic-like virus','Virus with 2x severity of 2009 pandemic','2x severity with 10-week delay in school openings','location','NE')%,h2,h3
grid on
grid minor
box on