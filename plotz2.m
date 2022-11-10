function f=plotz2(cellZ2,real,NNsums)
%{
Plot cumulative admissions in fall wave for each state
cellZ2 - cell array of cumulative admissions by age, 
sampled from MCMC output
real - array of true cumulative admissions by age fall wave
%}

numplots=10;
numages=5;
%cell in alphabetical order of 2-letter state code
if numplots==8
    names={'CA','CO','CT','GA','MD','NY','OR','TN'};%Not MN/NM/US
elseif numplots==6
    names={'CO','CT','GA','MD','NY','OR','TN'};%Not CA/MN/NM/US
elseif numplots==10
    names={'California','Colorado','Connecticut','Georgia','Maryland','Minnesota','New Mexico','New York','Oregon','Tennessee'};%Not CA/MN/NM/US
end
hospUnder=1;
hospUnderMin=1.7/2.74;
hospUnderMax=4.5/2.74;

legString={'0-4','-17','-49','-64','65+'};
cmap=lines(7);
transp=.1;

lc=length(cellZ2);
fs=10; lw=2;
figure('DefaultAxesPosition', [0.1, 0.1, 0.875, 0.875])
tiledlayout(lc/2,2,'TileSpacing','compact')
%
X=zeros(numplots,numages,3);
%
for i=1:numplots
    nexttile
    
    ci=cellZ2{i};
    %
    %Comment out for Z1
    ci=ci';
    mu=[393.0256 999.3054 702.6867 406.0680 177.1958];
    sig=[129.6984 329.7708 231.8866 134.0024 58.4746];
    repmu=repmat(mu',1,size(ci,2));
    repsig=repmat(sig',1,size(ci,2));
    ci=ci./normrnd(repmu,repsig);
    ci=ci';
    %}
    ci=ci/NNsums(i)*1000;
    reali=real(:,i);
    maxy=max(max(ci));
    maxy=max(maxy,max(reali));
    [lp,nbar]=size(ci);
    onesi=ones(lp,1);
    hold on
    violinplot(ci);
    %
    for j=1:nbar
        plot(j+.1,hospUnder*reali(j),'kx','markersize',10,'linewidth',lw)
        plot([j,j]+.1,reali(j)*[hospUnderMin,hospUnderMax],'k-','linewidth',lw)
    end
    hold off
    if i>numplots-2
        xlabel('Age group')
    end
    if i==1 || rem(i,2)==1
        ylabel('Hosp./1000')
    else
        yticks(0:.5:1)
        yticklabels({'',''})
    end
    title(names{i})
    set(gca,'fontsize',fs)
    xticks(1:nbar)
    xticklabels(legString)
    axis ([0,nbar+1,0,.6])
    grid on
    grid minor
    box on
end