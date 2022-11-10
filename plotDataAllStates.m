function f=plotDataAllStates(cellNN,cellx,cellyx)
%Plot fluSurvNet data for all 10 states
%Inputs are cell arrays of data for each state:
%Age group populatios; time points; hospital admissions

names={'California','Connecticut','Colorado','Georgia','Maryland','Minnesota','New Mexico','New York','Oregon','Tennessee'};
lc=length(cellNN);
fs=10; lw=2;
figure;
tiledlayout(lc/2,2,'TileSpacing','compact')
for i=1:lc
    NNi=cellNN{i};
    xi=cellx{i};
    yi=cellyx{i};
    b=size(yi,1);
    yi=yi/nansum(NNi)*1000;
    %
    nexttile
    %
    hold on
    h=plot(xi,yi,'linewidth',lw);
    plot(35*[1,1],[0,100],'k--','linewidth',lw)
    hold off
    set(gca,'fontsize',fs)
    title(names{i})
    axis ([17,52,0,ceil(nanmax(nanmax(yi))*100)/100])
    grid on
    grid minor
    box on
    if i>lc-2
        xlabel('MMWR week');
    end
    if i==5
        ylabel('FluSurv-NET weekly hospitalisations/1000');
    end
end

lg=legend(nexttile(2),'0-4','5-19','20-49','50-64','65+');
lg.Location='northeastoutside';