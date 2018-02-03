if (PLOT_MODE == 'LOG')
    semilogx(sim1.lamSim, sim1.MEASURES.Acc-sim1.MEASURES.Tr,'o','MarkerSize',simMarker,  'color', 'black'); hold on
else
    plot(sim1.lamSim, sim1.MEASURES.Acc-sim1.MEASURES.Tr,'o','MarkerSize',simMarker,  'color', 'black'); hold on
end
%plot(sim1.lamSim, sim1.MEASURES.Acc-sim1.MEASURES.Tr,'o','MarkerSize',simMarker,  'color', 'black'); hold on
handle = title(strcat('N = ',int2str(plTerm)),'FontSize', fontSize,'Interpreter','latex')
set(handle,'Position',titlePos); 
switch nPlots
    case 1
        plot(lam(:,1), RES(1).DIS(:,plTerm), 'LineWidth', lw);
    case 2
        plot(lam(:,1), RES(1).DIS(:,plTerm), LineStyle1, ...
                    lam(:,2), RES(2).DIS(:,plTerm), LineStyle2, ...
                    'LineWidth', lw);   
    case 3
        plot(lam(DIAP{1},1), RES(1).DIS(DIAP{1},plTerm),...
                    lam(DIAP{2},2), RES(2).DIS(DIAP{2},plTerm),...
                    lam(DIAP{3},3), RES(3).DIS(DIAP{3},plTerm),...
                    'LineWidth', lw);          
%         plot(lam(:,1), RES(1).DIS(:,plTerm),...
%                     lam(:,2), RES(2).DIS(:,plTerm),...
%                     lam(:,3), RES(3).DIS(:,plTerm),...
%                     'LineWidth', lw);
end
           
plot([1],[0.054], 'p', 'MarkerSize',15, 'MarkerFaceColor', 'black')
plot([10],[0.12], 'p', 'MarkerSize',15, 'MarkerFaceColor', 'black')

ylim([0, maxDiscard]); %grid on;
xlim([minLam,maxLam]); 
xlabel('$\lambda$, Mbit/s','FontSize', fontSize,'Interpreter','latex'); %ylabel('average dircart rate, mbit/s', 'FontSize', fontSize);
ylabel('$\gamma_{dis}$, Mbit/s', 'FontSize', fontSize, 'Interpreter','latex');