if (PLOT_MODE == 'LOG')
    semilogx(sim1.lamSim, sim1.MEASURES.ED,'o', 'MarkerSize',simMarker,'color', 'black'); hold on
else    
    plot(sim1.lamSim, sim1.MEASURES.ED,'o', 'MarkerSize',simMarker,'color', 'black'); hold on
end    
handle = title(strcat('N = ',int2str(plTerm)),'FontSize', fontSize,'Interpreter','latex')
set(handle,'Position',titlePos); 
%plot(sim1.lamSim, sim1.MEASURES.EDs,'o', 'MarkerSize',simMarker,'color', 'black'); hold on
switch nPlots
    case 1
        plot(lam(:,1), RES(1).EDs(:,plTerm),  'LineWidth', lw);    
    case 2
        plot(lam(:,1), RES(1).ED(:,plTerm), LineStyle1, ...
                    lam(:,2), RES(2).ED(:,plTerm), LineStyle2, ...
                    'LineWidth', lw);
    case 3
        plot(lam(DIAP{1},1), RES(1).ED(DIAP{1},plTerm),...
                    lam(DIAP{2},2), RES(2).ED(DIAP{2},plTerm),...
                    lam(DIAP{3},3), RES(3).ED(DIAP{3},plTerm),...
                    'LineWidth', lw);            
end
plot([1],[0.96], 'p', 'MarkerSize',15, 'MarkerFaceColor', 'black')        
plot([10],[1.3], 'p', 'MarkerSize',15, 'MarkerFaceColor', 'black')

ylim([0 maxDelay]); %grid on;
xlim([minLam,maxLam]); 
xlabel('$\lambda$, Mbit/s','FontSize', fontSize,'Interpreter','latex');% ylabel('mean packet delay, s', 'FontSize', fontSize);
ylabel('$E[D]$, s', 'FontSize', fontSize,'Interpreter','latex'); % set(h,'Interpreter','latex')
