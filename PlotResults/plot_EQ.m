if (PLOT_MODE == 'LOG')
    semilogx(sim1.lamSim, sim1.MEASURES.EQ/EG,'o', 'MarkerSize',simMarker, 'color', 'black'); hold on
else
    plot(sim1.lamSim, sim1.MEASURES.EQ/EG,'o', 'MarkerSize',simMarker, 'color', 'black'); hold on
end    
%% Averege single queue content
switch nPlots
    case 1
        plot(lam(:,1), RES(1).Q(:,plTerm)/EG,  'LineWidth', lw);
    case 2
        plot(lam(:,1), RES(1).Q(:,plTerm)/EG, LineStyle1, ...
                    lam(:,2), RES(2).Q(:,plTerm)/EG, LineStyle2, ...
                    'LineWidth', lw);
    case 3
        plot(lam(DIAP{1},1), RES(1).Q(DIAP{1},plTerm)/EG,...
                    lam(DIAP{2},2), RES(2).Q(DIAP{2},plTerm)/EG,...
                    lam(DIAP{3},3), RES(3).Q(DIAP{3},plTerm)/EG,...
                    'LineWidth', lw);         
%         plot(lam(:,1), RES(1).Q(:,plTerm)/EG,...
%                     lam(:,2), RES(2).Q(:,plTerm)/EG,...
%                     lam(:,3), RES(3).Q(:,plTerm)/EG,...
%                     'LineWidth', lw);                
end        
plot([1],[0.846], 'p', 'MarkerSize',markerSize, 'MarkerFaceColor', 'black')          
plot([10.5],[2.77], 'p', 'MarkerSize',markerSize , 'MarkerFaceColor', 'black')

ylim([0 maxEQ]); %grid on;
xlim([minLam,maxLam]); 
xlabel('$\lambda$, Mbit/s','FontSize', fontSize,'Interpreter','latex'); 
ylabel('$E[Q]$', 'FontSize', fontSize, 'Interpreter','latex');
handle = title(strcat('N= ',int2str(plTerm)),'FontSize', fontSize,'Interpreter','latex')
set(handle,'Position',titlePos);   
