if (PLOT_MODE == 'LOG')
    semilogx(sim1.lamSim, sim1.lamSim-sim1.MEASURES.Acc,'o','MarkerSize',simMarker,  'color', 'black'); hold on
else
    plot(sim1.lamSim, sim1.lamSim-sim1.MEASURES.Acc,'o','MarkerSize',simMarker,  'color', 'black'); hold on
end

%plot(sim1.lamSim, sim1.lamSim-sim1.MEASURES.Acc,'o','MarkerSize',simMarker,  'color', 'black'); hold on
handle = title(strcat('N = ',int2str(plTerm)),'FontSize', fontSize,'Interpreter','latex')
set(handle,'Position',titlePos); 
switch nPlots
    case 1
        plot(lam(:,1), RES(1).REJ(:,plTerm), 'LineWidth', lw);    
    case 2
        plot(lam(:,1), RES(1).REJ(:,plTerm), LineStyle1, ...
                    lam(:,2), RES(2).REJ(:,plTerm), LineStyle2, ...
                    'LineWidth', lw);        
    case 3
        plot(lam(DIAP{1},1), RES(1).REJ(DIAP{1},plTerm),...
                    lam(DIAP{2},2), RES(2).REJ(DIAP{2},plTerm),...
                    lam(DIAP{3},3), RES(3).REJ(DIAP{3},plTerm),...
                    'LineWidth', lw);             
%         plot(lam(:,1), RES(1).REJ(:,plTerm),...
%                     lam(:,2), RES(2).REJ(:,plTerm),...
%                     lam(:,3), RES(3).REJ(:,plTerm),...
%                     'LineWidth', lw);
end

plot([1.01],[0.11], 'p', 'MarkerSize',15, 'MarkerFaceColor', 'black')
plot([10],[7.8] , 'p', 'MarkerSize',15, 'MarkerFaceColor', 'black')

xlim([minLam,maxLam]); 
xlabel('$\lambda$, Mbit/s','FontSize', fontSize,'Interpreter','latex'); %ylabel('average rejection rate mbit/s', 'FontSize', fontSize);
ylabel('$\gamma_{rej}$, Mbit/s', 'FontSize', fontSize,'Interpreter','latex');
ylim([0 maxReject]); %grid on; 