function [fgh] = plotZtemplateByColocationPnon(DataTable,mdl,roiName)
%to do sdd xlim 
fgh = figure;

mu_Pnonc = mean(DataTable.pNonc);
pNonc = linspace(0,6,144);

subplot(1,2,1);
hold on;
[~,y,low,upp] = getXEUL(mdl,...
    1,pNonc-mu_Pnonc,-1);
shadedErrorBar(pNonc',y',[upp'-y';y'-low'],...
    'lineprops',{'Color',[0.9686 0.5216 0.0392]});
s = DataTable.colocation==-1;
scatter(DataTable.pNonc(s),DataTable.zTemplate(s),30,'k','filled');
line([0,6],[0,0],'Color','k','LineStyle','--');
ylim([-0.6,0.6]);
axis square;
title('Colocation = -1');

subplot(1,2,2);
hold on;
[~,y,low,upp] = getXEUL(mdl,...
    1,pNonc-mu_Pnonc,+1);
shadedErrorBar(pNonc',y',[upp'-y';y'-low'],...
    'lineprops',{'Color',[0.0392 0.7216 0.9686]});
s = DataTable.colocation==+1;
scatter(DataTable.pNonc(s),DataTable.zTemplate(s),30,'k','filled');
line([0,6],[0,0],'Color','k','LineStyle','--');
ylim([-0.6,0.6]);
axis square;
title('Colocation = +1');
hold off;

annotation(fgh, 'textbox', ...
    [0.25 0.9 0.5 0.05], ...   % [x y width height] in normalized figure units
    'String', roiName, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');

return
