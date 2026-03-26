function [fgh] = plotRoiEffect(cDT,mdl)
mu = mean(cDT.pNonc);
x = linspace(0,6,144);
fgh = figure;
subplot(1,2,1);
[~,y,low,upp] = getXEUL(mdl,...
    1,x-mu,-1);
shadedErrorBar(x',y',[upp'-y';y'-low'],...
    'lineprops',{'Color',[0.4,0.2,0.8]});
hold on;
s = cDT.colocation==-1;
scatter(cDT.pNonc(s),cDT.zTemplate(s),30,'k','filled');
line([0,6],[0,0],'Color','k','LineStyle','--');
ylim([-0.6,0.6]);
axis square;
title('Colocation = -1');

subplot(1,2,2);
[~,y,low,upp] = getXEUL(mdl,...
    1,x-mu,+1);
shadedErrorBar(x',y',[upp'-y';y'-low'],...
    'lineprops',{'Color',[0.2,0.8,0.4]});
hold on;
s = cDT.colocation==+1;
scatter(cDT.pNonc(s),cDT.zTemplate(s),30,'k','filled');
line([0,6],[0,0],'Color','k','LineStyle','--');
ylim([-0.6,0.6]);
axis square;
title('Colocation = +1');
hold off;
return
