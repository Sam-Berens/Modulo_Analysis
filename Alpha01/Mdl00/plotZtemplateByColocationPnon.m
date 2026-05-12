function [] = plotZtemplateByColocationPnon(DT,mdl,yLims)
mu = mean(DT.pNonc);
x = linspace(0,6,144);
figure;
subplot(1,2,1);
[~,y,low,upp] = getXEUL(mdl,...
    1,x-mu,0,-1);
shadedErrorBar(x',y',[upp'-y';y'-low'],...
    'lineprops',{'Color',[1 0.5412 0.2039]});
hold on;
s = DT.colocation==-1;
scatter(DT.pNonc(s),DT.zTemplate(s),30,'k','filled');
line([0,6],[0,0],'Color','k','LineStyle','--');
ylim(yLims);
axis square;
title('Colocation = -1');

subplot(1,2,2);
[~,y,low,upp] = getXEUL(mdl,...
    1,x-mu,0,+1);
shadedErrorBar(x',y',[upp'-y';y'-low'],...
    'lineprops',{'Color',[0.02 0.78 0.79]});
hold on;
s = DT.colocation==+1;
scatter(DT.pNonc(s),DT.zTemplate(s),30,'k','filled');
line([0,6],[0,0],'Color','k','LineStyle','--');
ylim(yLims);
axis square;
title('Colocation = +1');
return