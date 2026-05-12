function [] = plotZtemplateByHemiColocationPnon(DT,mdl,yLims)
mu = mean(DT.pNonc);
x = linspace(0,6,144);
figure;

%%
subplot(2,2,1);
colocation=-1;
hemisphere=-1;
[~,y,low,upp] = getXEUL(mdl,...
    1,x-mu,hemisphere,colocation);
shadedErrorBar(x',y',[upp'-y';y'-low'],...
    'lineprops',{'Color',[0.02 0.78 0.79]});
hold on;
s = DT.colocation==colocation & DT.hemisphere==hemisphere;
scatter(DT.pNonc(s),DT.zTemplate(s),30,'k','filled');
line([0,6],[0,0],'Color','k','LineStyle','--');
ylim(yLims);
axis square;
title('Left hippocampus, Colocation = -1');

%%
subplot(2,2,3);
colocation = +1;
hemisphere = -1;
[~,y,low,upp] = getXEUL(mdl,...
    1,x-mu,hemisphere,colocation);
shadedErrorBar(x',y',[upp'-y';y'-low'],...
    'lineprops',{'Color',[0.02 0.78 0.79]});
hold on;
s = DT.colocation==colocation & DT.hemisphere==hemisphere;
scatter(DT.pNonc(s),DT.zTemplate(s),30,'k','filled');
line([0,6],[0,0],'Color','k','LineStyle','--');
ylim(yLims);
axis square;
title('Left hippocampus, Colocation = +1');

%%
subplot(2,2,2);
colocation=-1;
hemisphere=+1;
[~,y,low,upp] = getXEUL(mdl,...
    1,x-mu,hemisphere,colocation);
shadedErrorBar(x',y',[upp'-y';y'-low'],...
    'lineprops',{'Color',[1 0.5412 0.2039]});
hold on;
s = DT.colocation==colocation & DT.hemisphere==hemisphere;
scatter(DT.pNonc(s),DT.zTemplate(s),30,'k','filled');
line([0,6],[0,0],'Color','k','LineStyle','--');
ylim(yLims);
axis square;
title('Right hippocampus, Colocation = -1');

%%
subplot(2,2,4);
colocation=+1;
hemisphere=+1;
[~,y,low,upp] = getXEUL(mdl,...
    1,x-mu,hemisphere,colocation);
shadedErrorBar(x',y',[upp'-y';y'-low'],...
    'lineprops',{'Color',[1 0.5412 0.2039]});
hold on;
s = DT.colocation==colocation & DT.hemisphere==hemisphere;
scatter(DT.pNonc(s),DT.zTemplate(s),30,'k','filled');
line([0,6],[0,0],'Color','k','LineStyle','--');
ylim(yLims);
axis square;
title('Right hippocampus, Colocation = +1');
return