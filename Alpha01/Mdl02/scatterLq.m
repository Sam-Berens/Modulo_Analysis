%% Get MNI Hippocampal masks
[dt02b,mmY] = getDataTable02b();

% Scatter plots
%coloc - 1, hem -1, weak
lables = {'l','r'};
gLables = {'Weak','Strong'};
subjectIds = getSubjectIds('G1');
sIds_S = getSubjectIds('G2');
sIds_W = subjectIds(~ismember(subjectIds,sIds_S));
panels = [1,2;3,4];
conditions = {'coloc=-1,lHipp','coloc=-1,rHipp';...
    'coloc=+1,lHipp','coloc=+1,rHipp'};
colours = [0.0392 0.7216 0.9686; 0.9686 0.5216 0.0392];

f1 = figure;
for iColoc = -1:2:1
    colocIdx = (iColoc >0) + 1;
    cPanles = panels(colocIdx,:);
    for iHemi=-1:2:1
        hemIdx = ((iHemi >0) + 1);
        cLabel = lables{hemIdx};
        y = mmY.(cLabel);
        s = dt02b.colocation == iColoc & dt02b.hemisphere == iHemi;
        nSubs = numel(sIds_W);
        subplot(2,2,cPanles(hemIdx))
        hold on
        for ii = 1:numel(sIds_W)
            cs = s & dt02b.subjectId ==sIds_W(ii);
            x = cell2mat(dt02b.lQ(cs));
            scatter(x, y, 8,...
                'Marker','x',...
                'MarkerEdgeAlpha',0.2,...
                'MarkerEdgeColor', colours(colocIdx,:));
        end
        for ii = 1:numel(sIds_W)
            cs = s & dt02b.subjectId ==sIds_W(ii);
            x = cell2mat(dt02b.lQ(cs));
            valid = ~isnan(x) & ~isnan(y);
            xv = x(valid);
            yv = y(valid);

            if numel(xv) < 2
                continue
            end

            p = polyfit(xv, yv, 1);
            xx = linspace(min(xv), max(xv), 100);
            yy = polyval(p, xx);
            plot(xx, yy, 'Color', colours(colocIdx,:).*0.6);
        end
        title(conditions{colocIdx,hemIdx});
        ylim([-45,10]);
        xlim([-8,3]);
        hold off;
    end
end

%% 

f2 = figure;
for iColoc = -1:2:1
    colocIdx = (iColoc >0) + 1;
    cPanles = panels(colocIdx,:);
    for iHemi=-1:2:1
        hemIdx = ((iHemi >0) + 1);
        cLabel = lables{hemIdx};
        y = mmY.(cLabel);
        s = dt02b.colocation == iColoc & dt02b.hemisphere == iHemi;
        nSubs = numel(sIds_S);
        subplot(2,2,cPanles(hemIdx))
        hold on
        for ii = 1:numel(sIds_S)
            cs = s & dt02b.subjectId ==sIds_S(ii);
            x = cell2mat(dt02b.lQ(cs));
            scatter(x, y, 8,...
                'Marker','x',...
                'MarkerEdgeAlpha',0.2,...
                'MarkerEdgeColor', colours(colocIdx,:));
        end 
        for ii = 1:numel(sIds_S)
            cs = s & dt02b.subjectId ==sIds_S(ii);
            x = cell2mat(dt02b.lQ(cs));
            valid = ~isnan(x) & ~isnan(y);
            xv = x(valid);
            yv = y(valid);

            if numel(xv) < 2
                continue
            end

            p = polyfit(xv, yv, 1);
            xx = linspace(min(xv), max(xv), 100);
            yy = polyval(p, xx);
            plot(xx, yy, 'Color', colours(colocIdx,:).*0.6);
        end
        title(conditions{colocIdx,hemIdx});
        ylim([-45,10]);
        xlim([-8,3])
        hold off;
    end
end
