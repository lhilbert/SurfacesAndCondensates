%% --- load analysis results

clear all
load('sortedResultsBundle')

numUniqConds = numel(uniqueCondNames);

numSurf = zeros(1,numUniqConds);
numSurf_CI = zeros(2,numUniqConds);

numCond = zeros(1,numUniqConds);
numCond_CI = zeros(2,numUniqConds);


SurfVol = zeros(1,numUniqConds);
SurfVol_CI = zeros(2,numUniqConds);

SurfSol = zeros(1,numUniqConds);
SurfSol_CI = zeros(2,numUniqConds);

CondVol = zeros(1,numUniqConds);
CondVol_CI = zeros(2,numUniqConds);

CondSol = zeros(1,numUniqConds);
CondSol_CI = zeros(2,numUniqConds);


SurfSurf_Int = zeros(1,numUniqConds);
SurfSurf_Int_CI = zeros(2,numUniqConds);

SurfCond_Int = zeros(1,numUniqConds);
SurfCond_Int_CI = zeros(2,numUniqConds);

CondSurf_Int = zeros(1,numUniqConds);
CondSurf_Int_CI = zeros(2,numUniqConds);

CondCond_Int = zeros(1,numUniqConds);
CondCond_Int_CI = zeros(2,numUniqConds);

SurfCond_dist = zeros(1,numUniqConds);
SurfCond_dist_CI = zeros(2,numUniqConds);

CondSurf_dist = zeros(1,numUniqConds);
CondSurf_dist_CI = zeros(2,numUniqConds);

n_boot = 500;

for cc = 1:numUniqConds

    thisCondInd = cc;

    numSurf(cc) = mean(sortedSurfaceNumCell{cc});
    numSurf_CI(:,cc) = bootci(...
        n_boot,@mean,sortedSurfaceNumCell{cc});

    numCond(cc) = mean(sortedDropletNumCell{cc});
    numCond_CI(:,cc) = bootci(...
        n_boot,@mean,sortedDropletNumCell{cc});

    % --- object volume and solidity calculation, start
    SurfVol(cc) = mean(sortedSurfaceVolCell{cc});
    SurfVol_CI(:,cc) = bootci(...
        n_boot,@mean,sortedSurfaceVolCell{cc});

    SurfSol(cc) = mean(sortedSurfaceSolCell{cc});
    SurfSol_CI(:,cc) = bootci(...
        n_boot,@mean,sortedSurfaceSolCell{cc});

    CondVol(cc) = mean(sortedDropletVolCell{cc});
    CondVol_CI(:,cc) = bootci(...
        n_boot,@mean,sortedDropletVolCell{cc});

    CondSol(cc) = mean(sortedDropletSolCell{cc});
    CondSol_CI(:,cc) = bootci(...
        n_boot,@mean,sortedDropletSolCell{cc});
    % --- object volume and solidity calculation, end

    % --- object intensity calculation, start
    SurfSurf_Int(cc) = mean(sortedSurfaceIntCell{1}{cc});
    try
        SurfSurf_Int_CI(:,cc) = bootci(...
        n_boot,@mean,sortedSurfaceIntCell{1}{cc});
    catch
        SurfSurf_Int_CI(:,cc) = SurfSurf_Int(cc);
    end

    SurfCond_Int(cc) = mean(sortedSurfaceIntCell{2}{cc});
    try
        SurfCond_Int_CI(:,cc) = bootci(...
        n_boot,@mean,sortedSurfaceIntCell{2}{cc});
    catch
        SurfCond_Int_CI(:,cc) = SurfCond_Int(cc);
    end

    CondSurf_Int(cc) = mean(sortedDropletIntCell{1}{cc});
    CondSurf_Int_CI(:,cc) = bootci(...
        n_boot,@mean,sortedDropletIntCell{1}{cc});
    CondCond_Int(cc) = mean(sortedDropletIntCell{2}{cc});
    CondCond_Int_CI(:,cc) = bootci(...
        n_boot,@mean,sortedDropletIntCell{2}{cc});
    % --- object intensity calculation, end
    
    SurfCond_dist(cc) = mean(sortedSurfDropDistCell{cc});
    SurfCond_dist_CI(:,cc) = bootci(...
        n_boot,@mean,sortedSurfDropDistCell{cc});

    CondSurf_dist(cc) = mean(sortedDropSurfDistCell{cc});
    try
        CondSurf_dist_CI(:,cc) = bootci(...
            n_boot,@mean,sortedDropSurfDistCell{cc});
    catch
        CondSurf_dist_CI(:,cc) = CondSurf_dist(cc);
    end
        

end


%% -- Overview figure

figure(1)
clf



subplot(2,6,1)
errorbar(1:numUniqConds, ...
    numSurf,...
    numSurf_CI(2,:)-numSurf,...
    numSurf-numSurf_CI(1,:),...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
set(gca,'XTick',1:numUniqConds,'XTickLabel',sortedCondNames,...
    'XTickLabelRotation',45,'XLim',[0.5,numUniqConds+0.5])
ylabel('Number surfaces')
title('Surfaces')

subplot(2,6,2)
errorbar(1:numUniqConds, ...
    SurfVol,...
    SurfVol_CI(2,:)-SurfVol,...
    SurfVol-SurfVol_CI(1,:),...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
set(gca,'XTick',1:numUniqConds,'XTickLabel',sortedCondNames,...
    'XTickLabelRotation',45,'XLim',[0.5,numUniqConds+0.5])
ylabel('Surf Vol [\mum^3]')

subplot(2,6,3)
errorbar(1:numUniqConds, ...
    SurfSol,...
    SurfSol_CI(2,:)-SurfSol,...
    SurfSol-SurfSol_CI(1,:),...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
set(gca,'XTick',1:numUniqConds,'XTickLabel',sortedCondNames,...
    'XTickLabelRotation',45,'XLim',[0.5,numUniqConds+0.5])
ylabel('Surf Solidity')



subplot(2,6,7)
errorbar(1:numUniqConds, ...
    SurfSurf_Int,...
    SurfSurf_Int_CI(2,:)-SurfSurf_Int,...
    SurfSurf_Int-SurfSurf_Int_CI(1,:),...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
set(gca,'XTick',1:numUniqConds,'XTickLabel',sortedCondNames,...
    'XTickLabelRotation',45,'XLim',[0.5,numUniqConds+0.5])
ylabel('Surf-Surf Int.')

subplot(2,6,8)
errorbar(1:numUniqConds, ...
    SurfCond_Int,...
    SurfCond_Int_CI(2,:)-SurfCond_Int,...
    SurfCond_Int-SurfCond_Int_CI(1,:),...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
set(gca,'XTick',1:numUniqConds,'XTickLabel',sortedCondNames,...
    'XTickLabelRotation',45,'XLim',[0.5,numUniqConds+0.5])
ylabel('Surf-Cond Int.')


subplot(2,6,9)
errorbar(1:numUniqConds, ...
    SurfCond_dist,...
    SurfCond_dist_CI(2,:)-SurfCond_dist,...
    SurfCond_dist-SurfCond_dist_CI(1,:),...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
set(gca,'XTick',1:numUniqConds,'XTickLabel',sortedCondNames,...
    'XTickLabelRotation',45,'XLim',[0.5,numUniqConds+0.5])
ylabel('Surf-Cond distance [\mum]')




subplot(2,6,4)
errorbar(1:numUniqConds, ...
    numCond,...
    numCond_CI(2,:)-numCond,...
    numCond-numCond_CI(1,:),...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
set(gca,'XTick',1:numUniqConds,'XTickLabel',sortedCondNames,...
    'XTickLabelRotation',45,'XLim',[0.5,numUniqConds+0.5])
ylabel('Number condensates')
title('Condensates')

subplot(2,6,5)
errorbar(1:numUniqConds, ...
    CondVol,...
    CondVol_CI(2,:)-CondVol,...
    CondVol-CondVol_CI(1,:),...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
set(gca,'XTick',1:numUniqConds,'XTickLabel',sortedCondNames,...
    'XTickLabelRotation',45,'XLim',[0.5,numUniqConds+0.5])
ylabel('Cond Vol [\mum^3]')

subplot(2,6,6)
errorbar(1:numUniqConds, ...
    CondSol,...
    CondSol_CI(2,:)-CondSol,...
    CondSol-CondSol_CI(1,:),...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
set(gca,'XTick',1:numUniqConds,'XTickLabel',sortedCondNames,...
    'XTickLabelRotation',45,'XLim',[0.5,numUniqConds+0.5])
ylabel('Cond Solidity')



subplot(2,6,10)
errorbar(1:numUniqConds, ...
    CondSurf_Int,...
    CondSurf_Int_CI(2,:)-CondSurf_Int,...
    CondSurf_Int-CondSurf_Int_CI(1,:),...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
set(gca,'XTick',1:numUniqConds,'XTickLabel',sortedCondNames,...
    'XTickLabelRotation',45,'XLim',[0.5,numUniqConds+0.5])
ylabel('Cond-Surf Int.')

subplot(2,6,11)
errorbar(1:numUniqConds, ...
    CondCond_Int,...
    CondCond_Int_CI(2,:)-CondCond_Int,...
    CondCond_Int-CondCond_Int_CI(1,:),...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
set(gca,'XTick',1:numUniqConds,'XTickLabel',sortedCondNames,...
    'XTickLabelRotation',45,'XLim',[0.5,numUniqConds+0.5])
ylabel('Cond-Cond Int.')




subplot(2,6,12)
errorbar(1:numUniqConds, ...
    CondSurf_dist,...
    CondSurf_dist_CI(2,:)-CondSurf_dist,...
    CondSurf_dist-CondSurf_dist_CI(1,:),...
    'ko-','Color',[0.0,0.0,0.0],'LineWidth',1,...
    'MarkerFaceColor',[0,0,0])
set(gca,'XTick',1:numUniqConds,'XTickLabel',sortedCondNames,...
    'XTickLabelRotation',45,'XLim',[0.5,numUniqConds+0.5])
ylabel('Cond-Surf distance [\mum]')


%% --- Scatter plots of intensities and solidity

figure(2)
clf

for cc = 1:numUniqConds

    thisCondInd = cc;
    


    % --- object volume and solidity calculation, start
    SurfVol = sortedSurfaceVolCell{cc};
    SurfSol = sortedSurfaceSolCell{cc};
    CondVol = sortedDropletVolCell{cc};
    CondSol = sortedDropletSolCell{cc};
    % --- object volume and solidity calculation, end

    % --- object intensity calculation, start
    SurfSurf_Int = sortedSurfaceIntCell{1}{cc};
    SurfCond_Int = sortedSurfaceIntCell{2}{cc};
    CondSurf_Int = sortedDropletIntCell{1}{cc};
    CondCond_Int = sortedDropletIntCell{2}{cc};
    % --- object intensity calculation, end
    
    subplot(2,numUniqConds,cc)
    scatter(SurfVol,SurfCond_Int,...
        5,SurfSol')
    colormap(parula)
    set(gca,'CLim',[0,1],'Box','on')
    set(gca,'XLim',[0,50],'YLim',[0,20000])
    a=colorbar;
    a.Label.String = 'Surface Solidity';
    title(sortedCondNames(cc))
    xlabel('Surf. Vol. [\mum^3]')
    ylabel('Surf. Cond. Int.')

    subplot(2,numUniqConds,numUniqConds+cc)
    scatter(CondSurf_Int,CondCond_Int,...
        5,CondVol')
    colormap(parula)
    set(gca,'CLim',[0,5],'Box','on')
    set(gca,'XLim',[0,20000],'YLim',[0,20000])
    a=colorbar;
    a.Label.String = 'Cond. Volume [\mum^3]';
    title(sortedCondNames(cc))
    xlabel('Int. Surface')
    ylabel('Int. X-Motif')
    

end

%% --- make plot of full and control condition combined

figure(3)

clf

plotConds = [3,4,5,6,1];

for_count = 0;

for cc = plotConds

    thisCondInd = cc;
    for_count = for_count + 1;
    
    % --- object volume and solidity calculation, start
    SurfVol = sortedSurfaceVolCell{cc};
    SurfSol = sortedSurfaceSolCell{cc};
    CondVol = sortedDropletVolCell{cc};
    CondSol = sortedDropletSolCell{cc};
    % --- object volume and solidity calculation, end

    % --- object intensity calculation, start
    SurfSurf_Int = sortedSurfaceIntCell{1}{cc};
    SurfCond_Int = sortedSurfaceIntCell{2}{cc};
    CondSurf_Int = sortedDropletIntCell{1}{cc};
    CondCond_Int = sortedDropletIntCell{2}{cc};
    % --- object intensity calculation, end
    
    scatter(CondSurf_Int,CondCond_Int,...
        CondVol./5,'filled')
    colormap(flipud(parula))
    set(gca,'CLim',[0,800],'Box','on')
    set(gca,'XLim',[0,15000],'YLim',[0,25000])
%     a=colorbar;
%     a.Label.String = 'Cond. Volume [\mum^3]';
    xlabel('Int. Surface')
    ylabel('Int. X-Motif')

    hold on

    sortCondVol = sort(CondVol,'descend');
    sprintf('%s, top 1-percentile: %2.2f',...
        sortedCondNames{cc},prctile(sortCondVol,99))

end

legend(sortedCondNames(plotConds))