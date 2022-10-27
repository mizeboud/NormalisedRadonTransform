%% Plot effect of normalisation for different example windows
% M. Izeboud, TU Delft, 09/2021
% Script to create Figure 1 and Figure 2 in the paper "Damage Detection on 
% Antarctic Ice Shelves using the Normalised Radon Transform" by Izeboud et
% al. (2022, in review at Remote Sensing of Environment). 
%
% This script requires the romaO and batlow colormap from Crameri, F. (2018)
%   Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862

%% Get colormaps
addpath('./matlab/functions/')

romaO = load('./files/romao.mat').romaO;
batlow = load('./files/batlow.mat').batlow;

%% PLOT EXAMPLE WINDOW WITH original VS norrmalised RADON TRANSFORM
% path2save = './example/';

RT_mark = '*';
RT_line = ':';
RTN_mark = 'o';
RTN_line = '-';
fontSize = 15;
RT_color = batlow(15,:);
RTN_color = batlow(200,:);

close all

CASES = ["uniform";"gray1";"gray2";"complex1";"complex2"];

% select single example:
CASES = "gray2";

for C = 1:length(CASES)
   
switch CASES(C)
    case "uniform"
        window = ones(9,9);
    case "gray1"
        window = zeros(9,9);
        window(4,:) = 0.5;
    case "gray2"
        window = zeros(9,9);
        window(4,:) = 0.5;
        window = window+0.5;
        window(window>1) = 1;
    case "complex1"
        window = zeros(9,9);
        window(6,:)=1; 
        xy = [1,7;2,8;3,9;];
        for i = 1:3
            window(xy(i,1),xy(i,2))=0.5;
        end
    case "complex2"
        window = zeros(9,9);
        window(6,:)=0.5;
        xy = [1,7;2,8;3,9;];
        for i = 1:3
            window(xy(i,1),xy(i,2))=1;
        end
        
end
    
RT = radon(window);
RTnorm = radoNorm(window); 


figPos = [441   368   534   429];

%% -- plot example window/IMG    
figure('Position',figPos) 
imagesc(window), axis equal, axis tight, axis off
    ax1 = gca; 
    colormap(ax1,'gray');caxis([0 1])
    ax1.XTick = 1:1:size(window,1);
    ax1.YTick = 1:1:size(window,2);
    ax2 = axes('Position',ax1.Position,...
      'Color','none', 'Ylim',ax1.YLim,'XLim',ax1.XLim,...
      'TickLength',[0 0],...
      'YTick', get(ax1,'ytick')-0.5,... 
      'XTick', get(ax1,'xtick')-0.5,  ...
      'YTickLabel', [],'XTickLabel', [], ...
      'gridcolor','r','linewidth',2);
    axis equal, axis tight 
    ax2.YAxis.Visible = 'off'; ax2.XAxis.Visible = 'off';
    ax2.XLim = ax2.XLim+0.5; ax2.YLim = ax2.YLim+0.5;
    grid on
    linkaxes([ax1 ax2],'xy')
    
%     saveas(gcf,[path2save 'normalised_radon_example_' char(CASES(C)) '_plt1'],'epsc');


  
%% -- plot full RT feature space

figure('position',figPos);
CLim = [0 round( max( max(RT,[],'all','omitnan'), max(RTnorm,[],'all','omitnan') ) ,1)];
clim_RT = 12;
cticks_RT = [0 6 12]; 
cticks_RTN= [0 0.5 1];

ax1= subplot(211);
    imagesc(RT,'AlphaData',~isnan(RT))
    colormap(batlow),caxis([0 12]);
    set(ax1,'fontsize',fontSize), axis(ax1,'tight')
    % Yaxis
    labels = nan*ones(1,size(RT,1)); labels(5) = 5; labels(10)=10;
    labels = strrep(strsplit(num2str(labels),' '),'NaN',' ');
    ax1.YTick = 1:size(RT,1); ax1.YTickLabel = labels;
    ax1.YLabel.String = 'Projection axis \rho';
    % Xaxis
    ax1.XTick = [0:30:180];
    % Grid
    ax1.XGrid='on';ax1.YGrid='on'; 
    ax1.GridAlpha = 0.5; ax1.GridColor = [0.5 0.5 0.5];
    ax1.XAxis.Visible = 'off';% hide x-axis (only show on bottom subplot)
    % Colorbar
    cbar=colorbar(); cbar.Ticks = cticks_RT; 
    cbar.Label.String = 'value of line integral';
    title('Radon Transform R(\rho,\theta)')

ax2 = subplot(212); 
    imagesc(RTnorm,'AlphaData',~isnan(RTnorm))
    colormap(batlow),caxis([0 1]);
    set(ax2,'fontsize',fontSize), axis(ax2,'tight');
    % Y axis
    labels = nan*ones(1,size(RTnorm,1)); labels(5) = 5; labels(10)=10;
    labels = strrep(strsplit(num2str(labels),' '),'NaN',' ');
    ax2.YTick = 1:size(RTnorm,1); ax2.YTickLabel = labels;
    ax2.YLabel.String = 'Projection axis \rho';
    % Grid
    ax2.XGrid='on'; ax2.YGrid='on';
    ax2.GridAlpha = 0.5; ax2.GridColor = [0.5 0.5 0.5];
    % X axis
    ax2.XTick = [0:30:180]; % or start fform 1 : [1, 20:20:180]
    ax2.XTickLabel = strsplit(num2str(ax2.XTick),' ');
    ax2.XLabel.String = 'Projection angle \theta';
    % Colorbar
    cbar = colorbar(); cbar.Ticks = cticks_RTN; 
    cbar.Label.String = 'value of line integral';
    title('Normalised Radon Transform Rnorm(\rho,\theta)')
    
    % Place axes [X Y W H]
    ax2.Position = [0.1300-0.03    0.1100+0.03    0.7750-0.02    0.3412];
    ax1.Position = [0.1300-0.03    0.5835-0.03    0.7750-0.02    0.3412];

    % -- save
%     saveas(gcf,[path2save 'normalised_radon_example_' char(CASES(C)) '_plt2b'],'epsc');

%% -- plot response signals and dominant signal extraction    
figure('position',[441   368   534   268]); 

    % -- calculate signal response 
    stdRad = std(RT,0,'omitnan');
    stdRad_norm = std(RTnorm,0,'omitnan'); 
    medRad = median(cat(1, circshift(stdRad,[0 -1]), stdRad, circshift(stdRad,[0 1])),1);  % running median filter
    medRad_norm = median(cat(1, circshift(stdRad_norm,[0 -1]), stdRad_norm, circshift(stdRad_norm,[0 1])),1);  % running median filter
    % -- dominant signal
    [pks,loc] = findpeaks(medRad, 'minpeakdistance', 10, 'sortstr', 'descend');
    [pks_norm,loc_norm] = findpeaks(medRad_norm, 'minpeakdistance', 10, 'sortstr', 'descend');
    
    % -- plot response for every theta
    yyaxis left
        p5 = plot(1:180,medRad,'color',RT_color,'linewidth',3,'displayname', ['\sigma(RT)']);
        hold on
        r2 = scatter(loc(1),pks(1),40,'k','o','linewidth',3,'displayname', '\sigma_{max} (RT)');
        % Y axis: even spacing
        ax = gca;
        ax.YTick = linspace(min(ax.YTick),max(ax.YTick),4);
        ax.YAxis(1).TickLabelFormat='%.1f';
        ax.YLabel.String = 'value of line integral';
    yyaxis right
        if ~isempty(pks_norm) 
            p6 = plot(1:180,medRad_norm,'color',RTN_color,'linewidth',3,'displayname', ['\sigma(RTnorm)']);
            hold on
            r3 = scatter(loc_norm(1),pks_norm(1),40,'k',RTN_mark,'linewidth',3,'displayname', '\sigma_{max}(RTnorm)');
        else 
            r3 = scatter(90,medRad_norm(90),1,'k',RTN_mark,'linewidth',3,'displayname', ['\sigma_{max}(RTnorm)'])    ;
            hold on
            p6 = plot(1:180,medRad_norm,'-','color',RTN_color,'linewidth',3,'displayname', ['RTnorm']);
        end
        legend_handle = legend([p5, p6, r2, r3],'location','bestoutside');
        % Y axis: even spacing
        ax = gca;
        ax.YTick = linspace(min(ax.YTick),max(ax.YTick),4);
        ax.YAxis(2).TickLabelFormat='%.2f';
        
    grid on
    
    % Y axis colors
    ax = gca;
    ax.YAxis(1).Color = RT_color;
    ax.YAxis(2).Color = RTN_color;

    % X axis
    ax.XLim = [0 180];
    ax.XTick = [0:30:180];
    ax.XTickLabel = strsplit(num2str(ax.XTick),' ');
    ax.XLabel.String = 'Projection angle \theta';
    set(ax,'fontsize',fontSize); set(legend_handle,'fontsize',fontSize);
    
    % Place axes [X Y W H] same width as fig2 
%     ax.Position = [0.1300-0.03    0.1100+0.03    0.7750-0.02    0.8150];
    
%     saveas(gcf,[path2save 'normalised_radon_example_' char(CASES(C)) '_plt3'],'epsc');
    
end



%% EXAMPLE SCHEMATIC: example window and two example projections

path2save = './example/';
TH = [45,90];

RT_color = batlow(15,:);
RTN_color = batlow(200,:);

CASES = ["uniform";"gray1";"gray2";"complex1";"complex2"];

CASES="gray2"; 

for C = 1:length(CASES)
   
switch CASES(C)
    case "uniform"
        window = ones(9,9);
    case "gray1"
        window = zeros(9,9);
        window(4,:) = 0.5;
    case "gray2"
        window = zeros(9,9);
        window(4,:) = 0.5;
        window = window+0.5;
        window(window>1) = 1;
    case "complex1"
        window = zeros(9,9);
        window(6,:)=1; 
        xy = [1,7;2,8;3,9;];
        for i = 1:3
            window(xy(i,1),xy(i,2))=0.5;
        end
    case "complex2"
        window = zeros(9,9);
        window(6,:)=0.5;
        xy = [1,7;2,8;3,9;];
        for i = 1:3
            window(xy(i,1),xy(i,2))=1;
        end
        
end

% -- calc RT
RT = radon(window);
RTnorm = radoNorm(window); 
    

figName = ['normalised_radon_schematic_' char(CASES(C))];

% -- IMG
figure('position',[ 441   377   560   420])
    imagesc(window), axis equal, axis tight, axis off
    colormap('gray'), caxis([0,1])
    ax1 = gca; 
    ax1.XTick = 1:1:size(window,1);
    ax1.YTick = 1:1:size(window,2);
    ax2 = axes('Position',ax1.Position,...
      'Color','none', 'Ylim',ax1.YLim,'XLim',ax1.XLim,...
      'TickLength',[0 0],...
      'YTick', get(ax1,'ytick')-0.5,... 
      'XTick', get(ax1,'xtick')-0.5,  ...
      'YTickLabel', [],'XTickLabel', [], ...
      'gridcolor','r','linewidth',2);
    axis equal, axis tight 
    ax2.YAxis.Visible = 'off'; ax2.XAxis.Visible = 'off';
    ax2.XLim = ax2.XLim+0.5; ax2.YLim = ax2.YLim+0.5;
    grid on
    linkaxes([ax1 ax2],'xy')
%     saveas(gcf,[path2save figName '_window'],'epsc');
%     close(gcf)


% --- RT FOR SELECTED THETA
    for theta = TH
    

        % shift RTnorm to be better aligned to RT when plotted

        RTnorm_plt = RTnorm;
        while size(RT,1) > size(RTnorm_plt,1)+1 
            dx = size(RT,1) - size(RTnorm,1) - 1;
            padd = nan*ones(1,180);
            if theta == 90
                RTnorm_plt = cat(1,padd,RTnorm_plt);
            elseif theta == 45
                RTnorm_plt = cat(1,RTnorm_plt,padd);
                RTnorm_plt = cat(1,padd,RTnorm_plt);
            end
            
        end


    figure('position',[ 441   527   560   270])
    yyaxis left %  RT
        p5 = plot(RT(:,theta),'color',RT_color,'linewidth',3,'displayname', ['RT(' num2str(theta) ')']);
        % Y axis: even spacing
        ax = gca;
        ax.YTick = linspace(min(ax.YTick),max(ax.YTick),3);
        ax.YAxis(1).TickLabelFormat='%.1f';
        ax.YAxis(1).FontSize=17;
        ax.Color = 'none';
        box off
        ax.FontSize = 14;
        ax.XLim = [1 size(RT,1)];
        ax.YLabel.String = 'R';
    yyaxis right % RT norm
        p6 = plot(RTnorm_plt(:,theta),'color',RTN_color,'linewidth',3,'displayname', ['RTnorm(' num2str(theta) ')']);
        ax = gca;
        ax.YTick = linspace(min(ax.YTick),max(ax.YTick),3);
        ax.YAxis(2).TickLabelFormat='%.2f';
        ax.YAxis(2).FontSize=17;
        ax.Color = 'none';
        box off
        ax.FontSize = 14;
        ax.XLim = [1 size(RT,1)];
        ax.XLabel.String = 'Projection axis \rho';
        ax.YLabel.String = 'Rnorm';
        
    title(['(Normalised) Radon Transform for theta=' num2str(theta)])
        
%         saveas(gcf,[path2save figName '_RT' num2str(theta)],'epsc');
        
    end

end 
