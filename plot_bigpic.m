function plot_bigpic(subj_name, results_saving_path, ...
                     cortex, bigpic_saving_path)
% -------------------------------------------------------------------------
% Plot all cluster all detection types
% -------------------------------------------------------------------------
% INPUTS:
% results_saving_path -
% bigpic_saving_path
%
% OUTPUTS:
% _______________________________________________________
%
%% 1. Export everything from brainstorm

close all
figure('visible','off')

%% 2. Spike detection

for spikes_detection = 1:3
    switch spikes_detection
        case 1 % visual markings
            spikes_extraction = 'visual';
        case 2 % ICA based
            spikes_extraction = 'ICA_based';
        case 3 % Spiking circus based
            spikes_extraction = 'SpyCir_based';
    end
    
    
    for channel_type_loop = [1 2]
        switch channel_type_loop % channels you want to analyse ('grad' or 'mag')
            case 1, channel_type = 'mag'; row_plot = 1;
            case 2, channel_type = 'grad';row_plot = 2;
        end
        
        if exist([results_saving_path spikes_extraction '_' channel_type '.mat'],'file')
            
            load([results_saving_path spikes_extraction '_' channel_type '.mat'],...
                'cluster','parameters')
            
            what = (spikes_detection-1)*6+(channel_type_loop-1)*3;
            switch spikes_detection
                case 1,what = 2;
                case 2,what = 1;
                case 3,what = 0;
            end
            
            cortex_lr = cortex;
            
            if spikes_detection == 1
                c = repmat([0 0 1],100,1);
            else
                
                c = hsv(length(cluster));
            end
            if ((length(cluster) > 1)|(spikes_detection == 1))
                data_lr = ones(length(cortex_lr.Vertices),1);
                mask_lr = zeros(size(data_lr));
                
                %     axes(ha(what+1))%subplot(2,6,what)
                axes('Position',[0+row_plot/8-0.1 0.0+what*.33   .30 .30])
                plot_brain_cmap2(cortex_lr, cortex_lr, [], data_lr, ...
                    mask_lr, 0.05)
                hold on
                for i = 1:length(cluster)
                    ind = cluster{1,i}(1,:);
                    scatter3(cortex.Vertices(ind,1), cortex.Vertices(ind,2), ...
                        cortex.Vertices(ind,3), 20, 'filled', 'MarkerEdgeColor','k',...
                        'MarkerFaceColor',c(i,:));
                end
                title([ spikes_extraction ' ' channel_type ])
                view(270, 90)
                
                
                axes('Position',[0.27+row_plot/10 0.0+what*.33   .30 .30])
                plot_brain_cmap2(cortex_lr, cortex_lr, [], data_lr, ...
                    mask_lr, 0.05)
                hold on
                for i = 1:length(cluster)
                    ind = cluster{1,i}(1,:);
                    scatter3(cortex.Vertices(ind,1), cortex.Vertices(ind,2), ...
                        cortex.Vertices(ind,3), 20, 'filled', 'MarkerEdgeColor','k',...
                        'MarkerFaceColor',c(i,:));
                end
                %title([ spikes_extraction ' ' channel_type ])
                
                view(0, 0)
                
                
                axes('Position',[0.55+row_plot/10  0+what*.33 .30 .30])
                plot_brain_cmap2(cortex_lr, cortex_lr, [], data_lr, ...
                    mask_lr, 0.05)
                hold on
                for i = 1:length(cluster)
                    ind = cluster{1,i}(1,:);
                    scatter3(cortex.Vertices(ind,1), cortex.Vertices(ind,2), ...
                        cortex.Vertices(ind,3), 20, 'filled', 'MarkerEdgeColor','k',...
                        'MarkerFaceColor',c(i,:));
                end
                %title([ spikes_extraction ' ' channel_type ])
                
                view(270, 0)
                
                hold on
                eval(['cluster' num2str(row_plot) ' = cluster;'])
            end
            
        end
        
        
    end
end
%%
set(gcf, 'Position', get(0, 'Screensize'));
xlsTAB = subj_name;
backsalsh = find(xlsTAB=='\');
if backsalsh
    xlsTAB(backsalsh)='_';
end
saveas(gcf, bigpic_saving_path)

end
