function plot_clusters(Data, channel_type, spikes_extraction, ...
    f_low, f_high, cortex, spike_trials, maxamp_spike_av, spike_ts, cluster, ...
    channels, G3, MRI, corr_thresh, save_clust, save_path, mute)
%
%  spike_ind, picked_components, picked_comp_top, spike_sources, ...
%     bf_ts, corr_thresh, hemi)

% -------------------------------------------------------------------------
% Visualization of automatic spike detection procedure
% -------------------------------------------------------------------------
% FORMAT:
%   epi_plot(cluster)
% INPUTS:
%   cluster -- structure with all detected clusters
%
%
% NOTE:
%
% _________________________________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com


% 1. Data is filtered for vizualization
Fs = 1/(Data.Time(2)-Data.Time(1)); % sampling frequency
R = G3.GridLoc; % locations of sources

% channel indices
if strcmp(channel_type, 'grad') == 1
    grad_idx = setdiff(1:306, 3:3:306);
    channel_idx = grad_idx(Data.ChannelFlag(grad_idx)~=-1);
elseif strcmp(channel_type, 'mag') == 1
    magn_idx = 3:3:306;
    channel_idx = magn_idx(Data.ChannelFlag(magn_idx)~=-1);
end

% channel names
k = 1;
for i = 1:length(channel_idx)
    namechan{i} = channels.Channel(channel_idx(i)).Name;
    k = k+1;
end

[b,a] = butter(4, [f_low f_high]/(Fs/2)); % filtering for visualization
Ff = filtfilt(b, a, Data.F(channel_idx,:)')';

% MAIN FIGURE
if mute
    h = figure('visible','off');
else
    h = figure('visible','on');
end
% top brain view
subplot(3,4,[5:6,9:10])
cortex_lr = cortex;
cortex_hr = cortex;
data_lr = ones(length(cortex_lr.Vertices),1);
mask_lr = zeros(size(data_lr));
plot_brain_cmap2(cortex_lr, cortex_lr, [], data_lr, ...
    mask_lr, 0.05)
axis equal
grid off
axis off
hold on
view(270, 90)

c = hsv(size(cluster,2));
for i = 1:length(cluster)
    ind = cluster{1,i}(1,:);
    scatter3(cortex.Vertices(ind,1), cortex.Vertices(ind,2), ...
        cortex.Vertices(ind,3), 150, 'filled', 'MarkerEdgeColor','k',...
        'MarkerFaceColor',c(i,:));
    %lgd = legend;
    %if size(cluster,2)>35
    %    %lgd.NumColumns = 2; %round(size(cluster,2)/35)
    %end
end

% bars for clicking
subplot(3,4,7)
        %scatter(1*ones(size(cluster,2),1),[1:size(cluster,2)],50,c,'filled')
        %text(1+0.1, 1, '1');
        %text(1+0.1, size(cluster,2), num2str(size(cluster,2)));
set(groot,'defaultAxesColorOrder',c) %colormap(c);
x = ones(2, size(cluster,2));
ff = area(x);
axis off
grid off
ylim([1, size(cluster,2)])
set(ff,'ButtonDownFcn', @cluststatistics, 'HitTest','on')

% plots for statistics
h1 = subplot(3,4,1);
title('Distribution of subcorrs')
axis equal
grid off
set(gca,'fontsize', 14)
h2 = subplot(3,4,2);
title('Distribution of events in time')
axis equal
grid off
set(gca,'fontsize', 14)
h3 = subplot(3,4,3);
title('Events timeseries (sensors)')
axis equal
grid off
set(gca,'fontsize', 14)
h4 = subplot(3,4,4);
title('Events timeseries (sources)')
axis equal
grid off
set(gca,'fontsize', 14)
subplot(3,4,7)
axis off
h5 = subplot(3,4,8);
axis equal
grid off
h6 = subplot(3,4,11);
axis equal
grid off
h7 = subplot(3,4,12);
axis equal
grid off


% button = uicontrol('Style', 'pushbutton',...
%     'String', 'Show timeseries',...
%     'Position', [300 15 310 30],...
%     'Callback', @(source,event)plotMEG(h),...
%     'FontSize', 14);

    function cluststatistics(source, event)
        
        loc = event.IntersectionPoint;
        clust_num = ceil(loc(2));

        % Distribution of subcorrs inside the cluster
        cla(h1)
        h1 = subplot(3,4,1);
        histogram(cluster{1,clust_num}(3,:), 'EdgeColor', 'k', 'FaceColor', ...
            c(clust_num,:), 'BinWidth', 0.001)
        xlim([corr_thresh-0.01 1])
        title('Distribution of subcorrs')
        set(gca,'fontsize', 14)
        
        % Distribution of this cluster events in time
        cla(h2)
        h2 = subplot(3,4,2);
        stem(cluster{1,clust_num}(2,:), ones(1, length(cluster{1,clust_num}(2,:))), ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', c(clust_num,:))
        title('Distribution of events in time')
        xlim([0 size(Ff, 2)])
        set(gca,'fontsize', 14)
        
        % Average spike on sensors
        cla(h3)
        h3 = subplot(3,4,3);
        plot(-40:80, maxamp_spike_av{clust_num}', 'Color', c(clust_num,:), 'LineWidth', 2)
        hold on
        xlim([-40 80])
        title('Average spike on 5 top amplitude sensors')
        set(gca,'fontsize', 14)
        
        % Source activations timeseries
        cla(h4)
        h4 = subplot(3,4,4);
        plot(-40:80, mean(spike_ts{clust_num}, 1), 'Color', c(clust_num, :), 'LineWidth', 2)
        title('Events timeseries (sources)')
        xlim([-40 80])
        setappdata(h, 'cluster', clust_num)
        set(gca,'fontsize', 14)
        
        
        coord_scs = R(cluster{clust_num}(1,:),:);
        coord_mri = cs_convert(MRI, 'scs', 'voxel', coord_scs);
        mnslice = round(mean(coord_mri, 1));
        coord_mri = round(coord_mri);
        
        cla(h5)
        h5 = subplot(3,4,8);
        mri = flipud(MRI.Cube(:,:,mnslice(3))');
        imagesc(mri)
        colormap('gray')
        axis equal
        grid off
        axis off
        hold on
        idx = coord_mri(:,1:2);
        scatter(idx(:,1), 257-idx(:,2),  50, 'filled', 'MarkerFaceColor', ...
            c(clust_num,:), 'MarkerEdgeColor', 'k');
        
        cla(h6)
        h6 = subplot(3,4,11);
        mri = flipud(squeeze(MRI.Cube(:,mnslice(2),:))');
        imagesc(mri)
        colormap('gray')
        axis equal
        grid off
        axis off
        hold on
        idx = coord_mri(:,[1,3]);
        scatter(idx(:,1), 257-idx(:,2), 50, 'filled', 'MarkerFaceColor', ...
            c(clust_num,:), 'MarkerEdgeColor', 'k');
        
        cla(h7)
        h7 = subplot(3,4,12);
        mri = flipud(fliplr(squeeze(MRI.Cube(mnslice(1),:,:))'));
        imagesc(mri)
        colormap('gray')
        hold on
        idx = coord_mri(:,[2,3]);
        scatter(257-idx(:,1), 257-idx(:,2), 50, 'filled', 'MarkerFaceColor', ...
            c(clust_num,:), 'MarkerEdgeColor', 'k');
        axis equal
        grid off
        axis off
        
        
    end



if save_clust
    for     clust_num = 1:length(cluster)%ceil(loc(2));
        
        % Distribution of subcorrs inside the cluster
        cla(h1)
        h1 = subplot(3,4,1);
        histogram(cluster{1,clust_num}(3,:), 'EdgeColor', 'k', 'FaceColor', ...
            c(clust_num,:), 'BinWidth', 0.001)
        xlim([corr_thresh-0.01 1])
        title('Distribution of subcorrs')
        set(gca,'fontsize', 14)
        
        % Distribution of this cluster events in time
        cla(h2)
        h2 = subplot(3,4,2);
        stem(cluster{1,clust_num}(2,:), ones(1, length(cluster{1,clust_num}(2,:))), ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', c(clust_num,:))
        title('Distribution of events in time')
        xlim([0 size(Ff, 2)])
        set(gca,'fontsize', 14)
        
        % Average spike on sensors
        cla(h3)
        h3 = subplot(3,4,3);
        plot(-40:80, maxamp_spike_av{clust_num}', 'Color', c(clust_num,:), 'LineWidth', 2)
        hold on
        xlim([-40 80])
        title('Average spike on 5 top amplitude sensors')
        set(gca,'fontsize', 14)
        
        % Source activations timeseries
        cla(h4)
        h4 = subplot(3,4,4);
        plot(-40:80, mean(spike_ts{clust_num}, 1), 'Color', c(clust_num, :), 'LineWidth', 2)
        title('Events timeseries (sources)')
        xlim([-40 80])
        setappdata(h, 'cluster', clust_num)
        set(gca,'fontsize', 14)
        
        
        coord_scs = R(cluster{clust_num}(1,:),:);
        coord_mri = cs_convert(MRI, 'scs', 'voxel', coord_scs);
        mnslice = round(mean(coord_mri, 1));
        coord_mri = round(coord_mri);
        
        cla(h5)
        h5 = subplot(3,4,8);
        mri = flipud(MRI.Cube(:,:,mnslice(3))');
        imagesc(mri)
        colormap('gray')
        axis equal
        grid off
        axis off
        hold on
        idx = coord_mri(:,1:2);
        scatter(idx(:,1), 257-idx(:,2),  50, 'filled', 'MarkerFaceColor', ...
            c(clust_num,:), 'MarkerEdgeColor', 'k');
        
        cla(h6)
        h6 = subplot(3,4,11);
        mri = flipud(squeeze(MRI.Cube(:,mnslice(2),:))');
        imagesc(mri)
        colormap('gray')
        axis equal
        grid off
        axis off
        hold on
        idx = coord_mri(:,[1,3]);
        scatter(idx(:,1), 257-idx(:,2), 50, 'filled', 'MarkerFaceColor', ...
            c(clust_num,:), 'MarkerEdgeColor', 'k');
        
        cla(h7)
        h7 = subplot(3,4,12);
        mri = flipud(fliplr(squeeze(MRI.Cube(mnslice(1),:,:))'));
        imagesc(mri)
        colormap('gray')
        hold on
        idx = coord_mri(:,[2,3]);
        scatter(257-idx(:,1), 257-idx(:,2), 50, 'filled', 'MarkerFaceColor', ...
            c(clust_num,:), 'MarkerEdgeColor', 'k');
        axis equal
        grid off
        axis off
        
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,[save_path 'Cluster_' spikes_extraction '_' channel_type '_' num2str(clust_num) '.bmp'])
        
        
    end
end

end