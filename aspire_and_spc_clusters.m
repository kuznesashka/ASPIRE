function aspire_and_spc_clusters(path_overlap, saving_path, ...
    channel_type, f_low, f_high, Data, channels, threshold)
%
% -------------------------------------------------------------------------
% Computat a matrix ASPIRE vs Spyking Circus clusters
% -------------------------------------------------------------------------
% INPUTS:
% threshold -- number of spikes in overlap that will be plotted 
% cluster_out_overlap_grad.mat
% cluster_out_overlap_mag.mat
%
%   Data -- brainstorm structure with artifact corrected maxfiltered MEG data
%   channel_type -- channels used ('grad' or 'mag')
%   f_low, f_high -- the bands for prefiltering before fitting
%   channels -- channels info
%_______________________________________________________
%

%% Compute overlap
overlap = load(path_overlap);
overlap = overlap.overlap;
%SPC_mag = load([path_cluster_out 'SpyCir_based_mag.csv']);

overlap(:,8) = overlap(:,8)+1; % Spyking Circus has template with the number 0

aspire_clust = unique(overlap(:,2));
spc_clust = unique(overlap(:,8));

overlap_clusters = zeros(length(aspire_clust),length(spc_clust));

for i = 1:length(overlap(:,2))
    
    ind1 = find(aspire_clust == overlap(i,2));
    ind2 = find(spc_clust == overlap(i,8));
    
    overlap_clusters(ind1, ind2) = overlap_clusters(ind1, ind2) +1;
end

%% Plot heatmap of overlap between clusters
fig = figure('visible','off');
subplot(1,1,1)
h = heatmap(overlap_clusters,'CellLabelColor','none');

h.Title = 'Overlap between Spyking Circus and ASPIRE clusters';
h.XLabel = 'Spyking Circus clusters';
h.XDisplayLabels = spc_clust-1;
h.YDisplayLabels = aspire_clust;
h.YLabel = 'ASPIRE clusters';

saveas(gcf, [saving_path channel_type],'png');


%% Plot topomaps for the higest overlaps
avgFC = load('avgFC.mat');

if strcmp(channel_type, 'grad') == 1
    
    % channel_idx
    grad_idx = setdiff(1:306, 3:3:306);
    cfg = [];
    cfg.layout_mag = 'neuromag306planar.lay';
    layout_plan = 'neuromag306planar.lay';
    channel_idx = grad_idx(Data.ChannelFlag(grad_idx)~=-1);
    
elseif strcmp(channel_type, 'mag') == 1
    magn_idx = 3:3:306;
    cfg = [];
    layout_mag = 'neuromag306mag.lay';
    channel_idx = magn_idx(Data.ChannelFlag(magn_idx)~=-1);
end

Fs = 1/(Data.Time(2)-Data.Time(1));
[b,a] = butter(4, [f_low f_high]/(Fs/2)); % butterworth filter before ICA
Ff = filtfilt(b, a, Data.F(channel_idx,:)')';

k = 1;
for i = 1:length(channel_idx)
    namechan{i} = channels.Channel(channel_idx(i)).Name;
    k = k+1;
end

% avarege the largest overlap
[rows,cols,~] = find(overlap_clusters>=threshold);

if ~isempty(rows)
    for i = 1:length(rows)
        aspire = aspire_clust(rows(i));
        spc = spc_clust(cols(i));
        spike_time = overlap((overlap(:,2)==aspire) & (overlap(:,8)==spc),1);
        
        spike = mean(Ff(:, spike_time),2);
        
        cfg = [];
        if strcmp(channel_type, 'grad') == 1
            cfg.layout = layout_plan;
        else
            cfg.layout = layout_mag;
        end
        plot_topo(spike, aspire, spc, saving_path, channel_type, ...
            namechan, avgFC, cfg, size(spike_time, 1))
    end
    
    
end
end


function plot_topo(spike, aspire, spc,  saving_path, channel_type, ...
                   namechan, avgFC, cfg, n_spikes)
fig = figure('visible','off');
data = avgFC.avgFC;
data.avg = spike;
data.var = spike;
data.dof = repmat(72, size(data.avg));
data.label = namechan';
data.time = 0;

cfg.xlim = [0 0];
cfg.parameter = 'avg';
data.dimord = 'chan_time';
cfg.comment = 'xlim';
ft_topoplotER(cfg,data);

title(['ASPIRE cluster ' num2str(aspire) '; Spyking Circus cluster: ' num2str(spc), ...
       '; Number of spikes: ' num2str(n_spikes)]);

saveas(gcf, [saving_path channel_type '_ASPIRE_' num2str(aspire) '_SPC_' num2str(spc)],'png');

end



