function TOPO = plot_spikes_ER_TOPO(spikes_fitted, path, Data, channel_type, f_low, f_high,  G3,  channels )

 

%% -------------------------------------------------------------------------
% Plot separate spikes

% assuming to have loaded brainstorm and fieldtrip
% brainstorm_190129\brainstorm3';
%  fieldtrip-20190203';
% -------------------------------------------------------------------------
% INPUTS:
%   spikes_fitted -- [N spikes x 5] [rec_time_seconds, detection sample, subspace corr,  leadfield_vertex_index, cluster]
%   Data -- brainstorm structure with artifact corrected maxfiltered MEG data
%   G3 -- brainstorm structure with forward operator
%   channel_type -- channels used ('grad' or 'mag')
%   f_low, f_high -- the bands for prefiltering before fitting
%   channels -- channels info
% _______________________________________________________
% Tommaso Fedele fedele.tm@gmail.com

spike_sec = spikes_fitted(:,1); % spikes time stamps in seconds in the recording time
spike_time = spikes_fitted(:,2); % spikes time stamps in milliseconds in the file time (sample/fsampe*1000)
gof             = spikes_fitted(:,3); % subspace correlation from MUSIC
cl                = spikes_fitted(end,end); % cluster

% prepare the data struct in fieldtrip format
data              = load('data_megin'); %default struct
data              = data.data_megin;
grad_idx      = setdiff(1:306, 3:3:306);
magn_idx    = 3:3:306;
layout_all     = 'neuromag306all.lay';
layout_cmb = 'neuromag306cmb.lay';
layout_mag = 'neuromag306mag.lay';

TOPO = []; % collect all topgraphies as columns

% channels selection
if strcmp(channel_type, 'grad') == 1
    channel_idx = grad_idx(Data.ChannelFlag(grad_idx)~=-1);
elseif strcmp(channel_type, 'mag') == 1
    channel_idx = magn_idx(Data.ChannelFlag(magn_idx)~=-1);
end

% channels: gradiometers or magnetometers
if strcmp(channel_type, 'grad') == 1
    grad_idx = setdiff(1:306, 3:3:306);
    channel_idx = grad_idx(Data.ChannelFlag(grad_idx)~=-1);
elseif strcmp(channel_type, 'mag') == 1
    magn_idx = 3:3:306;
    channel_idx = magn_idx(Data.ChannelFlag(magn_idx)~=-1);
end



% filtering
Fs = 1/(Data.Time(2)-Data.Time(1));
[b,a] = butter(4, [f_low f_high]/(Fs/2)); % butterworth filter before ICA
if 1 %for debug mode  == 0
    Ff = filtfilt(b, a, Data.F(channel_idx,:)')';
else
    Ff =    Data.F(channel_idx,:) ;
end

% svd forward model to 2 dimension
[G2, ~] = G3toG2(G3, channel_idx);

% labels for ERplot and TOPOplot
clear namechan*
k = 1;
for i = 1:length(channel_idx)
    %     namechan{i} = channels.Channel(channel_idx(i)).Name;
    namechanER{i} = channels.Channel(channel_idx(i)).Name;
    k = k+1;
end

k = 1;
for i = 1:length(magn_idx)
    %     namechan{i} = channels.Channel(channel_idx(i)).Name;
    namechanTOPO{i} = channels.Channel(magn_idx(i)).Name;
    k = k+1;
end


for i = 1:length(spike_time)
    
    %%  Calculate spike_ts
    spike = Ff(:,(spike_time(i)-20):(spike_time(i)+30)); % 50 ms for
    [U,S,V] = svd(spike);
    h = cumsum(diag(S)/sum(diag(S)));
    n = find(h>=0.95);
    corr = MUSIC_scan(G2, U(:,1:n(1)));
    [~, dip_ind] = max(corr); % vertex of max correlation
    G3.GridLoc(dip_ind,:)
    spike = Ff(:, (spike_time(i)-100):(spike_time(i)+100));
    g = G2(:,(dip_ind*2-1):dip_ind*2);
    [U,S, ~] = svd(spike);
    h = cumsum(diag(S)/sum(diag(S)));
    n = min(find(h>=0.95));
    [u, s, v] = svd(U(:,1:n)'*g);
    g_fixed = g*v(1,:)';
    spike_ts = spike'*g_fixed;
    
    %%  Plot
    figure('visible','off');
    
    subplot(233)
    
    plot(-100:100, spike_ts, 'LineWidth', 2);
    title(['Time:' num2str(spike_sec(i)) 's, corr: ' num2str(gof(i))]);
    ylim([min(spike_ts), max(spike_ts)])
    xlim([-100, 100])
    
    % ft_multiplotER
    subplot(2,3,[1:2,4:5])
    cfg = [];
    if strcmp(channel_type, 'grad') == 1
        cfg.layout = 'neuromag306all.lay';
        
    else
        cfg.layout =  'neuromag306mag.lay';
    end
    
    data.avg = spike;
    data.var = spike;
    data.dof = repmat(72, size(data.avg));
    data.label = namechanER';
    data.time = -100:100;
    
    cfg.parameter = 'avg';
    ft_multiplotER(cfg, data)
    
    
    subplot(2,3,6)
    vector =  data.avg(:,101);
    plot_topography(vector, channel_type, data.label)
    
    
    % collect all topographies
    TOPO = [TOPO data.avg(:,101)];
    
    
    colormap('jet')
    
    % screen settings
    %     set(gcf, 'Position', get(0, 'Screensize'),'PaperOrientation','portrait');
    set(gcf, 'Position', [3         779        2876        1634],'PaperOrientation','portrait');
    
    %%   save plot
    
    saveas(gcf,[path 'spike_' num2str(spike_time(i)/1000) '_s.png']);
    close(gcf);
end


figure,
plot_topography(mean(TOPO,2), channel_type, data.label)

set(gcf, 'Position', [3         779        2876        1634],'PaperOrientation','portrait');
saveas(gcf,[path 'average pattern spike.png']);


    
end


% %% Find nearest spike function
% % time - column of timestamps
% % value - one time point
% function [nearest, index] = find_nearest(stamps, value)
% [~,I] = min(abs(stamps-value));
% nearest = stamps(I);
% index = I;
% end



%%
% 
%  cfg = [];
%     if strcmp(channel_type, 'grad') == 1
%         % ft_combineplanar
%         cfg                 = [];
%         cfg.feedback        = 'no';
%         cfg.method          = 'template';
%         cfg.layout = layout_plan;
%         cfg.neighbours      = ft_prepare_neighbours(cfg, data);
%         cfg.planarmethod    = 'sincos';
%         data_planar        = ft_megplanar(cfg, data);
%         
%         cfg                 = [];
%         cfg.combinemethod  = 'sum';
%         data_planar_comb = ft_combineplanar(cfg,data_planar);
%         data = data_planar_comb;
%         % layout
%         cfg.layout = layout_plan;
%     else
%         % layout
%         cfg.layout = layout_mag;
%     end
%     cfg.xlim = [0 0];
%     cfg.parameter = 'avg';
%     cfg.comment = 'xlim';
%     data.dimord = 'chan_time';
%     ft_topoplotER(cfg, data); %