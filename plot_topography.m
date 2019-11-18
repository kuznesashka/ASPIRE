

function plot_topography(vector, channel_type, labels)


% input 
% vector: 306 coefficients for each sensor
% channel_type: 1 = MAG, 2 = GRAD, 
 
%%
% fieldtrip_path = 'E:\_USZserver\Matlab Toolbox\FIELDTRIP\fieldtrip-20190202';
% addpath(genpath(fieldtrip_path))
% ft_defaults

data              = load('data_megin');
data              = data.data_megin;
grad_idx      = setdiff(1:306, 3:3:306);
magn_idx    = 3:3:306;
layout_all     = 'neuromag306all.lay';
layout_cmb = 'neuromag306cmb.lay';
layout_mag = 'neuromag306mag.lay';
 
% spike_time = spikes_fitted(:,2)
% i = 1
% spike = Data.F(1:306,(spike_time(i)-20):(spike_time(i)+30)); % 50 ms for 
%  [U,S,V] = svd(spike);
%  spike = repmat(U(:,2), 1 , 10);
%%

%  % 
%  switch channel_type
%      case 1, layout        = 'neuromag306mag.lay'
%          channel_idx      = setdiff(1:306, 3:3:306);
%      case 2, layout        = 'neuromag306mag.lay'
%          channel_idx    = 3:3:306;          
%  end
 
 channel_idx = 1:306;
 k = 1;
for i = 1:length(channel_idx)
%     namechan{i} = channels.Channel(channel_idx(i)).Name;
     namechan{i} = channels.Channel(channel_idx(i)).Name;
    k = k+1;
end 
 
 
 data.avg = spike; % 50 ms for ,1,100);
 data.var  = spike;
 data.dof = repmat(72, size(data.avg));
 data.label = namechan';
 data.dimord = 'chan_time';
 data.time = 1:size(data.avg,2);
 
 
 %% planar
 cfg                 = [];
% cfg.feedback        = 'yes';
cfg.layout             = 'neuromag306all.lay';
 
 cfg.method = 'svd';
cfg.neighbours = 'neuromag306cmb_neighb.mat'
% avgFICplanar        = ft_megplanar(cfg, data); 


dataComb = ft_combineplanar(cfg,data); 
 

%%
%  cfg = [];
% cfg.layout = 'neuromag306mag.lay';
% cfg.layout = ft_prepare_layout(cfg);
% reading layout from file neuromag306all.lay
% the call to "ft_prepare_layout" took 0 seconds and required the additional allocation of an estimated 7 MB
% figure; ft_plot_layout(cfg.layout);
    close all

figure

zlim = 0.2*[-1 1]
  
subplot(221)
  cfg = []; 
  cfg.xlim = [1 1];
  cfg.zlim = zlim;
 cfg.layout  = 'neuromag306all.lay';  
 cfg.layout = ft_prepare_layout(cfg)
 cfg.channel = 'MEG*2'
 obj = ft_topoplotER(cfg, data);
   text(-1,.51,'MEG*2 grad vert')
  
subplot(222)
  cfg = []; 
  cfg.xlim = [1 1];
   cfg.zlim = zlim;
cfg.layout  = 'neuromag306all.lay';  
 cfg.layout = ft_prepare_layout(cfg)
  cfg.channel = 'MEG*3'
  ft_topoplotER(cfg, data );
   text(-1,.51,'MEG*3 grad hor')

  
subplot(223)
  cfg = []; 
  cfg.xlim = [1 1];
 cfg.layout  = 'neuromag306mag.lay';  
 cfg.layout = ft_prepare_layout(cfg)
  cfg.channel = 'MEG*1'
ft_topoplotER(cfg, data);
  text(-1,.51,'MEG*1 mag')

 
 subplot(224)
  cfg = []; 
  cfg.xlim = [1 1];
     cfg.zlim = zlim;
 cfg.layout  = 'neuromag306cmb.lay';  
 cfg.layout = ft_prepare_layout(cfg)
 ft_topoplotER(cfg, dataComb );
  text(-1,.51,'comb grad SVD')

 
 
 
 
%   figure
% 
%   cfg = [];
% cfg.showlabels = 'no';
%   cfg.layout        = 'neuromag306mag.lay'
%  ft_multiplotER(cfg, data);
 
end