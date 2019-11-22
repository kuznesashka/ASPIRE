function plot_topography(vector, channel_type, labels)



%% assumes fieldtrip-20190202 or later in the path

% input
% vector: Nch coefficients for each sensor
% channel_type: 1 = MAG, 2 = GRAD,
% labels  : Nch channel labels

% Tommaso Fedele, fedele.tm@gmail.com
% thanks to Lau Andersen http://megnord.org/2018/workshop_material/preprocess_MEG.html

data              = load('data_megin');
data              = data.data_megin;

layout_all     = 'neuromag306all.lay';
layout_cmb = 'neuromag306cmb.lay';
layout_mag = 'neuromag306mag.lay';


data.avg = vector; % 50 ms for ,1,100);
data.var  = vector;
data.dof = repmat(72, size(data.avg));
data.label = labels;
data.dimord = 'chan_time';
data.time = 1:size(data.avg,2);


if strcmp(channel_type,'grad')
    
    cfg                  = [];
    cfg.layout       = 'neuromag306all.lay';
    cfg.method    = 'svd';
    cfg.neighbours = 'neuromag306cmb_neighb.mat'
    dataComb         = ft_combineplanar(cfg,data);
    cfg             = [];
    cfg.layout  = 'neuromag306cmb.lay';
    cfg.layout  = ft_prepare_layout(cfg)
    ft_topoplotER(cfg, dataComb );
    
elseif strcmp(channel_type,'mag')
    
    cfg = [];
    cfg.layout  = 'neuromag306mag.lay';
    cfg.layout = ft_prepare_layout(cfg)
    cfg.channel = 'MEG*1'
    ft_topoplotER(cfg, data);
    
end
end