function plot_topography(vector, channel_type, labels)

% requisite: Fieldtrip 1902+

% INPUT
% vector: 306 coefficients for each sensor
% channel_type: 'mag' = MAG, 'grad' = GRAD,
% OUTPUT
% Figure of topographies in 2x2 plot
% MEG2 grad       MEG3 grad
% MEG1 mag      MEG2+3 SVD Combined

% Tommaso Fedele, fedele.tm@gmail.com
% thanks to Lau Andersen http://megnord.org/2018/workshop_material/preprocess_MEG.html

%%

%prepare 
data              = load('data_megin'); %default struct
data              = data.data_megin;
grad_idx      = setdiff(1:306, 3:3:306);
magn_idx    = 3:3:306;
layout_all     = 'neuromag306all.lay';
layout_cmb = 'neuromag306cmb.lay';
layout_mag = 'neuromag306mag.lay';


 
 switch channel_type
    case 'grad', layout        = 'neuromag306mag.lay'
        channel_idx      = setdiff(1:306, 3:3:306);
    case 'mag', layout        = 'neuromag306mag.lay'
        channel_idx    = 3:3:306;
 end
 
 

if size(vector,1) == 1, % must be a column
    vector = vector';
end

data.avg       = repmat(vector,1,10); % 50 ms for ,1,100);
data.var        = data.avg;
data.dof       = repmat(72, size(data.avg));
data.label     = labels;
data.dimord = 'chan_time';
data.time     = 1:size(data.avg,2);

switch channel_type
    case 'grad',
        cfg                 = [];
        cfg.layout      = 'neuromag306all.lay';
        cfg.method = 'svd';
        cfg.neighbours = 'neuromag306cmb_neighb.mat'
        dataComb = ft_combineplanar(cfg,data);
        
        cfg = [];
        cfg.xlim = [1 1];
        cfg.layout  = 'neuromag306cmb.lay';
        cfg.layout = ft_prepare_layout(cfg)
        ft_topoplotER(cfg, dataComb );
         
    case 'mag'
        cfg = [];
        cfg.xlim = [1 1];
        cfg.layout  = 'neuromag306mag.lay';
        cfg.layout = ft_prepare_layout(cfg)
        cfg.channel = 'MEG*1'
        ft_topoplotER(cfg, data);
         
end


end