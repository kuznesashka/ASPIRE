function [IndMax, ValMax, ind_m, spikeind] = dip_fit(spike_ind, Data, G3, channels, channel_idx, t1, t2, t3, MRI, cortex)

Voxels = cs_convert(MRI, 'scs', 'voxel', cortex.Vertices);
IndMax      = zeros(spike_ind)
ValMax      = zeros(spike_ind)
ind_m       = []
spikeind    = zeros(spike_ind)
% Convert head model
ftHeadmodel = out_fieldtrip_headmodel(G3, channels, channel_idx);

% Convert data file
ftData = out_fieldtrip_data(Data, channels, channel_idx, 0);

% Generate rough grid for first estimation
% ??? GridOptions = bst_get('GridOptions_dipfit');
% ??? GridLoc = bst_sourcegrid(GridOptions, G3.SurfaceFile);
DipoleModel = 'moving';
NumDipoles = 1;
TimeWindow = [1]; % samples (ms)
SymmetryConstraint = [];
% Prepare FieldTrip cfg structure
cfg = [];
cfg.channel     = {channels.Channel(channel_idx).Name};
cfg.headmodel   = ftHeadmodel;
cfg.latency     = TimeWindow;
cfg.numdipoles  = NumDipoles;
cfg.model       = DipoleModel;
cfg.nonlinear   = 'yes';
% cfg.grid.pos    = GridLoc;
% cfg.grid.inside = ones(size(GridLoc,1),1);
% cfg.grid.unit   = 'm';
cfg.grid.resolution = 1;
cfg.grid.unit   = 'cm';
cfg.symmetry    = SymmetryConstraint;
cfg.feedback    = 'textbar';
cfg.gridsearch  = 'yes';
cfg.senstype    = 'MEG';

if exist('fminunc', 'file')
        cfg.dipfit.optimfun = 'fminunc';
else 
        cfg.dipfit.optimfun = 'fminsearch';
end

for i = 1:length(spike_ind)
    t = spike_ind(i)

    if t1 == 0
        timelockeds.time   = [1];
        timelockeds.avg    = ftData.trial{1}(:, t);
    else
        avr1 = t - (t3-t1);
        avr2 = t - (t3-t2);
        timelockeds.time   = [1];
        timelockeds.avg    = mean(ftData.trial{1}(:, avr1:avr2),2);
    end

    timelockeds.label  = ftData.label;
    timelockeds.dimord = 'chan_time';
    timelockeds        = copyfields(ftData, timelockeds, {'grad', 'elec', 'opto', 'cfg', 'trialinfo', 'topo', 'topodimord', 'topolabel', 'unmixing', 'unmixingdimord'});
        

    ftDipole = ft_dipolefitting(cfg, timelockeds);
    nTime = length(ftDipole.time);
    switch (DipoleModel)
    case 'moving'
                dipPos = cat(1, ftDipole.dip.pos)';
                dipMom = reshape(cat(2, ftDipole.dip.mom), 3, []);
                dipRv  = cat(2, ftDipole.dip.rv);
    case 'regional'
                dipPos = repmat(ftDipole.dip.pos, nTime, 1)';
                dipMom = reshape(ftDipole.dip.mom, 3, []);
                dipRv  = ftDipole.dip.rv;
    end
    Goodness  = 1 - dipRv;
    mri_cord = cs_convert(MRI, 'scs', 'voxel', dipPos);
    %compute Euclidean distances:
    distances =  sum((Voxels - mri_cord) .^ 2, 2);
    %find the smallest distance:
    IndMax(i)   = find(distances == min(distances));
    ValMax(i)   = Goodness
    spikeind(i) = int64(mean([avr1 avr2]))
 
end 

if corr_thresh ~= 0.0
%     quant = 0.95;
    quant = corr_thresh;
else
    quant = 0.0;
end

% visual detection - corr_thresh
% ICA and SPC detection - prctile(ValMax,corr_thresh_prctile);
corr_thresh = quantile(ValMax, quant);
ind_m = find((ValMax > corr_thresh));
% channel_type_loop
disp(['Subcorr threshold: ', num2str(corr_thresh), ' Number of spike found: ', ...
    num2str(length(ind_m))]);


