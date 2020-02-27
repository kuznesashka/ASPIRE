function dip_fit_average(Data, evoked_data, G3, channels, channel_idx, ...
    			t1, ~, ~, t4, MRI, cortex, save_evoked)

% -------------------------------------------------------------------------
% Dipole fitting
% -------------------------------------------------------------------------
% INPUTS:
%     Data : brainstorm structure with artifact corrected maxfiltered MEG data
%     G3 : brainstorm structure with forward operator
%     channels : channels info
%     channel_idx : indices of selected sensors
%     t1 : begin of early period in samples(ms)
%     t2 : end of early period in samples(ms)
%     t3 : max peak (ms)
%     MRI : structure from bst
%     cortes : 
%
% OUTPUTS:
%
%     IndMax : locations of fitted dipoles (indices of cortex.Vertices array)
%     ValMax : values of subspace correlation (Goodness)
%     ind_m : indices of spikes survived the subcorr threshold
%     spikeind : time indices. If t1==0 equivalent of spike_ind
%              If t1 ~= 0 spikeind = spike_ind - (t3-t1)
% _______________________________________________________
%
IndMax   = [];
ValMax   = [];
spikeind = [];
Voxels = cs_convert(MRI, 'scs', 'voxel', cortex.Vertices);
% Convert head model
ftHeadmodel = out_fieldtrip_headmodel(G3, channels, channel_idx);

% Convert data file
ftData = out_fieldtrip_data(Data, channels, channel_idx, 0);

% Generate rough grid for first estimation
% ??? GridOptions = bst_get('GridOptions_dipfit');
% ??? GridLoc = bst_sourcegrid(GridOptions, G3.SurfaceFile);

% Prepare FieldTrip cfg structure
DipoleModel = 'moving';
cfg = [];
cfg.channel     = {channels.Channel(channel_idx).Name};
cfg.headmodel   = ftHeadmodel;
cfg.latency     = 1;
cfg.numdipoles  = 1;
cfg.model       = DipoleModel;
cfg.nonlinear   = 'yes';
% cfg.grid.pos    = GridLoc;
% cfg.grid.inside = ones(size(GridLoc,1),1);
% cfg.grid.unit   = 'm';
% cfg.grid.resolution = 2;
cfg.grid.unit   = 'cm';
cfg.symmetry    = [];
cfg.feedback    = 'textbar';
cfg.gridsearch  = 'yes';
cfg.senstype    = 'MEG';

if exist('fminunc', 'file')
        cfg.dipfit.optimfun = 'fminunc';
else 
        cfg.dipfit.optimfun = 'fminsearch';
end

for i = 1:ceil((t4-t1)/10)

    timelockeds.time   = 1;
    timelockeds.avg    = mean(evoked_data.evoked(:,t1+(i-1)*10:t1+i*10),2);
    timelockeds.label  = ftData.label;
    timelockeds.dimord = 'chan_time';
    timelockeds        = copyfields(ftData, timelockeds, {'grad', 'elec', 'opto', 'cfg', 'trialinfo', 'topo', 'topodimord', 'topolabel', 'unmixing', 'unmixingdimord'});
   
    ftDipole = ft_dipolefitting(cfg, timelockeds);
    nTime = length(ftDipole.time);
    switch (DipoleModel)
    case 'moving'
                dipPos = cat(1, ftDipole.dip.pos)';
                % dipMom = reshape(cat(2, ftDipole.dip.mom), 3, []);
                dipRv  = cat(2, ftDipole.dip.rv);
    case 'regional'
                dipPos = repmat(ftDipole.dip.pos, nTime, 1)';
                % dipMom = reshape(ftDipole.dip.mom, 3, []);
                % dipRv  = ftDipole.dip.rv;
    end
    Goodness  = 1 - dipRv;
    mri_cord = cs_convert(MRI, 'scs', 'voxel', dipPos);
    %compute Euclidean distances:
    distances =  sum((Voxels - mri_cord) .^ 2, 2);
    %find the smallest distance:
    IndMax(i)   = find(distances == min(distances));
    ValMax(i)   = Goodness;
    spikeind(i) = t1+(i-1)*10+5;
end

save(save_evoked, 'IndMax','ValMax','spikeind')




