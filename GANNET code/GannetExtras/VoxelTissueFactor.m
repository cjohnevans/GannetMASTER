function voxelfactor = VoxelTissueFactor(fcsf, fgm, fwm)
% function voxelfactor = VoxelTissueFactor(voxel_composition)
%   Calculate tissue correction factor for a voxel with composition given by voxel_composition
%   to correct Gannet GABAIU values for csf,gm,wm composition of voxel.
%   voxelfactor = Correction value for GABAIU (multiply GABAIU by this number)
%   fcsf = CSF volume fraction
%   fgm = grey matter volume fraction
%   fwm = white matter volume fraction

beta_avg = 0.65; % average water visibility used in Gannet 
beta_csf = 1;
beta_gm = 0.78;
beta_wm = 0.65;

ftot = fcsf + fgm + fwm;

%rescale in case of rounding errors leading to ftot <> 1.0

fcsf = fcsf ./ ftot;
fgm  = fgm ./ ftot;
fwm = fwm ./ ftot;

% if (sum(ftot < 0.99) > 0 )
%     disp('fcsf + fgm + fwm is less than 0.99')
%     return
% elseif (sum(ftot > 1.01 ) > 0)
%     disp('fcsf + fgm + fwm is greater than 1.01')
%     return
% end
%     
% water visibility correction factor
tissuebeta = (beta_csf*fcsf + beta_gm*fgm + beta_wm*fwm);
tissuebeta_corr = tissuebeta ./ beta_avg;

% GABA csf correction
GABAcsf_corr = 1 ./ (fgm + fwm);


voxelfactor = tissuebeta_corr .* GABAcsf_corr;


