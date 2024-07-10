function modulation = matRad_modulateCT(ct, modulation, varargin)
% matRad dose calculation wrapper bypassing dij calculation
% 
% call
%   modulation = matrad_modulateCT(ct, modulation, num_repetitions, Pmod)
%   modulation = matrad_modulateCT(ct, modulation, 100, pln.propOpt.ModulationPower)    
% input
%   ct:             ct cube          
%   modulation: modulation cube created in matRad_calcDoseInit when
%               modulation calculation is enabled
%  num_repetitions: number of CTs created equals number of repetitions for
%               modulation, default value is 100
%  Pmod         modulation power, currently only 3 Modpowers are supported
%               (250, 450, 800) 
% output
%   modualtion:  modulation cube with randomized CT cubes
%                additionally: some info of the statistics in .metrics
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% run precalculated modulation tables 
matRad_cfg =  MatRad_Config.instance();
load pMod_tables_250_450_800.mat

poisson_dens(2,:) =  1:size(poisson_dens,2);
poisson_dens = fliplr(poisson_dens');

%% hier muss die richtige tabelle ausgewählt werden
if ~nargin
    matRad_cfg.dispInfo(['using default vaules for CT modulation (Pmod = 250 & 100 CTs)', newline]);
    varargin{1,1} = 100;
    varargin{1,2} = 250;
end
if varargin{1,2} == 250
    interp_poission_pMod = interp_poission_pMod250;
elseif varargin{1,2} == 450
    interp_poission_pMod = interp_poission_pMod450;
elseif varargin{1,2} == 800
    interp_poission_pMod = interp_poission_pMod800;
else 
    matRad_cfg.dispInfo(['current modulation power not supported for randomization using default valaue of Pmod = 250\mum', newline]);
end

if varargin{1,1} < 100
    matRad_cfg.dispInfo(['Cave: not tested for repetitions < 100' newline]);
end

% replace minimum value from 0 density to 1e-6 to make sure TOPAS
% simulation ist run properly! 
% is taken care of in precompute_Densfnc just make sure again!
interp_poission_pMod(interp_poission_pMod == 0) = 0.000001;

%% we continue with modulation instead of ct. 
% cube{2} = all voxels inside lung with its final SPR assigned all
% other voxels = 0
% cube{1} = logical mask with 1 inside lung

% define the number of randomizations
num_repetitions = varargin{1,1};
% applying an addition filter for the HU values we are going to modulate =>
% add in varargin later
HU_schwelle = [-900 -100];
modulation.cube{1}(modulation.cubeHU{1}<HU_schwelle(1) | modulation.cubeHU{1}>HU_schwelle(2)) = 0;

% checking for HU values outside the range of [-1000,2995] for the density
% calculation
modulation.cubeHU{2} = ct.cubeHU{1};
modulation.cubeHU{2}(modulation.cubeHU{2} < -1000) = -1000;
modulation.cubeHU{2}(modulation.cubeHU{2} > 2995) = 2995;


%% modulating the ct here
%%%%%%%%%%%%%%% 1 %%%%%%%%%%%%%%%%%%
% convert HU-Value to density according to MC-Calculation
% Implementation of the following formula:
% needed contants are loaded from file. 
%iv:Ge/Patient/SchneiderHounsfieldUnitSections = 9 -1000 -98 15 23 101 2001 2995 2996 3090
%uv:Ge/Patient/SchneiderDensityOffset = 8 0.00121 1.018 1.03 1.003 1.017 2.201 4.54 1
%uv:Ge/Patient/SchneiderDensityFactor = 8 0.001029700665188 0.000893 0.0 0.001169 0.000592 0.0005 0.0 0.0
% uv:Ge/Patient/SchneiderDensityFactorOffset = 8 1000. 0. 1000. 0. 0. -2000. 0. 0.0
% Formula: Density = (Offset + (Factor*(FactorOffset + HU[-1000,2995] ))) * DensityCorrection
% lookup tables for the density calculation from Schneider et a. 
% load C:\Users\Matth\OneDrive\Matlab\_matRad\lungmodulation\density_schneider.mat
% density correction from Schneider et al., values taken from the TOPAS
% calculation input.
% load C:\Users\Matth\OneDrive\Matlab\_matRad\lungmodulation\DensityCorrection.mat 

%% for testing purpose
% testcube = randi([-900 -100], 5,5,3);
% columnwise_CTcubeHU(:,1)  = testcube(:);
%% preallocate
columnwise_CTcubeHU = zeros(numel(modulation.cubeHU{2}),3);
% we use the ct / modulation columnwise to make it faster:
columnwise_CTcubeHU(:,1) = modulation.cubeHU{2}(:);

% use density calculated via HLUT instead of Schneider implementation
columnwise_CTcubeHU(:,2) = ct.cube{1}(:);

% calculation of Density = (Offset + (Factor*(FactorOffset + HU[-1000,2995] ))) * DensityCorrection
% columnwise_CTcubeHU(:,2) = (((interp1(density_schneider(:,1), density_schneider(:,4), columnwise_CTcubeHU(:,1),'previous')...
%     + columnwise_CTcubeHU(:,1) ).*interp1(density_schneider(:,1), density_schneider(:,3), columnwise_CTcubeHU(:,1),'previous') )...
%     + interp1(density_schneider(:,1), density_schneider(:,2), columnwise_CTcubeHU(:,1), 'previous'))...
%     .*DensityCorrection(columnwise_CTcubeHU(:,1)+1001);
% single calcluation steps for check
% Offset(:,1) = interp1(density_schneider(:,1), density_schneider(:,2), HU, 'previous')
% Offset = interp1(density_schneider(:,1), density_schneider(:,2), HU, 'previous')
% Factor(:,1) = interp1(density_schneider(:,1), density_schneider(:,3), HU,'previous')
% Factor = interp1(density_schneider(:,1), density_schneider(:,3), HU,'previous')
% FactorOffset(:,1) = interp1(density_schneider(:,1), density_schneider(:,4), HU,'previous')
% FactorOffset = interp1(density_schneider(:,1), density_schneider(:,4), HU,'previous')
% % Density = (((FactorOffset + HU).*Factor)+Offset).*DensityCorrection;
% ct.dens{1} = reshape(matRad_interp1(hlut(:,1),hlut(:,2),double(ct.cubeHU{1}(:))),ct.cubeDim);
% Density = (((interp1(density_schneider(:,1), density_schneider(:,4), HU,'previous') + HU) .*interp1(density_schneider(:,1), density_schneider(:,3), HU,'previous') )+interp1(density_schneider(:,1), density_schneider(:,2), HU, 'previous'))*DensityCorrection(HU_scaled);

% introduce the scaling factor for the water equivalent Voxel length based
% on the paper from ringbaek et al. 
% d = (ct.resolution.x^(3/2)+ct.resolution.y^(3/2)+ct.resolution.z^(3/2))^-(2/3);
if numel(unique([ct.resolution.x, ct.resolution.y, ct.resolution.z])) ~= 1
   matRad_cfg.dispInfo(['Cave: assymetric Voxelsize, Pmod scaling should be considered' newline]);
end
% scale the density by the Voxelsize
columnwise_CTcubeHU(:,2) = columnwise_CTcubeHU(:,2).*ct.resolution.x;

% get index position for look up in poission distr. 
columnwise_CTcubeHU(:,3) = interp1(poisson_dens(:,2),poisson_dens(:,1),columnwise_CTcubeHU(:,2), 'nearest', 'extrap');

% use index for the column and and random number for the row 
% get a a new density from the resulting poisson dens table 
% columnwise_CTcubeHU(:,4) = interp_poission_pMod(sub2ind(size(interp_poission_pMod),randi(1001,size(columnwise_CTcubeHU,1),1),columnwise_CTcubeHU(:,3)));

%% creating random numbers, masking, reshape 
% % add original Density values outside the lung
modulation.metrics.origCTdens = mean(ct.cube{1}(modulation.cube{1} == 1));
tic
% create a random integer to draw a new density for a given density based
% on the poisson density distribution
% preallocate
columnwise_modulationCube = zeros(numel(modulation.cubeHU{2}),num_repetitions);
for i = 1: num_repetitions
    columnwise_modulationCube(:,i) = interp_poission_pMod(sub2ind(size(interp_poission_pMod),randi(1001,size(columnwise_CTcubeHU,1),1),columnwise_CTcubeHU(:,3)));
end
% due to statistics add a scaling factor if resulting density is not
% equivalent to original density
columnwise_modulationCube= columnwise_modulationCube.*columnwise_CTcubeHU(:,2)./mean(columnwise_modulationCube,2);
% now scale back the density and the voxelsize  
columnwise_modulationCube = columnwise_modulationCube./ct.resolution.x;


% reshape data back to ct cube
for i = 1: num_repetitions
    modulation.modCube{i,1} = reshape(columnwise_modulationCube(:,i),ct.cubeDim);
    % add original Density values outside the lung
    modulation.modCube{i,1}(modulation.cube{1} == 0) = ct.cube{1}(modulation.cube{1} == 0);
    % compute some metrics to check wether the calculations are meaninful
    % modulation.metrics.meandata(:,i) = mean(modulation.modCube{i,1}(modulation.cube{1} == 1));
    modulation.metrics.meandata(:,i) = mean(modulation.modCube{i,1}(modulation.cube{1} == 1));
end
toc

% modulation.modCube{i,1} = reshape(
modulation.metrics.meandens = mean(modulation.metrics.meandata);
modulation.metrics.num_repetitions = num_repetitions;
modulation.metrics.poission = interp_poission_pMod;
% modulation.metrics.Pmod = 250;

