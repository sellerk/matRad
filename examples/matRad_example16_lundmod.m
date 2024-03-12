clear, clc

%% minimal beispiel
load \\david.lse.thm.de\Temp\jessica\phantom\PL666931029189007_60_depth.mat


% machine
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'Generic';
% beam geometry settings
pln.numOfFractions  = 1;
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [270]; % [?]
pln.propStf.couchAngles     = [0]; % [?]

pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below





%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% calc dose
dij = matRad_calcParticleDose(ct,stf,pln,cst);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);


%% recalcs with lungmodulation
%% lungmodulation
pln.propOpt.lungModulation = true;
pln.propOpt.ModulationPower = 800;
resultGUI_lungmod = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
resultGUI.physicalDose_conv800 = resultGUI_lungmod .physicalDose;

matRad_compareDose(resultGUI.physicalDose,resultGUI_lungmod.physicalDose,ct,cst,[1 1 0]);




%% lungmodulation with CT randomization

if ~exist('modulation', "var")
    disp('perform matRad_calcDoseInit to obtain modulation struct')
    return
else
    pln.propDoseCalc.useGivenEqDensityCube = true;
    ct = matRad_calcWaterEqD(ct, pln);
    modulate_CT
end

tempstore_ct = ct;

pln.propDoseCalc.useGivenEqDensityCube = true;
pln.propOpt.lungModulation = false;
dose_sum = zeros(size(resultGUI.physicalDose));


for i = 1:size(modulation.modCube,1)
    disp(i)
    ct.cube{1} = modulation.modCube{i};
    resultGUI_recalc = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
    dose_sum = dose_sum + resultGUI_recalc.physicalDose;
end

dose_sum_div = dose_sum./size(modulation.modCube,1);
% resultGUI.physicalDose_sum = dose_sum;
resultGUI.physicalDose_CTrand_800 = dose_sum_div;
ct = tempstore_ct;
clearvars tempstore_ct

matRadGUI


matRad_compareDose(resultGUI.physicalDose_CTrand_800,resultGUI_lungmod.physicalDose,ct,cst,[1 1 0]);
matRad_compareDose(resultGUI.physicalDose_CTrand_800,resultGUI.physicalDose,ct,cst,[1 1 0]);