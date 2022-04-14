function [stf, pln, beamweight] = matrad_PBP2stf(ct, pln, PhysicalBeamPlan)
% matRad function to import a matRad stf struct from dicom RTPLAN data
% 
% call
%   [stf, pln] = matrad_PBP2stf(ct, pln, PhysicalBeamPlan)
%
% input
%   ct:             ct imported by the matRad_importDicomCt function
%   pln:            matRad pln struct with meta information
%   PhysicalBeamPlan:   	PBP as imported by simpleLoadXMLPlan()
%
% output
%   stf             matRad stf struct
%   pln:            matRad pln struct. 
%                   Note: pln is input and output since pln.bixelWidth is 
%                   determined here.
%   beamweight      additional output of the beamweights (assign to
%                   resultGUI.w for recalculation purposes)
%
% References

matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispInfo('matRad: Generating PBP Export...');


fileName = [pln.radiationMode '_' pln.machine];
try
   load([matRad_cfg.matRadRoot filesep 'basedata' filesep fileName]);
   SAD = machine.meta.SAD;
catch
   matRad_cfg.dispError('Could not find the following machine file: %s',fileName); 
end

% Preallocate stf

stf.gantryAngle   = pln.propStf.gantryAngles;
stf.couchAngle    = pln.propStf.couchAngles;
stf.bixelWidth    = pln.propStf.bixelWidth;
stf.radiationMode = pln.radiationMode ;
stf.SAD = SAD;
stf.isoCenter     = pln.propStf.isoCenter;
stf.numOfRays = [];
stf.ray = [];
stf.sourcePoint_bev = [0 -SAD 0];
stf.sourcePoint = [];
stf.numOfBixelsPerRay = [];
stf.totalNumOfBixels = [];

rotMat_vectors_T = transpose(matRad_getRotationMatrix(pln.propStf.gantryAngles,pln.propStf.couchAngles));
stf.sourcePoint = stf.sourcePoint_bev*rotMat_vectors_T;

[uniquePositions, ~, ic]  =unique(PhysicalBeamPlan.Allpoints(:,1:2),'rows');
a_counts = accumarray(ic,1);
numOfuniquePositions = size(uniquePositions,1);
beamweight =[];
% test = (PhysicalBeamPlan.Allpoints(PhysicalBeamPlan.Allpoints(uniquePositions(1,1),1),PhysicalBeamPlan.Allpoints(uniquePositions(1,2),2) ))
for rayloop = 1:numOfuniquePositions
%     disp([uniquePositions(rayloop,1), uniquePositions(rayloop,2)])
    stf.ray(rayloop).rayPos_bev = [uniquePositions(rayloop,1), 0, uniquePositions(rayloop,2)]; % X-, Z-, Y-Position
    stf.ray(rayloop).targetPoint_bev = [2*uniquePositions(rayloop,1), stf.SAD, 2*uniquePositions(rayloop,2)];
    stf.ray(rayloop).rayPos = stf.ray(rayloop).rayPos_bev*rotMat_vectors_T;
    stf.ray(rayloop).targetPoint = stf.ray(rayloop).targetPoint_bev*rotMat_vectors_T;
    stf.ray(rayloop).energy = PhysicalBeamPlan.Allpoints(PhysicalBeamPlan.Allpoints(:,1) == uniquePositions(rayloop,1) & PhysicalBeamPlan.Allpoints(:,2) == uniquePositions(rayloop,2),4)';
    stf.ray(rayloop).focusFWHM = PhysicalBeamPlan.Allpoints(PhysicalBeamPlan.Allpoints(:,1) == uniquePositions(rayloop,1) & PhysicalBeamPlan.Allpoints(:,2) == uniquePositions(rayloop,2),5)';
    stf.ray(rayloop).focusIx = ones(1,size(stf.ray(rayloop).focusFWHM,2));
    stf.ray(rayloop).beamweight = PhysicalBeamPlan.Allpoints(PhysicalBeamPlan.Allpoints(:,1) == uniquePositions(rayloop,1) & PhysicalBeamPlan.Allpoints(:,2) == uniquePositions(rayloop,2),3)'./1e6;
    beamweight = [beamweight; PhysicalBeamPlan.Allpoints(PhysicalBeamPlan.Allpoints(:,1) == uniquePositions(rayloop,1) & PhysicalBeamPlan.Allpoints(:,2) == uniquePositions(rayloop,2),3)];
end

beamweight = beamweight./1e6;
% rotMat_vectors_T = transpose(matRad_getRotationMatrix(pln.propStf.gantryAngles,pln.propStf.couchAngles));
% % 
% 

% 
% % Save ray and target position in lps system.
% stf.ray(j).rayPos      = stf(i).ray(j).rayPos_bev*rotMat_vectors_T;
% stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev*rotMat_vectors_T;

stf.numOfRays = numOfuniquePositions;
stf.numOfBixelsPerRay = [a_counts'];
stf.totalNumOfBixels = [size(PhysicalBeamPlan.Allpoints,1)];


% include rangeshifter data if not yet available
if strcmp(pln.radiationMode, 'protons') || strcmp(pln.radiationMode, 'carbon')
    for j = 1:stf.numOfRays
        for k = 1:numel(stf.ray(j).energy)
%             stf.ray(j).rangeShifter(k).ID = 1;
%             stf.ray(j).rangeShifter(k).eqThickness = 50; %23.3 %(MIT)
%             stf.ray(j).rangeShifter(k).sourceRashiDistance = 7884-300;
            stf.ray(j).rangeShifter(k).ID = 0;
            stf.ray(j).rangeShifter(k).eqThickness = 0;
            stf.ray(j).rangeShifter(k).sourceRashiDistance = 0;
        end
    end
end



matRad_cfg.dispInfo('done \n');


end

