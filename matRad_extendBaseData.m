function [Z_adapted, depths_adapted] = matRad_extendBaseData(baseData)
% matRad prolongation of basedata.Z into negative values to account for
% convolution artifacts in the edge areas
% 
% call
%  machine.data(ix) = matRad_extendBaseData(machine.data(ix));
%
% input
%   baseData:           matRad basedata of format machine.data(energyIx)
%   modulationDepth:    current modulation depth in mm
%   Pmod:               Modulation Power in µm
%
% output
%   baseData:           matRad basedata of format machine.data(energyIx)
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%e

% always double the lenth of the basedata into the negative in this way it
% is always clear how many values were taken independent of basedata
% generation
addlength = size(baseData.depths,1);

adddepths = flipud((baseData.depths(2:addlength)*-1));
depths = [adddepths; baseData.depths];
% addZ = (ones(addlength-1,1).*baseData.Z.profileORG(1)); % temp adapt 
addZ = (ones(addlength-1,1).*baseData.Z(1)); 
% Z = [addZ;baseData.Z.profileORG];% temp adapt 
Z = [addZ;baseData.Z];

Z_adapted = Z;
depths_adapted = depths;
% Z_adapted_fft = fft(Z_adapted);
% baseData;

end

