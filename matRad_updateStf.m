function stf = matRad_updateStf(stf, pln)
% matRad steering information update if FWHM is not added
% 
% call
%   stf = matRad_updateStf(stf, plan)
%
% input
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%
% output
%   stf:        updated matRad steering information struct
%
% References
%   -
%

fileName = [pln.radiationMode '_' pln.machine];
try
   load([fileparts(mfilename('fullpath')) filesep 'basedata' filesep fileName]);
catch
   matRad_cfg.dispError('Could not find the following machine file: %s',fileName); 
end

for i  = size(stf,2)
    for j = stf(i).numOfRays:-1:1
        
        [~, vEnergyIx] = min(abs(bsxfun(@minus,[machine.data.energy]',...
            repmat(stf(i).ray(j).energy,length([machine.data]),1))));
        for k = 1:stf(1).numOfBixelsPerRay(j)
            stf(i).ray(j).focusFWHM(k) = machine.data(vEnergyIx(k)).initFocus.SisFWHMAtIso;
        end
    end
end