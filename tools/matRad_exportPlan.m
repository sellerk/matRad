function PBP = matRad_exportPlan(stf, pln, beamweight, exportbool, exportformat, varargin)
% call
%   PBP = matRad_exportPlan(stf, pln, resultGUI.w, exportbool, exportformat, filename, pathname)
%
% input
%   stf:           matRad steering information struct
%   pln:           matRad pln struct with meta information
%   beamweight     Voxelweight information (e.g resultGUI.w)
%   exportbool     specify wether output should be created or not
%   exportformat   string specify wether exportformat is '.rst' or '.xml'
%   filename       varargin(1): output filename as string
%   outputdir      varargin(2): output Pathname as string
%
% output
%   PBP:           sorted Struct used to create phys beam plan file
%
%   ToDo:
%   Check if x and Y directions need rotation witch rotation matrix ro be
%   correct
%   Check if multiple Beam export works
%
%% check varargin
if nargin <= 5 || isempty(varargin(1))
    exportname = 'Plan_dummyname_';
elseif ~ischar(varargin{1})
    disp('wrong export pathname format, provide string')
    return
else
    exportname = varargin{1};
end

if ~exportbool || nargin < 5 || (nargin == 6 && isempty(varargin(2))) || ~((size(varargin,2) >0) && ~isfolder(varargin(2)))
    exportpath = cd;
else
    exportpath = varargin{2};
end

% check if weight vector is available, either in function call or in stf - otherwise dose calculation not possible
if isempty(beamweight) && ~isfield([stf.ray],'weight')
     error('No weight vector available. Please provide w or add info to stf')
end

% copy bixel weight vector into stf struct
if ~isempty(beamweight)
    if sum([stf.totalNumOfBixels]) ~= numel(beamweight)
        error('weighting does not match steering information')
    end
    counter = 0;
    for i = 1:size(stf,2)
        for j = 1:stf(i).numOfRays
            for k = 1:stf(i).numOfBixelsPerRay(j)
                counter = counter + 1;
                stf(i).ray(j).weight(k) = beamweight(counter);
            end
        end
    end
else % weights need to be in stf!
    beamweight = NaN*ones(sum([stf.totalNumOfBixels]),1);
    counter = 0;
    for i = 1:size(stf,2)
        for j = 1:stf(i).numOfRays
            for k = 1:stf(i).numOfBixelsPerRay(j)
                counter = counter + 1;
                beamweight(counter) = stf(i).ray(j).weight(k);
            end
        end
    end    
end


% Looping over matRad stf Struct needed because Physical Beamplan xml
% needs energy-sorting
for beamcounter = 1 : size(stf,2)
    loopcounter = 1;
    for rayloop = 1 : size(stf(beamcounter).ray,2)
        for bixelloop = 1 : (stf(beamcounter).numOfBixelsPerRay(rayloop))
            PBP_data{beamcounter}(loopcounter,1) = stf(beamcounter).ray(rayloop).rayPos_bev(1); % X-Position
            PBP_data{beamcounter}(loopcounter,2) = stf(beamcounter).ray(rayloop).rayPos_bev(3); % Y-Position
            PBP_data{beamcounter}(loopcounter,3) = beamweight(loopcounter).*1e6; % beamweight
            PBP_data{beamcounter}(loopcounter,4) = stf(beamcounter).ray(rayloop).energy(bixelloop); % energy
            PBP_data{beamcounter}(loopcounter,5) = stf(beamcounter).ray(rayloop).focusFWHM(bixelloop); % FWHM

            loopcounter = loopcounter+1;
        end
    end
    PBP_help(beamcounter).NumVoxel = loopcounter;
end



%% reshape the data into printable format
for beamcounter = 1 : size(stf,2)
    [PBP_help(beamcounter).Energies, ia]= unique(PBP_data{1,beamcounter}(:,4));
    PBP_help(beamcounter).foci = PBP_data{1,beamcounter}(ia,5);
    for IESloop = 1 : numel(PBP_help(beamcounter).Energies)
        PBP.IES(IESloop).energy =  PBP_help(beamcounter).Energies(IESloop);
        PBP.IES(IESloop).focus =   PBP_help(beamcounter).foci(IESloop);
        PBP.IES(IESloop).data{:,1} = PBP_data{1,beamcounter}((PBP_data{1,beamcounter}(:,4) ==  PBP_help(beamcounter).Energies(IESloop)),1); % X-Position
        PBP.IES(IESloop).data{:,2} = PBP_data{1,beamcounter}((PBP_data{1,beamcounter}(:,4) ==  PBP_help(beamcounter).Energies(IESloop)),2); % Y-Position
        PBP.IES(IESloop).data{:,3} = PBP_data{1,beamcounter}((PBP_data{1,beamcounter}(:,4) ==  PBP_help(beamcounter).Energies(IESloop)),3); % beamweight
        PBP_help(beamcounter).PBP_x{IESloop}(:,1) = PBP_data{1,beamcounter}((PBP_data{1,beamcounter}(:,4) ==  PBP_help(beamcounter).Energies(IESloop)),1);
        PBP_help(beamcounter).PBP_y{IESloop}(:,1) = PBP_data{1,beamcounter}((PBP_data{1,beamcounter}(:,4) ==  PBP_help(beamcounter).Energies(IESloop)),2);
        PBP_help(beamcounter).PBP_partcount{IESloop}(:,1) = PBP_data{1,beamcounter}((PBP_data{1,beamcounter}(:,4) ==  PBP_help(beamcounter).Energies(IESloop)),3);
        PBP.IES(IESloop).data =  cell2mat(PBP.IES(IESloop).data);
    end
    PBP_help(beamcounter).NumIES = IESloop;


    %% header needed for phys Beamplan
    if exportbool
        switch exportformat
            case 'xml'
                PBP_header_1 = ['<?xml version="1.0" encoding="utf-8"?>','\n',...
                    '<PTTxPlan>', '\n',...
                    '\t', '<Beam uid="5dca1dcc-790b-4ae0-a162-e620f7d118f4">', '\n',...
                    '\t','\t','<RstFormat>PT_2004</RstFormat>', '\n',...
                    '\t','\t','<Patient id="PT-2004-01" name="Ronny Tester" sex="M" birthDate="2004-01-01" />', '\n',...
                    '\t','\t','<TxInitiation therapist="USER NAME" dateTime="2007-01-23T13:52:27.2343750+01:00" />', '\n'];
                proton_header = [...
                    '\t','\t','<TxRoom name="Room2" projectile="PROTON" charge="" mass="" atomicNumber="" />', '\n',...
                    '\t','\t','<BAMS rippleFilter="254" rangeShifter="3" rangeShifterDistance="20.0" />', '\n'];
                PBP_header_2 = [...
                    '\t','\t','<TxTable lateral="150.0" longitudinal="-400.0" vertical="-100.5" roll=".0" pitch=".0" isocentricAngle=".0" />', '\n',...
                    '\t','\t','<Gantry angle="90.0" />', '\n'];

                PBP_header = [PBP_header_1, proton_header, PBP_header_2];

                %% print to file

                beamplan_name = [exportname , num2str(beamcounter)];
                fid_PBP = fopen([exportpath filesep beamplan_name, '.xml'], 'w');
                fprintf(fid_PBP, PBP_header);
                for IESloop = 1: numel(PBP_help(beamcounter).Energies)
                    fprintf(fid_PBP, ['\t','\t','<IES number="', num2str(IESloop), '" energy="', num2str(PBP_help(beamcounter).Energies(IESloop),'%.3f') '" focus="', num2str(PBP_help(beamcounter).foci(IESloop),'%.3f'), '">', '\n']);
                    for voxelloop = 1:numel(PBP_help(beamcounter).PBP_x{1,IESloop})
                        fprintf(fid_PBP, ['\t','\t','\t','<Voxel x="',...
                            num2str(PBP_help(beamcounter).PBP_x{1,IESloop}(voxelloop,1),'%.3f'),...
                            '" y="', num2str(PBP_help(beamcounter).PBP_y{1,IESloop}(voxelloop,1),'%.3f'),...
                            '" particles="', num2str(PBP_help(beamcounter).PBP_partcount{1,IESloop}(voxelloop,1),'%.2f'),...
                            '"/>', '\n']);
                    end
                    fprintf(fid_PBP, ['\t','\t','</IES>', '\n']);
                end
                fclose(fid_PBP);

            case 'rst'
                LibP = loadLib_P('C:\Users\Matth\OneDrive\Matlab\ParticleTools\lib\witt\pFocGaussFit\', 'LIBC_p_MIT.csv');
                gantryangle = mod(stf.gantryAngle-90,360);
                couchangle =  mod(stf.couchAngle-90,360);
                minpart_all = min(PBP_data{1,beamcounter}(:,3));
                maxpart_all = max(PBP_data{1,beamcounter}(:,3));
                sumpart_all = sum(PBP_data{1,beamcounter}(:,3));
                PatientID = 'generic-PatientID';
                rst_header = ['rstfile 20030630' newline ...
                    'sistable MIT_proton_2018' newline ...
                    'patient_id ' PatientID, newline ...
                    'machine# 0' newline ...
                    'projectile 1H' newline ...
                    'charge 1' newline ...
                    'mass 1' newline ...
                    'gantryangle ' num2str(gantryangle), newline ...
                    'couchangle ' num2str(couchangle), newline ...
                    'stereotacticcoordinates' newline ...
                    'bolus 0' newline ...
                    'ripplefilter 0', newline ...
                    '#submachines ' num2str(IESloop), newline ...
                    '#particles ' num2str(minpart_all, '%.10E'), ' ', num2str(maxpart_all, '%.10E'), ' ', num2str(sumpart_all, '%.10E'), newline ...
                    ];

                %% print to file
                beamplan_name = [exportname , num2str(beamcounter)];
                fid_rst = fopen([exportpath filesep beamplan_name, '.rst'], 'w');
                fprintf(fid_rst, rst_header);
                for IESloop = 1: numel(PBP_help(beamcounter).Energies)
                    %simple names for Export
                    PBP_help.PBP_partcount{1, IESloop} = PBP_help.PBP_partcount{1, IESloop};
                    [iE, ~] = find(LibP == PBP_help(beamcounter).Energies(IESloop)); % energy index according to SIS-file
                    energyMEV = PBP_help(beamcounter).Energies(IESloop);
                    [~, iFoc] = find(LibP(iE,3:end) == PBP_help(beamcounter).foci(IESloop)); % Focus index according to SIS-file
                    FWHM = PBP_help(beamcounter).foci(IESloop); % FWHM according to SIS-file
                    minpart_IES = min(PBP_help.PBP_partcount{1, IESloop});
                    maxpart_IES = max(PBP_help.PBP_partcount{1, IESloop});
                    sumpart_IES = sum(PBP_help.PBP_partcount{1, IESloop});
                    numpoints_IES = size(PBP_help.PBP_partcount{1, IESloop},1);
                    submachine = ['submachine# ', num2str(iE), ' ', num2str(energyMEV), ' ', num2str(iFoc), ' ',  num2str(FWHM), newline ...
                        '#particles ', num2str(minpart_IES, '%.10E'), ' ', num2str(maxpart_IES, '%.10E'), ' ', num2str(sumpart_IES, '%.10E') , newline ...
                        'stepsize ' num2str(pln.propStf.bixelWidth) ,' ',  num2str(pln.propStf.bixelWidth), newline ...
                        '#points ', num2str(numpoints_IES), newline];

                    fprintf(fid_rst, submachine);

                    for voxelloop = 1:numel(PBP_help(beamcounter).PBP_x{1,IESloop})
                        fprintf(fid_rst, [num2str(PBP_help(beamcounter).PBP_y{1,IESloop}(voxelloop,1),'%.3f'), ' ', ...
                            num2str((PBP_help(beamcounter).PBP_x{1,IESloop}(voxelloop,1)*-1),'%.3f'), ' ', ...
                            num2str(PBP_help(beamcounter).PBP_partcount{1,IESloop}(voxelloop,1),'%.2f'),...
                            newline]);
                    end
                end
                fclose(fid_rst);

            otherwise
                disp('wrong format, provide string with "rst" or "xml"')
        end
    else
        return
    end
end
PBP(beamcounter).Allpoints = cell2mat(PBP_data);
end

