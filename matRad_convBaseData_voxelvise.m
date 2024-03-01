function [Zmod_vox] = matRad_convBaseData_voxelvise(baseData, radiologicalDepth, modulationDepth, Pmod, varargin)
% matRad particle dose calculation wrapper
% 
% call
%   conv = matRad_calcParticleDose(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   baseData:           matRad basedata of format machine.data(energyIx)
%   radiologicalDepth:  computation depth in which the dose is needed
%   modulationDepth:    current modulation depth in mm
%   Pmod:               Modulation Power in µm
%
% output
%   Zmod:              modulated DDD for given lungdepth 
%
% References
%   [1] https://www.thm.de/lse/images/user/KZink-105/Abschlussarbeiten/Masterarbeit_Matthias_Witt_2014.pdf
%   [2] http://archiv.ub.uni-marburg.de/diss/z2020/0261
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% skip if current modulation depth is 0 => no modulating material in
% beampath => no convolution neccessary 
% this is already skipped in the call in calcParticleDoseBixel and stays
% here for safety
if modulationDepth == 0
    Zmod_vox = matRad_interp1(baseData.depths, baseData.Z, radiologicalDepth);
    return
else
    %  P_Mod = sigma^2/ mu
    mu = modulationDepth;
    sigma = sqrt((Pmod/1000)*mu);
    
    % Abfrage ob sigma 0 ist => NaN
    if sigma == 0
        Zmod_vox = matRad_interp1(baseData.depths, baseData.Z, radiologicalDepth);
        return
    end
    
    % adapt mu to the center of the prolonged base data
    % in that way no shift of the bragg curve is introduced
    % based on the way the basedata is prolonged in matRad_extendBaseData
    % this value should always be 0!
    mu_adapt = baseData.depths_adapted(ceil(size(baseData.depths_adapted,1)/2));
    
    % calculate normaldistribution
%     normvert(:,1) = baseData.depths_adapted;
%     normvert(:,2) = ( 1 / sqrt(2*pi()*sigma.^2) ) * exp(- ((mu_adapt - normvert(:,1)).^2)/(2*sigma.^2) );
    normvert = ( 1 / sqrt(2*pi()*sigma.^2) ) * exp(- ((mu_adapt - baseData.depths_adapted).^2)/(2*sigma.^2) );
%     normvert(:,3) = ( 1 / sqrt(2*pi()*sigma.^2) ) * exp(- ((mu - normvert(:,1)).^2)/(2*sigma.^2) );
    % Convolution
    try 
        conv_res = convnfft(baseData.Z_adapted, normvert, 'same');
    catch
        conv_res = conv(baseData.Z_adapted, normvert, 'same');
    end
    
    % shorten to relevant length (negative values omitted)
    
    conv_res = conv_res(ceil(size(baseData.depths_adapted,1)/2):end);
    % introduce a scaling in a way that the area under of modulated and
    % unmodulated bragg curves are the same => Dosiserhaltung!
    % scaling_conv =  sum(conv_res) / sum(baseData.Z.profileORG); %temp MaW
    scaling_conv =  sum(conv_res) / sum(baseData.Z); %temp MaW
    Zmod = conv_res./scaling_conv;
end



Zmod_vox = matRad_interp1(baseData.depths, Zmod, radiologicalDepth);



% % %% uncomment to see the basedata convolution
% if ~isempty(varargin)
%     figure
%     axes1 = axes('Parent',gcf,'YColor',[0.2 0.2 0.2],'XColor',[0.2 0.2 0.2],...
%         'TickDir','in',...
%         'TickLength',[0.01 0.01],...
%         'Position',[0.13 0.40990371389271 0.775 0.56009628610729],...
%         'LineWidth',3,...
%         'FontWeight','bold',...
%         'FontSize',14,...
%         'XTickLabel', []);
%     box(axes1,'on');
%     hold(axes1,'all');
%     min_ref = min(baseData.depths);
%     max_ref = max(baseData.depths);
% 
%     plot(baseData.depths, baseData.Z,'LineWidth', 3, 'Color', 'black')
%     plot(baseData.depths, Zmod,'LineWidth', 3, 'Color', 'blue');
%     updateAxesProperties(gca, 'LineWidth', 3, 'xlim', [min_ref max_ref], 'FontSize', 20)
%     ylimits = get(gca, 'ylim');
%     set(gca, 'ylim', [0 ylimits(2)])
%     xlabel('Range in water /mm')
%     ylabel('MeV cm^2/g per primary')
%     legend({'BaseData', 'Modulated BaseData'}, 'location', 'NorthWest');
% 
% 
%     axes1 = axes('Parent',gcf,'YColor',[0.2 0.2 0.2],'XColor',[0.2 0.2 0.2],...
%         'TickDir','in',...
%         'TickLength',[0.01 0.01],...
%         'Position',[0.13 0.0839064649243466 0.775 0.325997248968363],...
%         'LineWidth',3,...
%         'FontWeight','bold',...
%         'FontSize',14);
%     box(axes1,'on');
%     hold(axes1,'all');
%     plot(baseData.depths_adapted +mu , normvert,'LineWidth', 3);
%     updateAxesProperties(gca, 'LineWidth', 3,  'FontSize', 20, 'xlim', [min_ref max_ref])%, 'ylim', [0 0.08])
%     xlabel('range in water /mm')
%     ylabel('frequency')
%     legend({'modulationfunction'}, 'Location', 'NorthWest');
% 
%     % Create textbox
%     annotation(gcf,'textbox',...
%         [0.162823361823361 0.183521492945849 0.132763532763533 0.105360443622921],...
%         'String',{['used parameter:', newline, '\sigma: ', num2str(sigma,3), ' mm', newline , '\mu: ', num2str(mu,3), ' mm', newline , 'P_{mod}: ', num2str(Pmod,3), '\mum']},..., sprintf('\n'), 'Amplitude: ', num2str(w(3)), sprintf('\n'), 'Mod_Pow: ' num2str(abs(w(1))/sqrt(abs(w(2))))]},...
%         'FitBoxToText','on', ...
%         'FontSize', 15, ...
%         'FontWeight', 'normal',...
%         'Interpreter', 'tex');
% 
%     % %     figure, hold on
%     % %     plot(baseData.depths, baseData.Z, '-')
%     % %     plot(baseData.depths, conv_res, '--')
%     % %     line([baseData.peakPos baseData.peakPos], [0 50], 'Color', 'black')
%     % %     figure, hold on
%     % %     plot(normvert(:,1),normvert(:,2))
% % end

end
