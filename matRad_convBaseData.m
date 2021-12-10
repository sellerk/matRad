function [Zmod, normvert, mu, sigma, Pmod] = matRad_convBaseData(baseData, modulationDepth, Pmod)
% matRad convolution of base data 
% 
% call
%   [Zmod, normvert, mu, sigma, Pmod] = matRad_convBaseData(baseData, modulationDepth, Pmod)
%
% input
%   baseData:           matRad basedata of format machine.data(energyIx)
%   modulationDepth:    current modulation depth in mm
%   Pmod:               Modulation Power in µm
%
% output
%   Zmod:              modulated DDD for given lungdepth
%   normvert:          calculated normaldistribution for modulation
%   mu:                used mu for calculation of normaldistribution
%   sigma:             used sigma for calculation of normaldistribution
%   Pmod:              used Pmod for calculation of normaldistribution    
%
% References
%   [1] https://www.thm.de/lse/images/user/KZink-105/Abschlussarbeiten/Masterarbeit_Matthias_Witt_2014.pdf
%   [2] http://archiv.ub.uni-marburg.de/diss/z2020/0261
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  P_Mod = sigma^2/ mu
mu = modulationDepth;
sigma = sqrt((Pmod/1000)*mu);

% Abfrage ob sigma 0 ist => NaN 
if sigma == 0
    Zmod = baseData.Z;
    return
end

% prolong basedata into the negative to bypass convolution effects at the edge
% use 5 sigmas to make sure nothing is missing
addlength = ceil(5*sigma)*100;
adddepths = flipud((baseData.depths(1:addlength)*-1));
depths = [adddepths; baseData.depths];
addZ = (ones(addlength,1).*baseData.Z(1));
Z = [addZ;baseData.Z];

% adapt mu to the center of the prolonged base data
% in that way no shift of the bragg curve is introduced
mu_atapt = depths((size(depths,1)/2)+1);

% calculate normaldistribution
normvert(:,1) = depths;
normvert(:,2) = ( 1 / sqrt(2*pi()*sigma.^2) ) * exp(- ((mu_atapt - normvert(:,1)).^2)/(2*sigma.^2) );
normvert(:,3) = ( 1 / sqrt(2*pi()*sigma.^2) ) * exp(- ((mu - normvert(:,1)).^2)/(2*sigma.^2) );
% Concolution
conv_res = conv(Z, normvert(:,2), 'same');
% shorten to relevant length (negative values omitted)
conv_res = conv_res(addlength+1:end);
% introduce a scaling in a way that the area under of modulated and
% unmodulated bragg curves are the same => Dosiserhaltung!
scaling_conv =  sum(conv_res) / sum(baseData.Z);
Zmod = conv_res./scaling_conv;

%% uncomment to see the basedata convolution
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
%     plot(normvert(:,1), normvert(:,3),'LineWidth', 3);
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
% %     figure, hold on
% %     plot(baseData.depths, baseData.Z, '-')
% %     plot(baseData.depths, conv_res, '--')
% %     line([baseData.peakPos baseData.peakPos], [0 50], 'Color', 'black')
% %     figure, hold on
% %     plot(normvert(:,1),normvert(:,2))
% end
%%
end
