% This matlab script combines 7 cubes of GASS II observational data,
% integrates temperature values along line-of-sight velocity axis and 
% visualises the field shown in Figure 1 of the paper Henderson et al
% "Statistical Topology and the Random Interstellar Medium"

% Each observational data cube is a gas temperature as a function of 3 variables:
% longitute (l), latitude (b) and line-of-sight velocitiy (v_los).
% For downloading the cubes from the server see Instructions.pdf

  % Cube centers at longitudes = 178, 216, 254, 292, 330, 8, 46 degrees;
% latitudes = -50 degrees and v_los = 20 km/s to42 km/s


%% Uncompress GNU zip files in the current directory
gunzip('*.gz');

%%
clc; clear; close all;

Longitude = [178, 216, 254, 292, 330, 8, 46];

filename = dir(['gass_', num2str(Longitude(1)), '*.fits']);
dataCube = flipud(fitsread(filename.name));
% Uncomment next line to display metadata for the Header/Data Unit (HDU) found in the FITS file 
% fitsdisp(filename.name)

LRef = Longitude(1);
LRefPx = 157;
LIncrem = -1.244579217373E-01;
BRef = 0;
BRefPx = 782;
BIncrem = 8.000000000000E-02;
VRef = 1.978800010681E+04; 
VRefPx = 1;
VIncrem = 8.245000243187E+02;

n_L = 313;
n_V = 28;

% To find longitude (L) and  latitude (B) values for the current data cube,
% based on the reference points.
for i = 1 : n_L
    L(i) = LRef + (i - LRefPx) * LIncrem;
    B(i) = BRef + (i - BRefPx) * BIncrem;
end 
B = fliplr(B);

% To find LINE-OF-SIGHT VELOCITIES values for the current data cube
for i = 1 : n_V
    V(i) = VRef + (i - VRefPx) * VIncrem;
end 
V = V/1000; % to convert velocity to [km/sec]

% Every loop add one more data cube on the left side of the previous portion
% of the data
for lg = 2 : 7
    % Left data cube
    filenameL = dir(['gass_', num2str(Longitude(lg)), '*.fits']);
    % filenameL.name
    LeftPart = flipud(fitsread(filenameL.name));
    % fitsdisp(filename)

    LRefL = Longitude(lg); 

    n_LLeft = 313;
    for i = 1 : n_LLeft
        LLeft(i) = LRefL + (i - LRefPx) * LIncrem;
        if (LLeft(i) < 0)
            LLeft(i) = 360 - abs(LLeft(i));
        end
    end 

    Odds = abs(L - LLeft(n_LLeft));
    OddsMin = min(abs(Odds(:)));
    ind = find(Odds == OddsMin);

    % Horizontal concatenation: LeftPart + RightPart
    dataCube = [LeftPart dataCube(:, ind + 1 : n_L, :)];

    % L for the joined cube: LeftPart + RightPart
    L = [LLeft L(ind + 1 : n_L)];
    n_L = length(L);

end;    

data = trapz(dataCube(1 : 256, 62 : 1907, 2 : 26), 3);  % integrates along the dimension 3 

clearvars -except B L data 

%% To visualize the field as in Figure 1

LTicks = [145, 305, 465, 626, 787, 947, 1108, 1268, 1429, 1589, 1750];
% To create the new tick labels
LLabels = sprintf('%3.0f\n', L(LTicks));

BTicks = [32, 95, 157, 220];
% create the new tick labels
BLabels = sprintf('%3.0f\n', B(BTicks));

figure  
image(wcodemat(data, 190))
axis on, colormap pink, colorbar, axis image
set(gca, 'FontName', 'Times New Roman', 'XTick', LTicks, 'XTickLabel', LLabels, 'YTick', BTicks, 'YTickLabel', BLabels,...
    'TickDir', 'out', 'Box', 'on', 'TickLength',[0.008 0.008], 'FontSize', 14) 
%title('GASS II:  brightness temperature {\it T}_B , [K]', 'FontSize', 13, 'FontWeight', 'normal')
xlabel('Longitude {\it l}, [\circ]', 'FontSize', 15)
ylabel('Latitude {\it b}, [\circ]', 'FontSize', 15)

set(gcf, 'Position', get(0, 'Screensize'));

% To add white frames for 3 regions
hold on
rectangle('Position',[3, 3, 252, 252], 'LineWidth', 1.8, 'LineStyle','--', 'EdgeColor', 'w')
daspect([1,1,1])
rectangle('Position',[796, 3, 252, 252], 'LineWidth', 1.8, 'LineStyle','--', 'EdgeColor', 'w')
daspect([1,1,1])
rectangle('Position',[1591, 3, 252, 252], 'LineWidth', 1.8, 'LineStyle','--', 'EdgeColor', 'w')
daspect([1,1,1])
set(gcf, 'Color', 'w');

clearvars -except B L data











