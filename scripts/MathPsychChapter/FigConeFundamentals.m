function varargout = FigTrichromMetam(varargin)
%
% Show human cone sensitivities and a metamer.
%
% Used to produce panels for Figure 1 of the 2015 Annual Review Vision
% Science paper.
%   Panel A: FigTrichromCones.pdf
%   Panel B: FigTrichromMetam.pdf
%
% See also FigDichromMetam, FigSpatioCHromaticAliasing.
%
% 3/23/15   dhb  Wrote it.
% 4/20/15   dhb  Bring into the validation fold
% 4/29/15   dhb  Added some extra explanatory plots
% 07/18/20  dhb  Moved some calculations quickly to physical units.

% This needs an old version of ISETBio and to be cleaned up so everything
% makes sense in physical units.  Better to express spectra in radiance,
% but currently in retinal irradiance because I am lazy.
    
%% Clear
clear; close all;

%% Get standard LMS cone spectral sensitivities in quantal units
%
% This gets us the standard CIE fundamentals, for a 2 degree field.
% We could adjust for observer age, if we wanted.  32 years old is
% the standard default.
S = [380 1 401];
wls = SToWls(S);
coneParams = DefaultConeParams('cie_asano');
coneParams.ageYears = 32;
coneParams.fieldSizeDegrees = 2;
[~,T_energy,T_quanta] = ComputeObserverFundamentals(coneParams,S);

%% Figure parameters
curDir = pwd;
masterFigParamsDir = getpref('bfScripts','masterFigParamsDir');
cd(masterFigParamsDir);
figParams = MasterFigParams;
cd(curDir);
if (exist('../SecondaryFigParams','file'))
    cd ..
    figParams = SecondaryFigParams(figParams);
    cd(curDir);
end


%% Make the cone sensitivity figure
figParams.figName = 'FigTrichromCones_2';
figParams.xLimLow = 380;
figParams.xLimHigh = 720;
figParams.xTicks = [400 500 600 700];
figParams.xTickLabels = {'400' '500' '600' '700'};
figParams.yLimLow = -0.01;
figParams.yLimHigh = 0.76;
figParams.yTicks = [0.0 0.25 0.5 0.75];
figParams.yTickLabels = {' 0.00 ' ' 0.25 ' ' 0.50 ' ' 0.75 '};
theFig = figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
plot(wls,T_quanta(1,:)','r','LineWidth',figParams.lineWidth);
plot(wls,T_quanta(2,:)','g','LineWidth',figParams.lineWidth);
plot(wls,T_quanta(3,:)','b','LineWidth',figParams.lineWidth);
xlim([figParams.xLimLow figParams.xLimHigh]);
set(gca,'XTick',figParams.xTicks);
set(gca,'XTickLabel',figParams.xTickLabels);
xlabel('Wavelength (nm)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylim([figParams.yLimLow figParams.yLimHigh]);
ylabel('Quantal Efficiency','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
set(gca,'YTick',figParams.yTicks);
set(gca,'YTickLabel',figParams.yTickLabels);
legend({' L cones ' ' M cones ' ' S cones '},'Location','NorthEast','FontSize',figParams.legendFontSize);
axis('square');
title('Cone Sensitivities (2°)')
%set(gca,'XMinorTick','on');
FigureSave(fullfile(pwd,[mfilename '_' figParams.figName]),theFig,figParams.figType);

coneParams.fieldSizeDegrees = 10;
[~,T_energy,T_quanta] = ComputeObserverFundamentals(coneParams,S);
figParams.figName = 'FigTrichromCones_10';
figParams.xLimLow = 380;
figParams.xLimHigh = 720;
figParams.xTicks = [400 500 600 700];
figParams.xTickLabels = {'400' '500' '600' '700'};
figParams.yLimLow = -0.01;
figParams.yLimHigh = 0.76;
figParams.yTicks = [0.0 0.25 0.5 0.75];
figParams.yTickLabels = {' 0.00 ' ' 0.25 ' ' 0.50 ' ' 0.75 '};
theFig = figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
plot(wls,T_quanta(1,:)','r','LineWidth',figParams.lineWidth);
plot(wls,T_quanta(2,:)','g','LineWidth',figParams.lineWidth);
plot(wls,T_quanta(3,:)','b','LineWidth',figParams.lineWidth);
xlim([figParams.xLimLow figParams.xLimHigh]);
set(gca,'XTick',figParams.xTicks);
set(gca,'XTickLabel',figParams.xTickLabels);
xlabel('Wavelength (nm)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylim([figParams.yLimLow figParams.yLimHigh]);
ylabel('Quantal Efficiency','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
set(gca,'YTick',figParams.yTicks);
set(gca,'YTickLabel',figParams.yTickLabels);
legend({' L cones ' ' M cones ' ' S cones '},'Location','NorthEast','FontSize',figParams.legendFontSize);
axis('square');
title('Cone Sensitivities (10°)')
%set(gca,'XMinorTick','on');
FigureSave(fullfile(pwd,[mfilename '_' figParams.figName]),theFig,figParams.figType);

end