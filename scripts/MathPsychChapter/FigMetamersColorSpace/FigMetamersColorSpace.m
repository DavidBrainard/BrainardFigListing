function varargout = FigMetamersColorSpace(varargin)
%
% Compute some metamers and plot in tristimulus space
%
% Used for PNAS News and Views
%
% See also FigDichromMetam, FigSpatioChromaticAliasing.
%
% 05/15/22  dhb  Wrote it.
    
%% Clear
clear; close all;

%% Get standard XYZ tristimulus cmfs
%
% This gets us the standard CIE fundamentals, for a 2 degree field.
% We could adjust for observer age, if we wanted.  32 years old is
% the standard default.
S = [400 1 301];
wls = SToWls(S);
load T_xyz1931.mat
T = SplineCmf(S_xyz1931,T_xyz1931,S);

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

% Load in a monitor
load B_monitor.mat
B_metamer = SplineSpd(S_monitor,B_monitor(:,1:3),S);
B_metamer = MakeGaussBasis(wls,[430 520 650],8*[120 120 120]);

% Get a spectrum
spd1 = GenerateBlackBody(4000,wls);
spd1 = spd1/max(spd1);
spd1 =spd1/2;
XYZ1 = T*spd1;
metamer1 = B_metamer*inv(T*B_metamer)*XYZ1;
XYZ1Check = T*metamer1;
if (max(abs(XYZ1-XYZ1Check)./XYZ1)  > 1e-10)
    error('Did not compute an actual metamer');
end

% Make metamer figure 1
figParams.figName = 'FigMetamer1';
figParams.xLimLow = 380;
figParams.xLimHigh = 720;
figParams.xTicks = [400 500 600 700];
figParams.xTickLabels = {'400' '500' '600' '700'};
figParams.yLimLow = 0.0;
figParams.yLimHigh = 1.5;
figParams.yTicks = [0.0 0.5 1 1.5];
figParams.yTickLabels = {' 0.00 ' ' 0.50 ' ' 1.00 ' ' 1.50 '};
theFig1 = figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
plot(wls,spd1,'Color',[0 0 0],'LineStyle','-','LineWidth',figParams.lineWidth);
plot(wls,metamer1,'Color',[1 0 0],'LineStyle','-','LineWidth',figParams.lineWidth);
xlim([figParams.xLimLow figParams.xLimHigh]);
set(gca,'XTick',figParams.xTicks);
set(gca,'XTickLabel',figParams.xTickLabels);
xlabel('Wavelength (nm)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylim([figParams.yLimLow figParams.yLimHigh]);
ylabel('Power (arb units)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
set(gca,'YTick',figParams.yTicks);
set(gca,'YTickLabel',figParams.yTickLabels);
%legend({' L cones ' ' M cones ' ' S cones '},'Location','NorthEast','FontSize',figParams.legendFontSize);
axis('square');
%title('Cone Fundamentals (2°)')
%set(gca,'XMinorTick','on');
% FigureSave(fullfile(pwd,[mfilename '_' figParams.figName]),theFig,figParams.figType);

% Get a spectrum
spd2 = GenerateBlackBody(20000,wls);
spd2 = spd2/max(spd2);
XYZ2 = T*spd2;
metamer2 = B_metamer*inv(T*B_metamer)*XYZ2;
XYZ2Check = T*metamer2;
if (max(abs(XYZ2-XYZ2Check)./XYZ2)  > 1e-10)
    error('Did not compute an actual metamer');
end

% Make metamer figure 2
figParams.figName = 'FigMetamer2';
figParams.xLimLow = 380;
figParams.xLimHigh = 720;
figParams.xTicks = [400 500 600 700];
figParams.xTickLabels = {'400' '500' '600' '700'};
figParams.yLimLow = 0;
figParams.yLimHigh = 1.5;
figParams.yTicks = [0.0 0.5 1 1.5];
figParams.yTickLabels = {' 0.00 ' ' 0.50 ' ' 1.00 ' ' 1.50 '};
theFig2 = figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
plot(wls,spd2,'Color',[0 0 0],'LineStyle','-','LineWidth',figParams.lineWidth);
plot(wls,metamer2,'Color',[0 0 1],'LineStyle','-','LineWidth',figParams.lineWidth);
xlim([figParams.xLimLow figParams.xLimHigh]);
set(gca,'XTick',figParams.xTicks);
set(gca,'XTickLabel',figParams.xTickLabels);
xlabel('Wavelength (nm)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylim([figParams.yLimLow figParams.yLimHigh]);
ylabel('Power (arb units)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
set(gca,'YTick',figParams.yTicks);
set(gca,'YTickLabel',figParams.yTickLabels);
%legend({' L cones ' ' M cones ' ' S cones '},'Location','NorthEast','FontSize',figParams.legendFontSize);
axis('square');
%title('Cone Fundamentals (2°)')
%set(gca,'XMinorTick','on');
% FigureSave(fullfile(pwd,[mfilename '_' figParams.figName]),theFig,figParams.figType);

% Make KXYZ figure 3
figParams.figName = 'FigXYZ';
figParams.xLimLow = 380;
figParams.xLimHigh = 720;
figParams.xTicks = [400 500 600 700];
figParams.xTickLabels = {'400' '500' '600' '700'};
figParams.yLimLow = 0;
figParams.yLimHigh = 1;
figParams.yTicks = [0.0 0.2 0.4 0.6 0.8 1.0];
figParams.yTickLabels = {' 0 ' ' 50 ' ' 100 '};
theFig2 = figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
plot3(XYZ1(1),XYZ1(2),XYZ1(3),'o','Color',[1 0 0],'MarkerFaceColor',[1 0 0],'MarkerSize',20);
plot3(XYZ2(1),XYZ2(2),XYZ2(3),'o','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerSize',20);
view(30, 15);
xlim([0 100]);
ylim([0 100]);
zlim([0 100]);
xlabel('X','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylabel('Y','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
%set(gca,'XTickLabel',figParams.yTickLabels);
set(gca,'YTickLabel',figParams.yTickLabels);
set(gca,'ZTickLabel',figParams.yTickLabels);
zlabel('Z','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
axis('square');
grid('on');
ax = gca; 
ax.YAxisLocation = 'origin';


end