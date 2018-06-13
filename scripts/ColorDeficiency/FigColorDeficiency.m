function varargout = FigColorDeficiency(varargin)
%
% Show confusions produced by color deficiency before/after "treatment"
%
% 6/13/17   dhb  Wrote it

%{
    tbUse({'isetbio', 'Psychtoolbox-3'});
%}

    
%% Clear
ieInit;

%% Hello
mfiledir = fileparts(mfilename('fullpath'));
curDir = pwd;
cd(mfiledir);
outputDir = fullfile(mfiledir,sprintf('%s_Output',mfilename));
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

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

%% Cone sensitivities.
%
S = [380 4 101]; wls = SToWls(S);
ptbPhotoreceptorsSS2 = DefaultPhotoreceptors('CIE2Deg');
ptbPhotoreceptorsSS2.nomogram.S = S;
ptbPhotoreceptorsSS2.nomogram.source = 'StockmanSharpe';
ptbPhotoreceptorsSS2.nomogram.lambdaMax = [558.9 530 420.7]';
ptbPhotoreceptorsSS2 = rmfield(ptbPhotoreceptorsSS2,'absorbance');
ptbPhotoreceptorsSS2 = FillInPhotoreceptors(ptbPhotoreceptorsSS2);
T_conesQETrichrom = [ptbPhotoreceptorsSS2.isomerizationAbsorptance(1,:) ; ptbPhotoreceptorsSS2.isomerizationAbsorptance(2,:) ; ptbPhotoreceptorsSS2.isomerizationAbsorptance(3,:)];

%% We also need some anomolous cones
%
ptbPhotoreceptorsAnom = DefaultPhotoreceptors('CIE2Deg');
ptbPhotoreceptorsAnom.nomogram.S = S;
ptbPhotoreceptorsAnom.nomogram.source = 'StockmanSharpe';
ptbPhotoreceptorsAnom.nomogram.lambdaMax = [558.9 557 420.7]';
ptbPhotoreceptorsAnom = rmfield(ptbPhotoreceptorsAnom,'absorbance');
ptbPhotoreceptorsAnom = FillInPhotoreceptors(ptbPhotoreceptorsAnom);
T_conesQEAnom = [ptbPhotoreceptorsAnom.isomerizationAbsorptance(1,:) ; ptbPhotoreceptorsAnom.isomerizationAbsorptance(2,:) ; ptbPhotoreceptorsAnom.isomerizationAbsorptance(3,:)];

% Energy sensitivities.  We need these to convert the isomerizations
% into the right relative scaling for SRGB rendering.
T_conesTrichrom = EnergyToQuanta(S,T_conesQETrichrom')';
T_conesTrichrom = T_conesTrichrom/max(T_conesTrichrom(:));
T_conesAnom = EnergyToQuanta(S,T_conesQEAnom')';
T_conesAnom = T_conesAnom/max(T_conesAnom(:));
for ii = 1:3
    renderingScaleFactors(ii) = 1/max(T_conesTrichrom(ii,:));
end

% Generate hues around a hue circle
munsellValue = 5;
munsellChroma = 5;
nHues = 85;
hueAngles = linspace(0,360,nHues+1);
hueAngles = hueAngles(1:end-1);
munsellData = MunsellPreprocessTable;
[~,Xx,trix,vx,Xy,triy,vy,XY,triY,vY] = MunsellGetxyY(0,munsellValue,munsellChroma,munsellData);
for hh = 1:nHues
    hueCircle_xyY(:,hh) = MunsellGetxyY(hueAngles(hh),munsellValue,munsellChroma,[],Xx,trix,vx,Xy,triy,vy,XY,triY,vY);
end
hueCircle_XYZ = xyYToXYZ(hueCircle_xyY);

%% Plot the xy values of the hue circle
%
% This is a sanity check
theFig = figure; clf; hold on
plot(hueCircle_xyY(1,:),hueCircle_xyY(2,:),'ro','MarkerFaceColor','r','MarkerSize',8);

%% Convert Munsell XYZ to reflectance spectra
load spd_CIEC.mat
spd_CIEC = SplineSpd(S_CIEC,spd_CIEC,S);
load T_xyz1931.mat
T_xyz = SplineCmf(S_xyz1931,T_xyz1931,S);
load B_nickerson
B_sur = SplineSrf(S_nickerson,B_nickerson(:,1:3),S);
M_wgtsToXYZ = T_xyz*diag(spd_CIEC)*B_sur;
M_XYZToSur = B_sur*inv(M_wgtsToXYZ);
hueCircle_sur = M_XYZToSur*hueCircle_XYZ;

%% Convert surfaces to LMS for normal and anom
hueCircle_LMSTrichrom = T_conesTrichrom*diag(spd_CIEC)*hueCircle_sur;
hueCircle_LMSAnom = T_conesAnom*diag(spd_CIEC)*hueCircle_sur;

figure; clf; hold on
plot3(hueCircle_LMSTrichrom(1,:),hueCircle_LMSTrichrom(2,:),hueCircle_LMSTrichrom(3,:), ...
    'ro','MarkerFaceColor','r','MarkerSize',8);
plot3(hueCircle_LMSTrichrom(1,:),hueCircle_LMSTrichrom(2,:),hueCircle_LMSTrichrom(3,:), ...
    'r');
xlabel('L'); ylabel('M'); zlabel('S');
title('Normal Trichrom');

figure; clf; hold on
plot3(hueCircle_LMSAnom(1,:),hueCircle_LMSAnom(2,:),hueCircle_LMSAnom(3,:), ...
    'ro','MarkerFaceColor','r','MarkerSize',8);
plot3(hueCircle_LMSAnom(1,:),hueCircle_LMSAnom(2,:),hueCircle_LMSAnom(3,:), ...
    'r');
xlabel('L'); ylabel('M'); zlabel('S');
title('Anomolous');

% %% Make a bar plot of cone responses to the two lights
% mid1 = 1; mid2 = 2.5; mid3 = 4.0; sep = 0.25;
% figParams.figName = 'FigTrichromHistoTrichomMetamers';
% figParams.barWidth = 0.4;
% figParams.barLineWidth = 5;
% figure; clf; hold on
% set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
% set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
% xlim([sep mid3+sep+figParams.barWidth]);
% h1 = bar([mid1 - sep]',[theSpdCones(1)],figParams.barWidth,'FaceColor','r','LineWidth',figParams.barLineWidth);
% h2 = bar([mid1 + sep]',[theMetamerCones(1)],figParams.barWidth,'FaceColor','r','LineWidth',figParams.barLineWidth,'LineStyle',':');
% h3 = bar([mid2 - sep]',[theSpdCones(2)],figParams.barWidth,'FaceColor','g','LineWidth',figParams.barLineWidth);
% h4 = bar([mid2 + sep]',[theMetamerCones(2)],figParams.barWidth,'FaceColor','g','LineWidth',figParams.barLineWidth,'LineStyle',':');
% h5 = bar([mid3 - sep]',[theSpdCones(3)],figParams.barWidth,'FaceColor','b','LineWidth',figParams.barLineWidth);
% h6 = bar([mid3 + sep]',[theMetamerCones(3)],figParams.barWidth,'FaceColor','b','LineWidth',figParams.barLineWidth,'LineStyle',':');
% set(gca,'XTick',[mid1 mid2 mid3]);
% set(gca,'XTickLabel',{'L' 'Ms' 'S'},'FontName',figParams.fontName,'FontSize',figParams.labelFontSize+4);
% ylim([0 1]);
% set(gca,'YTick',[0 0.5 1]);
% set(gca,'YTickLabel',{' 0.00 ', ' 0.50 ', ' 1.00 '},'FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
% ylabel('Isomerization Rate (arbitrary units)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
% FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),gcf,figParams.figType);
% 
% % and just to the first light (for explanatory purposes)
% figParams.figName = 'FigTrichromHistoTrichomMetamersOne';
% figParams.barWidth = 0.4;
% figParams.barLineWidth = 5;
% figure; clf; hold on
% set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
% set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
% xlim([sep mid3+sep+figParams.barWidth]);
% h1 = bar([mid1 - sep]',[theSpdCones(1)],figParams.barWidth,'FaceColor','r','LineWidth',figParams.barLineWidth);
% h3 = bar([mid2 - sep]',[theSpdCones(2)],figParams.barWidth,'FaceColor','g','LineWidth',figParams.barLineWidth);
% h5 = bar([mid3 - sep]',[theSpdCones(3)],figParams.barWidth,'FaceColor','b','LineWidth',figParams.barLineWidth);
% set(gca,'XTick',[mid1 mid2 mid3]);
% set(gca,'XTickLabel',{'L' 'M' 'S'},'FontName',figParams.fontName,'FontSize',figParams.labelFontSize+4);
% ylim([0 1]);
% set(gca,'YTick',[0 0.5 1]);
% set(gca,'YTickLabel',{' 0.00 ', ' 0.50 ', ' 1.00 '},'FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
% ylabel('Isomerization Rate (arbitrary units)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
% FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),gcf,figParams.figType);
% 
% % and just L cone response to first light  (again, for explanatory purposes)
% figParams.figName = 'FigTrichromHistoTrichomMetamersOneLOnly';
% figParams.barWidth = 0.4;
% figParams.barLineWidth = 5;
% figure; clf; hold on
% set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
% set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
% xlim([sep mid3+sep+figParams.barWidth]);
% h1 = bar([mid1 - sep]',[theSpdCones(1)],figParams.barWidth,'FaceColor','r','LineWidth',figParams.barLineWidth);
% set(gca,'XTick',[mid1 mid2 mid3]);
% set(gca,'XTickLabel',{'L' 'M' 'S'},'FontName',figParams.fontName,'FontSize',figParams.labelFontSize+4);
% ylim([0 1]);
% set(gca,'YTick',[0 0.5 1]);
% set(gca,'YTickLabel',{' 0.00 ', ' 0.50 ', ' 1.00 '},'FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
% ylabel('Isomerization Rate (arbitrary units)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
% FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),gcf,figParams.figType);

%% Make the cone sensitivity figure
figParams.figName = 'FigTrichromCones';
figParams.xLimLow = 380;
figParams.xLimHigh = 720;
figParams.xTicks = [400 500 600 700];
figParams.xTickLabels = {'400' '500' '600' '700'};
figParams.yLimLow = -0.01;
figParams.yLimHigh = 0.51;
figParams.yTicks = [0.0 0.25 0.5];
figParams.yTickLabels = {' 0.00 ' ' 0.25 ' ' 0.50 '};
figParams.lineWidth = 5;
theFig = figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
plot(wls,T_conesQETrichrom(1,:)','r','LineWidth',figParams.lineWidth);
plot(wls,T_conesQETrichrom(2,:)','g','LineWidth',figParams.lineWidth);
plot(wls,T_conesQETrichrom(3,:)','b','LineWidth',figParams.lineWidth);
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
%set(gca,'XMinorTick','on');
FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),theFig,figParams.figType);

% and an anomalous version for explanatory purposes
figParams.figName = 'FigTrichromAnomCones';
figParams.xLimLow = 380;
figParams.xLimHigh = 720;
figParams.xTicks = [400 500 600 700];
figParams.xTickLabels = {'400' '500' '600' '700'};
figParams.yLimLow = -0.01;
figParams.yLimHigh = 0.51;
figParams.yTicks = [0.0 0.25 0.5];
figParams.yTickLabels = {' 0.00 ' ' 0.25 ' ' 0.50 '};
theFig = figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
plot(wls,T_conesQEAnom(1,:)','r','LineWidth',figParams.lineWidth);
plot(wls,T_conesQEAnom(2,:)','y','LineWidth',figParams.lineWidth);
plot(wls,T_conesQEAnom(3,:)','b','LineWidth',figParams.lineWidth);
xlim([figParams.xLimLow figParams.xLimHigh]);
set(gca,'XTick',figParams.xTicks);
set(gca,'XTickLabel',figParams.xTickLabels);
xlabel('Wavelength (nm)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylim([figParams.yLimLow figParams.yLimHigh]);
ylabel('Quantal Efficiency','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
set(gca,'YTick',figParams.yTicks);
set(gca,'YTickLabel',figParams.yTickLabels);
legend({' L cones ' ' M'' cones ' ' S cones '},'Location','NorthEast','FontSize',figParams.legendFontSize);
axis('square');
%set(gca,'XMinorTick','on');
FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),theFig,figParams.figType);

cd(curDir);


end


