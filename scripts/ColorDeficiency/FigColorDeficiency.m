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
outputDir = sprintf('%s_Output',mfilename);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

%% Figure parameters
curDir = pwd;
masterFigParamsDir = fullfile(theProjectRootDir,'figparams');
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
% Let's be pros and show
% unnormalized quantal sensitivities.  And, use
% isetbio to generate the data for the plots.
sensor = sensorCreate('human');
wave   = sensorGet(sensor,'wave');
S_in = WlsToS(wave);
T_conesQE_in = sensorGet(sensor,'spectral qe')';
S = [390 5 65];
wls = SToWls(S);
T_conesQETrichrom = SplineCmf(S_in,T_conesQE_in(2:4,:),S);

%% We also need some anomolous cones
%
% To get these, use PTB routines.  See v_Cones for some validation of how
% these compare to what is returned by isetbio, but it's pretty good.  The
% absorbance based agreement is better than the nomogram based agreement,
% because the SS nomogram doesn't perfectly match up with the tabulated SS
% absorbance.
%
% I didn't actually add any validation checks on these here, because I am
% too lazy to regenerate the validation data for this function.
ptbPhotoreceptorsSS2 = ptb.StockmanSharpePhotoreceptors(SToWls(S));
ptbPhotoreceptorsAnom = DefaultPhotoreceptors('CIE2Deg');
ptbPhotoreceptorsAnom.nomogram.S = S;
ptbPhotoreceptorsAnom.nomogram.source = 'StockmanSharpe';
ptbPhotoreceptorsAnom.nomogram.lambdaMax = [558.9 545 420.7]';
ptbPhotoreceptorsAnom = rmfield(ptbPhotoreceptorsAnom,'absorbance');
ptbPhotoreceptorsAnom = FillInPhotoreceptors(ptbPhotoreceptorsAnom);
T_conesQETrichromAnom = [T_conesQETrichrom(1,:) ; ptbPhotoreceptorsAnom.isomerizationAbsorptance(2,:) ; T_conesQETrichrom(3,:)];

% Energy sensitivities.  We need these to convert the isomerizations
% into the right relative scaling for SRGB rendering.
T_conesTrichrom = EnergyToQuanta(S,T_conesQETrichrom')';
T_conesTrichrom = T_conesTrichrom/max(T_conesTrichrom(:));
for ii = 1:3
    renderingScaleFactors(ii) = 1/max(T_conesTrichrom(ii,:));
end

% Generate hues around a hue circle
% Precompute, and simple test
munsellValue = 5;
munsellChroma = 5;
nHues = 85;
hueAngles = linspace(0,360,nHues+1);
hueAngles = hueAngles(1:end-1);
munsellData = MunsellPreprocessTable;
[~,Xx,trix,vx,Xy,triy,vy,XY,triY,vY] = MunsellGetxyY(0,munsellValue,munsellChroma,munsellData);
[hueCircle_xyY] = MunsellGetxyY(hueAngles,munsellValue,munsellChroma,[],Xx,trix,vx,Xy,triy,vy,XY,triY,vY);

% Compute monitor metamer for spd
weights = inv(T_conesQETrichrom*B)*T_conesQETrichrom*theSpd;
theMetamer = B*weights;
theMetamerCones = T_conesQETrichrom*theMetamer;

%% Make a bar plot of cone responses to the two lights
mid1 = 1; mid2 = 2.5; mid3 = 4.0; sep = 0.25;
figParams.figName = 'FigTrichromHistoTrichomMetamers';
figParams.barWidth = 0.4;
figParams.barLineWidth = 5;
figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
xlim([sep mid3+sep+figParams.barWidth]);
h1 = bar([mid1 - sep]',[theSpdCones(1)],figParams.barWidth,'FaceColor','r','LineWidth',figParams.barLineWidth);
h2 = bar([mid1 + sep]',[theMetamerCones(1)],figParams.barWidth,'FaceColor','r','LineWidth',figParams.barLineWidth,'LineStyle',':');
h3 = bar([mid2 - sep]',[theSpdCones(2)],figParams.barWidth,'FaceColor','g','LineWidth',figParams.barLineWidth);
h4 = bar([mid2 + sep]',[theMetamerCones(2)],figParams.barWidth,'FaceColor','g','LineWidth',figParams.barLineWidth,'LineStyle',':');
h5 = bar([mid3 - sep]',[theSpdCones(3)],figParams.barWidth,'FaceColor','b','LineWidth',figParams.barLineWidth);
h6 = bar([mid3 + sep]',[theMetamerCones(3)],figParams.barWidth,'FaceColor','b','LineWidth',figParams.barLineWidth,'LineStyle',':');
set(gca,'XTick',[mid1 mid2 mid3]);
set(gca,'XTickLabel',{'L' 'Ms' 'S'},'FontName',figParams.fontName,'FontSize',figParams.labelFontSize+4);
ylim([0 1]);
set(gca,'YTick',[0 0.5 1]);
set(gca,'YTickLabel',{' 0.00 ', ' 0.50 ', ' 1.00 '},'FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylabel('Isomerization Rate (arbitrary units)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),gcf,figParams.figType);

% and just to the first light (for explanatory purposes)
figParams.figName = 'FigTrichromHistoTrichomMetamersOne';
figParams.barWidth = 0.4;
figParams.barLineWidth = 5;
figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
xlim([sep mid3+sep+figParams.barWidth]);
h1 = bar([mid1 - sep]',[theSpdCones(1)],figParams.barWidth,'FaceColor','r','LineWidth',figParams.barLineWidth);
h3 = bar([mid2 - sep]',[theSpdCones(2)],figParams.barWidth,'FaceColor','g','LineWidth',figParams.barLineWidth);
h5 = bar([mid3 - sep]',[theSpdCones(3)],figParams.barWidth,'FaceColor','b','LineWidth',figParams.barLineWidth);
set(gca,'XTick',[mid1 mid2 mid3]);
set(gca,'XTickLabel',{'L' 'M' 'S'},'FontName',figParams.fontName,'FontSize',figParams.labelFontSize+4);
ylim([0 1]);
set(gca,'YTick',[0 0.5 1]);
set(gca,'YTickLabel',{' 0.00 ', ' 0.50 ', ' 1.00 '},'FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylabel('Isomerization Rate (arbitrary units)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),gcf,figParams.figType);

% and just L cone response to first light  (again, for explanatory purposes)
figParams.figName = 'FigTrichromHistoTrichomMetamersOneLOnly';
figParams.barWidth = 0.4;
figParams.barLineWidth = 5;
figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
xlim([sep mid3+sep+figParams.barWidth]);
h1 = bar([mid1 - sep]',[theSpdCones(1)],figParams.barWidth,'FaceColor','r','LineWidth',figParams.barLineWidth);
set(gca,'XTick',[mid1 mid2 mid3]);
set(gca,'XTickLabel',{'L' 'M' 'S'},'FontName',figParams.fontName,'FontSize',figParams.labelFontSize+4);
ylim([0 1]);
set(gca,'YTick',[0 0.5 1]);
set(gca,'YTickLabel',{' 0.00 ', ' 0.50 ', ' 1.00 '},'FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylabel('Isomerization Rate (arbitrary units)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),gcf,figParams.figType);

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
plot(wls,T_conesQETrichromAnom(1,:)','r','LineWidth',figParams.lineWidth);
plot(wls,T_conesQETrichromAnom(2,:)','y','LineWidth',figParams.lineWidth);
plot(wls,T_conesQETrichromAnom(3,:)','b','LineWidth',figParams.lineWidth);
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

% and an tetrachomatic version for explanatory purposes
figParams.figName = 'FigTetraCones';
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
plot(wls,T_conesQETetra(1,:)','r','LineWidth',figParams.lineWidth);
plot(wls,T_conesQETetra(4,:)','y','LineWidth',figParams.lineWidth);
plot(wls,T_conesQETetra(2,:)','g','LineWidth',figParams.lineWidth);
plot(wls,T_conesQETetra(3,:)','b','LineWidth',figParams.lineWidth);
xlim([figParams.xLimLow figParams.xLimHigh]);
set(gca,'XTick',figParams.xTicks);
set(gca,'XTickLabel',figParams.xTickLabels);
xlabel('Wavelength (nm)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylim([figParams.yLimLow figParams.yLimHigh]);
ylabel('Quantal Efficiency','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
set(gca,'YTick',figParams.yTicks);
set(gca,'YTickLabel',figParams.yTickLabels);
legend({' L cones ' ' M'' cones ' ' M cones ' 'S cones'},'Location','NorthEast','FontSize',figParams.legendFontSize);
axis('square');
%set(gca,'XMinorTick','on');
FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),theFig,figParams.figType);

%% Make the metamer figure
figParams.figName = 'FigTrichromMetam';
figParams.xLimLow = 380;
figParams.xLimHigh = 720;
figParams.xTicks = [400 500 600 700];
figParams.xTickLabels = {'400' '500' '600' '700'};
figParams.yLimLow = -0.01;
figParams.yLimHigh = 1.01;
figParams.yTicks = [0.0 0.25 0.5 0.75 1.0];
figParams.yTickLabels = {' 0.00 ' ' 0.25 ' ' 0.50 ' ' 0.75 ' ' 1.00 '};
theFig = figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
plot(wls,theSpd,'k','LineWidth',figParams.lineWidth);
plot(wls,theMetamer,'k:','LineWidth',figParams.lineWidth);
xlim([figParams.xLimLow figParams.xLimHigh]);
set(gca,'XTick',figParams.xTicks);
set(gca,'XTickLabel',figParams.xTickLabels);
xlabel('Wavelength (nm)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylim([figParams.yLimLow figParams.yLimHigh]);
ylabel('Power (arbitrary quantal units)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
set(gca,'YTick',figParams.yTicks);
set(gca,'YTickLabel',figParams.yTickLabels);
axis('square');
%set(gca,'XMinorTick','on');
FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),theFig,figParams.figType);

% a version with just the mixture
figParams.figName = 'FigTrichromMetamTwo';
figParams.xLimLow = 380;
figParams.xLimHigh = 720;
figParams.xTicks = [400 500 600 700];
figParams.xTickLabels = {'400' '500' '600' '700'};
figParams.yLimLow = -0.01;
figParams.yLimHigh = 1.01;
figParams.yTicks = [0.0 0.25 0.5 0.75 1.0];
figParams.yTickLabels = {' 0.00 ' ' 0.25 ' ' 0.50 ' ' 0.75 ' ' 1.00 '};
theFig = figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
plot(wls,theMetamer,'k:','LineWidth',figParams.lineWidth);
xlim([figParams.xLimLow figParams.xLimHigh]);
set(gca,'XTick',figParams.xTicks);
set(gca,'XTickLabel',figParams.xTickLabels);
xlabel('Wavelength (nm)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylim([figParams.yLimLow figParams.yLimHigh]);
ylabel('Power (arbitrary quantal units)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
set(gca,'YTick',figParams.yTicks);
set(gca,'YTickLabel',figParams.yTickLabels);
axis('square');
%set(gca,'XMinorTick','on');
FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),theFig,figParams.figType);

% and a version with just the first spectrum
figParams.figName = 'FigTrichromMetamOne';
figParams.xLimLow = 380;
figParams.xLimHigh = 720;
figParams.xTicks = [400 500 600 700];
figParams.xTickLabels = {'400' '500' '600' '700'};
figParams.yLimLow = -0.01;
figParams.yLimHigh = 1.01;
figParams.yTicks = [0.0 0.25 0.5 0.75 1.0];
figParams.yTickLabels = {' 0.00 ' ' 0.25 ' ' 0.50 ' ' 0.75 ' ' 1.00 '};
theFig = figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
plot(wls,theSpd,'k','LineWidth',figParams.lineWidth);
xlim([figParams.xLimLow figParams.xLimHigh]);
set(gca,'XTick',figParams.xTicks);
set(gca,'XTickLabel',figParams.xTickLabels);
xlabel('Wavelength (nm)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylim([figParams.yLimLow figParams.yLimHigh]);
ylabel('Power (arbitrary quantal units)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
set(gca,'YTick',figParams.yTicks);
set(gca,'YTickLabel',figParams.yTickLabels);
axis('square');
%set(gca,'XMinorTick','on');
FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),theFig,figParams.figType);

% a version with metamer components as r, g, and b
figParams.figName = 'FigTrichromMetamComponents';
figParams.xLimLow = 380;
figParams.xLimHigh = 720;
figParams.xTicks = [400 500 600 700];
figParams.xTickLabels = {'400' '500' '600' '700'};
figParams.yLimLow = -0.01;
figParams.yLimHigh = 1.01;
figParams.yTicks = [0.0 0.25 0.5 0.75 1.0];
figParams.yTickLabels = {' 0.00 ' ' 0.25 ' ' 0.50 ' ' 0.75 ' ' 1.00 '};
theFig = figure; clf; hold on
set(gcf,'Position',[100 100 figParams.sqSize figParams.sqSize]);
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
plot(wls,B(:,1)*weights(1)/(max(B(:,1)*weights(1))),'r','LineWidth',figParams.lineWidth);
plot(wls,0.5*B(:,2)*weights(2)/(max(B(:,2)*weights(2))),'g','LineWidth',figParams.lineWidth);
plot(wls,0.5*B(:,3)*weights(3)/(max(B(:,3)*weights(3))),'b','LineWidth',figParams.lineWidth);
xlim([figParams.xLimLow figParams.xLimHigh]);
set(gca,'XTick',figParams.xTicks);
set(gca,'XTickLabel',figParams.xTickLabels);
xlabel('Wavelength (nm)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylim([figParams.yLimLow figParams.yLimHigh]);
ylabel('Power (arbitrary quantal units)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
set(gca,'YTick',figParams.yTicks);
set(gca,'YTickLabel',figParams.yTickLabels);
axis('square');
%set(gca,'XMinorTick','on');
FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),theFig,figParams.figType);

end


