function varargout = FigColorDeficiency(varargin)
% Illustrate confusions produced by color deficiency before/after "treatment"
%
% Description:
%    Generate a Munsell hue circle (a la D100 test) and then use an
%    opponent space similarity measure to predict confusions.  Show how
%    these vary as a function of L and M spectral sensitivities.
%
% 6/13/17   dhb  Wrote it

%{
    tbUse({'BrainardLabBase'});
%}
   
%% Clear
clear; close all;

%% Parameters

% Anomolous lambda-max
%anomLambdaMax = [560 559 420.7]';
anomLambdaMax = [535 530 420.7]';
%anomLambdaMax = [540 530 420.7]';
%anomLambdaMax = [545 530 420.7]';
%anomLambdaMax = [560 530 420.7]';

% Hue circle parameters
munsellValue = 6;
munsellChroma = 5;
nHues = 60;

% How l, s and Lum are weighted in distance calc
lsLumWeights = [1 0.70 0.08];

% Do notch filter?
doNotch = 0;

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

%% Normal trichromatic cone sensitivities.
S = [380 4 101]; wls = SToWls(S);
ptbPhotoreceptorsTrichrom = DefaultPhotoreceptors('CIE2Deg');
ptbPhotoreceptorsTrichrom.nomogram.S = S;
ptbPhotoreceptorsTrichrom.nomogram.source = 'StockmanSharpe';
ptbPhotoreceptorsTrichrom.nomogram.lambdaMax = [560 530 420.7]';
ptbPhotoreceptorsTrichrom = rmfield(ptbPhotoreceptorsTrichrom,'absorbance');
ptbPhotoreceptorsTrichrom = FillInPhotoreceptors(ptbPhotoreceptorsTrichrom);
T_conesQETrichrom = [ptbPhotoreceptorsTrichrom.isomerizationAbsorptance(1,:) ; ptbPhotoreceptorsTrichrom.isomerizationAbsorptance(2,:) ; ptbPhotoreceptorsTrichrom.isomerizationAbsorptance(3,:)];

%% We also need some anomolous cones
ptbPhotoreceptorsAnom = DefaultPhotoreceptors('CIE2Deg');
ptbPhotoreceptorsAnom.nomogram.S = S;
ptbPhotoreceptorsAnom.nomogram.source = 'StockmanSharpe';
ptbPhotoreceptorsAnom.nomogram.lambdaMax = anomLambdaMax;
ptbPhotoreceptorsAnom = rmfield(ptbPhotoreceptorsAnom,'absorbance');
ptbPhotoreceptorsAnom = FillInPhotoreceptors(ptbPhotoreceptorsAnom);
T_conesQEAnom = [ptbPhotoreceptorsAnom.isomerizationAbsorptance(1,:) ; ptbPhotoreceptorsAnom.isomerizationAbsorptance(2,:) ; ptbPhotoreceptorsAnom.isomerizationAbsorptance(3,:)];

%% Energy sensitivities.
%
% We need these to actually work with our energy unit inputs.
T_conesTrichrom = EnergyToQuanta(S,T_conesQETrichrom')';
T_conesTrichrom = T_conesTrichrom/max(T_conesTrichrom(:));
T_conesAnom = EnergyToQuanta(S,T_conesQEAnom')';
T_conesAnom = T_conesAnom/max(T_conesAnom(:));

% Generate hues around a hue circle
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
%{
    theFig = figure; clf; hold on
    plot(hueCircle_xyY(1,:),hueCircle_xyY(2,:),'ro','MarkerFaceColor','r','MarkerSize',8);
%}

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

%% Estimates from normal cones to surface reflectance
M_wgtsToLMSTrichrom = T_conesTrichrom*diag(spd_CIEC)*B_sur;
M_LMSTrichromToSur = B_sur*inv(M_wgtsToLMSTrichrom);

%% Conversions from cone to CIE XYZ.
%
% This isn't exact because we're using nomogram-based cone
% fundamentals.
load T_xyzCIEPhys2
T_xyzPhys = SplineCmf(S_xyzCIEPhys2,T_xyzCIEPhys2,S);
M_ConesTrichromToXYZ = ((T_conesTrichrom')\(T_xyzPhys'))';
%{
    T_check = M_ConesTrichromToXYZ*T_conesTrichrom;
    figure; clf; hold on
    plot(T_check','r');
    plot(T_xyzPhys','g');
%}

%% Convert surfaces to LMS etc
%
% The ls and DKL spaces are not exactly according to standard conventions,
% but what I want here is an opponent + lum representation and this seems
% good enough.
white_XYZ = T_xyz*spd_CIEC;
hueCircle_LMSTrichrom = T_conesTrichrom*diag(spd_CIEC)*hueCircle_sur;
hueCircle_lsTrichrom = LMSToMacBoyn(hueCircle_LMSTrichrom);
hueCircle_LumTrichrom = 2*hueCircle_LMSTrichrom(1,:)+hueCircle_LMSTrichrom(2,:);
hueCircle_DKLTrichrom = [hueCircle_lsTrichrom ; hueCircle_LumTrichrom];
hueCircle_XYZTrichrom = M_ConesTrichromToXYZ*hueCircle_LMSTrichrom;
hueCircle_LabTrichrom = XYZToLab(hueCircle_XYZTrichrom,white_XYZ);
hueCircle_RGBTrichrom = double(SRGBGammaCorrect(XYZToSRGBPrimary(hueCircle_XYZTrichrom)))/255;

hueCircle_LMSAnom = T_conesAnom*diag(spd_CIEC)*hueCircle_sur;
hueCircle_lsAnom = LMSToMacBoyn(hueCircle_LMSAnom);
hueCircle_LumAnom = 2*hueCircle_LMSAnom(1,:)+hueCircle_LMSAnom(2,:);
hueCircle_DKLAnom = [hueCircle_lsAnom ; hueCircle_LumAnom];
hueCircle_XYZAnom = M_ConesTrichromToXYZ*hueCircle_LMSAnom;
hueCircle_LabAnom = XYZToLab(hueCircle_XYZAnom,white_XYZ);
hueCircle_RGBAnom = double(SRGBGammaCorrect(XYZToSRGBPrimary(hueCircle_XYZAnom)))/255;

%% Stimulus distances in weighted DKL
%
% Compute dissimilarity matrix. Start with interpoint distances.
stimulusDistances_Trichrom = zeros(nHues,nHues);
stimulusDistances_Anom = zeros(nHues,nHues);
for i = 1:nHues
    for j = 1:nHues
         stimulusDistances_Trichrom(i,j) = norm(diag(lsLumWeights)*(hueCircle_DKLTrichrom(:,i)-hueCircle_DKLTrichrom(:,j)));
         stimulusDistances_Anom(i,j) = norm(diag(lsLumWeights)*(hueCircle_DKLAnom(:,i)-hueCircle_DKLAnom(:,j)));
    end
end

%% LS plane plot of hue circle, with lines to closet neighbors for anom cones.
theFig = figure; clf; hold on
set(gcf,'Position',[100 100 1100 650]);
markerSize = 8;
subplot(1,2,2); hold on
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
for hh = 1:nHues
    [~,index] = sort(stimulusDistances_Anom(hh,:));
    closestIndex = index(2);
    plot([hueCircle_lsTrichrom(1,hh) hueCircle_lsTrichrom(1,closestIndex)],[hueCircle_lsTrichrom(2,hh) hueCircle_lsTrichrom(2,closestIndex)],'k','LineWidth',2);
    closestIndex = index(3);
    plot([hueCircle_lsTrichrom(1,hh) hueCircle_lsTrichrom(1,closestIndex)],[hueCircle_lsTrichrom(2,hh) hueCircle_lsTrichrom(2,closestIndex)],'k','LineWidth',2); 
end
for hh = 1:nHues
    plot(hueCircle_lsTrichrom(1,hh),hueCircle_lsTrichrom(2,hh),'o','Color',hueCircle_RGBTrichrom(:,hh),'MarkerFaceColor',hueCircle_RGBTrichrom(:,hh),'MarkerSize',markerSize);
end
xlim([0.60 0.75]); ylim([0 0.2]);
set(gca,'XTick',[0.6 0.65 0.7 0.75]); set(gca,'XTickLabel',{'0.60' '0.65' '0.70' '0.75'});
set(gca,'YTick',[0 0.05 0.1 0.15 0.2]); set(gca,'YTickLabel',{'0.00' '0.05' '0.10' '0.15' '0.20'});
xlabel('red-green','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylabel('blue-yellow','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
title({'Distance-Based' ; 'Confusions'},'FontName',figParams.fontName,'FontSize',figParams.titleFontSize-2);
axis('square');

subplot(1,2,1); hold on
set(gca,'FontName',figParams.fontName,'FontSize',figParams.axisFontSize,'LineWidth',figParams.axisLineWidth);
for hh = 1:nHues
    plot(hueCircle_lsAnom(1,hh),hueCircle_lsAnom(2,hh),'o','Color',hueCircle_RGBTrichrom(:,hh),'MarkerFaceColor',hueCircle_RGBTrichrom(:,hh),'MarkerSize',markerSize);
end
xlim([0.60 0.75]); ylim([0 0.2]);
set(gca,'XTick',[0.6 0.65 0.7 0.75]); set(gca,'XTickLabel',{'0.60' '0.65' '0.70' '0.75'});
set(gca,'YTick',[0 0.05 0.1 0.15 0.2]); set(gca,'YTickLabel',{'0.00' '0.05' '0.10' '0.15' '0.20'});
xlabel('red-green','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylabel('blue-yellow','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
title({'Cone-Based' ; 'Chromaticity'},'FontName',figParams.fontName,'FontSize',figParams.titleFontSize-2);
axis('square');
figParams.figName = sprintf('LS_%0.0f_%0.0f_%d_%d',ptbPhotoreceptorsAnom.nomogram.lambdaMax(1),ptbPhotoreceptorsAnom.nomogram.lambdaMax(2),nHues,doNotch);
FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),theFig,figParams.figType);

%% Make the cone sensitivity figure
figParams.figName = 'TrichromCones';
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
FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),theFig,figParams.figType);

%% And an anomalous version for explanatory purposes
figParams.figName = sprintf('AnomCones_%0.0f_%0.0f_%d',ptbPhotoreceptorsAnom.nomogram.lambdaMax(1),ptbPhotoreceptorsAnom.nomogram.lambdaMax(2),doNotch);
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
if (anomLambdaMax(1) == 560)
    plot(wls,T_conesQEAnom(1,:)','r','LineWidth',figParams.lineWidth);
    plot(wls,T_conesQEAnom(2,:)','y','LineWidth',figParams.lineWidth);
    plot(wls,T_conesQEAnom(3,:)','b','LineWidth',figParams.lineWidth);
else
    plot(wls,T_conesQEAnom(1,:)','y','LineWidth',figParams.lineWidth);
    plot(wls,T_conesQEAnom(2,:)','g','LineWidth',figParams.lineWidth);
    plot(wls,T_conesQEAnom(3,:)','b','LineWidth',figParams.lineWidth);
end
xlim([figParams.xLimLow figParams.xLimHigh]);
set(gca,'XTick',figParams.xTicks);
set(gca,'XTickLabel',figParams.xTickLabels);
xlabel('Wavelength (nm)','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
ylim([figParams.yLimLow figParams.yLimHigh]);
ylabel('Quantal Efficiency','FontName',figParams.fontName,'FontSize',figParams.labelFontSize);
set(gca,'YTick',figParams.yTicks);
set(gca,'YTickLabel',figParams.yTickLabels);
if (anomLambdaMax(1) == 560)
    legend({' L' ' L''' ' S'},'Location','NorthEast','FontSize',figParams.legendFontSize);
else
    legend({' M''' ' M' ' S'},'Location','NorthEast','FontSize',figParams.legendFontSize);
end
axis('square');
FigureSave(fullfile(outputDir,[mfilename '_' figParams.figName]),theFig,figParams.figType);

% %% 3D plots in DKL
% figure; clf; hold on
% for hh = 1:nHues
%     plot3(hueCircle_DKLTrichrom(1,hh),hueCircle_DKLTrichrom(2,hh),hueCircle_DKLTrichrom(3,hh), ...
%         'o','Color',hueCircle_RGBTrichrom(:,hh),'MarkerFaceColor',hueCircle_RGBTrichrom(:,hh),'MarkerSize',8);
% end
% xlabel('l'); ylabel('s'); zlabel('Lum');
% title('Normal Trichrom DKL');
% 
% figure; clf; hold on
% for hh = 1:nHues
%     plot3(hueCircle_DKLAnom(1,hh),hueCircle_DKLAnom(2,hh),hueCircle_DKLAnom(3,hh), ...
%         'o','Color',hueCircle_RGBAnom(:,hh),'MarkerFaceColor',hueCircle_RGBAnom(:,hh),'MarkerSize',8);
% end
% xlabel('l'); ylabel('s'); zlabel('Lum');
% title('Anomolous DKL');
% 
% figure; clf; hold on
% for hh = 1:nHues
%     plot3(hueCircle_LMSTrichrom(1,hh),hueCircle_LMSTrichrom(2,hh),hueCircle_LMSTrichrom(3,hh), ...
%         'o','Color',hueCircle_RGBTrichrom(:,hh),'MarkerFaceColor',hueCircle_RGBTrichrom(:,hh),'MarkerSize',8);
% end
% xlabel('L'); ylabel('M'); zlabel('S');
% title('Normal Trichrom LMS');
% 
% figure; clf; hold on
% for hh = 1:nHues
%     plot3(hueCircle_LMSAnom(1,hh),hueCircle_LMSAnom(2,hh),hueCircle_LMSAnom(3,hh), ...
%         'o','Color',hueCircle_RGBAnom(:,hh),'MarkerFaceColor',hueCircle_RGBAnom(:,hh),'MarkerSize',8);
% end
% xlabel('L'); ylabel('M'); zlabel('S');
% title('Anomolous LMS');

cd(curDir);


end


