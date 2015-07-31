function FigBasicAliasing
% FigBasicAliasing
%
% Show hyperspectral image and resulting cone planes
%
% 3/23/15   dhb  Wrote it.

%% Clear, define, etc.
close all; ieInit;

%% Figure parameters
curDir = pwd;
cd ..
figParams = FigParams;
cd(curDir);
figParams.imType = 'tiff';
figParams.resizeScale = 1;

%% Compute parameters
%
% The isetbio cone default is for 200 linear cones per degree in the fovea.
% This seems high.  Right now I don't have time to figure out how to change
% this, or even know whether I should.  So instead, we'll fix it by
% changing the grating frequency in proportion
isetBioConesPerDegree = 200;
desiredConesPerDegree = 120;
gratingMultiplier = isetBioConesPerDegree/desiredConesPerDegree;
sceneDegrees = 2;
extractDegrees = 1;
DELTAFCNOPTICS = 0;
SONLY = 1;
gratingCpd = 12;
if (DELTAFCNOPTICS)
    opticsStr = 'noblur';
else
    opticsStr = 'blur'
end
if (SONLY)
    sConeStr = 'sConeMod';
else
    sConeStr = 'achromMod';
end
outputSuffix = sprintf('_%dcpd_%s_%s',gratingCpd,opticsStr,sConeStr);

%% Color matching functions
S = [400 10 31];
wls = SToWls(S);
load('T_xyz1931');
T_xyz = SplineCmf(S_xyz1931,T_xyz1931,S);

%% Set up sensor
%
% Do this early to get cone spectral sensitivities
params.rgbDensities = [0.0 0.625 0.325 .05];
sensor = sensorCreate('human',[],params);
sensor = sensorSet(sensor, 'noise flag', 0);
sensor = sensorSet(sensor,'exp time',0.050);
sensor = sensorSet(sensor,'rows',128);
sensor = sensorSet(sensor,'cols',128);
T_conesQE = sensorGet(sensor,'spectral qe')';
T_conesQE = T_conesQE(2:4,:);
T_cones = EnergyToQuanta(wls,T_conesQE')';

%% Create the scene
%
% Choices are:
%   'freqorient'
%   'hyperspectral';
sceneType = 'harmonic';
switch (sceneType)
    case 'freqorient'
        parms.angles = linspace(0,pi/2,5);
        parms.freqs  =  [1,2,4,8,16];
        parms.blockSize = 64;
        parms.contrast  = .8;
        scene = sceneCreate('frequency orientation',parms);
        scene = sceneSet(scene,'wave',wls');
        scene = sceneSet(scene,'fov',sceneDegrees);
        
    case 'harmonic'
        parms.freq = round(gratingMultiplier*gratingCpd*sceneDegrees); parms.contrast = 1; parms.ph = 0;
        parms.ang= 0; parms.row = 400; parms.col = 400; parms.GaborFlag=0;
        scene = sceneCreate('harmonic',parms);
        scene = sceneSet(scene,'wave',wls');
        scene = sceneSet(scene,'fov',sceneDegrees);
        
        % Use a hyperpectral image to make the scene
    case 'hyperspectral'
        
        % Load a hyperspectral image.
        %
        % This should be an m by n by k matrix.
        
        % Get calibration factors
        calFactorsAll = importdata(fullfile('BearFruitGrayB','calibration.txt'));
        calFactors = calFactorsAll(:,2);
        
        imRawSizePixels = 2020;
        for i = 1:length(wls)
            imageNameIn = fullfile('BearFruitGrayB',num2str(wls(i)));
            fid=fopen(imageNameIn, 'r', 'b'); tmpImage = fread(fid, [imRawSizePixels,imRawSizePixels], 'ushort'); fclose(fid);
            tmpImage = imrotate(tmpImage,-90);
            %tmpImage = imresize(tmpImage,[imSizePixels imSizePixels]);
            if (i == 1)
                theHyperspectralImage = zeros(imRawSizePixels, imRawSizePixels,length(wls));
            end
            tmpImg = calFactors(i)*tmpImage;
            theHyperspectralImage(:,:,i) = tmpImg;
            clear tmpImg
        end
        
        % Crop/resize etc
        theHyperspectralImage = theHyperspectralImage(1:1350,1150:end,:);
        rowCropLow = 52; rowCropHigh = 52+383;
        colCropLow = 142; colCropHigh = 142+484;
        theHyperspectralImage = theHyperspectralImage(rowCropLow:rowCropHigh,colCropLow:colCropHigh,:);
        [hyperspectralRows,hyperspectralCols] = size(theHyperspectralImage(:,:,1));
        
        %% Convert energy to quantal units
        theHyperspectralQuantal = zeros(size(theHyperspectralImage));
        for w = 1:length(wls)
            temp = theHyperspectralImage(:,:,w);
            temp = temp(:);
            tempQuantal = EnergyToQuanta(wls(w),temp);
            theHyperspectralQuantal(:,:,w) = reshape(tempQuantal,hyperspectralRows,hyperspectralCols);
        end
        
        %% Make the isetbio multispectral scene
        %
        % We have not got the absolute units correct
        scene = sceneCreate('multispectral');
        scene = sceneSet(scene,'wave',wls');
        scene = sceneSet(scene, 'photons', theHyperspectralQuantal);
        scene = sceneSet(scene,'fov',sceneDegrees);
        vcAddAndSelectObject(scene); sceneWindow;
        
    otherwise
        error('Unknown scene type');
end

%% Muck with modulation direction
% Optional.  Muck with photons to make a cone isolationg modulation
if (SONLY)
    % Get photon image from scene
    photons = sceneGet(scene,'photons');
    
    % Get mean level of the pattern, so that we can extract what
    % the modulation around this is.
    for w = 1:length(wls)
        temp = photons(:,:,w);
        backgroundPhotons(w) = mean(temp(:));
    end
    backgroundPhotons = backgroundPhotons';
    
    % Define an identity basis and expresss the background with
    % respect to it
    B_primary = 3*max(backgroundPhotons(:))*eye(S(3));
    backgroundPrimary = B_primary\backgroundPhotons;
    
    % Use silent substitution toolbox machinery to find us a
    % spectral modulation that isolates the S cones.
    whichPrimariesToPin = [1 size(B_primary,1)];
    primaryHeadRoom = 0.05;
    ambientSpd = zeros(size(B_primary,2),1);
    maxPowerDiff = 10000*max(backgroundPhotons(:));
    whichReceptorsToTarget = [3];
    whichReceptorsToIgnore = [];
    whichReceptorsToMinimize = [];
    desiredContrast = 0.75;
    modulationPrimary = ReceptorIsolate(T_conesQE,whichReceptorsToTarget, whichReceptorsToIgnore, whichReceptorsToMinimize, ...
        B_primary, backgroundPrimary, backgroundPrimary, whichPrimariesToPin,...
        primaryHeadRoom, maxPowerDiff, desiredContrast, ambientSpd);
    
    % Check that we to a sensible modulation with desired
    % properties
    backgroundReceptors = T_conesQE*(B_primary*backgroundPrimary + ambientSpd);
    modulationReceptors = T_conesQE*B_primary*(modulationPrimary - backgroundPrimary);
    contrastReceptors = modulationReceptors ./ backgroundReceptors;
    for n = 1:size(T_conesQE,1)
        fprintf('\t%d: contrast = %0.4f\n',n,contrastReceptors(n));
    end
    figure; clf; hold on
    plot(wls,modulationPrimary,'r');
    plot(wls,backgroundPrimary,'k');
    ylim([0 1]);
    
    % Muck with the photons in the scene to make it the desired
    % spectral modulation.
    modspec = B_primary*(modulationPrimary - backgroundPrimary);
    for w = 1:length(wls)
        temp = photons(:,:,w);
        imageMean = mean(temp(:));
        modulation = temp-imageMean;
        modulation = modulation/max(modulation(:));
        modulation = modulation*modspec(w);
        temp = imageMean+modulation;
        photons(:,:,w) = temp;
    end
    scene = sceneSet(scene,'photons',photons);
end

%% Look at scene
vcAddAndSelectObject(scene); sceneWindow;
        
%% Render the scene as sRGB
tempHyperspectralQuantal = scene.data.photons;
[mScene,nScene,~] = size(tempHyperspectralQuantal);
rowBorder = round(mScene*extractDegrees/(2*gratingMultiplier*sceneDegrees));
colBorder = round(nScene*extractDegrees/(2*gratingMultiplier*sceneDegrees));
tempHyperspectralEnergy = zeros(size(tempHyperspectralQuantal));
for w = 1:length(wls)
    temp = tempHyperspectralQuantal(:,:,w);
    temp = temp(:);
    tempQuantal = QuantaToEnergy(wls(w),temp);
    tempHyperspectralEnergy(:,:,w) = reshape(tempQuantal,mScene,nScene);
end
xyzImage = zeros(mScene,nScene,3);
for i = 1:length(wls)
    xyzImage(:,:,1) = xyzImage(:,:,1) + T_xyz(1,i)*tempHyperspectralEnergy(:,:,i);
    xyzImage(:,:,2) = xyzImage(:,:,2) + T_xyz(2,i)*tempHyperspectralEnergy(:,:,i);
    xyzImage(:,:,3) = xyzImage(:,:,3) + T_xyz(3,i)*tempHyperspectralEnergy(:,:,i);
end
[xyzCal,nX,nY] = ImageToCalFormat(xyzImage); clear xyzImage;
srgbPrimaryCal = XYZToSRGBPrimary(xyzCal); clear xyzCal;
sRGBCal = SRGBGammaCorrect(srgbPrimaryCal); clear srgbPrimaryCal;
sRGBImage = uint8(CalFormatToImage(sRGBCal,nX,nY)); clear SRGBCal;
figure; clf;
imshow(sRGBImage);
imwrite(sRGBImage(rowBorder:end-rowBorder,colBorder:end-colBorder,:),['sRGBScene' outputSuffix '.' figParams.imType],figParams.imType);

%% Make optical image
DELTAFCNOPTICS = true;
oi = oiCreate('human');
if (DELTAFCNOPTICS)
    % Kluge.  Replace OTF with ones to get delta function optics
    oi.optics.OTF.OTF = ones(size(oi.optics.OTF.OTF));
end
oi = oiCompute(oi,scene);
vcAddAndSelectObject(oi); oiWindow;

%% Render the oi as sRGB
tempHyperspectralQuantal = oi.data.photons;
[mOi,nOi,~] = size(tempHyperspectralQuantal);
rowExtraBorder = round((mOi-mScene)/2);
colExtraBorder = round((nOi-nScene)/2);
tempHyperspectralEnergy = zeros(size(tempHyperspectralQuantal));
for w = 1:length(wls)
    temp = tempHyperspectralQuantal(:,:,w);
    temp = temp(:);
    tempQuantal = QuantaToEnergy(wls(w),temp);
    tempHyperspectralEnergy(:,:,w) = reshape(tempQuantal,mOi,nOi);
end

xyzImage = zeros(mOi,nOi,3);
for i = 1:length(wls)
    xyzImage(:,:,1) = xyzImage(:,:,1) + T_xyz(1,i)*tempHyperspectralEnergy(:,:,i);
    xyzImage(:,:,2) = xyzImage(:,:,2) + T_xyz(2,i)*tempHyperspectralEnergy(:,:,i);
    xyzImage(:,:,3) = xyzImage(:,:,3) + T_xyz(3,i)*tempHyperspectralEnergy(:,:,i);
end
[xyzCal,nX,nY] = ImageToCalFormat(xyzImage); clear xyzImage;
srgbPrimaryCal = XYZToSRGBPrimary(xyzCal); clear xyzCal;
sRGBCal = SRGBGammaCorrect(srgbPrimaryCal); clear srgbPrimaryCal;
sRGBImage = uint8(CalFormatToImage(sRGBCal,nX,nY)); clear SRGBCal;
figure; clf;
imshow(sRGBImage);
imwrite(sRGBImage(rowExtraBorder+rowBorder:end-rowBorder-rowExtraBorder,colExtraBorder+colBorder:end-colBorder-colExtraBorder,:),['sRGBOptics' outputSuffix '.' figParams.imType],figParams.imType);


%% Make sensor image
[sensor, ~] = sensorSetSizeToFOV(sensor,sceneDegrees,scene,oi);
sensor = sensorCompute(sensor,oi);
vcAddAndSelectObject(sensor); sensorWindow('scale',1);

%% Demosaic
cfa = sensorGet(sensor,'cfa');
mosaicImage = sensorGet(sensor,'volts');
[imageHeight,imageWidth] = size(mosaicImage);
mask = cfa.pattern;

% Use griddata to do separate submosaic interpolation
%
% Be sure to deal with isetbio's indexing that 1 -> black.
method = 'linear';
[X Y] = meshgrid(1:imageWidth,1:imageHeight);
nCones = zeros(size(T_conesQE,1),1);
for n = 1:size(T_conesQE,1)
    index = find(cfa.pattern == n+1);
    nCones(n) = length(index(:));
    x1 = X(index);
    y1 = Y(index);
    z1 = mosaicImage(index);
    interpImageLMS(:,:,n) = griddata(x1,y1,z1,X,Y,method);
end
for n = 1:size(T_conesQE,1)
    fprintf('Fraction cone class %d = %0.2f\n',n,nCones(n)/sum(nCones));
end

% Render by converting from LMS to sRGB
M_LMSToXYZ = ((T_cones)'\(T_xyz'))';
T_xyzCheck = M_LMSToXYZ*T_cones;
figure; clf; hold on
plot(wls,T_xyz','k');
plot(wls,T_xyzCheck','r');

[interpImageLMSCal,m,n] = ImageToCalFormat(interpImageLMS);
interpImageXYZCal = M_LMSToXYZ*interpImageLMSCal;
interpImageSRGBPrimayCal = XYZToSRGBPrimary(interpImageXYZCal);
interpImageSRGBCal = SRGBGammaCorrect(interpImageSRGBPrimayCal);
interpImageSRGB = uint8(CalFormatToImage(interpImageSRGBCal,m,n));
figure; clf;
imshow(interpImageSRGB);
imwrite(interpImageSRGB(rowBorder:end-rowBorder,colBorder:end-colBorder,:),['sRGBInterp' outputSuffix '.' figParams.imType],figParams.imType);


