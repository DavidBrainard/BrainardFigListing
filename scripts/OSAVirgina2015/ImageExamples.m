% IllumExampleImages
%
% Get some example images for use in talks/posters.
%
% This isn't really done, because we want sRGB and not
% the settings for our lab monitors.  And maybe some
% scaling.
%
% 6/7/15  dhb  Wrote it.

%% Clear
clear; close all;

%% Start grabbing them
outputDir = 'xImageExamples';
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

%% Get base dir
imageBaseDir = '/Volumes/ColorShare1/Users/Shared/Matlab/Analysis/BLIlluminationDiscriminationCalcs/Data/ImageData';

%% Get size to pull out
yRegion = 392:390+480;
xRegion = 490:490+480;
sizeFactor = 0.25;

%% Load up target image
theImageData = load(fullfile(imageBaseDir,'Standard','TestImage0.mat'));
theImage = theImageData.sensorImageLeftRGB(yRegion,xRegion,:);
imwrite(imresize(theImage,sizeFactor),fullfile(outputDir,'targetImage.tiff'),'tiff');

%% Blue series
theImageData = load(fullfile(imageBaseDir,'BlueIllumination','blue10L-RGB.mat'));
theImage = theImageData.sensorImageLeftRGB(yRegion,xRegion,:);
imwrite(imresize(theImage,sizeFactor),fullfile(outputDir,'blueImage10.tiff'),'tiff');

theImageData = load(fullfile(imageBaseDir,'BlueIllumination','blue20L-RGB.mat'));
theImage = theImageData.sensorImageLeftRGB(yRegion,xRegion,:);
imwrite(imresize(theImage,sizeFactor),fullfile(outputDir,'blueImage20.tiff'),'tiff');

theImageData = load(fullfile(imageBaseDir,'BlueIllumination','blue30L-RGB.mat'));
theImage = theImageData.sensorImageLeftRGB(yRegion,xRegion,:);
imwrite(imresize(theImage,sizeFactor),fullfile(outputDir,'blueImage30.tiff'),'tiff');

theImageData = load(fullfile(imageBaseDir,'BlueIllumination','blue40L-RGB.mat'));
theImage = theImageData.sensorImageLeftRGB(yRegion,xRegion,:);
imwrite(imresize(theImage,sizeFactor),fullfile(outputDir,'blueImage40.tiff'),'tiff');

theImageData = load(fullfile(imageBaseDir,'BlueIllumination','blue50L-RGB.mat'));
theImage = theImageData.sensorImageLeftRGB(yRegion,xRegion,:);
imwrite(imresize(theImage,sizeFactor),fullfile(outputDir,'blueImage50.tiff'),'tiff');