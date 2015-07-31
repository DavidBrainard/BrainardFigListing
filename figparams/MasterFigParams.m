function figParams = FigParams
% figParams = FigParams
%
% Set figure parameters.
%
% 4/29/14  dhb  Wrote it.

figParams.fontName = 'Helvetica';
figParams.markerSize = 26;
figParams.lineWidth = 8;
figParams.axisLineWidth = 3;
figParams.dashedLineWidth = 8;
figParams.axisFontSize = 32;
figParams.labelFontSize = 36;
figParams.titleFontSize = 36;
figParams.legendFontSize = 28;
figParams.figType = {'pdf'};
figParams.figDir = 'FigureOutput';
figParams.sqSize = 675;
figParams.size = 800;
figParams.imType = 'tiff';

% These parameters have to do with generating the same spectrum across
% multiple scripts.
figParams.whichMCCSquare = 22;
figParams.whichMCCSquare2 = 19;
figParams.whichMCCSquare3 = 21;
figParams.mccSquareScale = 0.6*0.75;
figParams.mccSquareScale2 = 0.65*0.25;
figParams.mccSquareScale3 = -0.65*0.25;
figParams.spectralConstant = 0.05;

% These have to do with iamges for looking at
figParams.imageExtractPixels = 20;
figParams.imageRepFactor = 15;
figParams.mosaicImageSytle = 'standard';
