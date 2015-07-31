% SSToolboxBuglet
%
% This was written to reproduce a crash I got when desired contrast of 0.50
% was passed to the the routine that finds the silent substitution
% stimulus.  But it doesn't crash anymore, either here or in the bigger
% program.
%
% 3/23/15   dhb  Wrote it.

%% Clear, define, etc.
close all; ieInit;

%% Compute parameters
desiredSContrast = 0.5;

%% Set up sensor.  These happen to be quantal sensitivities.
S = [400 10 31];
wls = SToWls(S);
T_conesQE = [ ...
    0.0013    0.0046    0.0096    0.0143    0.0200    0.0242    0.0308    0.0463    0.0638    0.0854    0.1262    0.1898    0.2634    0.3167    0.3553    0.3723    0.3816    0.3820  0.3640    0.3424    0.3027    0.2519    0.1947    0.1385    0.0904    0.0553    0.0307    0.0162    0.0081    0.0038    0.0018 ;
    0.0011    0.0041    0.0099    0.0176    0.0283    0.0372    0.0486    0.0720    0.0944    0.1191    0.1643    0.2316    0.3010    0.3382    0.3530    0.3403    0.3139    0.2734    0.2157    0.1599    0.1068    0.0645    0.0360    0.0189    0.0094    0.0046    0.0021    0.0010    0.0005    0.0002    0.0001 ;
    0.0053    0.0211    0.0481    0.0695    0.0839    0.0791    0.0637    0.0512    0.0303    0.0161    0.0091    0.0044    0.0021    0.0009    0.0003    0.0001    0.0000    0.0000     0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000
    ];

%% Define a background spectrum
% The numbers are big because these are in quantal units.
backgroundSpectrum = 1.0e+15 * [ ...
    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433 3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433    3.8433 ...
    ]';
  
% Define an identity basis and expresss the background with
% respect to it
B_primary = 3*max(backgroundSpectrum(:))*eye(S(3));
backgroundPrimary = B_primary\backgroundSpectrum;
    
% Use silent substitution toolbox machinery to find us a
% spectral modulation that isolates the S cones.
whichPrimariesToPin = [1 size(B_primary,1)];
primaryHeadRoom = 0.05;
ambientSpd = zeros(size(B_primary,2),1);
maxPowerDiff = 10000*max(backgroundSpectrum(:));
whichReceptorsToTarget = [3];
whichReceptorsToIgnore = [];
whichReceptorsToMinimize = [];
modulationPrimary = ReceptorIsolate(T_conesQE,whichReceptorsToTarget, whichReceptorsToIgnore, whichReceptorsToMinimize, ...
B_primary, backgroundPrimary, backgroundPrimary, whichPrimariesToPin,...
primaryHeadRoom, maxPowerDiff, desiredSContrast, ambientSpd);

% Check that we got a sensible modulation with desired properties
backgroundReceptors = T_conesQE*(B_primary*backgroundPrimary + ambientSpd);
modulationReceptors = T_conesQE*B_primary*(modulationPrimary - backgroundPrimary);
contrastReceptors = modulationReceptors ./ backgroundReceptors;
fprintf('Maximized contrasts\n');
for n = 1:size(T_conesQE,1)
    fprintf('\tPre-blur unmodified spectra, cone class %d contrast = %0.4f\n',n,contrastReceptors(n));
end
modSpecFig = figure; clf; hold on
plot(wls,modulationPrimary,'r');
plot(wls,backgroundPrimary,'k');
ylim([0 1]);
    
  