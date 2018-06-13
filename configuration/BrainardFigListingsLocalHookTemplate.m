function BrainardFigListingsLocalHook
% BrainardFigListingsLocalHook
%
% Configure things for working on the BrainardFigListings project.
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by defalut,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUseProject('BrainardFigListings') to set up for
% this project.  You then edit your local copy to match your configuration.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

%% Say hello
fprintf('Running BrainardFigListings local hook\n');

%% Specify project name and location
projectName = 'BrainardFigListings';
projectBaseDir = tbLocateProject('BrainardFigListings');
if (ispref(projectName))
    rmpref(projectName);
end

% Specify project-specific preferences
p = struct(...
    'projectName',           'bfScripts', ...                                                                                       % The project's name (also the preferences group name)
    'validationRootDir',     fullfile(projectBaseDir), ...                                                                          % Directory location where the 'scripts' subdirectory resides.
    'alternateFastDataDir',  '',  ...                                                                                               % Alternate FAST (hash) data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/fast
    'alternateFullDataDir',  '/Volumes/Users1/Shared/Matlab/Analysis/BrainardFigListings/data/full', ...                            % Alternate FULL data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/full
    'clonedWikiLocation',    '/Users/Shared/GitWebSites/BrainardFigListings.wiki', ...                                              % Local path to the directory where the wiki is cloned. Only relevant for publishing tutorials.
    'clonedGhPagesLocation', '/Users/Shared/GitWebSites/BrainardFigListings', ...                                                   % Local path to the directory where the gh-pages repository is cloned. Only relevant for publishing tutorials.
    'githubRepoURL',         'https://davidbrainard.github.io/BrainardFigListings', ...                                             % Github URL for the project. This is only used for publishing tutorials.
    'generateGroundTruthDataIfNotFound', true, ...                                                                                  % Flag indicating whether to generate ground truth if one is not found
    'listingScript',         'bfValidateListAllValidationDirs', ...                                                                 % Routine that lists validation dirs
    'masterFigParamsDir',    fullfile(projectBaseDir,'figparams') ...                                                               ß% Specific for this project.
    );

generatePreferenceGroup(p);

UnitTest.usePreferencesForProject(p.projectName);

end

function generatePreferenceGroup(p)
% remove any existing preferences for this project
if ispref(p.projectName)
    rmpref(p.projectName);
end

% generate and save the project-specific preferences
setpref(p.projectName, 'projectSpecificPreferences', p);
fprintf('Generated and saved preferences specific to the ''%s'' project.\n', p.projectName);
end