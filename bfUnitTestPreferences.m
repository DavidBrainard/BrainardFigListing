% Method to set project-specific preferences. Generally, this script should
% be run once only. For different projects, copy this file to the projects' directory
% and adapt the p-struct according to that project's specifics.
%
% You can just run the distributed (Template) version to accept the
% defaults, or you can make a copy outside of the distribution and modify
% for your own system.
%
% To adapt this file to your own project, replace 'theProject', with your project's name
% and similarly for any other fields where 'theProject' appears

function bfUnitTestPreferences

    % Specify root directory
    theProjectRootDir = '/Users/Shared/Matlab/Analysis/BrainardFigListings';
    
    % I (Nicolas) added the following 4 lines for debugging purposes
    userName =  char(java.lang.System.getProperty('user.name'));
    if (strcmp(userName, 'nicolas'))
        theProjectRootDir = '/Users1/DropboxUPenn/Dropbox/xBrainardFigs';
    end
    
    % Specify project-specific preferences
    p = struct(...
            'projectName',           'bfScripts', ...                                                                                       % The project's name (also the preferences group name)
            'validationRootDir',     fullfile(theProjectRootDir), ...                                                                       % Directory location where the 'scripts' subdirectory resides.
            'alternateFastDataDir',  '',  ...                                                                                               % Alternate FAST (hash) data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/fast
            'alternateFullDataDir',  '/Volumes/Users1/Shared/Matlab/Analysis/BrainardFigListings/data/full', ...                            % Alternate FULL data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/full
            'clonedWikiLocation',    '/Users/Shared/GitWebSites/BrainardFigListings.wiki', ...                                             % Local path to the directory where the wiki is cloned. Only relevant for publishing tutorials.
            'clonedGhPagesLocation', '/Users/Shared/GitWebSites/BrainardFigListings/BrainardFigListings', ...                               % Local path to the directory where the gh-pages repository is cloned. Only relevant for publishing tutorials.
            'githubRepoURL',         'https://github.com/DavidBrainard/BrainardFigListings', ...                                            % Github URL for the project. This is only used for publishing tutorials.
            'generateGroundTruthDataIfNotFound', true, ...                                                                                  % Flag indicating whether to generate ground truth if one is not found
            'listingScript',         'bfValidateListAllValidationDirs', ...                                                                  % Routine that lists validation dirs
            'masterFigParamsDir',    fullfile(theProjectRootDir,'figparams') ...                                                            % Specific for this project.
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