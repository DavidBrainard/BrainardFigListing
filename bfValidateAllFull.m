function bfValidateAllFull

    % Use preferences for the 'theProject' project
    UnitTest.usePreferencesForProject('bfScripts', 'reset');

    %% Change some preferences:
    %% Run time error behavior
    % valid options are: 'rethrowExceptionAndAbort', 'catchExceptionAndContinue'
    %UnitTest.setPref('onRunTimeErrorBehavior', 'catchExceptionAndContinue');
    UnitTest.setPref('onRunTimeErrorBehavior', 'rethrowExceptionAndAbort');

    %% Plot generation
    UnitTest.setPref('generatePlots',  false); 
    UnitTest.setPref('closeFigsOnInit', true);
    
    %% Verbosity Level
    % valid options are: 'none', min', 'low', 'med', 'high', 'max'
    UnitTest.setPref('verbosity', 'med');
    
    %% Numeric tolerance for comparison to ground truth data
    UnitTest.setPref('numericTolerance', 500*eps);
    
    %% Whether to plot data that do not agree with the ground truth
    UnitTest.setPref('graphMismatchedData', true);
    
    %% Print current values of the project's prefs
    UnitTest.listPrefs();
    
    %% What to validate
    validateAllDirs = true;
    if (validateAllDirs)
        % List of script directories to validate. Each entry contains a cell array with 
        % with a validation script directory and an optional struct with
        % prefs that override the corresponding project's prefs.
        % At the moment only the 'generatePlots' pref can be overriden.
        
        % Get the validation rootDir
        rootDir = UnitTest.getPref('validationRootDir');

        % List of script directories to validate
        listingScript = UnitTest.getPref('listingScript');
        vScriptsList = eval(listingScript);
        
    else
        % Alternatively, you can provide a list of scripts to validate. 
        % In this case each entry contains a cell array with 
        % with a script name and an optional struct with
        % prefs that override the corresponding isetbioValidation prefs.
        % At the moment only the generatePlots pref can be overriden.
        
        % Get rootDir
        rootDir = UnitTest.getPref('validationRootDir');
        
        vScriptsList = {...
                {fullfile(rootDir, 'scripts', 'AnnReviewColor2015')} ... 
            };
    end
    
    %% How to validate. Uncomment one of the following five options

    % Run a RUN_TIME_ERRORS_ONLY validation session
    % UnitTest.runValidationSession(vScriptsList, 'RUN_TIME_ERRORS_ONLY')
    
    % Run a FAST validation session (comparing SHA-256 hash keys of the data)
    %UnitTest.runValidationSession(vScriptsList, 'FAST');
    
    % Run a FULL validation session (comparing actual data)
    UnitTest.runValidationSession(vScriptsList, 'FULL');
    
    % Run a PUBLISH validation session (comparing actual data and update github wiki)
    % UnitTest.runValidationSession(vScriptsList, 'PUBLISH);
    
    % Run a validation session without a specified mode. You will be
    % promped to select one of the available modes.
    % UnitTest.runValidationSession(vScriptsList);
end