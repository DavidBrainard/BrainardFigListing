function vScriptsList = bfValidateListAllValidationDirs
% ListBrainardFigsScripts
%
% List the script directories to be validated.
  
% Get the validation rootDir
rootDir = UnitTest.getPref('validationRootDir');
        
vScriptsList = {...
                {fullfile(rootDir, 'scripts', 'AnnReviewColor2015')} ... 
                {fullfile(rootDir, 'scripts', 'OSAVirgina2015')} ... 
            };