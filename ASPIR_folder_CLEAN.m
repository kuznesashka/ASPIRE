load cases_files_191022.mat


resultsdir_root = [hdisk 'Valerii\45_cases\'];
for case_n =  [1:4 7:length(cases_files.file_names_short)] 
for case_n =  [18:length(cases_files.file_names_short)] 
    
% Specify the folder where the files live.
subj_name = cases_files.cases{case_n};
myFolder = [resultsdir_root, subj_name, '\ASPIRE\'];
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
   
else
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.fig'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now deleting %s\n', fullFileName);
  delete(fullFileName);
end
end
end