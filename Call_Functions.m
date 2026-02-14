%% MaQiT v.0.1               Author: James M. Lea, 8/2/2018 
% Four methods can be used in this toolbox to calculate glacier terminus
% change. These can be used "out of the box" by simply using the GUI, or
% as functions for use in other scripts.
% The methods that are included in this toolbox are as follows, with
% descriptions of methods 1 and 2 given in Lea et al., 2014,
% 'Evaluation of new and existing methods of tracking glacier terminus
% change', Journal of Glaciology.
%
%       method = 1 (Centreline method)
%       method = 2 (Curvilinear box method)
%       method = 3 (Variable box method)
%       method = 4 (Multi-centreline method)


%% Variables to be defined by the user:
tic
clearvars -except handles hObject
%   1) Terminus change method to be used (see above for options)
    method=handles.method;
    
%       1a) If method 2 is chosen, need to define box width
            box_width=handles.BoxWidth;
%       1b) If method 4 is chose, need to choose a window of days for
%       change to be calculated over. If want all observations included,
%       set change_min to 0.
            change_min=handles.MinGap;
            change_max=handles.MaxGap;


%   2) Path to directory where glacier termini or terminus shapefile is stored
%      (remember to end directory name with a /).
%           NB: File and path names MUST be bounded by single inverted commas '',
%           so that the text appears pink/purple. The script will not work
%           otherwise.
    terminus_path=handles.TerminiShapefilePath;
    terminus_directory=strcat(fileparts(terminus_path),'/');
    
%       2a) Only needs to be defined terminus_format=2 ...
%           Name of terminus shapefile (including .shp file suffix)
            [~,terminus_filename]=fileparts(terminus_path);
            
            terminus_format=2;
    
%   3) Path to directory where centreline shapefile is stored (remember to end
%      directory name with a /)
    centreline_path=handles.CentrelinePath;
    centreline_directory=strcat(fileparts(centreline_path),'/');
    
%	4) Name of centreline shapefile (including .shp file suffix)
    [~,centreline_name]=fileparts(centreline_path);
    
%   6) Write to csv file (1=yes, 0=no)
    csv_path=handles.OutputPath
    if ~isempty(csv_path)
%       6a) Directory where csv file will be written
        csv_directory=strcat(fileparts(csv_path),'/');
%       6b) Name of csv file 
        [~,csv_filename]=fileparts(csv_path);
    end

            
%   7) Plot results on diagram with labels showing meters difference from
%   the most recent observation (1=yes, 0=no)
    plot_results=handles.PlotOutput;

            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Commands to run script. Below this line does not need to be modified.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function reads terminus and shapefile data
[centreline,termini,termini_date]=readGlacierData(centreline_directory,...
    centreline_name,terminus_directory,terminus_filename,terminus_format);

% Calculation of terminus change
if method==1
    Results=CentrelineMethod(centreline,termini,termini_date);
elseif method==2
    Results=CurvilinearBox(centreline,termini,termini_date,box_width);
elseif method==3
    Results=VariableBox(centreline,termini,termini_date);
elseif method==4
    Results=MultiCentrelineMethod(centreline,termini,termini_date,plot_results,change_min,change_max);
end

if ~isempty(csv_directory) % && method~=4
    writeOutputData(Results,csv_directory,csv_filename);
end

if plot_results==1 && method~=4
    plotTermini(Results)
end
    assignin('base','Results',Results);
msgbox('Calculation complete.')
toc
    
    