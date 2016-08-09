% This file demonstrates how to load a dataset and show the position/force
% visualizer. The dataset (which is commented out) is from March 2016.
%
% To get rid of the .mov: vlc_wrapper('cleanup')
% regexp(data{19},'["[a-z]+','match') (for finding the material)
% regexp(data{22},'[1-9]+/[1-9]','match') (for finding the tooling ball
% size)

function go_visforce(clearstuff)

%If this file is run initially, then clear everything; if called by
%go_visforce_compare (switching to other format) then don't (no sense in
%making a new figure window)
if ~exist('clearstuff','var')
    clear;clc;close all
else blah = clearstuff; %Just to get rid of that highlight, not actually used
end

%Making sure this path is added for future functions
addpath('vlc-matlab')

%Various data needed
mass = 0.1503;
com = [-0.0029   -0.0030    0.0348 ]';
% % % % % % % Old transformation matrices
% % % % H_vic2bod = [ 0.9912   -0.0236   -0.1303         0
% % % %               0.0162    0.9982   -0.0571   36.0685
% % % %               0.1314    0.0545    0.9898 -511.6330
% % % %                    0         0         0    1.0000 ];
% % % % H_bal2imu = [ 1.0000         0         0  254.3402   %Really end effector to body frame
% % % %                    0    1.0000         0         0
% % % %                    0         0    1.0000         0
% % % %                    0         0         0    1.0000 ];
% % % % % % % % % % % % % % % % % % % % % %

% % % % % % % % New transformation matrices
H_m402bod = [      0         0    1.0000  108.9900  %Really Mini 40 to IMU
              1.0000         0         0    0.5300
                   0    1.0000         0   -2.9800
                   0         0         0    1.0000 ];
% H_bal2imu = [ 1.0000         0         0  258.5469   %New Sphere calibration
%                     0    1.0000         0         0
%                     0         0    1.0000         0
%                     0         0         0    1.0000 ];
% H_vic2bod = [ 0.9912   0.0162   0.1314         0
%               -0.0236    0.9982   0.0545   36.4497
%               -0.1303    -0.0571    0.9898 -510.5435
%                    0         0         0    1.0000 ];
H_bal2imu = [ 1.0000         0         0  258.2339   %New Sphere calibration
                    0    1.0000         0         0
                    0         0    1.0000         0
                    0         0         0    1.0000 ];
 H_vic2bod = [ 0.9912   -.024   -.1303         0
              .0166    0.9982   -0.0577   25.9908
              0.1315    0.055    0.9898 -577.6098
                   0         0         0    1.0000 ];              
               
% % % % % % % % % % % % % % % % % % % % % % %

%Creating a window to choose which data set is desired
figPos = get(gcf,'Position');
set(gcf,'Position',[figPos(1:2), 500 400])

%Choosing the data set
uicontrol('Style','pushbutton', 'Position',[150 150 200 100],...
    'String','Choose Your Data Set Folder','Callback',@GetData)

%Textbox to show that folder has been selected
editDay = uicontrol('Style','edit','Position',[150 300 200 50],'String','Trial Day');
editTrial = uicontrol('Style','edit','Position',[150 250 200 50], 'String','Trial Type');

%Button to press once data has been selected
uicontrol('Style','pushbutton','Position',[150 50 200 100],'String','Begin Analysis',...
    'Callback',@StartAnalysis)

%Button to quit
uicontrol('Style','pushbutton','Position',[395 5 100 50],'String','Quit',...
    'Callback',@QuitFunc)

%Button to switch to go_visforce_compare (two trials)
uicontrol('Style','pushbutton','Position',[5 5 110 50],'String','Compare Two Trials',...
    'Callback',@SwitchFunc)

%Initializing some variables
foldername = [];
dateFormated = [];
trialName = [];

%Calling the data set
%[v, f, ~,~,~, a] = load_stick('/Users/alex/Documents/research/proton/code/nri/data/20160310/black1stick/');
%[v, f, ~,~,~, a] = load_stick('/Users/sarahallen/Documents/MATLAB/HapticsLab/mfi16-figures/Data/20160310/black1stick/');

%Processing data set
%[~,~,~,~,~,~, v, ~,~,~,~, f] = process_stick(v, f, a, mass, com, H_vic2bod, H_m402bod, H_bal2imu, -2.5887);
%visforce(f, v, '/Users/alex/Documents/research/proton/code/nri/data/20160310/black1stick/video_small.mov', 20, mass, com);
%visforce(f, v, '/Users/sarahallen/Documents/MATLAB/HapticsLab/mfi16-figures/Data/20160310/black1stick/video_small.mov', 20, mass, com);


    function GetData(~,~)
        %Selects the folder containing the data
        foldername = [uigetdir, '/'];
        
        %Making sure that a file was actually chosen; if not, the following
        %steps are skipped
        if ~strcmp(foldername,[0 '/'])
        
        %For the date: Finds where the is a string of 8 numbers, sets that
        %as the date. If there is no string of 8 numbers, then "Date
        %Unknown" is displayed
        numStart = regexp(foldername,'\d{8}','once'); %Finds where #s begin
        if ~isempty(numStart)
            dateRaw = foldername(numStart:numStart+7);
            dateFormated = [dateRaw(1:4),'-',dateRaw(5:6),'-',dateRaw(7:end)]; %Formatting
        else dateFormated = 'Date Unknown';
        end
        set(editDay,'String',dateFormated)
        
        %Putting in the name of the trial
            %Assumes that this "name" is the last part of the file
            %extension (TODO: Change this maybe? Make it a better name?
        slashes = find(foldername=='/');
        trialName = foldername(slashes(end-1)+1: end-1);
        
            %Checks to see if the chosen folder has the appropriate files
            %in it. If not, then "Choose folder with appropriate data
            %files" will appear
           if exist([foldername 'teensy.ft.csv'],'file') && exist([foldername 'teensy.acc.csv'],'file')...
                   && exist([foldername 'teensy.gyro.csv'],'file') && exist([foldername 'teensy.mag.csv'],'file')
               set(editTrial,'String',trialName)
           else set(editTrial,'String','Choose folder with appropriate data files')
               foldername = [];
           end
        
        else foldername = [];
        end
    end

    function StartAnalysis(~,~)
        %Creates the graphs/analysis
        if ~isempty(foldername) %If a data set has been chosen
            
            %Run the analysis
            [v, f, ~,~,~, a] = load_stick(foldername);
            [~,~,~,~,~,~, v, ~,vibration,~,~, f] = ...  
                process_stick(v, f, a, mass, com, H_vic2bod, H_m402bod, H_bal2imu, -4.3);
            visforce(f, v, [foldername 'video_small.mov'], 20, mass, com, dateFormated, trialName,vibration);
        end
    end

    function QuitFunc(~,~)
        %Closes out of the window
        close all
    end

    function SwitchFunc(~,~)
    %Switches to the other function (One Trial <--> Compare Trials)
    clf
    go_visforce_compare(1)
end
end