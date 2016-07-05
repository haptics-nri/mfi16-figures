%This file allows two data sets to be compared
%Can be easily switched back and forth between this and the "normal"
%   go_visforce (i.e. just with one data set)
function go_visforce_compare(clearstuff)

%If this file is run initially, then clear everything; if called by
%go_visforce (switching to other format) then don't (no sense in making a
%new figure window)
if ~exist('clearstuff','var')
    clear;clc;close all
else blah = clearstuff; %Just to get rid of that highlight, not actually used
end

%Making sure this path is added for future functions
addpath('vlc-matlab')

%Various data needed
mass = 0.1503;
com = [-0.0029   -0.0030    0.0348 ]';
H_vic2bod = [ 0.9912   -0.0236   -0.1303         0
    0.0162    0.9982   -0.0571   36.0685
    0.1314    0.0545    0.9898 -511.6330
    0         0         0    1.0000 ];
H_m402bod = [      0         0    1.0000  108.9900
    1.0000         0         0    0.5300
    0    1.0000         0   -2.9800
    0         0         0    1.0000 ];
H_bal2imu = [ 1.0000         0         0  254.3402
    0    1.0000         0         0
    0         0    1.0000         0
    0         0         0    1.0000 ];

%Creating a window to choose which data set is desired
figPos = get(gcf,'Position');
set(gcf,'Position',[figPos(1:2), 500 400])

%Choosing the first data set
uicontrol('Style','pushbutton', 'Position',[20 150 200 100],...
    'String','Choose Your Data Set Folder','Callback',{@GetData,1})

%Choosing the second data set
uicontrol('Style','pushbutton', 'Position',[280 150 200 100],...
    'String','Choose Your Data Set Folder','Callback',{@GetData,2})

%Textbox to show that the first folder has been selected
editDay = uicontrol('Style','edit','Position',[20 300 200 50],'String','Trial Day');
editTrial = uicontrol('Style','edit','Position',[20 250 200 50], 'String','Trial Type');

%Textbox to show that the second folder has been selected
editDay2 = uicontrol('Style','edit','Position',[280 300 200 50],'String','Trial Day');
editTrial2 = uicontrol('Style','edit','Position',[280 250 200 50], 'String','Trial Type');

%Button to press once data has been selected
uicontrol('Style','pushbutton','Position',[150 50 200 100],'String','Begin Comparison',...
    'Callback',@StartAnalysis)

%Button to quit
uicontrol('Style','pushbutton','Position',[395 5 100 50],'String','Quit',...
    'Callback',@QuitFunc)

%Button to switch to go_visforce (one trial)
uicontrol('Style','pushbutton','Position',[5 5 100 50],'String','One Trial Only',...
    'Callback',@SwitchFunc)


%Initializing some variables
foldername = cell(1,2);
dateFormated = cell(1,2);
trialName = cell(1,2);

    function GetData(~,~,trial)
        %Selects the folder containing the data
        foldername{trial} = [uigetdir, '/'];
        
        %Making sure that a file was actually chosen; if not, the following
        %steps are skipped
        if ~strcmp(foldername{trial},[0 '/'])
            
            %For the date: Finds where the is a string of 8 numbers, sets that
            %as the date. If there is no string of 8 numbers, then "Date
            %Unknown" is displayed
            numStart = regexp(foldername{trial},'\d{8}','once'); %Finds where #s begin
            if ~isempty(numStart)
                dateRaw = foldername{trial}(numStart:numStart+7);
                dateFormated{trial} = [dateRaw(1:4),'-',dateRaw(5:6),'-',dateRaw(7:end)]; %Formatting
            else dateFormated{trial} = 'Date Unknown';
            end
            
            %Changing the date for the correct text box
            if trial == 1
                set(editDay,'String',dateFormated{trial})
            else set(editDay2,'String',dateFormated{trial})
            end
            
            %Putting in the name of the trial
            %Assumes that this "name" is the last part of the file
            %extension (TODO: Change this maybe? Make it a better name?
            slashes = find(foldername{trial}=='/');
            trialName{trial} = foldername{trial}(slashes(end-1)+1: end-1);
            
            %Checks to see if the chosen folder has the appropriate files
            %in it. If not, then "Choose folder with appropriate data
            %files" will appear
            if exist([foldername{trial} 'teensy.ft.csv'],'file') && exist([foldername{trial} 'teensy.acc.csv'],'file')...
                    && exist([foldername{trial} 'teensy.gyro.csv'],'file') && exist([foldername{trial} 'teensy.mag.csv'],'file')
                
                %Changing the trial name for the correct text box
                if trial == 1
                    set(editTrial,'String',trialName{trial})
                else set(editTrial2,'String',trialName{trial})
                end
                
            else
                %Putting the error message in the correct text box
                if trial == 1
                    set(editTrial,'String','Choose folder with appropriate data files')
                else set(editTrial2,'String','Choose folder with appropriate data files')
                end
                foldername{trial} = [];
            end
            
        else foldername{trial} = [];
        end
    end

    function StartAnalysis(~,~)
        %Creates the graphs/analysis
        if ~isempty(foldername{1}) && ~isempty(foldername{2}) %If a data set has been chosen for both steps
            
            %Run the comparison
            %For first trial
            fprintf('\tFor first trial:\n');
            [v1, f1, ~,~,~, a1] = load_stick(foldername{1});
            [~,~,~,~,~,~, v1, ~,vib1,~,~, f1] = ...
                process_stick(v1, f1, a1, mass, com, H_vic2bod, H_m402bod, H_bal2imu, -4.3);
            
            fprintf('\n\tFor second trial:\n');
            %For second trial
            [v2, f2, ~,~,~, a2] = load_stick(foldername{2});
            [~,~,~,~,~,~, v2, ~,vib2,~,~, f2] = ...
                process_stick(v2, f2, a2, mass, com, H_vic2bod, H_m402bod, H_bal2imu, -4.3);
            
            %Comparing
            visforce_compare(f1, v1, [foldername{1} 'video_small.mov'], 20, dateFormated{1}, trialName{1}, vib1,...
                f2, v2, [foldername{2} 'video_small.mov'], 20, dateFormated{2}, trialName{2}, vib2, ...
                mass, com);
            
        end
    end

    function QuitFunc(~,~)
        %Closes out of the window
        close all
    end


    function SwitchFunc(~,~)
        %Switches to the other function (One Trial <--> Compare Trials)
        clf
        go_visforce(1)
    end
end