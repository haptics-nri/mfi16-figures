function [nums, ids, centers, p1s, p2s, p3s, p4s] = load_april(filename)
% Loads information about April tags in a set of images from CSV (generated
% by the bluefox processor on the NUC).
%
% Inputs: path+filename to the CSV file.
% Outputs:
%   - frame numbers (vector)
%   - tag IDs (cell array vectors)
%   - tag centers (cell array of 2xM matrices)
%   - tag corners #1 (cell array of 2xM matrices)
%   - tag corners #2 (cell array of 2xM matrices)
%   - tag corners #3 (cell array of 2xM matrices)
%   - tag corners #4 (cell array of 2xM matrices)

    % parse CSV
    table = readtable(filename, ...
                      'FileType','text', ...
                      'Delimiter','comma', ...
                      'Format','%d%s%q%q%q%q%q'); % '%q' means quoted string (with escaped commas)
    
    nums = table.FrameNumber;
    ids = cell(length(nums),1);
    centers = cell(length(nums),1);
    p1s = cell(length(nums),1);
    p2s = cell(length(nums),1);
    p3s = cell(length(nums),1);
    p4s = cell(length(nums),1);
    
    for i=1:length(nums)
        % the file format is contrived so we can use eval() directly
        ids{i} = eval(['[' table.TagIDs{i} ']'])';
        centers{i} = eval(['[' table.TagCenters{i} ']'])';
        p1s{i} = eval(['[' table.TagP1s{i} ']'])';
        p2s{i} = eval(['[' table.TagP2s{i} ']'])';
        p3s{i} = eval(['[' table.TagP3s{i} ']'])';
        p4s{i} = eval(['[' table.TagP4s{i} ']'])';
    end
end
