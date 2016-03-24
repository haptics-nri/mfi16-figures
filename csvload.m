function out = csvload(filename, colnames, readtable_args)
    if nargin < 3
        readtable_args = {};
    end

    if exist(filename, 'file')
        data = readtable(filename, 'FileType','text', readtable_args{:});
        cols = zeros(size(colnames));
        for i=1:length(colnames)
            found = find(strcmp(colnames{i}, data.Properties.VariableNames));
            if isempty(found)
                error('Could not find column %s', colnames{i});
            end
            cols(i) = found;
        end

        out = zeros(size(data.(cols(1)), 1), length(cols));
        for i=1:length(cols)
            out(:,i) = data.(cols(i));
        end
    else
        out = [];
    end

end