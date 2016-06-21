function flow = parse_flow(filename)

    flow = struct('title', '', ...
                  'stamp', '', ...
                  'steps', struct('name',{}, 'started',{}, 'actions', {}), ...
                  'messages', struct('text',{}, 'stamp',{}), ...
                  'answers', containers.Map);
    
    fid = fopen(filename);
    line = fgetl(fid);
    while ischar(line)
        line = strtrim(line);
        
        if ~isempty(line)
            tokens = regexp(line, '^(.*?) \[([^\]]+)\]$', 'tokens');
            text = tokens{1}{1};
            stamp = datestr(str2double(tokens{1}{2})/86400 + datenum(1970,1,1));
            
            if isempty(flow.title) % haven't read title yet: must be here
                flow.title = text;
                flow.stamp = stamp;
            elseif text(1) == '-' % beginning of a step
                i = length(flow.steps) + 1;
                flow.steps(i).name = text(3:end);
                flow.steps(i).started = stamp;
                flow.steps(i).actions = struct('name',{}, 'started',{});
            elseif text(1) == '>' % it's an answered question
                tokens = regexp(text, '^> "([^"]+)" \["([^"]+)"\]$', 'tokens');
                flow.answers(tokens{1}{1}) = struct('text', tokens{1}{2}, 'stamp', stamp);
            elseif text(1) == '"' % it's a message
                i = length(flow.messages) + 1;
                flow.messages(i).text = text(2:end-1);
                flow.messages(i).stamp = stamp;
            else % must be an action
                i = length(flow.steps);
                j = length(flow.steps(i).actions) + 1;
                flow.steps(i).actions(j).name = text;
                flow.steps(i).actions(j).started = stamp;
            end
        end
        
        line = fgetl(fid);
    end
    fclose(fid);

end
