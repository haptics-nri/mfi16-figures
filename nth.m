function x = nth(m, n, f, varargin)

    s = '[';
    for i=1:m
        s = [s 'a{' num2str(i) '}'];
        if i ~= m
            s = [s ', '];
        end
    end
    s = [s '] = feval(f, varargin{:});'];
    
    eval(s);
    if numel(n) == 1
        x = a{n};
    else
        x = a(n);
    end

end
