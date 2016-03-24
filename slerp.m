% adapted from http://www.alecjacobson.com/weblog/?p=981

function vq = slerp(x, v, xq)

    vq = zeros(length(xq), size(v,2));
    
    for i = 1:length(xq)
        q = xq(i);
        if q < x(1)
            error('SLERP: query point below input range');
        elseif q > x(end)
            error('SLERP: query point above input range');
        else
            exact = find(q == x, 1);
            if ~isempty(exact)
                vq(i,:) = v(exact,:);
            else
                j = find(q <= x, 1);
                above = j;
                below = j-1;
                frac = (q - x(below))/(x(above) - x(below));
                vq(i,:) = slerp_one(v(below,:), v(above,:), frac);
            end
        end
    end

end

function c = slerp_one(a,b,t)
    angle = acos( dot(a,b)/(norm(a)*norm(b)) );
    
    % easy degenerate case
    if (0==angle)
        c = a;
    % hard case
    elseif (pi==angle) 
        error('SLERP_ONE: opposite vectors');
    else
        c = (sin((1.0-t)*angle)/sin(angle))*a + (sin(t*angle)/sin(angle))*b; 
    end
end
