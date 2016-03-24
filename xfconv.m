% DWIM converter for various transformation formats
% - rotation vector (3x1, 1x3) <=> rotation matrix (3x3)
% - rot/trans vector (6x1, 1x6) <=> homogeneous matrix (4x4)
function varargout = xfconv(in)

    if numel(in) == 3 && any(size(in) == 1)
        ret = axisangle2mat(in);
    elseif all(size(in) == [3 3])
        ret = mat2axisangle(in);
    elseif numel(in) == 6 && any(size(in) == 1)
        ret = [axisangle2mat(in(4:6)) [in(1) in(2) in(3)]'
               0 0 0                  1                   ];
    elseif all(size(in) == [4 4])
        ret = [in(1:3,4)' mat2axisangle(in(1:3,1:3))];
    else
        throw(MException('xfconv:WrongInputSize', 'Input matrix/vector is the wrong size.'));
    end
    
    if nargout == 1
        varargout{1} = ret;
    elseif nargout == 2 && all(size(ret) == [4 4])
        varargout{1} = ret(1:3,1:3);
        varargout{2} = ret(1:3,4);
    else
        throw(MException('xfconv:WrongOutputSize', 'Wrong number of output arguments.'));
    end

end

function M = axisangle2mat(a)

    theta = norm(a);
    a = a/norm(a);
    a_cross = [    0 -a(3)  a(2)
                a(3)     0 -a(1)
               -a(2)  a(1)     0];
    
    M = eye(3) + a_cross*sin(theta) + a_cross*a_cross*(1-cos(theta));
    
    % force it to be a rotation matrix
    [U,~,V] = svd(M);
    M = U*V';

end

function a = mat2axisangle(M)

    %a = vrrotmat2vec(M);
    %a = a(1:3) * a(4);
    
    %r = logm(M);
    %a = [r(3,2) r(1,3) r(2,1)];
    
    theta = acos(0.5*(trace(M) - 1));
    r = [M(3,2)-M(2,3) M(1,3)-M(3,1) M(2,1)-M(1,2)]/(2*sin(theta));
    a = r*theta;

end