function eul = quat2eul( q, varargin )
%QUAT2EUL Convert quaternion to Euler angles
%   EUL = QUAT2EUL(Q) converts a unit quaternion rotation into the corresponding 
%   Euler angles. The input, Q, is an N-by-4 matrix containing N quaternions. 
%   Each quaternion represents a 3D rotation and is of the form q = [w x y z], 
%   with a scalar number as the first value. Each element of Q must be a real number.
%   The output, EUL, is an N-by-3 array of Euler rotation angles with each 
%   row representing one Euler angle set. Rotation angles are in radians.
%
%   EUL = QUAT2EUL(Q, SEQ) converts unit quaternion into Euler angles.
%   The Euler angles are specified by the axis rotation sequence, SEQ.
%
%   The default rotation sequence is 'ZYX', where the order of rotation
%   angles is Z Axis Rotation, Y Axis Rotation, and X Axis Rotation.
%
%   The following rotation sequences, SEQ, are supported: 'ZYX' and 'ZYZ'.
%
%   Example:
%      % Calculates Euler angles for a quaternion
%      % By default, the ZYX axis order will be used.
%      q = [sqrt(2)/2 0 sqrt(2)/2 0];
%      eul = quat2eul(q)
%
%      % Calculate the Euler angles for a ZYZ rotation
%      eulZYZ = quat2eul(q, 'ZYZ')
%
%   See also eul2quat

%   Copyright 2014-2015 The MathWorks, Inc.

%#codegen

%robotics.internal.validation.validateNumericMatrix(q, 'quat2eul', 'q', ...
%    'ncols', 4);
%seq = robotics.internal.validation.validateEulerSequence(varargin{:});
if nargin == 2
    seq = varargin{1};
else
    seq = 'ZYX';
end

% Normalize the quaternions
q = normalizeRows(q);

qw = q(:,1);
qx = q(:,2);
qy = q(:,3);
qz = q(:,4);

% Pre-allocate output
eul = zeros(size(q,1), 3, 'like', q);

% The parsed sequence will be in all upper-case letters and validated
switch seq
    case 'ZYX'
        % Cap all inputs to asin to 1, since values >1 produce complex
        % results
        % Since the quaternion is of unit length, this should never happen,
        % but some code generation configuration seem to hit this edge case
        % under some circumstances.
        aSinInput = -2*(qx.*qz-qw.*qy);
        aSinInput(aSinInput > 1) = 1;
        
        eul = [ atan2( 2*(qx.*qy+qw.*qz), qw.^2 + qx.^2 - qy.^2 - qz.^2 ), ...
            asin( aSinInput ), ...
            atan2( 2*(qy.*qz+qw.*qx), qw.^2 - qx.^2 - qy.^2 + qz.^2 )];
        
    case 'ZYZ'
        % Need to convert to intermediate rotation matrix here to avoid
        % singularities
        R = quat2rotm(q);
        eul = rotm2eul(R, 'ZYZ');
end

% Check for complex numbers
if ~isreal(eul)
    eul = real(eul);
end

end
