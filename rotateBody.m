function newPoints = rotateBody(oldValues, position, angles, type)
%ROTATEBODY Finds the coordinates of the new, rotated body
% oldValues is the coordinates of the original body
% position is the position of the body ([x, y, z])
% angles the angles of the body (rotation vector)
% type is the type of the body
%   1: cylinder
%   2: sphere
%
% newX, newY, and newZ are the new x, y, and z coordinates of the body

    %Giving new names to some variables
    x = position(1);
    y = position(2);
    z = position(3);
    xcil = oldValues{1};
    ycil = oldValues{2};
    zcil = oldValues{3};

    %Finding the rotation matrix
    rotFull = xfconv([angles(1),angles(2),angles(3)]);

    %Making full transformation (i.e. adding translation + slight rotation)
    angle = 25 ;%68.53; %25
    angleOffset = ... %Rotation around x axis
        [1 0 0; 0 cosd(angle) -sind(angle); 0 sind(angle) cosd(angle)];
    transformation = angleOffset * rotFull;
    transformation(:,4) = [x y z]';
    transformation(4,:) = [0 0 0 1];

    %Applying transformation to each point
    oldVector = [xcil(:) ycil(:) zcil(:) ones(numel(xcil),1)]';
    newValues = transformation * oldVector;

    %Outputting the values through newPoints
    %Slightly different depending on whether it is a cylinder/sphere
    if type == 1 %if it is a cylinder
        newPoints{1} = [newValues(1,1:2:end); newValues(1,2:2:end)];
        newPoints{2} = [newValues(2,1:2:end); newValues(2,2:2:end)];
        newPoints{3} = [newValues(3,1:2:end); newValues(3,2:2:end)];
    else %if it is a sphere
        for j = 1:21     %Try to find a non-for loop way of doing this
            newPoints{1}(j,:) = newValues(1,j:21:end);
            newPoints{2}(j,:) = newValues(2,j:21:end);
            newPoints{3}(j,:) = newValues(3,j:21:end);
        end
    end


end


