% construct surface for a magnified+rotated+translated cylinder or cone
% returns two cylinders: x/y/z the original and bx/by/bz the magnified
%   profile: the radius/radii of the cylinder (see docs for CYLINDER)
%   height: height of the transformed cylinder
%   offset: move the origin along the cylinder's own z-axis before rotating
%   axis: vector which the rotated cylinder will be aligned with (need not be normalized)
%   pos: where to translate the cylinder after rotating
%   magnify: scale factor for the magnified cylinder
%   mpos: where to translate the magnified cylinder
function [x, y, z, bx, by, bz] = cylinder_surf(profile, height, offset, axis, pos, magnify, mpos)

    axis = axis/norm(axis);
    [x, y, z] = cylinder(profile);
    z = z*height - offset;
    c = [reshape(x, [1 42]); reshape(y, [1 42]); reshape(z, [1 42])];
    c = xfconv(cross([0 0 1], axis)/norm(cross([0 0 1], axis))*acos(dot([0 0 1], axis)))*c;
    bc = c*magnify;
    x = reshape(c(1,:), [2 21]) + pos(1);
    y = reshape(c(2,:), [2 21]) + pos(2);
    z = reshape(c(3,:), [2 21]) + pos(3);
    bx = reshape(bc(1,:), [2 21]) + mpos(1);
    by = reshape(bc(2,:), [2 21]) + mpos(2);
    bz = reshape(bc(3,:), [2 21]) + mpos(3);

end
