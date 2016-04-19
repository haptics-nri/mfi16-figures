% c: center (3 x 1)
% r: radius (scalar)
% pts [optional]: point sets to plot on top (cell array of Nx3 matrices)

function sphereplot(c, r, pts)

    if ~isempty(c)
        [x,y,z] = sphere;
        x = x*r + c(1);
        y = y*r + c(2);
        z = z*r + c(3);
        mesh(x, y, z, 'FaceColor','none');
    end
    hold on
    if nargin > 2
        for i=1:length(pts)
            plot3(pts{i}(:,1), pts{i}(:,2), pts{i}(:,3), '.');
        end
    end
    if ~isempty(c)
        plot3(c(1), c(2), c(3), 'k.', 'markersize',30);
    end
    hold off
    axis equal vis3d
    xlabel X
    ylabel Y
    zlabel Z

end
