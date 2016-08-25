function [pos, q] = estimate_pose_pnp(sensor)
%   ESTIMATE_POSE 6DOF pose estimator based on apriltags
%   This function returns the positon and orientation of IMU with respect
%   to the world frame at an instant
%   
%   INPUTS:
%   sensor - struct stored in provided dataset, fields include
%          
%          - id: 1xn ids of detected tags
%          - p0, p1, p2, p3, p4: 2xn pixel position of center and
%                                four corners of detected tags
%            Y
%            ^ P3 == P2
%            | || P0 ||
%            | P4 == P1
%            o---------> X
%   
%   OUTPUTS:
%   pos - 3x1 position of the IMU in world frame 
%   q   - 4x1 quaternion of the IMU [w, x, y, z] 
    
    persistent quat position

    space_x = 6.71; %% space between consecutive tags along x-axis (mm) 
    space_y = 9.75; %% space between consecutive tags along y-axis (mm)
    gap = 107.52;   %% space between 5th and 6th column tags along x-axis (mm)
    % gap = 109.21; % new frame
    side = 22.71;   %% side of the April tag (mm)

    %%% Inverse of the camera matrix
    Kinv = inv([1796.43667267279, 0, 0; 0, 1798.87345781109, 0; 754.455028480654, 675.072500116700, 1]');

    %%% Conversion between RGB camera and IMU
    angle2 = 14.41757172;
    R_rgb2imu =  [0, -cosd(angle2),  sind(angle2);
                  1,    0,       0;
                  0, sind(angle2), cosd(angle2)];

    T_imu = [109.94483315; 0; -393.4643504];

    cam2imu = R_rgb2imu;
    XYZ = -cam2imu'*T_imu;
    
    %%% April tag IDs    
    UniqueIDs = 1:1:80; %% IDs of tags which have to be used in the algorithm
    [R,C] = ind2sub([8,10],UniqueIDs);
    IDs = sensor.id +1;
    
    extraIDs = 81:86; %% The tags which are not supposed to be used and are on the frame
    [IDs,idx_useful] = setdiff(IDs,extraIDs);
    
    if isempty(extraIDs)
        idx_useful = 1:length(IDs);
    end
    
    if ~isempty(IDs)
        
        %%% Assigning the world frame coordinates to the observed tags
        [~,idx1] = ismember(IDs,UniqueIDs);
        R = R(idx1);
        C = C(idx1);
        Xw_P0 = space_x/2 + (C-1)*(space_x + side);
        Yw_P0 = space_y/2+ (R-1)*(space_y + side);

        idx2 = find(C>5);
        Xw_P0(idx2) = Xw_P0(idx2)+ gap;

        Xw_P1 = Xw_P0 - side/2;
        Yw_P1 = Yw_P0 - side/2;

        Xw_P2 = Xw_P0 - side/2;
        Yw_P2 = Yw_P0 + side/2;

        Xw_P3 = Xw_P0 + side/2;
        Yw_P3 = Yw_P0 + side/2;

        Xw_P4 = Xw_P0 + side/2;
        Yw_P4 = Yw_P0 - side/2;
        
        Pw = [Xw_P0, Xw_P1, Xw_P2, Xw_P3, Xw_P4; Yw_P0, Yw_P1, Yw_P2, Yw_P3, Yw_P4];
        
        %%%% Adding a rotation noise to the world frame coordinates
%         angle_z = -0.006;
%         R_align_z = [cos(angle_z), sin(angle_z), 0;
%                    -sin(angle_z), cos(angle_z), 0;
%                    0         ,          0, 1];
%         
%         angle_y = -0.0294;
%         R_align_y = [cos(angle_y), 0, -sin(angle_y);
%                      0, 1,  0;
%                     sin(angle_y),0, cos(angle_y)];
%             
%         angle_x = -0.016;
%         R_align_x = [1,            0,           0;
%                      0, cos(angle_x), sin(angle_x);
%                      0, -sin(angle_x), cos(angle_x)];  
%                  
%         Pw = R_align_z*R_align_x*R_align_y*[Pw;zeros(1,size(Pw,2))];
%         Pw = Pw(1:2,:);
        %%%%%%
                
        Pp = [sensor.p0(1,idx_useful), sensor.p1(1,idx_useful), sensor.p2(1,idx_useful), sensor.p3(1,idx_useful), sensor.p4(1,idx_useful);
              sensor.p0(2,idx_useful), sensor.p1(2,idx_useful), sensor.p2(2,idx_useful), sensor.p3(2,idx_useful), sensor.p4(2,idx_useful)];  
   

        Pw = Pw';
        Pp = Pp';
        num = size(Pw,1);
        
        %%% PNP Algorithm

        D1 = [-Pw, -ones(num,1), zeros(num,3), Pw(:,1).*Pp(:,1), Pw(:,2).*Pp(:,1), ones(num,1).*Pp(:,1)];
        D2 = [zeros(num,3), -Pw, -ones(num,1), Pw(:,1).*Pp(:,2), Pw(:,2).*Pp(:,2), ones(num,1).*Pp(:,2)];
        D = [D1;D2];

        [V,~] = eig(D'*D);
        A = V(:,1);
        
        % The alternative method of doing what is done by above two lines 
        
%         [U,S,V] = svd(D);        
%         A = V(:,end);
        H = [A(1),A(2),A(3);A(4),A(5),A(6);A(7),A(8),A(9)];
        h = Kinv*H;

        h2 = [h(:,1),h(:,2),cross(h(:,1),h(:,2))];
        [U2,~,V2] = svd(h2);

        R = U2*[1,0,0;0,1,0;0,0,det(U2*V2')]*V2';
        T = h(:,3)/norm(h(:,1));  
        
        pos = [R', -R'*T;0,0,0,1]*[XYZ;1];
        pos = pos(1:3);
        
        FinalR = R'*cam2imu';
        
        %%% Sometimes you may end up with a solution which gives your position of camera at negative z-axis   
        if pos(3)<0
            pos = [[-1,0,0;0,-1,0;0,0,1]*R', [-1,0,0;0,-1,0;0,0,1]*R'*T;0,0,0,1]*[XYZ;1];
            pos = pos(1:3);
            FinalR = [-1,0,0;0,-1,0;0,0,1]*FinalR;            
        end

        q = rotm2quat(FinalR);
        q = q';          
           
        quat = q;
        position = pos;
        
    else %% If no april tag is detected, return the position and orientation from previous iteration
                
        pos = position;
        q = quat;
    end

end