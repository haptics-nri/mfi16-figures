clear all
close all

%cd Data
SubFilename = 'absPlastic_Subject1.mat';
load(SubFilename);
%cd ..

Ts = 1/10000;
animate = true;

Aconvx = 10.3734*9.81; %conversion factor from V to m/s^2
Aconvy = 10.3896*9.81;
Aconvz = 10.3950*9.81;

aX = data(:,1); %get x, y, and z accelerations (need to get conversion to g's)
aX = (aX-mean(aX))*Aconvx;
aY = data(:,2);
aY = (aY-mean(aY))*Aconvy;
aZ = data(:,3);
aZ = (aZ-mean(aZ))*Aconvz;
time = nitime(1,:); %time vector for acceleration (Ts = 1/10000)
time = time-time(1);
timePos = trackedtime-trackedtime(1); %time vector for position
Tspos = max(timePos)/length(timePos);

vst = [119.65;-3.30;2.03];

A = unwrap(trackedorient(1,1:end)*pi/180);
E = unwrap(trackedorient(2,1:end)*pi/180);
R = unwrap(trackedorient(3,1:end)*pi/180);
posWorld = trackedpos;

FT=data(:,4:9)-data(:,10:15);
FT=FT.';
FT=FT-diag(FTbias)*ones(size(FT));
workmatrix= [-0.0148464530000000,-0.0528270810000000,-0.128986314000000,-3.30285491500000,0.260661584000000,3.21999616900000;0.0362472370000000,3.95397478000000,-0.0378257670000000,-1.96783033200000,-0.160302624000000,-1.83848942100000;3.68742806500000,0.0353190710000000,3.89915414400000,0.164660002000000,3.57345682600000,-0.0243610810000000;0.706659571000000,24.2776279500000,20.7647702100000,-11.1733273100000,-21.0068402500000,-11.1552407700000;-23.8273414000000,0.0191202520000000,13.7248583200000,21.0056583200000,9.70693341200000,-19.9421094100000;0.215633100000000,13.5789013200000,0.736703546000000,14.1776746400000,1.02066334000000,14.1765185000000;];
v=workmatrix*FT;
v=v.';
fX = interp1(time,v(:,3),timePos,'spline')';
fY = interp1(time,v(:,2),timePos,'spline')';
fZ = interp1(time,-v(:,1),timePos,'spline')';

force = [fX fY fZ];

V = zeros(3,length(A));
T = zeros(3,length(A));
F = zeros(length(A),3);

%Convert sensor position to tooltip position
for k = 1:length(A)
    R11 = cos(E(k))*cos(A(k));
    R12 = cos(E(k))*sin(A(k));
    R13 = -sin(E(k));
    R21 = -cos(R(k))*sin(A(k))+sin(R(k))*sin(E(k))*cos(A(k));
    R22 = cos(R(k))*cos(A(k))+sin(R(k))*sin(E(k))*sin(A(k));
    R23 = sin(R(k))*cos(E(k));
    R31 = sin(R(k))*sin(A(k))+cos(R(k))*sin(E(k))*cos(A(k));
    R32 = -sin(R(k))*cos(A(k))+cos(R(k))*sin(E(k))*sin(A(k));
    R33 = cos(R(k))*cos(E(k));
    
    Rot = [R11,R12,R13;R21,R22,R23;R31,R32,R33];
    
    V(:,k) = Rot'*vst;
    T(:,k) = V(:,k) + posWorld(:,k);  
    
    F(k,:) = Rot'*force(k,:)';
        
end
win = barthannwin(5);          % create a 50-point bartlett-hanning window 
win = win/sum(win);             % normalize the window

x = T(1,:); %tooltip x position
y = T(2,:); %tooltip y position
z = T(3,:); %tooltip y position
x = filtfilt(win,1,x);
y = filtfilt(win,1,y);
z = filtfilt(win,1,z);

forceZ = F(2:end,3); %Units are N
forceX = F(2:end,1);
forceY = F(2:end,2);

accel = proj321_OA(aX,aY,aZ); %combine 3 accel axes
speedX = filtfilt(win,1,diff(x)/Tspos); %calculate x and y speeds
speedY = filtfilt(win,1,diff(y)/Tspos);
speedZ = filtfilt(win,1,diff(z)/Tspos);
eDir = [speedX'*Tspos speedY'*Tspos speedZ'*Tspos];
len = sqrt(dot(eDir,eDir,2));
eDir(:,1) = eDir(:,1)./len;
eDir(:,2) = eDir(:,2)./len;
eDir(:,3) = eDir(:,3)./len;
speed = sqrt(speedX.^2 + speedY.^2 + speedZ.^2); %calculated speed

speed = filtfilt(win,1,speed);

b = regress(z',[ones(size(x')) x' y']);
x1fit = min(x):10:max(x);
x2fit = min(y):10:max(y);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;
v1 = [X1FIT(end,1)-X1FIT(1,1),X2FIT(end,1)-X2FIT(1,1),YFIT(end,1)-YFIT(1,1)];
v2 = [X1FIT(1,end)-X1FIT(1,1),X2FIT(1,end)-X2FIT(1,1),YFIT(1,end)-YFIT(1,1)];

vN = cross(v1,v2);
vN = vN/sqrt(dot(vN,vN));
eN = repmat(vN,length(forceX),1);
forceN = dot([forceX,forceY,forceZ],eN,2);
forceFriction = -dot([forceX,forceY,forceZ],eDir,2);
forceN = filtfilt(win,1,forceN);
forceFriction = filtfilt(win,1,forceFriction);

mu = median(abs(forceFriction./forceN));

GraphingTimeDelay = 0.02; % The length of time that Matlab should pause between positions when graphing, if at all, in seconds.
force_scale = 25; % Amount by which to scale forces for plot, in mm/N.

if animate
    for i = 1:length(x)-1
        if (i == 1)
            % Open figure 1.
            figure(2); clf;

            % Also plot the tip position of the haptic camera and the force vector,
            % using hold on and hold off, also keeping a handle to the plot so
            % we can update the data points later.
            hold on;
            htip = plot3(x(i),y(i),z(i),'-','color',[0 .7 0],'linewidth',2);
            hforceN = plot3(x(i) + force_scale*[0 forceN(i)*eN(i,1)], y(i) + force_scale*[0 forceN(i)*eN(i,2)], z(i) + force_scale*[0 forceN(i)*eN(i,3)], '-', 'color',[.7 0 .7],'linewidth',2);
            hforceF = plot3(x(i) + force_scale*[0 forceFriction(i)*eDir(i,1)], y(i) + [0 forceFriction(i)*eDir(i,2)], z(i) + [0 forceFriction(i)*eDir(i,3)], '-', 'color','red','linewidth',2);

            hold off;

            % Label the axes.
            xlabel('X (mm)');
            ylabel('Y (mm)');
            zlabel('Z (mm)');

            % Turn on the grid and the box.
            grid on;
            box on;

            % Set the axis limits.
            axis([floor(min(x)/50)*50 ceil(max(x)/50)*50 floor(min(y)/50)*50 ceil(max(y)/50)*50 0 250])

            % Set the view
            view(42, 30);

            % Set the axis properties for 3D visualization, which makes one
            % unit the same in every direction, and enables rotation.
            axis vis3d;

            title('Force vector over time');
        else
            % Once the animation has been set up, we don't need to reformat the
            % whole plot.  We just set the data to the correct new values for
            % the animation, the tip history, the force, and the text
            % showing the elapsed time.
            set(htip,'xdata',x(1:i),'ydata',y(1:i),'zdata',z(1:i))
            set(hforceN, 'xdata',x(i) + force_scale*[0 forceN(i)*eN(i,1)], 'ydata', y(i) + force_scale*[0 forceN(i)*eN(i,2)], 'zdata', z(i) + force_scale*[0 forceN(i)*eN(i,3)]);
            set(hforceF, 'xdata',x(i) + force_scale*[0 forceFriction(i)*eDir(i,1)], 'ydata', y(i) + [0 forceFriction(i)*eDir(i,2)], 'zdata', z(i) + [0 forceFriction(i)*eDir(i,3)]);

        end

        % Pause for a short duration so that the viewer can watch animation evolve over time.
        pause(GraphingTimeDelay)

    end
end
