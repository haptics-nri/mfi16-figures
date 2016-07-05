function [mu,R2] = FitFriction_animate
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
close all
clear all

plotF = true;
animate = true;

FricFilename = 'MDF_Subject9.mat';
load(FricFilename);

Ts = 1/10000;
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
timePos = trackedtime-trackedtime(1)+1/240; %time vector for position
Tspos = max(timePos)/length(timePos);

vst = [119.65;-3.30;2.03];
vcl = vst*0.9;

clippedTime = time(1:100000);
endPos = find(timePos>clippedTime(end),1,'first')-1;
if isempty(endPos)
    endPos = length(timePos);
end
timePos = timePos(1:endPos);

A = unwrap(trackedorient(1,1:endPos)*pi/180);
E = unwrap(trackedorient(2,1:endPos)*pi/180);
R = unwrap(trackedorient(3,1:endPos)*pi/180);
posWorld = trackedpos(:,1:endPos);

FT=data(1:100000,4:9)-data(1:100000,10:15);
FT=FT.';
FT=FT-diag(FTbias)*ones(size(FT));
workmatrix= [-0.0148464530000000,-0.0528270810000000,-0.128986314000000,-3.30285491500000,0.260661584000000,3.21999616900000;0.0362472370000000,3.95397478000000,-0.0378257670000000,-1.96783033200000,-0.160302624000000,-1.83848942100000;3.68742806500000,0.0353190710000000,3.89915414400000,0.164660002000000,3.57345682600000,-0.0243610810000000;0.706659571000000,24.2776279500000,20.7647702100000,-11.1733273100000,-21.0068402500000,-11.1552407700000;-23.8273414000000,0.0191202520000000,13.7248583200000,21.0056583200000,9.70693341200000,-19.9421094100000;0.215633100000000,13.5789013200000,0.736703546000000,14.1776746400000,1.02066334000000,14.1765185000000;];
v=workmatrix*FT;
v=v.';

fX = interp1(clippedTime,v(:,3),timePos,'cubic')';
fY = interp1(clippedTime,v(:,2),timePos,'cubic')';
fZ = interp1(clippedTime,-v(:,1),timePos,'cubic')';

force = [fX fY fZ];

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
    Vcl(:,k) = Rot'*vcl;
    T(:,k) = V(:,k) + posWorld(:,k);  
    Tcl(:,k) = Vcl(:,k) + posWorld(:,k);
    
    F(k,:) = Rot'*force(k,:)';
        
end
win = barthannwin(5);          % create a 50-point bartlett-hanning window 
win = win/sum(win);             % normalize the window

xend = posWorld(1,:);
yend = posWorld(2,:);
zend = posWorld(3,:);

xcl = Tcl(1,:);
ycl = Tcl(2,:);
zcl = Tcl(3,:);

x = T(1,:); %tooltip x position
y = T(2,:); %tooltip y position
z = T(3,:); %tooltip y position
z = ones(size(z))*mean(z);

forceZ = F(2:end,3); %Units are N
forceX = F(2:end,1);
forceY = F(2:end,2);

accel = proj321_OA(aX,aY,aZ); %combine 3 accel axes
speedX = diff(x)/Tspos; %calculate x and y speeds
speedY = diff(y)/Tspos;
speedZ = diff(z)/Tspos;


eDir = [speedX'*Tspos speedY'*Tspos speedZ'*Tspos];
len = sqrt(dot(eDir,eDir,2));
eDir(:,1) = eDir(:,1)./len;
eDir(:,2) = eDir(:,2)./len;
eDir(:,3) = eDir(:,3)./len;
nan_loc = find(isnan(eDir));
eDir(nan_loc) = eDir(nan_loc-1);

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

forceFriction = abs(forceFriction);
forceN = filtfilt(win,1,forceN);
forceFriction = filtfilt(win,1,forceFriction);

p= robustfit(forceN,forceFriction);
Y = p(1) + p(2)*forceN;

temp = corrcoef(Y,forceFriction);
R2 = temp(1,2)^2;

% mu = p(2);
mu = forceN\forceFriction;
Y = mu*forceN;

if plotF
    figure;
    plot(forceN,forceFriction,'.','Color',[148 138 84]/255)
    hold on;
%     plot(forceN,Y,'Color',[255 51 0]/255,'LineWidth',1.5);
    plot([0; forceN],[0; Y],'Color',0.2*[1 1 1],'LineWidth',1.5);
    xlabel('Normal Force (N)');
    ylabel('Friction Force (N)');
%     axis([0 0.5*ceil(max(forceN)/0.5) 0 0.5*ceil(max(forceFriction)/0.5)]) 

    R = corrcoef(forceFriction,Y);
    R2 = R(1,2)^2
    
    axis equal
%     axis([0 1.5 0 1])
%     set(gca,'XTick',[0 0.5 1 1.5]);
    axis([0 4.25 0 1])
    set(gca,'XTick',[0 1 2 3 4]);
    xstart = 0.637;
%     harr = annotation('textarrow',[xstart xstart-0.10], [0.5 0.5],'String',...
%         strcat(' $\mu_k = ',{' '},num2str(mu,2),'$'),'FontSize',12,'FontWeight','bold','LineWidth',...
%         1.5,'HeadLength',15,'Interpreter','LaTEX');
%     xbeg = 0.7;
%     xstop = xbeg+0.15;
%     yoff = 0.025;
    xbeg = 2.9;
    xstop = xbeg+0.25;
    yoff = 0.025;
% %     text(0.4,0.75+yoff/2,strcat(' $\mu_k = ',{' '},num2str(mu,'%6.2f'),'$'),'FontSize',12,'FontWeight','bold','Interpreter','LaTEX');
%     text(0.35,0.55,strcat(' $\mu_k = ',{' '},num2str(mu,'%6.2f'),'$'),'FontSize',12,'FontWeight','bold','Interpreter','LaTEX');
    text(2.25,0.75,strcat(' $\mu_k = ',{' '},num2str(mu,'%6.2f'),'$'),'FontSize',12,'FontWeight','bold','Interpreter','LaTEX');
    line([xbeg xstop],[mu*xstop+yoff mu*xstop+yoff],'Color','k','LineWidth',1.25)
    line([xbeg xbeg],[mu*xbeg+yoff mu*xstop+yoff],'Color','k','LineWidth',1.25);
    
    set(gcf,'PaperPosition',[1 1 5 2.5]);
    print('-depsc','frictionMDF.eps');
end

GraphingTimeDelay = 0.00; % The length of time that Matlab should pause between positions when graphing, if at all, in seconds.
force_scale = 100; % Amount by which to scale forces for plot, in mm/N.
if animate
    
    writerObj = VideoWriter('FrictionForces_slow.avi');
    writerObj.FrameRate = 24; %48
    opengl('software')
    open(writerObj);
    for i = 1:length(x)-1
        if (i == 1)
            % Open figure 1.
            hfig = figure;
            Coord=get(hfig,'Position');
            set(hfig,'Position',[Coord(1),Coord(2),1.5*Coord(3),Coord(4)]);
            set(hfig,'color','w');
            % Also plot the tip position of the robot and the force vector,
            % using hold on and hold off, also keeping a handle to the plot so
            % we can update the data points later.
            subplot(2,2,[1 3]);
            set(gca,'Zdir','reverse')
            hold on;
            x1fit = 100:10:200;
            x2fit = 200:10:350;
            [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
            YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;
            surface(X1FIT,X2FIT,YFIT,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.5 0.5 0.5]);
            alpha(0.5)
            
            [xcyl, ycyl, zcyl] = cylinder(5);
            [xcone ycone zcone] = cylinder([0 5]);
            hc = surf([xcyl(1,:)+xcl(i);xcyl(2,:)+xend(i)],[ycyl(1,:)+ycl(i);ycyl(2,:)+yend(i)],[zcyl(1,:)+zcl(i);zcyl(2,:)+zend(i)]);
            % set(hc,'FaceColor',[148 54 52]/255);
            set(hc,'FaceColor','none');
            set(hc,'LineWidth',1);
            set(hc,'EdgeAlpha',0.2);
            hcone = surf([xcone(1,:)+x(i);xcone(2,:)+xcl(i)],[ycone(1,:)+y(i);ycone(2,:)+ycl(i)],[zcone(1,:)+z(i);zcone(2,:)+zcl(i)]);
%             set(hc,'FaceColor',[148 54 42]/255);
            set(hcone,'FaceColor','none');
            set(hcone,'LineWidth',1);
            set(hcone,'EdgeAlpha',0.2);
            
%             htip = plot3(x(i),y(i),z(i),'-','color',[0.7 0.7 0.7],'linewidth',2);
%             hforce = plot3(x(i) + force_scale*[0 test(i)*eDir(i,1)], y(i) + force_scale*[0 test(i)*eDir(i,2)], z(i) + force_scale*[0 test(i)*eDir(i,3)], '-', 'color',[.7 0 .7],'linewidth',2);
            hforce = plot3(x(i) - force_scale*[0 forceFriction(i)*eDir(i,1)], y(i) - force_scale*[0 forceFriction(i)*eDir(i,2)], z(i) - force_scale*[0 forceFriction(i)*eDir(i,3)], '-', 'color',[153 102 51]/255,'linewidth',3);
            hfn = plot3(x(i) + force_scale*[0 forceN(i)*eN(i,1)], y(i) + force_scale*[0 forceN(i)*eN(i,2)], z(i) + force_scale*[0 forceN(i)*eN(i,3)], '-', 'color',[137 117 43]/255,'linewidth',3);
            hold off;

            % Label the axes.
            xlabel('X (mm)');
            ylabel('Y (mm)');
            zlabel('Z (mm)');

            % Turn on the grid and the box.
            grid on;
            box on;

            % Set the axis limits.
            axis([floor(min(x)/50)*50 ceil(max(x)/50)*50 floor(min(y)/50)*50 ceil(max(y)/50)*50 -150 100])

            % Set the view to see the robot easily.
            view(42, 30);

            % Set the axis properties for 3D visualization, which makes one
            % unit the same in every direction, and enables rotation.
            axis vis3d;

            subplot(2,2,2);
            hn = plot(timePos(i),forceN(i),'Color',[137 117 43]/255,'LineWidth',1.5);
            axis([0 10 0 1.5]);
            set(gca, 'XTick', []);
            ylabel('Normal Force (N)')
            subplot(2,2,4);
            hf = plot(timePos(i),forceFriction(i),'Color',[153 102 51]/255,'LineWidth',1.5);
            axis([0 10 0 1]);
            xlabel('Time (s)')
            ylabel('Friction Force (N)');
            
            frame = getframe(hfig);
            writeVideo(writerObj,frame);
        else%if mod(i,5)== 1
            % Once the animation has been set up, we don't need to reformat the
            % whole plot.  We just set the data to the correct new values for
            % the robot animation, the tip history, the force, and the text
            % showing the elapsed time.
%             set(htip,'xdata',x(1:i),'ydata',y(1:i),'zdata',z(1:i))
%             set(hforce, 'xdata',x(i) + force_scale*[0 test(i)*eDir(i,1)], 'ydata', y(i) + force_scale*[0 test(i)*eDir(i,2)], 'zdata', z(i) + force_scale*[0 test(i)*eDir(i,3)]);
            set(hc,'xdata',[xcyl(1,:)+xcl(i);xcyl(2,:)+xend(i)],'ydata',[ycyl(1,:)+ycl(i);ycyl(2,:)+yend(i)],'zdata',[zcyl(1,:)+zcl(i);zcyl(2,:)+zend(i)]);
            set(hcone,'xdata',[xcone(1,:)+x(i);xcone(2,:)+xcl(i)],'ydata',[ycone(1,:)+y(i);ycone(2,:)+ycl(i)],'zdata',[zcone(1,:)+z(i);zcone(2,:)+zcl(i)]);
            set(hforce, 'xdata',x(i) - force_scale*[0 forceFriction(i)*eDir(i,1)], 'ydata', y(i) - force_scale*[0 forceFriction(i)*eDir(i,2)], 'zdata', z(i) - force_scale*[0 forceFriction(i)*eDir(i,3)]);
            set(hfn, 'xdata',x(i) + force_scale*[0 forceN(i)*eN(i,1)], 'ydata', y(i) + force_scale*[0 forceN(i)*eN(i,2)], 'zdata', z(i) + force_scale*[0 forceN(i)*eN(i,3)]);
            
            set(hn,'xdata',timePos(1:i),'ydata',forceN(1:i));
            set(hf,'xdata',timePos(1:i),'ydata',forceFriction(1:i));    
            frame = getframe(hfig);
            writeVideo(writerObj,frame);
        end

        % Pause for a short duration so that the viewer can watch animation evolve over time.
        pause(GraphingTimeDelay)

    end
    close(writerObj);
end

end

