clc
clear
close all

load('way_point.mat')

animate = 1;            % animation on/off
plot_er = 0;            % plot error path following
plot_state = 1;         % plot state x-y
plot_attack = 1;        % plot sensor attack vs time

% feedback control options
feedback_mode = 3;      % 1 = ideal, 2 = meas, 3 = est AO, 4 = est AEKF

% param truck trailer path following
vt = 2;          % max speed
k1 = 0.5;        % koef path following
k2 = 1;          % koef path following
dt = 0.05;       % sampling time

% make waypoint
time_s = [0 10 20 30 40 50 60 70 80 90 100];
time_array = 0:1:max(time_s);      

way_px = spline(time_s,way_p_1(:,1),time_array);  
way_py = spline(time_s,way_p_1(:,2),time_array);

% inisial state simulation
xs(1) = -10;
ys(1) = 4;
vs(1) = 0;
yaws(1) = 1.484;

v_act = 0;
r_act = 0;

% matrix for observer
A = zeros(4,4);
C = [1 0 0 0; 0 1 0 0];
D = eye(2);
B = [0 0; 0 0; 1 0; 0 1];

As = 10*eye(2);

Ab = [A zeros(4,2); As*C -As];
Bb = [B; zeros(2,2)];
Cb = [zeros(2,4) eye(2)];
Dd = [zeros(4,2); As*D];

Ab = eye(6)+Ab*dt;

x           = [xs(1) ys(1) 0 yaws(1) 0 0]';
xbar        = [xs(1) ys(1) 0 yaws(1) 0 0]';
xhat        = [xs(1) ys(1) 0 yaws(1) 0 0]';
xreal       = [xs(1) ys(1) 0 yaws(1)]';
xArray      = [];
xrealArray  = [];
xbarArray   = [];
xhatArray   = [];
yreal       = [xs(1) ys(1)]';
yrealArray  = [];

theta           = [0;0];
thetabar        = [0;0];
thetahat        = [0;0];
thetaArray      = [];
thetabarArray   = [];
thetahatArray   = [];

% parameter observer
lambdax = 0.999;
lambdat = 0.8;

Rx = 10000*diag([100 100]);
Rt = 0.1*diag([10 10]);
Px = 0.001*eye(6);
Pt = 0.001*eye(2);
Gamma = [zeros(2,6)]';

% parameter AEKF
lambda = 0.93;
P  = 0.1*eye(6);
Q  = 0.1*diag([10 10 1 1 1 1]);
R  = 10000*diag([1 1]);
S  = 0.1*eye(2);
Kappa = [zeros(2,6)]';

i = 1;

while true
    
    if i*dt < 59
        theta = zeros(2,1);
    else
        if mod(i,1) == 0
            theta(1) = theta(1)-0.02;
        end
        if mod(i,1) == 0
            theta(2) = theta(2)-0.01;
        end
    end
    
    xrealArray  = [xrealArray xreal];
    xbarArray  = [xbarArray xbar];
    xhatArray  = [xhatArray xhat];
    xArray  = [xArray x];
    
    yrealArray  = [yrealArray yreal];
    
    thetaArray    = [thetaArray theta];
    thetabarArray = [thetabarArray thetabar];
    thetahatArray = [thetahatArray thetahat];
    
    if feedback_mode == 1
        statex = xs(i);
        statey = ys(i);
    elseif feedback_mode == 2
        statex = yreal(1);
        statey = yreal(2);
    elseif feedback_mode == 3
        statex = xbar(1);
        statey = xbar(2);
    elseif feedback_mode == 4
        statex = xhat(1);
        statey = xhat(2);
    end
    
    % calculate error position and orientation
    [valx,idx] = min(abs(sqrt((way_px-statex).^2+(way_py-statey).^2)));
    
    if idx == 1
        idx = 2;
    end
    
    way_yaw = atan2(way_py(idx)-way_py(idx-1),way_px(idx)-way_px(idx-1));
    
    er(i)   = way_yaw-yaws(i);
    
    dxw     = way_px(idx) - way_px(idx-1);
    dyw     = way_py(idx) - way_py(idx-1);
    curv    = dxw * way_py(idx-1) -  dyw * way_px(idx-1);
    ep(i)   = ((statex*dyw)+ curv - (statey*dxw))/sqrt(dxw^2+dyw^2) + 10^(-32);
    
    % calculate control signal
    v_cd    = vt/(2+sqrt(abs(ep(i))+abs(er(i))));
    
    if abs(er(i)) < 0.0001 
        r_cd = k1*er(i) + k2*ep(i)*v_cd;
    else
        r_cd = k1*er(i) + k2*ep(i)*v_cd*(sin(er(i))/er(i));
    end
    
    % pid actuator
    v_act = cs_long_control(v_act, v_cd);
    r_act = cs_lat_control(r_act, r_cd);
    
    % sensor speed and steer
    a(i) = (v_act - vs(i))/dt;
    
    v(i) = v_act;
    r(i) = r_act;
    
    u = [a(i); r(i)];
    
    % simulation
    xs(i+1)     = xs(i)+dt*vs(i)*cos(yaws(i));
    ys(i+1)     = ys(i)+dt*vs(i)*sin(yaws(i));
    vs(i+1)     = vs(i)+dt*a(i);
    yaws(i+1)   = yaws(i)+dt*r(i);
    
    xreal = [xs(i+1) ys(i+1) vs(i+1) yaws(i+1)]';
    yreal = C*xreal+D*theta;
    
    % simulation with sensor attack
    fx = [x(3)*cos(x(4)); x(3)*sin(x(4)); a(i); r(i);...
         0; 0];
    
    x = Ab*x + fx*dt + Dd*theta*dt;
    y = Cb*x;
      
    % Estimation using observer
    Kx      = Px*Cb'*inv(Cb*Px*Cb'+Rx);
    Kt      = Pt*Gamma'*Cb'*inv(Cb*Gamma*Pt*Gamma'*Cb'+Rt);
    Gamma   = (eye(6)-Kx*Cb)*Gamma;
    xbar    = xbar+(Kx+Gamma*Kt)*(y-Cb*xbar);
    thetabar = thetabar-Kt*(y-Cb*xbar);
    
    fxbar = [xbar(3)*cos(xbar(4)); xbar(3)*sin(xbar(4)); a(i); r(i);...
         0; 0];
    
    xbar = Ab*xbar+fxbar*dt+Dd*thetabar*dt;
    thetabar = thetabar;

    Fxbar = [0 0 cos(xbar(4)) -xbar(3)*sin(xbar(4));
        0 0 sin(xbar(4)) xbar(3)*cos(xbar(4));
        0 0 0 0;
        0 0 0 0];
    
    Fxb = [Fxbar zeros(4,2); zeros(2,6)];
    
    Px      = (1/lambdax)*(Ab+Fxb*dt)*(eye(6)-Kx*Cb)*Px*(Ab+Fxb*dt)';
    Pt      = (1/lambdat)*(eye(2)-Kt*Cb*Gamma)*Pt;
    Gamma   = (Ab+Fxb*dt)*Gamma-Dd*dt;
    
    % Estimation using adaptive EKF
    
    Fxhat = [0 0 cos(xhat(4)) -xhat(3)*sin(xhat(4));
        0 0 sin(xhat(4)) xhat(3)*cos(xhat(4));
        0 0 0 0;
        0 0 0 0];

    Fxh = [Fxhat zeros(4,2); zeros(2,6)];
    
    fxhat = [xhat(3)*cos(xhat(4)); xhat(3)*sin(xhat(4)); ...
        a(i); r(i); 0; 0];

    P = (Ab+Fxh*dt)*P*(Ab+Fxh*dt)' + Q;
    Sigma = Cb*P*Cb' + R;
    Kf = P*Cb'/(Sigma);
    P = (eye(6)-Kf*Cb)*P;
    
    Omega = Cb*(Ab+Fxh*dt)*Kappa+Cb*Dd*dt;
    Kappa = (eye(6)-Kf*Cb)*(Ab+Fxh*dt)*Kappa+(eye(6)-Kf*Cb)*Dd*dt;
    Lambda = inv(lambda*Sigma+Omega*S*Omega');
    Delta = S*Omega'*Lambda;
    S = (1\lambda)*S-(1\lambda)*S*Omega'*Lambda*Omega*S;

    ytilde = y-Cb*((Ab+Fxh*dt)*xhat+Bb*dt*u+Dd*dt*thetahat);

    thetahat_prev = thetahat;
    thetahat = thetahat + Delta*ytilde;
    xhat = Ab*xhat+fxhat*dt+Dd*dt*thetahat_prev...
        +Kf*dt*ytilde+Kappa*dt*(thetahat-thetahat_prev);
   
    
    if idx == length(way_px)
        stop = 1;
    else
        stop = 0;
    end
    
    i = i+1;
    
    if stop == 1 || i > 5000
        a(i) = a(i-1);
        v(i) = v(i-1);
        r(i) = r(i-1);
        ep(i) = ep(i-1);
        er(i) = er(i-1);
        xrealArray  = [xrealArray xreal];
        xbarArray  = [xbarArray xbar];
        xhatArray  = [xhatArray xhat];
        xArray  = [xArray x];
        yrealArray  = [yrealArray yreal];
        thetaArray    = [thetaArray theta];
        thetabarArray = [thetabarArray thetabar];
        thetahatArray = [thetahatArray thetahat];
        break
    end
    
end

if animate == 1 
    figure(2); 

    fh = figure(2);
    fh.WindowState = 'maximized';

    img = imread('background.png');
    image('CData',img,'XData',[-20 80],'YData',[40 -10])
    
    hold on
    
    plot(way_px,way_py,'g','LineWidth',3)
    
    axis equal
    axis([-20 80 -10 40])
    grid on
    xlabel('x(m)');ylabel('y(m)')


    for j = 1:size(xbarArray,2)-1 

        X   = xbarArray(1,j);
        Y   = xbarArray(2,j);
        Yaw = xbarArray(4,j);
        Xa  = yrealArray(1,j);
        Ya  = yrealArray(2,j);
        Yawa = atan2(yrealArray(2,j+1)-yrealArray(2,j),yrealArray(1,j+1)-yrealArray(1,j));
        
        if j < 3
            Yawa = Yaw;
        end        
        
        % plot x-y
        plot(yrealArray(1,1:j),yrealArray(2,1:j),'r','LineWidth',3);
        plot(xbarArray(1,1:j),xbarArray(2,1:j),':y','LineWidth',3);
        
        % draw the ship
        [X0a,Y0a,X1a,Y1a,X2a,Y2a,X3a,Y3a,X4a,Y4a] = draw_2D_ship(Xa,Ya,Yawa);
        [X0,Y0,X1,Y1,X2,Y2,X3,Y3,X4,Y4] = draw_2D_ship(X,Y,Yaw);

        if j==1
            
            D1a = patch('Faces',1:32,'Vertices',[X0a; Y0a]','Facecolor',[0.6 0 0.3],'FaceAlpha',1);
            D2a = patch('Faces',1:21,'Vertices',[X1a; Y1a]','Facecolor',[0.8 0 0.4],'FaceAlpha',1);
            D3a = patch('Faces',1:21,'Vertices',[X2a; Y2a]','Facecolor',[1 0 0.5],'FaceAlpha',1);
            D4a = patch('Faces',1:12,'Vertices',[X3a; Y3a]','Facecolor',[1 0 0.5],'FaceAlpha',1,'Edgecolor',[0.8 0 0.0]);
            D5a = patch('Faces',1:13,'Vertices',[X4a; Y4a]','Facecolor',[0.8 0 0.4],'FaceAlpha',1);

            D1 = patch('Faces',1:32,'Vertices',[X0; Y0]','Facecolor',[0 0.3 0.6]);
            D2 = patch('Faces',1:21,'Vertices',[X1; Y1]','Facecolor',[0 0.4 0.8]);
            D3 = patch('Faces',1:21,'Vertices',[X2; Y2]','Facecolor',[0 0.5 1]);
            D4 = patch('Faces',1:12,'Vertices',[X3; Y3]','Facecolor',[0 0.5 1],'Edgecolor',[0 0 0.8]);
            D5 = patch('Faces',1:13,'Vertices',[X4; Y4]','Facecolor',[0 0.4 0.8]);

        elseif j==2
            
            D1a = patch('Faces',1:32,'Vertices',[X0a; Y0a]','Facecolor',[0.6 0 0.3],'FaceAlpha',1);
            D2a = patch('Faces',1:21,'Vertices',[X1a; Y1a]','Facecolor',[0.8 0 0.4],'FaceAlpha',1);
            D3a = patch('Faces',1:21,'Vertices',[X2a; Y2a]','Facecolor',[1 0 0.5],'FaceAlpha',1);
            D4a = patch('Faces',1:12,'Vertices',[X3a; Y3a]','Facecolor',[1 0 0.5],'FaceAlpha',1,'Edgecolor',[0.8 0 0.0]);
            D5a = patch('Faces',1:13,'Vertices',[X4a; Y4a]','Facecolor',[0.8 0 0.4],'FaceAlpha',1);

            D1 = patch('Faces',1:32,'Vertices',[X0; Y0]','Facecolor',[0 0.3 0.6]);
            D2 = patch('Faces',1:21,'Vertices',[X1; Y1]','Facecolor',[0 0.4 0.8]);
            D3 = patch('Faces',1:21,'Vertices',[X2; Y2]','Facecolor',[0 0.5 1]);
            D4 = patch('Faces',1:12,'Vertices',[X3; Y3]','Facecolor',[0 0.5 1],'Edgecolor',[0 0 0.8]);
            D5 = patch('Faces',1:13,'Vertices',[X4; Y4]','Facecolor',[0 0.4 0.8]);
            
            pause(1);

        else
            
            set(D1a,'Faces',1:32,'Vertices',[X0a; Y0a]','Facecolor',[0.6 0 0.3],'FaceAlpha',1);
            set(D2a,'Faces',1:21,'Vertices',[X1a; Y1a]','Facecolor',[0.8 0 0.4],'FaceAlpha',1);
            set(D3a,'Faces',1:21,'Vertices',[X2a; Y2a]','Facecolor',[1 0 0.5],'FaceAlpha',1);
            set(D4a,'Faces',1:12,'Vertices',[X3a; Y3a]','Facecolor',[1 0 0.5],'FaceAlpha',1,'Edgecolor',[0.8 0 0.0]);
            set(D5a,'Faces',1:13,'Vertices',[X4a; Y4a]','Facecolor',[0.8 0 0.4],'FaceAlpha',1);

            set(D1,'Faces',1:32,'Vertices',[X0; Y0]','Facecolor',[0 0.3 0.6]);
            set(D2,'Faces',1:21,'Vertices',[X1; Y1]','Facecolor',[0 0.4 0.8]);
            set(D3,'Faces',1:21,'Vertices',[X2; Y2]','Facecolor',[0 0.5 1]);
            set(D4,'Faces',1:12,'Vertices',[X3; Y3]','Facecolor',[0 0.5 1],'Edgecolor',[0 0.0 0.8]);
            set(D5,'Faces',1:13,'Vertices',[X4; Y4]','Facecolor',[0 0.4 0.8]);

        end

        pause(0.01)   

    end
    
end


if plot_state == 1
    figure(3)
    fh = figure(3);
    fh.WindowState = 'maximized';
    img = imread('background.png');
    image('CData',img,'XData',[-20 80],'YData',[40 -10])
    hold on
    plot(way_px,way_py,'LineWidth',2,'Color','g')                                           % plot way point
%     plot(xrealArray(1,:),xrealArray(2,:),'r','LineWidth',4)                                 % plot real state
    plot(yrealArray(1,:),yrealArray(2,:),'r','LineWidth',4)                                 % plot real measurement
    plot(xbarArray(1,:),xbarArray(2,:),':','LineWidth',4,'Color','y')                       % plot estimated value using AO
%     plot(xhatArray(1,:),xhatArray(2,:),':','LineWidth',4,'Color','#0000FF')                 % plot estimated value using AEKF
    plot(way_px(1),way_py(1),'.b','MarkerSize',40)                                          % plot initial point
    grid on
    grid minor
    axis equal
    axis([-20 80 -10 40])
    legend('Desired path','Path under cyber-attacks','Estimate','Initial point')
    xlabel('x (m)')
    ylabel('y (m)')
    set(gca,'color','white','FontSize',20)
end

if plot_er == 1
    figure(1)
    subplot(2,1,1)
    plot(r)
    hold on
    plot(v)
    legend('delta','v')
    subplot(2,1,2)
    plot(er)
    hold on
    plot(ep)
    legend('er','ep')
end

if plot_attack == 1
    tt = [1:length(thetaArray(1,:))]*0.05;
    figure(4)
    fh = figure(4);
    fh.Position = [0 50 950 550];
    subplot(2,1,1)
    plot(tt,thetaArray(1,:),'k','LineWidth',4)
    hold on
    plot(tt,thetabarArray(1,:),':b','LineWidth',4)
    grid on
    grid minor
    xlim([0 90])
    set(gca,'color','white','FontSize',20)
    legend('\theta_{Actual}','\theta_{Estimate}','location','northwest')
    ylabel({'$\theta_1$'},'Interpreter','latex')
    subplot(2,1,2)
    plot(tt,thetaArray(2,:),'k','LineWidth',4)
    hold on
    plot(tt,thetabarArray(2,:),':b','LineWidth',4)
    grid on
    grid minor
    xlim([0 90])
    set(gca,'color','white','FontSize',20)
    ylabel({'$\theta_2$'},'Interpreter','latex')
    xlabel('Time (s)')
end

function cs_long = cs_long_control(cs_long_act, cs_long_cd)
    er_cs_long = cs_long_cd-cs_long_act;
    cs_long = cs_long_act + 0.02 * er_cs_long;
end

function cs_lat = cs_lat_control(cs_lat_act, cs_lat_cd)
    er_cs_lat = cs_lat_cd-cs_lat_act;
    cs_lat = cs_lat_act + 0.2 * er_cs_lat;
end