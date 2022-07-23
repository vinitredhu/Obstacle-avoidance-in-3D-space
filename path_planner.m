clc
clear
close all
figure
hold on
grid on
format long
xlabel('x')
ylabel('y')
zlabel('z')
view(33,14)
xlim([-50 50])
ylim([-50 50])
zlim([-50 50])
o_d = [5 5 5;5 10 -1;35 35 35;32 25 35;-25 19 44];
q_goal = [40 40 40];
plot3(q_goal(1),q_goal(2),q_goal(3),'b*')
d_s_g = 10;
q_o_i = 25;
zeta = 0.5;
eta = 5;
d0 = [0 0 0];
plot3(d0(1),d0(2),d0(3),'b*')
alpha = 0.05;



for i = 1:length(o_d)
    x = o_d(i,1);
    y = o_d(i,2);
    z = o_d(i,3);
    plot3(x,y,z,'mo')
end

k = 0;
while d0 ~= q_goal
    plot3(d0(1),d0(2),d0(3),'g.')
    pause(0.005)
    [uatt,x_out_a,y_out_a,z_out_a] = uattractive(d0,q_goal,d_s_g,zeta);
    [u_r,x_out_r,y_out_r,z_out_r] = urepulsive(d0,o_d,q_goal,q_o_i,eta);
    d0(1) = d0(1)- alpha*(x_out_r + x_out_a);
    d0(2) = d0(2)- alpha*(y_out_r + y_out_a);
    d0(3) = d0(3)- alpha*(z_out_r + z_out_a);
    d0;
    u = uatt + u_r;
end


function [uatt,x_out_a,y_out_a,z_out_a] = uattractive(d0,q_goal,d_s_g,zeta)
    x0 = d0(1);
    y0 = d0(2);
    z0 = d0(3);
    x = q_goal(1);
    y = q_goal(2);
    z = q_goal(3);
    datt = (((x-x0)^2+(y-y0)^2+(z-z0)^2)^0.5);
    if datt <= d_s_g
        uatt = 0.5*zeta*(datt^2);
        x_out_a = zeta*(x0-x);
        y_out_a = zeta*(y0-y);
        z_out_a = zeta*(z0-z);
    else
        uatt = (d_s_g*zeta*datt)-0.5*zeta*((d_s_g)^2);
        x_out_a = zeta*(x0-x)/datt;
        y_out_a = zeta*(y0-y)/datt;
        z_out_a = zeta*(z0-z)/datt;
    end
end

function [u_r,x_out_r,y_out_r,z_out_r] = urepulsive(d0,o_d,q_goal,q_o_i,eta)
    n = size(o_d(:,1));
    u_t = 0;
    x_out_t = 0;
    y_out_t = 0;
    z_out_t = 0;
    u = zeros(n);
    x_out = zeros(n);
    y_out = zeros(n);
    z_out = zeros(n);
    x0 = d0(1);
    y0 = d0(2);
    z0 = d0(3);
    x = q_goal(1);
    y = q_goal(2);
    z = q_goal(3);
    datt = (((x-x0)^2+(y-y0)^2+(z-z0)^2)^0.5);
    [d,dx,dy,dz] = dists(d0,o_d);
    for i = 1:length(d)
        if d(i)<= q_o_i
            u(i) = (0.5*eta*((1/d(i))-(1/q_o_i))^2);
            x_out(i) = eta*((-1/d(i))+(1/q_o_i))*(1/d(i)^2)*dx(i);
            y_out(i) = eta*((-1/d(i))+(1/q_o_i))*(1/d(i)^2)*dy(i);
            z_out(i) = eta*((-1/d(i))+(1/q_o_i))*(1/d(i)^2)*dz(i);
        else
            u(i) = 0;
            x_out(i) = 0;
            y_out(i) = 0;
            z_out(i) = 0;
        end
        u_t = u_t+u(i);
        x_out_t = x_out_t + x_out(i);
        y_out_t = y_out_t + y_out(i);
        z_out_t = z_out_t + z_out(i);
    end
    u_r = u_t;
    x_out_r = x_out_t;
    y_out_r = y_out_t;
    z_out_r = z_out_t;
end
    
function [d,dx,dy,dz] = dists(d0,o_d)
    x0 = d0(1);
    y0 = d0(2);
    z0 = d0(3);
    n = size(o_d(:,1));
    d = zeros(n);
    dx = zeros(n);
    dy = zeros(n);
    dz = zeros(n);
    for i = 1:length(o_d)
        x = o_d(i,1);
        y = o_d(i,2);
        z = o_d(i,3);
        d1 = (((x0-x)^2+(y0-y)^2+(z0-z)^2)^0.5);
        dx(i) = x0-x;
        dy(i) = y0-y;
        dz(i) = z0-z;
        d(i) = d1;
    end
end



