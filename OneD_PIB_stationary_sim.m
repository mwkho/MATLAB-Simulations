clear all;
%{
%%%%%%%%%%
Simulation for One - Dimensional Particle in a Box (stationary states)
Quadrature used: Trapezoidal Rule

with 
hbar =1
L =10 
mass =1
%%%%%%%%%%
%}

%%%%%stationary state animation%%%%%
% system parameters
n=1;
h_bar= 1;
m=1;
L =10;
num_points = 1000;

%grid
x = linspace(0,L,num_points);

%basis function and enrgy value
fx = sqrt(2/L)* sin((pi/L)*n*x);
En =((n*pi*h_bar)^2 / (2*m*L));

% simulation
figure,hold off
for t=0:0.05:2
    wf_real= fx*cos((En\h_bar)*t);
    wf_im = fx*sin(-(En\h_bar)*t);
    
    %plotting real, imaginary and probability density
    p1 = plot(x,wf_real ,'r');
    hold on
    p2 = plot(x,wf_im, 'b');
    p3 = plot(x, abs(wf_real+1i*wf_im).^2, 'g');
    hold off
    
    %labeling
    legend([p1,p2,p3], 'Re(\Psi(x,t))', 'Im(\Psi(x,t)', '|\Psi(x,t)|^2')
    xlabel('x')
    axis([0,L, -1,1])
    pause(0.01)
end
