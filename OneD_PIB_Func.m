function OneD_PIB_Func = OneD_PIB_Func(h_bar, mass,L, num_points,x, n_basis,gx, time)
close all;
%{
%%%%%%%%%%
Simulation for One - Dimensional Particle in a Box (non - stationary states)
Quadrature used: Trapezoidal Rule

with 
hbar =1
L = 10
mass = 1

Legend:
red - Real WF
blue - Imaginary WF
green - probability density of WF
%%%%%%%%%%
%}

%system parameters
%h_bar= 1;
%mass=1;
%L =10;
%num_points = 1000;
%x = linspace(0,L,num_points);
dx = x(2)-x(1);

% number of basis used
%n_basis = 5;
num_vec = 1:n_basis;

%using trapezoidal rule to find coefficients
w = ones([1,num_points]);
w(1) = 0.5; 
w(num_points) = 0.5;
w = dx* w;

% defining a non stationary function and normalizing it
%gx = sqrt(2/L)*(sin((pi/L)*2*x)+ sin((pi/L)*1*x)+sin((pi/L)*3*x)+sin((pi/L)*4*x));
A = 1/sqrt((gx.*w)*gx');
gx = A.*gx;

% matrix of functions
psix = sqrt(2/L)* sin((pi/L)*num_vec'*x);

%finding coefficient of each wavefuncton
cn = psix *(gx.*w)';

%animating real and complex components
En_vec = (num_vec).^2*(pi*h_bar)^2 / (2*mass*L);

% attach time dependence and do a final summation
figure
for t=0:0.05:time
    ht_real = cos((En_vec/h_bar) * t);
    ht_im = sin(-(En_vec/h_bar) * t);
    
    wf_real = (cn' .* ht_real)* psix;
    wf_im = (cn' .* ht_im)* psix;

    % plotting real, imaginary components, and the prob density
    plot(x,wf_real ,'r',x,wf_im, 'b',x, abs(wf_real+1i*wf_im).^2, 'g');
    
    %labelling graph
    %legend([p1,p2,p3], 'Re(\Psi(x,t))', 'Im(\Psi(x,t)', '|\Psi(x,t)|^2')
    xlabel('x')
    axis([0,L, -1,1])
    pause(0.01)
end
end