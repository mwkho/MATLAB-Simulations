clear all;
%{
%%%%%%%%%%
Animation for One - Dimensional Particle in a Box
By: Martin Ho

First part simulates nth stationary state in a box of length L.
Second part simulates a non - stationary state in a box of length L.
%%%%%%%%%%
%}

%%%%%stationary state animation%%%%%
% system parameters
n=1;
h_bar= 1;
m=1;
L =10;
num_points = 1000;

x = linspace(0,L,num_points);
dx = x(2)-x(1);

%basis function 
fx = sqrt(2/L)* sin((pi/L)*n*x);
En =((n*pi*h_bar)^2 / (2*m*L));

% animating real and complex components
%{
figure,hold off
for t=0:0.05:1000
    wf_real= fx*cos((En\h_bar)*t);
    wf_im = fx*sin(-(En\h_bar)*t);
    plot(x,wf_real ,'r', x,wf_im, 'b')
    axis([0,L, -1,1])
    pause(0.01)
end
%}

%%%%%non-stationary state animation%%%%%
%using trapezoidal rule to find coefficients
w = dx *ones([1,num_points]);
w(1) = 0.5; 
w(num_points) = 0.5;

% defining a non stationary function and normalizing it
gx = cos(x);
A = 1/sqrt((gx.*w)*gx');
gx = A.*gx;

% number of basis used
n_basis = 2;
num_vec = 1:n_basis;

% matrix of functions
psix = sqrt(2/L)* sin((pi/L)*num_vec'*x);

%finding coefficient of each wavefuncton
cn = psix *(gx.*w)';

%animating real and complex components
En_vec = (num_vec).^2*(pi*h_bar)^2 / (2*m*L);

% attach time dependence and do a final summation
figure,hold off
for t=0:0.05:100
    ht_real = cos((En_vec/h_bar) * t);
    ht_im = sin(-(En_vec/h_bar) * t);
    wf_real = (cn' .* ht_real) * psix;
    wf_im = (cn' .* ht_im) * psix;
    % plotting real, imaginary components, and the prob density
    plot(x,wf_real ,'r', x,wf_im, 'b', x, abs(wf_real+1i*wf_im).^2, 'g')
    axis([0,L, -1,1])
    pause(0.01)
end
