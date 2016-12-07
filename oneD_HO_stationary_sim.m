clear all;

%stationary state animation
%system parameters
n=0;
h_bar= 1;
mass=1;
freq=0.5;
xbound =20;
num_points=1000;
x= linspace(-xbound, xbound,num_points);

%making the x^2 
x_sqr = x.*x;

%making gaussian 
gaussian = exp(-mass*freq*x_sqr/(2*h_bar));

% coefficients
n_coeff = (factorial(n)*2^n)^-0.5*(freq*mass/(pi*h_bar))^0.25;

%making hermite polynomials of degree n
n_hermipol = hermite(n);
n_hermite = polyval(n_hermipol,sqrt(mass*freq/h_bar)*x);

n_wf= n_coeff*gaussian.*n_hermite;

%energy
En = freq*(n+0.5);

%animation
figure,hold off
for t=0:0.05:10
    % the real and the complex functions
    n_real = n_wf*cos(-t*En);
    n_im = n_wf*sin(-t*En);
    plot(x, n_real,'r',x, n_im,'b',x,abs(n_real+1i*n_im).^2, 'g')
    axis([-xbound,xbound, -1,1])
    pause(0.01)
end
