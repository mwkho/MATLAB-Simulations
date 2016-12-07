clear all;

%non stationary state
%system parameters
h_bar= 0.2;
mass=0.51;
freq=0.5;
xbound =70;
num_points=1000;
x= linspace(-xbound, xbound,num_points);
dx =x(2) -x(1);

eps = sqrt((mass * freq) / h_bar);
x_sqr =x.^2;

%number of basis to use
num_basis = 5;
num_vec =  0: num_basis-1;

root_vec=  2.^(-0.5*num_vec);
fact_vec= gamma(num_vec+1);

%defining basis functions
gaussian = exp(-0.5*eps^2 * x_sqr);
H = zeros(num_basis,num_points);

for i = 1:num_basis
    hermi = hermite(i-1);
    H(i,:)=polyval(hermi, x);
end

fx = sqrt(eps)*((root_vec.* fact_vec)' * gaussian).*H;

%using trapezoidal rule to integrate
w = ones([1,num_points]);
w(1) = 0.5;
w(num_points) = 0.5;
w= dx*w;

% defining non-stationary state and normalizing
gx =0.01*x.^2+x;
A =1/sqrt((gx.*gx)*w');
gx = A*gx;

% finding coefficients for expansion of basis functions
cn =fx*(gx.*w)';
%cn'*cn

%simulation
En_vec = freq*h_bar*(num_vec +0.5);

figure, hold off
for t=0:0.05:10
    wf_re = (cn'.*cos((En_vec/h_bar) * t)) *fx;
    wf_im = (cn'.*sin((-En_vec/h_bar) * t)) *fx;
    
    %plotting real, imaginary componenets of non stationary state and
    %probabilty density
    plot(x, wf_re, 'r', x, wf_im,'b',x, abs(wf_re+1i*wf_im),'g');
    
    %labelling graph - with shrunken 
    xlim([-xbound, xbound]);
    ylim([-1,1]);
    pause(0.01)
end
%close all;
