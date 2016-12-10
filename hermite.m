%Constructing a recursive function for the hermite polynomial. This wiill

%compute the coefficients of the hermite polynomial at some value "n".

function h=hermite(n)

if n==0, 

   h=1;

elseif n==1,

    h=[2 0];

else

    A=zeros(1,n+1);

    A(1:n)=2*hermite(n-1);

    B=zeros(1,n+1);

    B(3:end)= 2*(n-1)*hermite(n-2);

    h=A-B;

end;











