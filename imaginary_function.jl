function int_(x)
	z=convert(Int64,floor(x));
end

function f(x)
	z=sqrt(x) * exp(-2 * x * pi) * gamma(1 + -4im * x) * (4 * x) ^ (4im * x) * exp(1im * (3 - x)) / (x - 3);
end


#plot([1,2],[1,34]);


dx=0.01;
L=4;
N=int_(L/dx+1);

a=zeros(N,3);

for i=1:N
	x=(i-1)*dx;
	a[i,1]=x;
#	a[i,2]=gamma(x);
	a[i,2]=real(f(x));
	a[i,3]=imag(f(x));
end

plot(a[:,1],a[:,2]);
plot(a[:,1],a[:,3]);
