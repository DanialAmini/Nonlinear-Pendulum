#cd("D:/julia_files/pendulum")
#include("pendul_sav1.jl")

function rond(x)
    z=convert(Int64,floor(x));
end

#using Plots
using PyPlot
#pyplot()

T=1000;
dt=0.028;
N=rond(T/dt+1);
t=zeros(N,1);
theta=zeros(N,2);
theta_d=zeros(N,2);

theta[1,1]=pi/1.02;
theta[2,1]=theta[1,1];

for i=1:N
    t[i,1]=(i-1)*dt;
end

for i=1:N-1
    thetaS=theta[i,1]+dt*theta_d[i,1];
    theta_dS=theta_d[i,1]-dt*sin(theta[i,1]);

    theta[i+1,1]=theta[i,1]+dt*(theta_d[i,1]+theta_dS)/2;
    theta_d[i+1,1]=theta_d[i,1]-dt*(sin(theta[i,1])+sin(thetaS))/2;
end

res_=zeros(N,4);
for i=1:N
    res_[i,1]=t[i,1];
    res_[i,2]=theta[i,1];
    res_[i,3]=theta_d[i,1];
    res_[i,4]=1-cos(theta[i,1])+1/2*theta_d[i,1]^2;
end

#plot(t[:,1],theta[:,1]);
fff=open(string("res_",round(theta[1,1],digits=3),".txt"),"w");
for i=1:N
    for j=1:4
        write(fff,string(res_[i,j]));
        if j<4
            write(fff,"\t");
        end
    end
    write(fff,"\r\n");
end
close(fff);

PyPlot.subplot(3,1,3);
plot(res_[:,1],res_[:,2]);
PyPlot.xlabel("t");
PyPlot.ylabel("theta");
PyPlot.xlim((0,50));
PyPlot.title("first oscillations");

#PyPlot.figure();
PyPlot.subplot(3,1,1);
plot(res_[:,1],res_[:,2]);
PyPlot.xlabel("t");
PyPlot.ylabel("theta");
PyPlot.title("long time");

#PyPlot.figure();
PyPlot.subplot(3,1,2);
plot(res_[:,1],res_[:,4]);
PyPlot.xlabel("t");
PyPlot.ylabel("energy");
PyPlot.title("long time");

PyPlot.tight_layout();
#PyPlot.close();

