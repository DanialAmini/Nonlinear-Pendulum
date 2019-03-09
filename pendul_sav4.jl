#cd("D:/julia_files/pendulum")
#include("pendul_sav2.jl")

#run it with julia client, juno IDE or VS code IDE

#add analytical solution from this paper:
#Ref [1]
#A comprehensive analytical solution of the nonlinear pendulum
#European Journal of Physics

function rond(x)
    z=convert(Int64,floor(x));   #julia can't mix float type with integer type
end

function theta_analytical_noswing(tt)
    z=2*asin(sin(theta0/2)*Jacobi.sn(K(sin(theta0/2)^2)-w0*tt,sin(theta0/2)^2));
end


function file_writer(Nrow,Ncol,Matr,filename_,header_)   #writes array to file
    fff=open(filename_,"w");
    write(fff,header_,"\r\n")
    for i=1:Nrow
        for j=1:Ncol
            write(fff,string(Matr[i,j]));
            if j<Ncol
                write(fff,"\t");
            end
        end
        write(fff,"\r\n");
    end
    close(fff);
end

cd("D:/julia_files/pendulum")

#using Plots
 using PyPlot  #you can't use Plots & PyPlots together
 PyPlot.close();  #close all existing charts

import PyPlot: pygui;
PyPlot.pygui(true);

#for Elliptic integrals & Jacobi elliptic functions
using Elliptic;
import Elliptic;

m=1; #mass
g=9.81; #gravity
L=9.81; #length of string
w0=sqrt(g/L); #natural frequency

theta0=pi/1.02;  #initial conditions
theta0_d=0;

Ep=4*w0^2;
E0=theta0_d^2+Ep*(sin(theta0/2))^2;
k=sqrt(E0/Ep); #ifE0>Ep: swing-over takes place, if E0<Ep, oscillate back&forth
xG=1/k;



T=1000;
dt=0.002;           dt2=0.002;
N=rond(T/dt+1);     N2=rond(T/dt2+1);
t=zeros(N,1);       t2=zeros(N2,1);
theta=zeros(N,1);   theta2=zeros(N2,1);
theta_d=zeros(N,1); theta2_d=zeros(N2,1);

theta[1,1]=theta0;   #initialize
theta_d[1,1]=theta0_d;

theta2[1,1]=theta[1,1];
theta2_d[1,1]=theta_d[1,1];


T_period=4*Elliptic.K( (sin(theta0/2))^2 )/w0; #no swing-over
#Ref. [2]
#formula from Exact solution for the nonlinear pendulum
#by A. Beléndez1; C. Pascual; D.I. Méndez; T. Beléndez; C. Neipp
#web page: http://www.scielo.br/scielo.php?script=sci_arttext&pid=S1806-11172007000400024

dt_analyt=T_period/1000;
N_analyt=rond(T/dt_analyt+1);
t_analyt=zeros(N_analyt,1);
theta_analyt=zeros(N_analyt,1);
for i=1:N_analyt
    tt=(i-1)*dt_analyt;
    t_analyt[i,1]=tt;
    #if tt<=T_period
    theta_analyt[i,1]=theta_analytical_noswing(tt);

end


for i=1:N
    t[i,1]=(i-1)*dt;
end
for i=1:N2
    t2[i,1]=(i-1)*dt2;
end

for i=1:N-1
    thetaS=theta[i,1]+dt*theta_d[i,1];
    theta_dS=theta_d[i,1]-dt*w0^2*sin(theta[i,1]);

    theta[i+1,1]=theta[i,1]+dt*(theta_d[i,1]+theta_dS)/2;
    theta_d[i+1,1]=theta_d[i,1]-dt*w0^2*(sin(theta[i,1])+sin(thetaS))/2;
end

for i=1:N2-1
    for counter_=1:5
        theta2s=theta2[i,1];
        theta2[i + 1, 1] = (cos(theta2s) * dt2 ^ 2 * theta2s * w0 ^ 2 - sin(theta2[i, 1]) * dt2 ^ 2 * w0 ^ 2 - sin(theta2s) * dt2 ^ 2 * w0 ^ 2 + 4 * theta2_d[i, 1] * dt2 + 4 * theta2[i, 1]) / (cos(theta2s) * dt2 ^ 2 * w0 ^ 2 + 4);
        theta2_d[i + 1, 1] = -(cos(theta2s) * dt2 ^ 2 * w0 ^ 2 * theta2_d[i, 1] + 2 * cos(theta2s) * dt2 * w0 ^ 2 * theta2[i, 1] - 2 * cos(theta2s) * dt2 * theta2s * w0 ^ 2 + 2 * sin(theta2[i, 1]) * dt2 * w0 ^ 2 + 2 * sin(theta2s) * dt2 * w0 ^ 2 - 4 * theta2_d[i, 1]) / (cos(theta2s) * dt2 ^ 2 * w0 ^ 2 + 4);
        theta2s=theta2[i+1,1];
    end
end

res_=zeros(N,4);
for i=1:N
    res_[i,1]=t[i,1];
    res_[i,2]=theta[i,1];
    res_[i,3]=theta_d[i,1];
    #res_[i,4]=1-cos(theta[i,1])+1/2*theta_d[i,1]^2;
    res_[i,4]=theta_d[i,1]^2+Ep*(sin(theta[i,1]/2))^2;
end

res2_=zeros(N2,4);
for i=1:N2
    res2_[i,1]=t2[i,1];
    res2_[i,2]=theta2[i,1];
    res2_[i,3]=theta2_d[i,1];
    #res2_[i,4]=1-cos(theta2[i,1])+1/2*theta2_d[i,1]^2;
    res2_[i,4]=theta2_d[i,1]^2+Ep*(sin(theta2[i,1]/2))^2;
end

#plot(t[:,1],theta[:,1]);
header_=string("theta\"+sin(theta)=0",",","Corrector-Predictor,dt=",dt,"\r\n","t\ttheta\ttheta_d\tenergy");
filename_=string("res_ ",round(theta[1,1],digits=3),".txt");
file_writer(N,4,res_,filename_,header_);

PyPlot.subplot(3,1,3);
plot(t_analyt[:,1],theta_analyt[:,1],"k",label="analytical");
plot(res_[:,1],res_[:,2],"-b",label=string("corrector-predictor,dt=",dt));
plot(res2_[:,1],res2_[:,2],":r",label=string("mid-point,dt=",dt));


PyPlot.xlabel("t");
PyPlot.ylabel("theta");
PyPlot.xlim((0,50));
PyPlot.title("first oscillations");
PyPlot.legend();

#PyPlot.figure();
PyPlot.subplot(3,1,1);
plot(t_analyt[:,1],theta_analyt[:,1],"k",label="analytical");
plot(res_[:,1],res_[:,2],"-b",label=string("corrector-predictor,dt=",dt));
plot(res2_[:,1],res2_[:,2],":r",label=string("mid-point,dt=",dt));

PyPlot.xlabel("t");
PyPlot.ylabel("theta");
PyPlot.title("long time");
PyPlot.legend();

#PyPlot.figure();
PyPlot.subplot(3,1,2);
plot([0,T],[E0,E0],"k");#,label("analytical"));
plot(res_[:,1],res_[:,4],"-b");#,label(string("corrector-predictor,dt=",dt)));
plot(res2_[:,1],res2_[:,4],":r");#,label(string("mid-point,dt=",dt)));

PyPlot.xlabel("t");
PyPlot.ylabel("energy");
PyPlot.title("long time");

PyPlot.tight_layout();

savefig("fn.png");
#gcf();
#gcf();


#savefig("fn.pdf");
