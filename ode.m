size=256;
total=size^2;

start_infected=1;

days=500;

s(1)=(total-start_infected)/total;
i(1)=(start_infected)/total;
r(1)=0;
d(1)=0;

beta=.059;
gamma_r=.043;
gamma_d=.003225;

for k=2:days
  s(k)=s(k-1)-beta*s(k-1)*i(k-1);
  i(k)=i(k-1)+beta*s(k-1)*i(k-1)-(gamma_r+gamma_d)*i(k-1);
  r(k)=r(k-1)+gamma_r*i(k-1);
  d(k)=d(k-1)+gamma_d*i(k-1);
endfor

x=1:days;
clf;
hold on
semilogy(x,s(x)*total,'b','LineWidth',1)
semilogy(x,i(x)*total,'r','LineWidth',1)
semilogy(x,r(x)*total,'g','LineWidth',1)
semilogy(x,d(x)*total,'k','LineWidth',1)
set(gca,"ylim",[0 size^2])
set(gca,"xtick",0:100:days)
set(gca,"fontsize",15)
xlabel("Days", 'fontsize', 15)
ylabel("Number of people", 'fontsize', 15)
%legend("Susceptible","Infected","Recovered","Dead")
title({"ODE Model"})
hold off
