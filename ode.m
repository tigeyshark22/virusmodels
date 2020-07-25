size=256;
total=size^2;

start_infected=10;

days=300;

s(1)=(total-start_infected)/total;
i(1)=(start_infected)/total;
r(1)=0;
d(1)=0;

beta=.1;
gamma_r=.04;
gamma_d=.003;

for k=2:days
  s(k)=s(k-1)-beta*s(k-1)*i(k-1);
  i(k)=i(k-1)+beta*s(k-1)*i(k-1)-(gamma_r+gamma_d)*i(k-1);
  r(k)=r(k-1)+gamma_r*i(k-1);
  d(k)=d(k-1)+gamma_d*i(k-1);
endfor

x=1:days;
clf;
hold on
plot(x,s(x)*total,'b','LineWidth',1)
plot(x,i(x)*total,'r','LineWidth',1)
plot(x,r(x)*total,'g','LineWidth',1)
plot(x,d(x)*total,'k','LineWidth',1)
xlabel("Days")
ylabel("Number of people")
legend("Susceptible","Infected","Recovered","Dead")
title({"ODE Model"})
hold off