function ode_func (days2,beta2,r2,d2,legendon,exportgraph,exportlog,gname,logname,)

size=256;
total=size^2;

start_infected=1;

days=days2;

s(1)=(total-start_infected)/total;
i(1)=(start_infected)/total;
r(1)=0;
d(1)=0;

beta=beta2;
gamma_r=r2;
gamma_d=d2;

for k=2:days
  s(k)=s(k-1)-beta*s(k-1)*i(k-1);
  i(k)=i(k-1)+beta*s(k-1)*i(k-1)-(gamma_r+gamma_d)*i(k-1);
  r(k)=r(k-1)+gamma_r*i(k-1);
  d(k)=d(k-1)+gamma_d*i(k-1);
endfor

if exportgraph==1
  figure(1)
  x=1:days;
  clf;
  hold on
  plot(x,s(x)*total,'b','LineWidth',2)
  plot(x,i(x)*total,'r','LineWidth',2)
  plot(x,r(x)*total,'g','LineWidth',2)
  plot(x,d(x)*total,'k','LineWidth',2)
  set(gca,"ylim",[0 size^2])
  set(gca,"xtick",0:100:days)
  set(gca,"fontsize",15)
  xlabel("Days", 'fontsize', 15)
  ylabel("Number of people", 'fontsize', 15)
  if legendon==1
    legend("Susceptible","Infected","Recovered","Dead")
  endif
  title({"ODE Model"})
  hold off
  print (['C:\Users\Owen Yang\Dropbox\virusmodels\' gname], "-dpng")
endif

if exportlog==1
  figure(2)
  x=1:days;
  clf;
  hold on
  semilogy(x,s(x)*total,'b','LineWidth',2)
  semilogy(x,i(x)*total,'r','LineWidth',2)
  semilogy(x,r(x)*total,'g','LineWidth',2)
  semilogy(x,d(x)*total,'k','LineWidth',2)
  set(gca,"ylim",[0 size^2])
  set(gca,"xtick",0:100:days)
  set(gca,"fontsize",15)
  xlabel("Days", 'fontsize', 15)
  ylabel("Number of people", 'fontsize', 15)
  %legend("Susceptible","Infected","Recovered","Dead")
  title({"ODE Model"})
  hold off
  print (['C:\Users\Owen Yang\Dropbox\virusmodels\' logname], "-dpng")
endif

endfunction