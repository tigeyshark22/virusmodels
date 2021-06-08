%SIRD model
clear

size=256; %side length of square

days=150;

lattice_si=zeros(size,size);
lattice_rd=zeros(size,size); %0 is susceptible
lattice_vaccinated=zeros(size,size);
lattice_i_days=zeros(size,size);

initial_infections=1;
infection_rate=.05; %the rate that each additional neighbor multiplies the infection probability by
infection_radius=3; %how much taxicab distance away someone can be and still infect
infection_factor=2; %chance goes down by a factor of this for every further distance
death_chance=.0003; %chances are PER DAY
recovery_chance=.004;
%long_connections=5; %how many "longer distance" connections can infect people

mask_probability=0;
mask_effectiveness=.8; %how much of the time it prevents infection

elderly_proportion=0; %reflects current statistics for the world
elderly_size=floor(size*sqrt(elderly_proportion));
elderly_check=ones(elderly_size,elderly_size);
elderly_vulnerability=2; %how many times as vulnerable they are

vaccination_rate=0; %proportion of susceptible population vaccinated on a given day
vaccination_spread_rate=.5; %chance of spread from vaccinated people (proportion that the computed
                            %infection rate/factor result is multiplied by)
vaccination_threshold=.12; %proportion of people I/R/D before vaccination starts

for x=1:initial_infections
  %initial_x=floor(size*rand())+1; %first position of infected
  %initial_y=floor(size*rand())+1;
  initial_x=floor(size/2)+1;
  initial_y=floor(size/2)+1;
  lattice_si(initial_x,initial_y)=1; %1 is infected
endfor

i=rand(size);
mask_wearers=1-(mask_effectiveness*(i<mask_probability));

if elderly_proportion==0
  elderly_matrix=ones(size,size);
else
  initial_x=floor((size-elderly_size)*rand())+1+floor(elderly_size/2); %first position of infected
  initial_y=floor((size-elderly_size)*rand())+1+floor(elderly_size/2);
  elderly_matrix=zeros(size,size);
  elderly_matrix(initial_x,initial_y)=1;
  elderly_matrix=1+(elderly_vulnerability-1)*(conv2(elderly_matrix,elderly_check,"same")>0);
endif

%setting up infection mechanism
for radius=1:infection_radius
  for x_offset=0:radius
    infection_matrix(infection_radius+1+x_offset,infection_radius+1+radius-x_offset)=1/infection_factor^(radius-1); 
    infection_matrix(infection_radius+1-x_offset,infection_radius+1+radius-x_offset)=1/infection_factor^(radius-1);
    infection_matrix(infection_radius+1+x_offset,infection_radius+1-radius+x_offset)=1/infection_factor^(radius-1);
    infection_matrix(infection_radius+1-x_offset,infection_radius+1-radius+x_offset)=1/infection_factor^(radius-1);    
  endfor
endfor

for j=1:days
  i=rand(size);
  
  lattice_i_days=lattice_i_days+(lattice_si==1); %each day the chance of death/recovery goes up linearly

  lattice_neighbors=(infection_rate*conv2(lattice_si,infection_matrix,"same")).*mask_wearers.*elderly_matrix-i;
  lattice_si_temp=(((lattice_si==0) & not(lattice_vaccinated)) & (lattice_neighbors>0))+lattice_si;
  %lattice_si_temp=lattice_si | (lattice_neighbors>0);
  
  i=rand(size);
  lattice_rd=lattice_rd+3*(lattice_si.*(i>(1-(death_chance*lattice_i_days))))+2*(lattice_si.*(i<(recovery_chance*lattice_i_days)));
  lattice_si=lattice_si_temp & not(lattice_rd);
  
  end_s(j)=sum(sum((lattice_si+lattice_rd)==0 & lattice_vaccinated==0));
  end_i(j)=sum(sum((lattice_si+lattice_rd)==1));
  end_r(j)=sum(sum((lattice_si+lattice_rd)==2));
  end_d(j)=sum(sum((lattice_si+lattice_rd)==3));
  end_v(j)=sum(sum(lattice_vaccinated==1));
  
  i=rand(size);
  if j>1 && (end_i(j-1)+end_r(j-1)+end_d(j-1))>size^2*vaccination_threshold %checks if vaccination occurs
   lattice_vaccinated=lattice_vaccinated+((lattice_vaccinated==0).*((lattice_si+lattice_rd)==0).*(i<vaccination_rate));
  endif
  
  if any(any(lattice_si))==0
    break
  endif
endfor

%plot

figure(1)
clf;
x=1:size;
y=1:size;
hold on
contourf(x,y,lattice_si(x,y)+lattice_rd(x,y),0:4)
if elderly_proportion!=0
  plot(initial_x-elderly_size/2:initial_x+elderly_size/2,initial_y+elderly_size/2:initial_y+elderly_size/2,'r');
  plot(initial_x-elderly_size/2:initial_x+elderly_size/2,initial_y-elderly_size/2:initial_y-elderly_size/2,'r');
  plot(initial_x+elderly_size/2:initial_x+elderly_size/2,initial_y-elderly_size/2:initial_y+elderly_size/2,'r');
  plot(initial_x-elderly_size/2:initial_x-elderly_size/2,initial_y-elderly_size/2:initial_y+elderly_size/2,'r');
endif
hold off
title ({"Final lattice"});

figure(2)
clf;
x=1:j;
hold on
plot(x,end_s(x),'b','LineWidth',2)
plot(x,end_i(x),'r','LineWidth',2)
plot(x,end_r(x),'g','LineWidth',2)
plot(x,end_d(x),'k','LineWidth',2)
%plot(x,end_v(x),'y','LineWidth',2)
set(gca,"ylim",[0 size^2])
set(gca,"xtick",0:100:days)
set(gca,"fontsize",15)
xlabel("Days", 'fontsize', 15)
ylabel("Number of people", 'fontsize', 15)
legend("Susceptible (not vaccinated)","Infected","Recovered","Dead","Vaccinated")
title({"Stochastic Model"})

figure(3)
clf;
hold on
x=2:j;
plot(x, (end_s(x)-end_s(x-1))./(end_s(x-1).*end_i(x-1)),'k','LineWidth',1)
set(gca,"fontsize",15)
set(gca,"ylim",[-.000005 0])
title({"Change in S divided by SI"})
xlabel("Days")
hold off