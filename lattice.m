%SIRD model
clear

size=256; %side length of square

days=300;

lattice_si=zeros(size,size);
lattice_rd=zeros(size,size); %0 is susceptible
lattice_i_days=zeros(size,size);

initial_infections=1;
infection_rate=.05; %the rate that each additional neighbor multiplies the infection by
infection_radius=3; %how much taxicab distance away someone can be and still infect
infection_factor=2; %chance goes down by a factor of this for every further distance
death_chance=.0003; %chances are PER DAY
recovery_chance=.004;
%long_connections=5; %how many "longer distance" connections can infect people

mask_probability=0;
mask_effectiveness=.5; %how much of the time it prevents infection

elderly_probability=0; %reflects current statistics for the world
elderly_vulnerability=6; %how many times as vulnerable they are

for x=1:initial_infections
  initial_x=floor(size*rand())+1; %first position of infected
  initial_y=floor(size*rand())+1;
  initial_x=floor(size/2)+1;
  initial_y=floor(size/2)+1;
  lattice_si(initial_x,initial_y)=1; %1 is infected
endfor

i=rand(size);
mask_wearers=1-(mask_effectiveness*(i<mask_probability));

i=rand(size);
elderly=1+(elderly_vulnerability-1)*(i<elderly_probability);

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
  lattice_si_temp=lattice_si;
  
  lattice_i_days=lattice_i_days+(lattice_si==1); %each day the chance of death/recovery goes up linearly
    
  lattice_neighbors=(infection_rate*conv2(lattice_si,infection_matrix,"same")).*mask_wearers.*elderly-i;
  lattice_si_temp=lattice_si | (lattice_neighbors>0);
  lattice_rd=lattice_rd+3*(lattice_si.*(i>(1-(death_chance*lattice_i_days))))+2*(lattice_si.*(i<(recovery_chance*lattice_i_days)));
  lattice_si=lattice_si_temp & not(lattice_rd);
  
  end_s(j)=sum(sum((lattice_si+lattice_rd)==0));
  end_i(j)=sum(sum((lattice_si+lattice_rd)==1));
  end_r(j)=sum(sum((lattice_si+lattice_rd)==2));
  end_d(j)=sum(sum((lattice_si+lattice_rd)==3));
  
  if any(any(lattice_si))==0
    break
  endif
endfor

%plot

figure(1)
clf;
x=1:size;
y=1:size;
%contourf(x,y,lattice_si(x,y)+lattice_rd(x,y),0:4)
title ({"Final lattice"});

figure(2)
clf;
x=1:days;
hold on
plot(x,end_s(x),'b','LineWidth',1)
plot(x,end_i(x),'r','LineWidth',1)
plot(x,end_r(x),'g','LineWidth',1)
plot(x,end_d(x),'k','LineWidth',1)
xlabel("Days")
ylabel("Number of people")
legend("Susceptible","Infected","Recovered","Dead")
