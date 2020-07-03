%SIRD model
clear

size=50; %side length of square

days=30;

lattice_si=zeros(size,size);
lattice_rd=zeros(size,size); %0 is susceptible

initial_infections=1;

for x=1:initial_infections
  initial_x=floor(size*rand())+1; %first position of infected
  initial_y=floor(size*rand())+1;
  lattice_si(initial_x,initial_y)=1; %1 is infected
endfor

infection_rate=.1; %the rate that each additional neighbor multiplies the infection by
infection_radius=3; %how much taxicab distance away someone can be and still infect
infection_factor=2; %chance goes down by a factor of this for every further distance
death_chance=.01;
recovery_chance=.05;
long_connections=5; %how many "longer distance" connections can infect people

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
    
  lattice_neighbors=infection_rate*conv2(lattice_si,infection_matrix,"same")-i;
  lattice_si_temp=lattice_si | (lattice_neighbors>0);
  lattice_rd=lattice_rd+3*(lattice_si.*(i>(1-death_chance))) + 2*(lattice_si.*(i<recovery_chance));
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
colormap ("default");
x=1:size;
y=1:size;
contourf(x,y,lattice_si(x,y)+lattice_rd(x,y),0:4)
title ({"Final lattice"});

figure(2)
x=1:days;
hold on
plot(x,end_s(x),'b','LineWidth',1)
plot(x,end_i(x),'r','LineWidth',1)
plot(x,end_r(x),'g','LineWidth',1)
plot(x,end_d(x),'k','LineWidth',1)
