%SIRD model
clear

size=50; %side length of square

days=30;

infection_rate=.1; %the rate that each additional neighbor multiplies the infection by
infection_radius=3; %how much taxicab distance away someone can be and still infect
infection_factor=2; %chance goes down by a factor of this for every further distance
death_chance=.01;
recovery_chance=.05;
long_connections=5; %how many "longer distance" connections can infect people

times=100; %ensemble only

%setting up infection mechanism
for radius=1:infection_radius
  for x_offset=0:radius
    infection_matrix(infection_radius+1+x_offset,infection_radius+1+radius-x_offset)=1/infection_factor^(radius-1); 
    infection_matrix(infection_radius+1-x_offset,infection_radius+1+radius-x_offset)=1/infection_factor^(radius-1);
    infection_matrix(infection_radius+1+x_offset,infection_radius+1-radius+x_offset)=1/infection_factor^(radius-1);
    infection_matrix(infection_radius+1-x_offset,infection_radius+1-radius+x_offset)=1/infection_factor^(radius-1);    
  endfor
endfor

for iter=1:times
  lattice_si=zeros(size,size);
  lattice_rd=zeros(size,size); %0 is susceptible

  initial_infections=1;

  for x=1:initial_infections
    initial_x=floor(size*rand())+1; %first position of infected
    initial_y=floor(size*rand())+1;
    lattice_si(initial_x,initial_y)=1; %1 is infected
  endfor

  for j=1:days
    i=rand(size);
    lattice_si_temp=lattice_si;
      
    lattice_neighbors=infection_rate*conv2(lattice_si,infection_matrix,"same")-i;
    lattice_si_temp=lattice_si | (lattice_neighbors>0);
    lattice_rd=lattice_rd+3*(lattice_si.*(i>(1-death_chance))) + 2*(lattice_si.*(i<recovery_chance));
    lattice_si=lattice_si_temp & not(lattice_rd);
    
    end_s(j,iter)=sum(sum((lattice_si+lattice_rd)==0));
    end_i(j,iter)=sum(sum((lattice_si+lattice_rd)==1));
    end_r(j,iter)=sum(sum((lattice_si+lattice_rd)==2));
    end_d(j,iter)=sum(sum((lattice_si+lattice_rd)==3));
    
    if any(any(lattice_si))==0
      break
    endif
  endfor
endfor

%statistics
for j=1:days
  end_s_average(j)=sum(end_s(j,:))/times;
  end_i_average(j)=sum(end_i(j,:))/times;
  end_r_average(j)=sum(end_r(j,:))/times;
  end_d_average(j)=sum(end_d(j,:))/times;
  end_si_average(j)=sum(end_s(j,:).*end_i(j,:))/times;
  end_sr_average(j)=sum(end_s(j,:).*end_r(j,:))/times;
  end_sd_average(j)=sum(end_s(j,:).*end_d(j,:))/times;
  end_ir_average(j)=sum(end_i(j,:).*end_r(j,:))/times;
  end_id_average(j)=sum(end_i(j,:).*end_d(j,:))/times;
  end_rd_average(j)=sum(end_r(j,:).*end_d(j,:))/times;
endfor

figure(1)
x=1:days;
hold on
plot(x,end_s_average(x),'b','LineWidth',1)
plot(x,end_i_average(x),'r','LineWidth',1)
plot(x,end_r_average(x),'g','LineWidth',1)
plot(x,end_d_average(x),'k','LineWidth',1)

figure(2)
hold on
plot(x,(end_si_average(x)-end_s_average(x).*end_i_average(x))./end_si_average(x),'r','LineWidth',1)
plot(x,(end_sr_average(x)-end_s_average(x).*end_r_average(x))./end_sr_average(x),'g','LineWidth',1)
plot(x,(end_sd_average(x)-end_s_average(x).*end_d_average(x))./end_sd_average(x),'k','LineWidth',1)
plot(x,(end_ir_average(x)-end_i_average(x).*end_r_average(x))./end_ir_average(x),'c','LineWidth',1)
plot(x,(end_id_average(x)-end_i_average(x).*end_d_average(x))./end_id_average(x),'m','LineWidth',1)
plot(x,(end_rd_average(x)-end_r_average(x).*end_d_average(x))./end_rd_average(x),'b','LineWidth',1)
