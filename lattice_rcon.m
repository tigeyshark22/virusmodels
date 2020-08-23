
%SIRD model
clear

size=256; %side length of square

days=300;

infection_rate=.05; %the rate that each additional neighbor multiplies the infection by
infection_radius=3; %how much taxicab distance away someone can be and still infect
infection_factor=2; %chance goes down by a factor of this for every further distance
death_chance=.0003;
recovery_chance=.004;
long_connections=5; %how many "longer distance" connections can infect people
success_connect=.5;

times=20; %ensemble only

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
  lattice_si(:,:,iter)=zeros(size,size);
  lattice_rd(:,:,iter)=zeros(size,size); %0 is susceptible
  lattice_i_days=zeros(size,size);
  lattice_random(:,:,iter)=floor(rand(infection_radius*2+1)+success_connect);
  
  initial_infections=1;

  for x=1:initial_infections
%    initial_x=floor(size*rand())+1; %first position of infected
%    initial_y=floor(size*rand())+1;
    initial_x=floor(size/2)+1;
    initial_y=floor(size/2)+1;
    lattice_si(initial_x,initial_y,iter)=1;
  endfor

  for j=1:days
    i=rand(size);
    lattice_si_temp=lattice_si(:,:,iter);
    
    lattice_i_days=lattice_i_days+(lattice_si(:,:,iter)==1); %each day the chance of recovery/death goes up linearly
      
    lattice_neighbors=infection_rate*conv2(lattice_si(:,:,iter),infection_matrix.*lattice_random(:,:,iter),"same")-i;
    lattice_si_temp=lattice_si(:,:,iter) | (lattice_neighbors>0);
    lattice_rd(:,:,iter)=lattice_rd(:,:,iter)+3*(lattice_si(:,:,iter).*(i>(1-(death_chance*lattice_i_days))))+2*(lattice_si(:,:,iter).*(i<(recovery_chance*lattice_i_days)));
    lattice_si(:,:,iter)=lattice_si_temp & not(lattice_rd(:,:,iter));
    
    end_s(j,iter)=sum(sum((lattice_si(:,:,iter)+lattice_rd(:,:,iter))==0));
    end_i(j,iter)=sum(sum((lattice_si(:,:,iter)+lattice_rd(:,:,iter))==1));
    end_r(j,iter)=sum(sum((lattice_si(:,:,iter)+lattice_rd(:,:,iter))==2));
    end_d(j,iter)=sum(sum((lattice_si(:,:,iter)+lattice_rd(:,:,iter))==3));
    
    if any(any(lattice_si(:,:,iter)))==0
      end_s(:,iter)=end_s(:,iter-1);
      end_i(:,iter)=end_i(:,iter-1);
      end_r(:,iter)=end_r(:,iter-1);
      end_d(:,iter)=end_d(:,iter-1);
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
  end_ss_average(j)=sum(end_s(j,:).*end_s(j,:))/times;
  end_si_average(j)=sum(end_s(j,:).*end_i(j,:))/times;
  end_sr_average(j)=sum(end_s(j,:).*end_r(j,:))/times;
  end_sd_average(j)=sum(end_s(j,:).*end_d(j,:))/times;
  end_ii_average(j)=sum(end_i(j,:).*end_i(j,:))/times;
  end_ir_average(j)=sum(end_i(j,:).*end_r(j,:))/times;
  end_id_average(j)=sum(end_i(j,:).*end_d(j,:))/times;
  end_rr_average(j)=sum(end_r(j,:).*end_r(j,:))/times;
  end_rd_average(j)=sum(end_r(j,:).*end_d(j,:))/times;
  end_dd_average(j)=sum(end_d(j,:).*end_d(j,:))/times;
endfor

for x=1:size
  for y=1:size
    lattice_si_average(x,y)=sum(lattice_si(x,y,:))/times;
    lattice_rd_average(x,y)=sum(lattice_rd(x,y,:))/times;
  endfor
endfor
figure(1)
x=1:days;
clf;
hold on
plot(x,end_s_average(x),'b','LineWidth',1)
plot(x,end_i_average(x),'r','LineWidth',1)
plot(x,end_r_average(x),'g','LineWidth',1)
plot(x,end_d_average(x),'k','LineWidth',1)
set(gca,"ylim",[0 size^2])
set(gca,"xtick",0:100:days)
set(gca,"fontsize",15)
xlabel("Days", 'fontsize', 15)
ylabel("Number of people", 'fontsize', 15)
%legend("Susceptible","Infected","Recovered","Dead")
%title({"Stochastic Model"})
hold off

figure(2)
clf;
hold on
plot(x,(end_si_average(x)-end_s_average(x).*end_i_average(x))./((end_s_average(x)+1).*(end_i_average(x)+1)),'r','LineWidth',1)
plot(x,(end_sr_average(x)-end_s_average(x).*end_r_average(x))./((end_s_average(x)+1).*(end_r_average(x)+1)),'g','LineWidth',1)
plot(x,(end_sd_average(x)-end_s_average(x).*end_d_average(x))./((end_s_average(x)+1).*(end_d_average(x)+1)),'k','LineWidth',1)
plot(x,(end_ir_average(x)-end_i_average(x).*end_r_average(x))./((end_i_average(x)+1).*(end_r_average(x)+1)),'c','LineWidth',1)
plot(x,(end_id_average(x)-end_i_average(x).*end_d_average(x))./((end_i_average(x)+1).*(end_d_average(x)+1)),'m','LineWidth',1)
plot(x,(end_rd_average(x)-end_r_average(x).*end_d_average(x))./((end_r_average(x)+1).*(end_d_average(x)+1)),'b','LineWidth',1)
legend('SI','SR','SD','IR','ID','RD')
hold off

figure(3)
clf;
hold on
plot(x,(end_ss_average(x)-end_s_average(x).*end_s_average(x))./((end_s_average(x)+1).*(end_s_average(x)+1)),'b','LineWidth',1)
plot(x,(end_ii_average(x)-end_i_average(x).*end_i_average(x))./((end_i_average(x)+1).*(end_i_average(x)+1)),'r','LineWidth',1)
plot(x,(end_rr_average(x)-end_r_average(x).*end_r_average(x))./((end_r_average(x)+1).*(end_r_average(x)+1)),'g','LineWidth',1)
plot(x,(end_dd_average(x)-end_d_average(x).*end_d_average(x))./((end_d_average(x)+1).*(end_d_average(x)+1)),'k','LineWidth',1)
legend('SS','II','RR','DD')
hold off

figure(4)
clf;
hold on
x=2:days;
plot(x, (end_s_average(x)-end_s_average(x-1))./(end_s_average(x-1).*end_i_average(x-1)),'k','LineWidth',1)
title({"Change in S divided by SI"})
xlabel("Days")
hold off

figure(5)
clf;
hold on
x=1:size;
y=1:size;
%contourf(x,y,lattice_si_average(x,y)+lattice_rd_average(x,y),0:4)
hold off
