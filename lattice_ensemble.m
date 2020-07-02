%SIRD model
clear
pkg load statistics

size=50; %side length of square

days=30;

initial_x=floor(size*rand())+1; %first position of infected
initial_y=floor(size*rand())+1;

infection_rate=.2; %the rate that each additional neighbor multiplies the infection by
infection_radius=3; %how much taxicab distance away someone can be and still infect
infection_factor=2; %chance goes down by a factor of this for every further distance
death_chance=.01;
recovery_chance=.05;
long_connections=5; %how many "longer distance" connections can infect people

lattice_si=zeros(size,size);
lattice_rd=zeros(size,size);
%0 is susceptible
lattice_si(initial_x,initial_y)=1; %1 is infected

%setting up infection mechanism
for radius=1:infection_radius
  for x_offset=0:radius
    infection_matrix(infection_radius+1+x_offset,infection_radius+1+radius-x_offset)=1/infection_factor^(radius-1); 
    infection_matrix(infection_radius+1-x_offset,infection_radius+1+radius-x_offset)=1/infection_factor^(radius-1);
    infection_matrix(infection_radius+1+x_offset,infection_radius+1-radius+x_offset)=1/infection_factor^(radius-1);
    infection_matrix(infection_radius+1-x_offset,infection_radius+1-radius+x_offset)=1/infection_factor^(radius-1);    
  endfor
endfor

for times=1:100
  for j=1:days
    i=rand(size);
    lattice_si_temp=lattice_si;
      
    lattice_neighbors=infection_rate*conv2(lattice_si,infection_matrix,"same")-i;
    lattice_si_temp=lattice_si | (lattice_neighbors>0);
    lattice_rd=lattice_rd+3*(lattice_si.*(i>(1-death_chance))) + 2*(lattice_si.*(i<recovery_chance));
    lattice_si=lattice_si_temp & not(lattice_rd);
    
    coords=floor(random("uniform",1,size,[4, long_connections]));
    for connectioncount=1:long_connections
      if abs(coords(1,connectioncount)-coords(3,connectioncount))+abs(coords(2,connectioncount)-coords(4,connectioncount))>=infection_radius
        x=coords(1,connectioncount);
        y=coords(2,connectioncount);
        x2=coords(3,connectioncount);
        y2=coords(4,connectioncount);
        
        for jj=1:2
          if lattice_si(x,y)==0 && lattice_rd(x,y)==0
            infected_neighbors=0;
            infection_multiply=lattice_si(x.-infection_radius:x.+infection_radius,y.-infection_radius:y.+infection_radius);
            for row=1:2.*infection_radius.+1
              infected_neighbors+=infection_multiply(row,:)*infection_matrix(:,row);
            endfor
            
            if lattice_si(x2,y2)==1
              infected_neighbors+=1;
            endif
            
            if i(x,y)<=infection_rate*infected_neighbors %because i remains the same, this can be repeated in the main loop without conflict
              lattice_si_temp(x,y)=1;
            endif
            x=coords(3,connectioncount);
            y=coords(4,connectioncount);
            x2=coords(1,connectioncount);
            y2=coords(2,connectioncount);
          endif
        endfor
      endif
    endfor
    
    end_s(j)=sum(sum((lattice_si+lattice_rd)==0));
    end_i(j)=sum(sum((lattice_si+lattice_rd)==1));
    end_r(j)=sum(sum((lattice_si+lattice_rd)==2));
    end_d(j)=sum(sum((lattice_si+lattice_rd)==3));
    
    if any(any(lattice_si))==0
      break
    endif
  endfor
endfor

%statistics
