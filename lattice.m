%SIRD model
clear

pkg load statistics

size=50; %side length of square

days=365;

initial_x=6; %first position of infected
initial_y=8;

infection_rate=5; %for each additional neighbor the percent chance of infection goes up by this number
infection_radius=3; %how much taxicab distance away someone can be and still infect
infection_factor=2; %chance goes down by a factor of this for every further distance
death_chance=5;
recovery_chance=10;
r_infection_rate=0; %how much each recovered person counts by as a neighbor, should be a decimal
long_connections=5; %how many "longer distance" connections can infect people

lattice_si=zeros(size.+2.*infection_radius,size.+2.*infection_radius);
lattice_rd=zeros(size.+2.*infection_radius,size.+2.*infection_radius);
%0 is susceptible, the borders are just there to prevent errors
lattice_si(initial_x.+infection_radius,initial_y.+infection_radius)=1; %1 is infected

infection_matrix=zeros(2.*infection_radius.+1,2.*infection_radius.+1); %setting up infection mechanism
for radius=1:infection_radius
  for x_offset=0:radius
    infection_matrix(infection_radius.+1.+x_offset,infection_radius.+1.+radius.-x_offset)=1./infection_factor.^(radius.-1); 
    infection_matrix(infection_radius.+1.-x_offset,infection_radius.+1.+radius.-x_offset)=1./infection_factor.^(radius.-1);
    infection_matrix(infection_radius.+1.+x_offset,infection_radius.+1.-radius.+x_offset)=1./infection_factor.^(radius.-1);
    infection_matrix(infection_radius.+1.-x_offset,infection_radius.+1.-radius.+x_offset)=1./infection_factor.^(radius.-1);    
  endfor
endfor

for j=1:days
  i=rand(size+2.*infection_radius);
  i=floor(100*i); %the probabilities will be out of 100, numbers 0 to 99
  lattice_si_temp=lattice_si;
  infected_count=0;
  
  coords=floor(random("uniform", infection_radius.+1, size.+infection_radius, [4, long_connections]));
  
  for x=infection_radius.+1:size.+infection_radius
    for y=infection_radius.+1:size.+infection_radius
      if lattice_si(x,y)==0 %basic infection mechanism
        if lattice_rd(x,y)==0
          infected_neighbors=0;
          infection_multiply=lattice_si(x.-infection_radius:x.+infection_radius,y.-infection_radius:y.+infection_radius);
          for row=1:2.*infection_radius.+1
            infected_neighbors+=infection_multiply(row,:)*infection_matrix(:,row);
          endfor
          
          for connectioncount=1:long_connections
            if abs(coords(1,connectioncount)-coords(3,connectioncount))+abs(coords(2,connectioncount)-coords(4,connectioncount))>=infection_radius
              if coords(1,connectioncount)==x && coords(2, connectioncount)==y
                if lattice_si(coords(3,connectioncount),coords(4,connectioncount))==1
                  infected_neighbors+=1;
                endif
              endif
              if coords(3,connectioncount)==x && coords(4, connectioncount)==y
                if lattice_si(coords(1,connectioncount),coords(2,connectioncount))==1
                  infected_neighbors+=1;
                endif
              endif        
            endif
          endfor
          
          if i(x,y)<=infection_rate*infected_neighbors
            lattice_si_temp(x,y)=1;
          endif
        endif
      endif
      
      if lattice_si(x,y)==1 %determines what happens to infected
        infected_count+=1;
        if i(x,y)>=100-death_chance
          lattice_rd(x,y)=3; %3 is dead
          lattice_si_temp(x,y)=0;
        endif
        
        if i(x,y)<=recovery_chance-1
          lattice_rd(x,y)=2; %2 is recovered
          lattice_si_temp(x,y)=r_infection_rate;
        endif
      endif
      
    endfor
  endfor
  lattice_si=lattice_si_temp;
  if infected_count==0
    break
  endif
endfor
%plot
end_s=0;
end_i=0;
end_r=0;
end_d=0;
figure(2)
hold on
for x=infection_radius.+1:size.+infection_radius
  for y=infection_radius.+1:size.+infection_radius
    if lattice_rd(x,y)==2
      plot(x, y, 'og')
      end_r+=1;
    endif
    if lattice_rd(x,y)==3
      plot(x, y, 'ok')
      end_d+=1;
    endif
    if lattice_rd(x,y)==0
      if lattice_si(x,y)==0
        plot(x, y, 'ob')
        end_s+=1;
      endif
      if lattice_si(x,y)==1
        plot(x, y, 'or')
        end_i+=1;
      endif
    endif
  endfor
endfor
end_s
end_i
end_r
end_d