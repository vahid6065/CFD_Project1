function [ Delta_y ] = Error_Slope( Eror,H ,Grid)
e = zeros(1,1);
delta_e=zeros(1,1);
Delta_y = zeros(1,1);

for jj=1:Grid-1
     delta_e(jj)=(H(jj+1)/H(jj));
end
for i=1:4
    for j=1:Grid-2
        e(i,j)=(Eror(i,j+1)/Eror(i,j));

    end

    
    for j=1:Grid-2
          Delta_y(i,j)=abs(log10(e(i,j)/delta_e(j)));
    end

end

end