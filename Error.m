function [ error_3 ] = Error(U_Mean,Grid )
error_3=zeros(1,1);
    for i=1:4
        for j=1:Grid-1
            error_3(i,j)=U_Mean(i,j+1)-U_Mean(i,j);
        end
    end
end