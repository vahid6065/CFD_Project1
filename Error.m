function [ error ] = Error(U_Mean,Grid )
error=zeros(1,1);
%%%%Calculate the relative error value
    for i=1:4
        for j=1:Grid-1
            error(i,j)=U_Mean(i,j+1)-U_Mean(i,j);
        end
    end
end