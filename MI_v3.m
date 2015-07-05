function [Ixy] = MI_v3(neighbor,particle,OnMat)
%% MI_v3 This calculates the mutual information between two nuclei
%   We call this in MutualInformation_v3.m for two nuclei in the same bin a
%   given TDU apart. See MutualInformation_v3.m for more information. 
%%
X=OnMat(particle,1:end);
Y=OnMat(neighbor,1:end);
Px1=mean(X); Py1=mean(Y); %p(x=1) and p(y=1) 
Px0=1-Px1;Py0=1-Py1;
Px=[Px0,Px1]; Py=[Py0,Py1];
Ixy=0;%I(X,Y)
for ii=0:1
    for jj=0:1
        Pxy=length(find(X==ii & Y==jj))/length(X); %p(x,y)
        if Pxy~=0 %otherwise log2(0)*0=Nan
            Ixy=Ixy+(Pxy*log2(Pxy/(Px(ii+1)*Py(jj+1))));
        else
            Ixy=Ixy;
        end
    end
    
end
end

