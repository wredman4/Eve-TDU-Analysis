function [F]= MutualInfoMovie(time)
%% Mutual Info For Movie
%   Because the movie maker script requires some changes in the
%   MutualInformation_v3.m, I decided it would be easier just to make
%   another copy and make the changes here.
Folder = 'Eve'; 
Names={'2014-03-14-Eve2B','2014-03-19-Eve2A', '2014-03-20-Eve2A',...
    '2014-03-20-Eve2B', '2014-03-20-Eve2C','2014-03-20-Eve2D',...
    '2014-03-20-Eve2E'}; 
NC=14;maxTDU=5;
NumNames=length(Names);
BinSize=0.015; 
Bins=0.32:BinSize:0.48;
%% MI Calculation
MIbyNamebyBin=struct('Bin',{});
for ii=1:NumNames
    
     Name = Names{ii}; 
     Directory = [Folder '/' Name];
     load(['\\Client\C$\Users\wtredman\Desktop\GregorLab\' Directory '\_data_NC' num2str(NC) '.mat'])
     APpos=data.ParticleCenters_x_APpos;
     
     dataMI=struct('Bin',[],'MI',{});
     dataBin=struct('MI',{});
     %%
     X=zeros(length(APpos),floor(time*length(data.NC_Frames))); 
     for jj=1:length(APpos)
         frames=data.AllFrames{jj}{1};
         frames=1+frames-min(data.NC_Frames);
         frames(frames<1)=[]; frames(frames>floor(time*length(data.NC_Frames)))=[];
         X(jj,frames)=1;%X stores when in time a given nucleus is on
     end
     %% 
     for jj=1:length(Bins)-1
         partInbin=find(APpos>Bins(jj) & APpos<Bins(jj+1)); 
         if length(partInbin)~=1 && isempty(partInbin)~=1        
             for kk=1:length(partInbin)            
                 TDUMIarray=zeros(1,maxTDU);                 
                 for ll=1:maxTDU                 
                     neighbors=find(data.ParticleNucleiTopoDistances(partInbin(kk),:)==ll);                                         
                     if isempty(neighbors)~=1  
                        sameTDUarray=zeros(1,1);                        
                        for mm=1:length(neighbors)    
                            if isempty(intersect(neighbors(mm),partInbin))~=1
                                sameTDUarray(mm)=MI_v3(partInbin(kk),neighbors(mm),X);
                            end                  
                         end                         
                         TDUMIarray(ll)=mean(sameTDUarray);%averaging the MI of all neighbors for a given particle for a given TDU                    
                     else                         
                         TDUMIarray(ll)=nan;                     
                     end
                 end                 
                 dataMI(partInbin(kk)).MI=TDUMIarray;
                 dataMI(partInbin(kk)).Bin=jj;
             end
         end
     end     
     %%
     for jj=1:length(Bins)-1     
         ind=find([dataMI.Bin]==jj);         
         if isempty(ind)~=1         
             binMIvalues=[dataMI(ind).MI];
             binMIvaluesarray=zeros(1,maxTDU);             
             for kk=1:maxTDU             
                 TDUarray=zeros(1,length(binMIvalues)/maxTDU);
                 numberofdatasets=1:length(binMIvalues)/maxTDU;
                 TDUarray(numberofdatasets)=binMIvalues(kk+(numberofdatasets-1)*maxTDU);
                 TDUarray=TDUarray(~isnan(TDUarray));                
                 if isempty(TDUarray)~=1                 
                     binMIvaluesarray(kk)=mean(TDUarray);                 
                 else                     
                     binMIvaluesarray(kk)=nan;                
                 end
             end             
             dataBin(jj).MI=binMIvaluesarray;
         else
             dataBin(jj).MI=[nan,nan,nan,nan,nan];
         end
     end     
     MIbyNamebyBin(ii).Bin=dataBin;
end
%% Averaging MI by Bin
%AverageMI=struct('Bin',{});
AverageMI=zeros(maxTDU,length(Bins)-1);
for ii=1:length(Bins)-1
    for jj=1:maxTDU
        array2=zeros(1,NumNames);
        for kk=1:NumNames
            mi=MIbyNamebyBin(kk).Bin(ii).MI;
            if isempty(mi)==1 || isnan(mi(jj))==1
                array2(kk)=nan;
            else
                array2(kk)=mi(jj);
            end
        end
        AverageMI(jj,ii)=mean(array2(~isnan(array2)));
    end
end
%% Movie Frames
close all
figure
MeanBins= (Bins(2:end) + Bins(1:end-1))/2;
for ii=1:maxTDU 
    plot(MeanBins,AverageMI(ii,:),'-o','color',[(ii-1)/maxTDU, .1 1-(ii-1)/maxTDU],...
    'LineWidth', 1.5, 'MarkerFaceColor', [(ii-1)/maxTDU, .1 1-(ii-1)/maxTDU]), hold on
end
F=getframe;
xlabel('AP');
ylabel('Averaged Mutual Information');
%   Or if you want to see the 1 vs. 2-5 movie
% figure 
% %plot(MeanBins,AverageMI(1,:),'-o','color',[22 165 100]/255), hold on
% two2fiveTDU=mean(AverageMI(2:end, :), 1);
% %plot(MeanBins,two2fiveTDU,'-o','color',[237 164 17]/255);
% xlabel('AP');
% ylabel('Averaged Mutual Information')
% %F=getframe;



            
         
                 