%% Turn Off v3
%   The goal of this script is to look at how "sharp" the wave is. In
%   particular, we want to know how closely neighboring nuclei turn off as
%   opposed to non-neighboring nuclei. This gives figure 5 in the new
%   write-up. 
%%
Folder = 'Eve'; 
Names={'2014-03-14-Eve2B','2014-03-19-Eve2A', '2014-03-20-Eve2A',...
    '2014-03-20-Eve2B', '2014-03-20-Eve2C','2014-03-20-Eve2D',...
    '2014-03-20-Eve2E'}; 
NC=14;maxTDU=5;
NumNames=length(Names);
BinSize=0.02; 
Bins=0.32:BinSize:0.5;
OffTime={};
%%
for ii=1:NumNames
    
     Name = Names{ii}; 
     Directory = [Folder '/' Name];
     load(['\\Client\C$\Users\wtredman\Desktop\GregorLab\' Directory '\_data_NC' num2str(NC) '.mat'])
     APpos=data.ParticleCenters_x_APpos;
     
     OffTimeMat=zeros(maxTDU,length(Bins)-1);
     %This code is very similar to that of the mutual information stuff,
     %but a little bit more streamlined
     for jj=1:length(Bins)-1
         partInbin=find(APpos>Bins(jj) & APpos<Bins(jj+1)); 
         if length(partInbin)~=1 && isempty(partInbin)~=1        
             partInbinMat=zeros(maxTDU,length(partInbin));
             for kk=1:length(partInbin)     
                 partFinalTime=data.AllFrames{partInbin(kk)}{1}(end);
                 TDUdiffarray=zeros(1,maxTDU);                 
                 for ll=1:maxTDU                 
                     neighbors=find(data.ParticleNucleiTopoDistances(partInbin(kk),:)==ll);                                         
                     if isempty(neighbors)~=1  
                        diffarray=nan(1,1);                        
                        for mm=1:length(neighbors)    
                            if isempty(intersect(neighbors(mm),partInbin))~=1
                                neighFinalTime=data.AllFrames{neighbors(mm)}{1}(end);
                                diffarray=abs(neighFinalTime-partFinalTime);
                            end
                        end
                        TDUdiffarray(ll)=mean(diffarray);
                     else
                         TDUdiffarray(ll)=nan;
                     end
                 end
                 partInbinMat(:,kk)=TDUdiffarray;
             end
             OffTimeMat(:,jj)=nanmean(partInbinMat')';
         else
             OffTimeMat(:,jj)=[nan nan nan nan nan]';
         end
     end
     OffTime{ii}=OffTimeMat;
end
%%
AveragedOffMat=zeros(maxTDU,length(Bins)-1);
StdOffMat=zeros(maxTDU,length(Bins)-1);
for ii=1:length(Bins)-1
    for jj=1:maxTDU
        array=zeros(1,NumNames);
        for kk=1:NumNames
            array(kk)=OffTime{kk}(jj,ii);
        end
        AveragedOffMat(jj,ii)=nanmean(array);
        StdOffMat(jj,ii)=nanstd(array);
    end
end
%% Plotting
%%  Difference in Off Time all 5 TDU
figure
MeanBins=(Bins(1:end-1)+Bins(2:end))/2;
for ii=1:maxTDU
    errorbar(MeanBins,AveragedOffMat(ii,:),StdOffMat(ii,:),'-o','color',...
        [(ii-1)/maxTDU, .1 1-(ii-1)/maxTDU],'LineWidth', 1.5,...
        'MarkerFaceColor', [(ii-1)/maxTDU, .1 1-(ii-1)/maxTDU]), hold on
end
xlabel('AP');
ylabel('Difference in Final Turn Off Time');
%%  Difference in Off Time 1 TDU vs. 2-5 TDU
figure
errorbar(MeanBins,AveragedOffMat(1,:),StdOffMat(1,:),'-o','color',...
    [22 165 100]/255), hold on
two2fiveTDU=mean(AveragedOffMat(2:end,:),1);
errorbar(MeanBins,two2fiveTDU,mean(StdOffMat(2:end,:),1),'-o','color',[237 164 17]/255);
xlabel('AP');
ylabel('Difference in Final Turn Off Time');

        
    
    
                        
     