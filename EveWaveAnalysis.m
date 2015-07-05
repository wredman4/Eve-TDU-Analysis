%% Stripe Speed Script
%   This script looks at both the width of the stripe through time, and
%   then, after applying a smoother and the Shape Language Modelling toolbox 
%   (this is something you must have downloaded before:http://www.mathworks
%   .com/matlabcentral/fileexchange/24443-slm-shape-language-modeling),
%   differentiates and finds the velocity of the stripe. This is the
%   script that generated figure 6 in the revised write-up. 
%%
Folder = 'Eve'; 
Names={'2014-03-14-Eve2B','2014-03-19-Eve2A', '2014-03-20-Eve2A',...
    '2014-03-20-Eve2B', '2014-03-20-Eve2C','2014-03-20-Eve2D',...
    '2014-03-20-Eve2E'}; 
NC=14;maxTDU=5;
NumNames=length(Names);
PostMat=zeros(NumNames,90);AntMat=zeros(NumNames,90);
for ii=1:NumNames
    Name = Names{ii}; 
     Directory = [Folder '/' Name];
     load(['\\Client\C$\Users\wtredman\Desktop\GregorLab\' Directory '\_data_NC' num2str(NC) '.mat']);
     
     X=zeros(length(data.ParticleCenters_x_APpos),1);
     for jj=1:length(data.ParticleCenters_x_APpos)
         frames=data.AllFrames{jj}{1};
         frames=1+frames-min(data.NC_Frames);
         frames(frames<1)=[]; 
         X(jj,frames)=data.ParticleCenters_x_APpos(jj);
     end
     
     PostMat(ii,1:length(max(X)))=max(X); %most posterior active nuclei
     X(X==0)=1;%this allows us to use the same array but now find the min
     AntMat(ii,1:length(max(X)))=min(X);
end
%%
AntMat(AntMat==1)=nan;AntMat(AntMat==0)=nan;%this gets rid of the bias from earlier
PostMat(PostMat==0)=nan;
AverageAnt=nanmean(AntMat);
AveragePost=nanmean(PostMat);

AverageAnt=AverageAnt(1:80);%After 80, more than half of the data sets are done
AveragePost=AveragePost(1:80);

figure
plot((1:length(AverageAnt))/80,AverageAnt,'bo'); hold on
plot((1:length(AveragePost))/80,AveragePost,'ro');

%% Smoothing
tol=0.00025*(1/80); %this was determined by playing around. Seemed to get a good result
AntcurveFit=spaps((1:80)/80,AverageAnt,tol);%this is where the smoother is applied, which makes the slm fitting easier and more accurate
PostcurveFit=spaps((1:80)/80,AveragePost,tol);
fnplt(spaps((1:80)/80,AverageAnt,tol),'b',2);
fnplt(spaps((1:80)/80,AveragePost,tol),'r',2);
slmAnt=slmengine((1:length(AntcurveFit.coefs))/(length(AntcurveFit.coefs),AntcurveFit.coefs,'plot','on','knots',40,'increasing',[0.012,0.037;0.126,0.48;0.521,1],'decreasing',[0,0.125;0.481,0.52]);
slmPost=slmengine((1:length(PostcurveFit.coefs))/length(PostcurveFit.coefs),PostcurveFit.coefs,'plot','on','knots',20,'increasing',[0,0.1625],'decreasing',[0.1626,1]);

AntVelocity=diff(slmAnt.y);%slmeval(slmAnt.y,slmAnt);
PostVelocity=diff(slmPost.y);%slmeval(slmPost.y,slmPost);

%% Plotting
%To make all the velocity less noisy, averaging every 4 data points
AntVel=zeros(1,20);PostVel=zeros(1,20);
every4=1:4:length(AntVelocity);
for ii=1:(length(AntVelocity)/4)-1
    AntVel(ii)=mean(AntVelocity(every4(ii):every4(ii+1)));
    PostVel(ii)=mean(PostVelocity(every4(ii):every4(ii+1)));
end
%Total time velocity plot
figure
plot(1:81,AntVelocity,'bo-'); hold on
plot(1:81,PostVelocity,'ro-');

%Every 4 averaged velocity
figure
plot(1:20,AntVel,'bo-');hold on
plot(1:20,PostVel,'ro-')
xlabel('Percent Time');
ylabel('Velocity (AP/Time)');


