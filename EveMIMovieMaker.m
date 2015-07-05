%% Eve MI Movie Maker
%   This is part of the time analysis that we are going to do to make the
%   Eve stuff more concrete. Also, to see if we can reconcile its time
%   scale with that of Hunchback. 
%%
percenttime=0.02; %time is the percent of total time we want to look at 
M = struct('cdata',[],'colormap',[]);
for ii=1:1/percenttime
    time=percenttime*ii;
    [F]=MutualInfo_v3(time);
    M(ii).cdata=F.cdata;
    M(ii).colormap=F.colormap;
end
movie(M,1)