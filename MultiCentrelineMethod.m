function Results=MultiCentrelineMethod(centreline,termini,termini_date,plot_output,change_min,change_max)

% MultiCentrelineMethod
% Method takes initial centreline and creates equidistantly spaced
% centrelines either side at a given spacing. Terminus change along each of
% these centrelines is calculated, and a figure produced showing how
% termcentreline_nameinus position varies across the glacier front.

method_name='Multi-centreline method';
disp(strcat('Method: ',method_name))



%% Creates centrelines
box_precision=50;
centreline_spacing=10;
num_centrelines=2000;
disp('Processing:')
centreline.X(isnan(centreline.X))=[];
centreline.Y(isnan(centreline.Y))=[];


for n=1:round(num_centrelines/2,0)
    for m=1:length(centreline.X)-1
        ind=n;
        dist1=(centreline_spacing*round(num_centrelines/2,0))-(n-1)*centreline_spacing;
        angle=atan2d(centreline.Y(m)-centreline.Y(m+1),centreline.X(m)-centreline.X(m+1));
        x_diff=dist1*cosd(angle-90);
        y_diff=dist1*sind(angle-90);
        centrelinesX(ind,m)=centreline.X(m)+x_diff;
        centrelinesY(ind,m)=centreline.Y(m)+y_diff;
        dist2curve=nanmin(((centrelinesX(ind,m)-centreline.X).^2+(centrelinesY(ind,m)-centreline.Y).^2).^0.5);
        if dist2curve<dist1*0.99
            centrelinesX(ind,m)=NaN;
            centrelinesY(ind,m)=NaN;
        end
        if m==length(centreline.X)-1  && dist2curve>dist1*0.99
            centrelinesX(ind,m+1)=centreline.X(m+1)+x_diff;
            centrelinesY(ind,m+1)=centreline.Y(m+1)+y_diff;
        elseif m==length(centreline.X)-1 && dist2curve<dist1*0.99
            centrelinesX(ind,m+1)=NaN;
            centrelinesY(ind,m+1)=NaN;
        end
        
        
        ind=num_centrelines-(n-1)+1;
        angle=atan2d(centreline.Y(m)-centreline.Y(m+1),centreline.X(m)-centreline.X(m+1));
        x_diff=dist1*cosd(angle-90);
        y_diff=dist1*sind(angle-90);
        centrelinesX(ind,m)=centreline.X(m)-x_diff;
        centrelinesY(ind,m)=centreline.Y(m)-y_diff;
        dist2curve=nanmin(((centrelinesX(ind,m)-centreline.X).^2+(centrelinesY(ind,m)-centreline.Y).^2).^0.5);
        if dist2curve<dist1*0.99
            centrelinesX(ind,m)=NaN;
            centrelinesY(ind,m)=NaN;
        end
        if m==length(centreline.X)-1 && dist2curve>dist1*0.99
            centrelinesX(ind,m+1)=centreline.X(m+1)-x_diff;
            centrelinesY(ind,m+1)=centreline.Y(m+1)-y_diff;
        elseif m==length(centreline.X)-1 && dist2curve<dist1*0.99
            centrelinesX(ind,m+1)=NaN;
            centrelinesY(ind,m+1)=NaN;
        end
        
    end
    [~,dist2centreline]=distance2curve([centrelinesX(n,:),centrelinesY(n,:)],[centreline.X,centreline.Y]);
end
centrelinesX(round(num_centrelines/2,0)+1,:)=centreline.X;
centrelinesY(round(num_centrelines/2,0)+1,:)=centreline.Y;

centrelinesX(centrelinesX==0)=NaN;
centrelinesy(centrelinesY==0)=NaN;


disp('Centreline geometries calculated.')
%% calculates intercepts (or lack thereof) for each centreline
disp('Determining where terminus intercepts with centrelines...')
interceptX=nan(length(centrelinesX(:,1)),length(termini_date(:,1)));
interceptY=interceptX;

h=waitbar(0,'Calculating...')
for n=1:length(termini(:,1))
    
    term_intercept_dummyx=0;
    term_intercept_dummyx1=0;
    for m=1:round(num_centrelines/2,0)
        %to save calculation time, works from the middle of the terminus
        %outwards to edges of terminus
        ind=round(num_centrelines/2,0)+1-m;
        if ~isempty(term_intercept_dummyx)
            [term_intercept_dummyx,term_intercept_dummyy]=polyxpoly(centrelinesX(ind,:),...
                centrelinesY(ind,:),termini{n,1}.X,termini{n,1}.Y);
            if ~isempty(term_intercept_dummyx)
                interceptX(ind,n)=term_intercept_dummyx(1,1);
                interceptY(ind,n)=term_intercept_dummyy(1,1);
            end
        end
        if ~isempty(term_intercept_dummyx) && m==1  %gets centreline
            [term_intercept_dummyx,term_intercept_dummyy]=polyxpoly(centrelinesX(ind+1,:),...
                centrelinesY(ind+1,:),termini{n,1}.X,termini{n,1}.Y);
            if ~isempty(term_intercept_dummyx)
                interceptX(ind+1,n)=term_intercept_dummyx(1,1);
                interceptY(ind+1,n)=term_intercept_dummyy(1,1);
            end
        end
        ind=round(num_centrelines/2,0)+1+m;
        if ~isempty(term_intercept_dummyx1)
            [term_intercept_dummyx1,term_intercept_dummyy1]=polyxpoly(centrelinesX(ind,:),...
                centrelinesY(ind,:),termini{n,1}.X,termini{n,1}.Y);
            if ~isempty(term_intercept_dummyx1)
                interceptX(ind,n)=term_intercept_dummyx1(1,1);
                interceptY(ind,n)=term_intercept_dummyy1(1,1);
            end
        end             
    end
    waitbar(n/length(termini(:,1)))
end

centrelinesX(~any(~isnan(interceptX), 2),:) = [];    %removes centrelines where there are no intercepts
centrelinesY(~any(~isnan(interceptX), 2),:) = [];
interceptX(~any(~isnan(interceptX), 2),:) = [];  %removes rows where there are no intercepts
interceptY(~any(~isnan(interceptY), 2),:) = [];

%References each intercept to a date and distance _across_ terminus
[date_grid,dist_grid]=meshgrid(termini_date(:,4),[0:centreline_spacing:centreline_spacing*length(interceptX(:,1))-centreline_spacing]);

disp('Intercepts calculated.')

%% Finds common baseline for all centrelines
disp('Calculating distance of intercepts along centrelines')
for n=length(centrelinesX(1,:)):-1:1
    num_nans=sum(isnan(centrelinesX(:,n)));
    if num_nans>0
        centrelinesX(:,n)=[];
        centrelinesY(:,n)=[];
    end
end

%% calculates distance along centreline
for n=1:length(date_grid(1,:))
    for m=1:length(centrelinesX(:,1))
        cutoff=0;
        distance_raw(m,n)=0;
        for p=1:length(centrelinesX(1,:))-1
            if cutoff==0
                dist_to_point=((centrelinesX(m,p)-centrelinesX(m,p+1))^2+...
                    (centrelinesY(m,p)-centrelinesY(m,p+1))^2)^0.5;
                dist_to_intercept=((centrelinesX(m,p)-interceptX(m,n))^2+...
                    (centrelinesY(m,p)-interceptY(m,n))^2)^0.5;
                if dist_to_point<dist_to_intercept
                    distance_raw(m,n)=distance_raw(m,n)+dist_to_point;
                elseif dist_to_point>dist_to_intercept
                    distance_raw(m,n)=distance_raw(m,n)+dist_to_intercept;
                    cutoff=1;
                end
            end
        end
    end
end
distance_raw(distance_raw==0)=NaN;
Results1.Method=method_name;
Results1.DateAll=termini_date;
Results1.Date=termini_date;
Results1.DistanceRaw=distance_raw;
Results1.Distance=distance_raw-nanmin(distance_raw(:));
Results1.DistanceFullRes=Results1.Distance;
Results1.DistanceChange=nan(size(Results1.Distance));
Results1.DistanceChange(:,2:end)=diff(Results1.Distance,1,2);
Results1.RateChange=nan(size(Results1.Distance));
Results1.RateChange(:,2:end)=Results1.DistanceChange(:,2:end)./(diff(date_grid,1,2)./365);
Results1.Distance1D=nanmean(Results1.Distance,1)';
Results1.Distance1DFullRes=Results1.Distance1D;
Results1.DistanceChange1D=diff(Results1.Distance1D);
Results1.RateChange1D=Results1.DistanceChange1D./(diff(termini_date(:,4))./365);
Results1.DistAcrossAll(:,1)=0:centreline_spacing:(length(Results1.Distance(:,1))*centreline_spacing)-centreline_spacing;
Results1.DistAcross(:,1)=0:centreline_spacing:(length(Results1.Distance(:,1))*centreline_spacing)-centreline_spacing;

if change_min~=0
    a=0;
    ind_check=1;
    for n=1:length(Results1.Date(:,1))-1
        if n>=ind_check
            check=0;
            for m=n+1:length(Results1.Date(:,1))
                if check==0
                   if Results1.Date(m,4)-Results1.Date(n,4)>=change_min &&...
                            Results1.Date(m,4)-Results1.Date(n,4)<=change_max
                        a=a+1;
                        Results.Date(a:a+1,:)=[Results1.Date(n,:);Results1.Date(m,:)];
                        Results.DistanceRaw(:,a:a+1)=[Results1.DistanceRaw(:,n),Results1.DistanceRaw(:,m)];
                        ind_check=m;
                        check=1;
                    elseif Results1.Date(m,4)-Results1.Date(n,4)>change_max
                        a=a+1;
                        Results.Date(a,:)=[Results1.Date(n,1:3),Results1.Date(n,4)+0.0001];
                        Results.DistanceRaw(:,a)=nan(length(Results1.DistanceRaw(:,1)),1);
                        ind_check=m;
                        check=1;
                    end
                end
            end
        end
    end
    [date_grid1,dist_grid1]=meshgrid(Results.Date(:,4),...
        [0:centreline_spacing:centreline_spacing*length(interceptX(:,1))-centreline_spacing]);
    Results.Method=Results1.Method;
    Results.DateAll=Results1.Date;
    Results.Distance=Results.DistanceRaw-nanmin(Results.DistanceRaw(:));
    Results.DistanceFullRes=Results1.Distance;
    Results.DistanceChange=nan(size(Results.Distance));
    Results.DistanceChange(:,2:end)=diff(Results.Distance,1,2);
    Results.RateChange=nan(size(Results.Distance));
    Results.RateChange(:,2:end)=Results.DistanceChange(:,2:end)./(diff(date_grid1,1,2)./365);
    Results.Distance1D=nanmean(Results.Distance,1)';
    Results.Distance1DFullRes=Results1.Distance1D;
    Results.DistanceChange1D=diff(Results.Distance1D);
    Results.RateChange1D=Results.DistanceChange1D./(diff(Results.Date(:,4))./365);
    Results.DistAcrossAll=Results1.DistAcross;
    Results.DistAcross(:,1)=0:centreline_spacing:(length(Results.Distance(:,1))*centreline_spacing)-centreline_spacing;
    
    
%     for n=1:length(Results1.Date(:,4))
%         date_diff(:,n)=Results1.Date(:,4)-Results1.Date(n,4);
%         date_diff(n,n)=NaN; %removes comparison with the same observation
%     end
%     date_diff(date_diff<0)=NaN; %removes negative values from matrix
%     date_diff(date_diff<change_min)=-1; %identifies if there are subsequent obs that occur before date+change_min
%     date_diff(date_diff>change_min&date_diff<change_max)=-2;    %identifies those in obs window
%     date_diff(date_diff>0)=NaN;
%     a=1;
%     
%     prev_ind=1
%     for n=1:length(date_diff(1,:))-1
%         holder_ind=1;
% 
%         if sum(~isnan(date_diff(:,n)))>0
%             holder_ind=nanmin(find(date_diff(:,n)==-2))
%             if ~isempty(holder_ind) && n>=prev_ind
%                 Results.DateDistance(a:a+2,:)=[Results1.Date(n,:);Results1.Date(holder_ind,:);[Results1.Date(holder_ind,1:3),Results1.Date(holder_ind,4)+0.0001]];
%                 Results.DistanceRaw(:,a:a+2)=[Results1.DistanceRaw(:,n),Results1.DistanceRaw(:,holder_ind),nan(length(Results1.DistanceRaw(:,n)),1)];
% 
%                 a=a+3;
%                 prev_ind=holder_ind;
%             end
%         end
%         
%     end
%     Results.DistanceRaw(Results.DistanceRaw==0)=NaN;
%     Results.Distance=Results.DistanceRaw-nanmin(Results.DistanceRaw(:));
%     
%     a=0;
%     for n=1:length(Results.Distance(1,:))
%         if sum(~isnan(Results.Distance(:,n)))>0
%             a=a+1;
%             Results.DateChange(a,:)=Results.DateDistance(n,:);
%             Results.DistanceNoGap(:,a)=Results.Distance(:,n);
%         end
%     end
%     
%     a=0;
%     DistanceChange=nan(size(Results.Distance));
%     RateChange=nan(size(Results.Distance));
%     DistanceChange=diff(Results.DistanceNoGap,1,2);
%     DistanceChange(:,end)=NaN;
%     a=1;
%     a=0;
%     for n=1:length(DistanceChange(1,:))-1
%         a=a+1
%         n
%         Results.DateChange(a,:)=Results.DateDistance(n,:);
%         Results.DistanceChange(:,a)=DistanceChange(:,n);
%         
%         if sum(DistanceChange(:,n+1)~=0&~isnan(DistanceChange(:,n+1)))>5
%             a=a+1;
%             Results.DateChange(a,:)=[Results.DateDistance(n,1:3),Results.DateDistance(n,4)+0.0001];
%             Results.DistanceChange(:,a)=NaN;
%         end
%     end
    
elseif change_min==0
    date_grid1=date_grid;
    dist_grid1=dist_grid;
    Results=Results1;
end
       
   close(h)             
%% Plotting
if plot_output==1
    custom_colormap
    colormap(cmap);
    %Gets annual labels
    ticks=[datenum(Results.Date(1,1),1,1):365:datenum(Results.Date(end,1),1,1)+365.25];
    for n=1:length(ticks)
        leap=leapyear(year(ticks(n)));
        if leap==1
            ticks(n+1:end)=ticks(n+1:end)+1;
        end
    end
    
    disp('Plotting output...')
    figure(1)
    %Absolute location of terminus at full resolution
    subplot(4,3,1:2)
    pcolor(date_grid,dist_grid,Results.DistanceFullRes)
    shading flat
    xticks(ticks)
    datetick('x','yy','keepticks')
    xlabel('Year')
    ylabel('Distance across terminus (m)')
    colorbar()
%     xlim([datenum(2013,1,1),datenum(2018,1,1)])
    grid on
    grid minor
    
    %Absolute location at defined interval
    subplot(4,3,4:5)
    pcolor(date_grid1,dist_grid1,Results.Distance)
    shading flat
    xticks(ticks)
    datetick('x','yy','keepticks')
    xlabel('Year')
    ylabel('Distance across terminus (m)')
    colorbar()
%     xlim([datenum(2013,1,1),datenum(2018,1,1)])
    grid on
    grid minor
    
    %Distance change
    subplot(4,3,7:8)
    pcolor(date_grid1,dist_grid1,Results.DistanceChange);shading flat
    shading flat
    xticks(ticks)
    datetick('x','yy','keepticks')
    xlabel('Year')
    ylabel('Distance across terminus (m)')
    colorbar()
%     xlim([datenum(2013,1,1),datenum(2018,1,1)])
    grid on
    grid minor
    
    %Rate change

    subplot(4,3,10:11)
    pcolor(date_grid1,dist_grid1,Results.RateChange)
    shading flat
    xticks(ticks)
    datetick('x','yy','keepticks')
    
    xlabel('Year')
    ylabel('Distance across terminus (m)')
    colorbar()
%     xlim([datenum(2013,1,1),datenum(2018,1,1)])
    grid on
    grid minor
    
    %Terminus change 1D full res
    subplot(4,3,3)
    plot(date_grid(1,:),Results.Distance1DFullRes,'-xk')
    xticks(ticks)
    datetick('x','yy','keepticks')
    xlabel('Year')
    ylabel('Terminus position (m)')
%     xlim([datenum(2013,1,1),datenum(2018,1,1)])
    grid on
    grid minor
       
    %Terminus change 1D
    subplot(4,3,6)
    plot(Results.Date(:,4),Results.Distance1D,'-xk')
    xticks(ticks)
    datetick('x','yy','keepticks')
    xlabel('Year')
    ylabel('Terminus position (m)')
%     xlim([datenum(2013,1,1),datenum(2018,1,1)])
    grid on
    grid minor
    
    %Distance change 1D
    subplot(4,3,9)
    plot(Results.Date(2:end,4),Results.DistanceChange1D,'-xk')
    xticks(ticks)
    datetick('x','yy','keepticks')
    xlabel('Year')
    ylabel('Distance change from previous obs (m)')
%     xlim([datenum(2013,1,1),datenum(2018,1,1)])
    grid on
    grid minor
    
    %Rate change 1D
    subplot(4,3,12)
    plot(Results.Date(2:end,4),Results.RateChange1D,'-xk')
    xticks(ticks)
    datetick('x','yy','keepticks')
    xlabel('Year')
    ylabel('Rate of change (m yr^-^1)')
%     xlim([datenum(2013,1,1),datenum(2018,1,1)])
    grid on
    grid minor
    
    
end
disp('Done!')

end