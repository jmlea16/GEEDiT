function Results=VariableBox(centreline,termini,termini_date)

% VariableBox
% 'Centreline' = centreline shapefile data that have already been read in to
% Matlab; 'termini' = terminus data that have already been read into Matlab;
% 'termini_date' = data extracted from termini shapefiles by readGlacierData
% function.
% A box is defined around the centreline, with its width determined by the 
% straight line width of the glacier terminus. The terminus is used to
% split this box in two, with distance from the end of the centreline given
% by the upstream box area divided by the box width.

method_name='Variable Box Method';
disp(strcat('Method: ',method_name))
    
%% Creates boxes for each terminus

box_precision=50;
disp('Processing:')
h=waitbar(0,'Calculating')
for n=1:length(termini)
   disp(strcat('Terminus date:',num2str(termini_date(n,3)),'/',num2str(termini_date(n,2)),...
       '/',num2str(termini_date(n,1)),'; Observation:',num2str(n),'/',num2str(length(termini))))
if length(termini{n}.X)>2
% make points at [box_precision] m intervals along centreline
a=0;
for m=2:length(centreline.X)
    num_points=(((centreline.Y(m)-centreline.Y(m-1))^2+(centreline.X(m)-...
        centreline.X(m-1))^2)^0.5)/box_precision;   %linear distance between centreline nodes div by box_precision
    node_distances(m-1,1)=((centreline.Y(m)-centreline.Y(m-1))^2+(centreline.X(m)-...
        centreline.X(m-1))^2)^0.5;  %generates vector of distances between centreline nodes (go to line 39)
    for p=1:floor(num_points)
        a=a+1;
        if p==1
            centrelineX(a,1)=centreline.X(m-1); %gets node coordinate of initial centreline
            centrelineY(a,1)=centreline.Y(m-1);
        else
            centrelineX(a,1)=centrelineX(a-1)+((centreline.X(m)-...
                centreline.X(m-1))/num_points); %gets remainder of points by adding increments
            centrelineY(a,1)=centrelineY(a-1)+((centreline.Y(m)-...
                centreline.Y(m-1))/num_points);
        end
    end
end

% Gives warning to simplify centreline if mean distance between centreline
% nodes is < the precision of hte box
mean_node_dist=mean(node_distances);
if mean_node_dist<box_precision
    disp('WARNING: If analysis is slow, recommend simplifying centreline (i.e. fewer nodes)')
end

%% converts xy centreline coordinates to sn coordinates, where n=0 is equiv to the centreline
%Define distances from centreline and determine which way (i.e. NtoS/StoN
%etc) the terminus has been digitised
[s_coord,n_coord]=xy2sn(centrelineX,centrelineY,centreline.X(1:end),centreline.Y(1:end));
[s_term,n_term]=xy2sn(termini{n}.X,termini{n}.Y,centreline.X,centreline.Y);
radiusUp=abs(nanmax(n_term));   %gets distances from centreline to box in sn space
radiusDown=abs(nanmin(n_term));
n_coordUp=n_coord+radiusUp;   %gets parallel lines delineating edge of box in sn space
n_coordDown=n_coord-radiusDown;

[x_coordUp,y_coordUp]=sn2xy(s_coord(2:end),n_coordUp(2:end),centreline.X(1:end),centreline.Y(1:end)); %converts back from sn coords to xy coords
[x_coordDown,y_coordDown]=sn2xy(s_coord(2:end),n_coordDown(2:end),centreline.X(1:end),centreline.Y(1:end));

% %Calculates upstream corners of box
coeffs=polyfit(centreline.X(1:2),centreline.Y(1:2),1);
gradient=-1/coeffs(1);
icept=centreline.Y(1)-gradient*centreline.X(1);

angle1=atan2d(gradient,1);
deltaXUp=radiusUp*sind(angle1);
deltaYUp=radiusUp*cosd(angle1);
deltaXDown=radiusDown*sind(angle1);
deltaYDown=radiusDown*cosd(angle1);
% 
x_coordUp=[centreline.X(1)-deltaXUp;centreline.X(1)+deltaXUp;x_coordUp];    %append both possible points
y_coordUp=[centreline.Y(1)-deltaYUp;centreline.Y(1)+deltaYUp;y_coordUp];
x_coordDown=[centreline.X(1)-deltaXDown;centreline.X(1)+deltaXDown;x_coordDown];    %append both possible points
y_coordDown=[centreline.Y(1)-deltaYDown;centreline.Y(1)+deltaYDown;y_coordDown];

%gets distance of box edges to centreline, and removes all points where
%<radius
[~,dist2centrelineUp]=distance2curve([centrelineX,centrelineY],[x_coordUp(1:end),y_coordUp(1:end)]);
[~,dist2centrelineDown]=distance2curve([centrelineX,centrelineY],[x_coordDown(1:end),y_coordDown(1:end)]);

x_coordUp(dist2centrelineUp<radiusUp.*0.9995|dist2centrelineUp>radiusUp.*1.0005)=[];
y_coordUp(dist2centrelineUp<radiusUp.*0.9995|dist2centrelineUp>radiusUp.*1.0005)=[];
x_coordDown(dist2centrelineDown<radiusDown.*0.9995|dist2centrelineDown>radiusDown.*1.0005)=[];
y_coordDown(dist2centrelineDown<radiusDown.*0.9995|dist2centrelineDown>radiusDown.*1.0005)=[];



%checks for superfluous points where the angle between successive
%points is 0 deg. This helps speed up calculation of polyxpoly
angle1=[];
angle2=[];
angle1(:,1)=atand(diff(y_coordUp)./diff(x_coordUp));
angle2(:,1)=atand(diff(y_coordDown)./diff(x_coordDown));

x_coordUp=x_coordUp([find(diff(angle1(:))~=0)+1;length(x_coordUp)]);
y_coordUp=y_coordUp([find(diff(angle1(:))~=0)+1;length(y_coordUp)]);
x_coordDown=x_coordDown([find(diff(angle2(:))~=0)+1;length(x_coordDown)]);
y_coordDown=y_coordDown([find(diff(angle2(:))~=0)+1;length(y_coordDown)]);

%constructs the box from bottom left clockwise to bottom right
bufferxy=[[x_coordUp,y_coordUp];flipud([x_coordDown,y_coordDown])];
%gets rid of NaNs
bufferxy(~any(~isnan(bufferxy),2),:)=[];
%closes box at end
bufferxy=[bufferxy;bufferxy(1,1:2)];

    %Calculates distance between nodes on terminus for use in
    %determining number of points between nodes
    for m=1:length(termini{n}.X)-1  %-1 is to prevent loop going out of bounds
        node_dist(m,1)=((termini{n}.Y(m)-termini{n}.Y(m+1))^2+(termini{n}.X(m)-termini{n}.X(m+1))^2)^0.5;
    end
    %Calculate points (nominally) every 10m between nodes of terminus
    %polyline. In reality, nodes will be spaced marginally less than 10m
    %apart as number of points required is rounded up to the nearest
    %whole number.
    for m=1:length(termini{n}.X)-1
        num_points=ceil(node_dist(m,1))/10;   %number of points to be included
        xy=[linspace(termini{n}.X(m),termini{n}.X(m+1),num_points).'...
            linspace(termini{n}.Y(m),termini{n}.Y(m+1),num_points).'];
        if m==1
            T_xy_1m=xy;
        end
        T_xy_1m=[T_xy_1m;xy];
    end
    %determine which points lie within the box. These are the ones
    %that will be used to construct the upper bound of the box.
    intercept=[];
    term_poly_logical=inpolygon(T_xy_1m(:,1),T_xy_1m(:,2),bufferxy(:,1),bufferxy(:,2));
    term_poly=T_xy_1m(find(term_poly_logical==1),1:2);
    [intercept(:,1),intercept(:,2)]=polyxpoly(T_xy_1m(:,1),T_xy_1m(:,2),bufferxy(:,1),bufferxy(:,2));
    
    %as is unknown whether terminus has been digitised N to S/S to N/W to
    %E/E to W etc, need to get terminus coordinates in correct order
    %otherwise box geometry will be incorrect
    [~,check]=distance2curve([x_coordUp,y_coordUp],[term_poly(end,1),term_poly(end,2)]);
    if check>radiusUp
        term_poly=flipud(term_poly);
        T_xy_1m=flipud(T_xy_1m);
    end
    
    %check if both sides of the terminus intersect with the box. If
    %not, then the shortest distance between a terminus endpoint and
    %the initial box will be used to 'close' the box.
    int_check1=[];  %preallocate/wipe values from previous terminus
    int_check2=[];
    %find terminus index value where centreline intersects it
    [~,centreline_dist_to_term]=distance2curve([centreline.X',centreline.Y'],T_xy_1m);
    [~,ind_split]=nanmin(centreline_dist_to_term);
    [~,int_check1(1,:)]=polyxpoly(T_xy_1m(1:ind_split,1),...
        T_xy_1m(1:ind_split,2),bufferxy(:,1),bufferxy(:,2)); %1st half of terminus
    [~,int_check2(1,:)]=polyxpoly(T_xy_1m(ind_split:end,1),...
        T_xy_1m(ind_split:end,2),bufferxy(:,1),bufferxy(:,2)); %1st half of terminus

    intercept_extra=[];
    if sum(int_check1(1,:))==0 || sum(int_check2(1,:))==0
        if sum(int_check1(1,:))==0
            %find distance to box node that end of terminus is closest to
            [intercept_extra(1,1:2),extrap_dist(1,1)]=distance2curve(bufferxy(:,1:2),T_xy_1m(1,1:2));
            intercept=[intercept;intercept_extra(1,:)];
        end
        if sum(int_check2(1,:))==0
            [intercept_extra(2,1:2),extrap_dist(2,1)]=distance2curve(bufferxy(:,1:2),T_xy_1m(end,1:2));
            intercept=[intercept_extra(2,:);intercept];
        end
    end

    %check if there are multiple intercepts that are a tiny distance
    %from one another and remove
    for m=1:length(intercept(:,1))-1
        if ((intercept(m,1)-intercept(m+1,1))^2+(intercept(m,2)-intercept(m+1,2))^2)^0.5<20
            intercept(m,:)=NaN;
        end
    end
    intercept(~any(~isnan(intercept),2),:)=[]; %cleans up the NaNs

    %append intercepts at start and end of the terminus, checking at
    %which end each should go
    perp_dist(1,1)=((intercept(1,1)-term_poly(1,1))^2+(intercept(1,2)...
        -term_poly(1,2))^2)^0.5;
    perp_dist(1,2)=((intercept(1,1)-term_poly(end,1))^2+(intercept(1,2)...
        -term_poly(end,2))^2)^0.5;
    perp_dist(2,1)=((intercept(end,1)-term_poly(1,1))^2+(intercept(end,2)...
        -term_poly(1,2))^2)^0.5;
    perp_dist(2,2)=((intercept(end,1)-term_poly(end,1))^2+(intercept(end,2)...
        -term_poly(end,2))^2)^0.5;
    if perp_dist(1,1)>perp_dist(1,2)
        term_poly=[intercept(1,1:2);flipud(term_poly);intercept(end,1:2)];
    elseif perp_dist(2,1)>perp_dist(2,2)
        term_poly=[intercept(1,1:2);term_poly;intercept(end,1:2)];
    end

    %Remove part of box that does not cover where the glacier is. Achieves
    %this by checking where the intersects are located on the box and
    %only taking the nodes that are "upstream" of the terminus. The
    %terminus points are also appended.

   a=0;
   ind1=0;
   ind2=0;
   for m=1:length(bufferxy(:,1))-1
       if ind1==0
          [~,position_intersect(m,1)]=inpolygon(intercept(1,1),intercept(1,2),...
              bufferxy(m:m+1,1),bufferxy(m:m+1,2));
          if position_intersect(m,1)==1
              ind1=m;
          end
       end
       if ind2==0
          [~,position_intersect(m,2)]=inpolygon(intercept(end,1),intercept(end,2),...
          bufferxy(m:m+1,1),bufferxy(m:m+1,2));
          if position_intersect(m,2)==1
              ind2=m+1;
          end
       end
   end

    box_final=[bufferxy(1:ind1,1:2);term_poly;bufferxy(ind2:end,1:2)];

    for m=1:length(termini{n}.X)-1  %-3 is due to last value being NaN and the value before being same as the value previous to that
        node_dist(m,1)=((termini{n}.Y(m)-termini{n}.Y(m+1))^2+(termini{n}.X(m)-termini{n}.X(m+1))^2)^0.5;
    end
    terminus_length(n,1)=sum(node_dist(:,1));

       Results.Method=method_name;
        Results.Date(n,:)=termini_date(n,1:end);
        Results.BoxWidth(n,1)=radiusUp+radiusDown;
        %calculate straight line width of terminus
        Results.TerminusWidth(n,1)=((termini{n}.X(1)-termini{n}.X(end-1))^2+...
            (termini{n}.Y(1)-termini{n}.Y(end-1))^2)^0.5;
        %calculate the path length along the glacier terminus.
        Results.TerminusPathLength(n,1)=terminus_length(n,1);
        Results.BoxArea(n,1)=polyarea(box_final(:,1),box_final(:,2));
        Results.BoxGeometry{n,1}=box_final;
        Results.RawDistance(n,1)=Results.BoxArea(n,1)/(radiusUp+radiusDown);
        Results.TerminusGeometry{n,1}=[termini{n}.X';termini{n}.Y'];
        [Results.CentrelineCut{n,1}(:,1),Results.CentrelineCut{n,1}(:,2)]=...
            polyxpoly(centreline.X,centreline.Y,termini{n}.X,termini{n}.Y);
        Results.Distance=Results.RawDistance-Results.RawDistance(end,1);
        Results.TerminusDetail(n,1)=terminus_length(n,1)/length(Results.TerminusGeometry{n,1}(1,:));
end

Results.Centreline(:,1)=centreline.X;
Results.Centreline(:,2)=centreline.Y;

distance(1,2:3)=NaN;
Results.TerminusChange(1,1)=nan;
Results.RateChange(1,1)=nan;    

waitbar(n/length(termini))
end
for m=2:length(Results.Distance(:,1))
    %calculate terminus change
    Results.TerminusChange(m,1)=Results.Distance(m,1)-Results.Distance(m-1,1);
    %calculate rate of terminus change in m/yr
    Results.RateChange(m,1)=Results.TerminusChange(m,1)/...
        ((Results.Date(m,4)-Results.Date(m-1,4))/365);
end

Results.BoxWidth(Results.BoxArea==0)=NaN;
%calculate straight line width of terminus
Results.TerminusWidth(Results.BoxArea==0)=NaN;
%calculate the path length along the glacier terminus.
Results.TerminusPathLength(Results.BoxArea==0)=NaN;
Results.RawDistance(Results.BoxArea==0)=NaN;


Results.Distance(Results.BoxArea==0)=NaN;
Results.Date(Results.BoxArea==0,:)=NaN;
Results.BoxArea(Results.BoxArea==0)=NaN;
    close(h)
end