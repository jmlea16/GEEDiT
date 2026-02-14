function Results=ExtrapCentreline(centreline,termini,termini_date,exponent)

% ExtrapCentreline
% 'Centreline' = centreline shapefile data that have already been read in to
% Matlab; 'termini' = terminus data that have already been read into Matlab;
% 'termini_date' = data extracted from termini shapefiles by readGlacierData
% function; 'exponent' is the power value in the inverse distance weighting
% calculation used to determine how the centreline is extrapolated.
% Function calculates the distance along a centreline a terminus is located.
% Achieves this by taking points located at 1m intervals along a terminus, 
% and for each point extrapolating the centreline distance value to that
% point. The overall terminus position is determined by taking the average
% value of all the terminus points.

    method_name='Extrapolated centreline method';

    if ~exist('exponent','var')
        exponent=2;
    end
    %interpolates values along centreline every 1m...
    %Step 1: calculate distancebetween nodes
    node_dist=[];
    num_points=[];
    xy_1m=[];
    for n=1:length(centreline.X)-3  %-3 is due to last value being NaN and the value before being same as the value previous to that
        node_dist(n,1)=((centreline.Y(n)-centreline.Y(n+1))^2+(centreline.X(n)-centreline.X(n+1))^2)^0.5;
    end
    if sum(node_dist<1)>=1
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        disp(strcat('WARNING: ', num2str(sum(node_dist<1)),'_of_',num2str(length(node_dist)),...
            ' nodes are <1 m apart. If this looks like too many,'))
        disp('        (1) check if your data are in UTM format, and (2) if they are,')
        disp('        consider reducing the number of nodes on your centreline')
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    end
    %Step 2: calculate coordinates and total distance along centreline at 
    %1m intervals
    for n=1:length(centreline.X)-3
        num_points=ceil(node_dist(n,1));   %number of points to be included
        xy=[linspace(centreline.X(n),centreline.X(n+1),num_points).'...
            linspace(centreline.Y(n),centreline.Y(n+1),num_points).'];
        if n==1
            xy_1m=xy;
        else
            xy_1m=[xy_1m;xy];
        end
        xy=[];
    end
    distance_calc=@(x,y) (x.^2+y.^2).^0.5;
    distance_centreline(:,1)=distance_calc(diff(xy_1m(:,1)),diff(xy_1m(:,2)));
    xy_1m(1,3)=0;
    xy_1m(2:end,3)=cumsum(distance_centreline);
    distance_centreline=[];
    
    %The above can be used for all the terminus observations, and therefore
    %does not have to be repeated within the script. The following steps
    %are therefore applied individually to each glacier terminus.
    for n=1:length(termini(:,1))
        n
        T_xy_1m=[];
        distance_tk=[];
    %Step 3: Calculate locations  every 1 m along a given terminus polyline.
        node_dist=[];
        xy=[];
        T_xy_1m=[];
        %Calculate distances between nodes on the terminus
        for m=1:length(termini{n}.X)-3  %-3 is due to last value being NaN and the value before being same as the value previous to that
            node_dist(m,1)=((termini{n}.Y(m)-termini{n}.Y(m+1))^2+(termini{n}.X(m)-termini{n}.X(m+1))^2)^0.5;
        end
        terminus_length(n,1)=sum(node_dist(:,1));
        %Calculate points every 1m between nodes of terminus polyline
        for m=1:length(termini{n}.X)-3
            num_points=ceil(node_dist(m,1));   %number of points to be included
            xy=[linspace(termini{n}.X(m),termini{n}.X(m+1),num_points).'...
                linspace(termini{n}.Y(m),termini{n}.Y(m+1),num_points).'];
            if m==1
                T_xy_1m=xy;
            end
            T_xy_1m=[T_xy_1m;xy];
        end
        %Equations 1-3 in Lea et al., 2014
        for m=1:length(T_xy_1m(:,1))
            distance_tk=[];
            %calculate distance from every point on centreline to t_k
            distance_tk(:,1)=((xy_1m(:,1)-T_xy_1m(m,1)).^2+(xy_1m(:,2)-T_xy_1m(m,2)).^2).^0.5;
            %calculate weighting (eqn. 1)
            weighting(:,1)=1./(distance_tk(:,1).^exponent);
            %normalise weighting (eqn. 2)
            weighting_norm(:,1)=weighting(:,1)./sum(weighting(:,1));
            %calculate position of t_k relative to the centreline
            position_tk(m,1)=sum(xy_1m(:,3).*weighting_norm(:,1));
        end
     %Step 4: Calculate overall terminus position (eqn. 5)
        Results.Method=method_name;
        Results.Date(n,:)=termini_date(n,1:end);
        %calculate straight line width of terminus
        Results.TerminusWidth(n,1)=((termini{n}.X(1)-termini{n}.X(end-1))^2+...
            (termini{n}.Y(1)-termini{n}.Y(end-1))^2)^0.5;
        %calculate the path length along the glacier terminus.
        Results.TerminusPathLength(n,1)=terminus_length(n,1);
        Results.RawDistance(n,1)=sum(position_tk(:,1))/length(position_tk);
        Results.TerminusGeometry{n,1}=[termini{n}.X;termini{n}.Y];
    end
    
    Results.Centreline(:,1)=centreline.X;
    Results.Centreline(:,2)=centreline.Y;
    Results.Distance=Results.RawDistance-Results.RawDistance(end,1);
    distance(1,2:3)=NaN;
    Results.TerminusChange(1,1)=nan;
    Results.RateChange(n,1)=nan;    
    for n=2:length(Results.Distance(:,1))
        %calculate terminus change
        Results.TerminusChange(n,1)=Results.Distance(n,1)-Results.Distance(n-1,1);
        %calculate rate of terminus change in m/yr
        Results.RateChange(n,1)=Results.TerminusChange(n,1)/...
            ((termini_date(n,4)-termini_date(n-1,4))/365);
    end
    