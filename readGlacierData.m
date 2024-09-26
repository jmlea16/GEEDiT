function [centreline,termini,termini_date] = readGlacierData(centreline_directory,...
    centreline_name,terminus_directory,terminus_filename,terminus_format)

% readGlacierData: Function takes centreline and terminus shapefile data,
% imports it into Matlab, and also extracts date information from termini.
% Tool also automatically sorts termini into chronological order from the
% earliest observation to the most recent.
% If using terminus_format = 1, terminus_filename can be replaced by []
    disp('Loading terminus data')
    if terminus_format==1
        terminus_filename=[];
    end

    %   Read centreline
    centreline=shaperead(strcat(centreline_directory,centreline_name));
    
    %gets rid of NaNs (if they exist)
    centreline.X(isnan(centreline.X))=[];
    centreline.Y(isnan(centreline.Y))=[];
    
    if centreline.X(1)>-180 && centreline.X(1)<180 && centreline.Y(1)>-90 ...
        && centreline.Y(1)<90
        [centreline.X,centreline.Y]=deg2utm(centreline.Y,centreline.X);
    end
       
    %gets rid of duplicate coordinates (if any)
    centreline1=unique([centreline.X;centreline.Y],'stable');
    
    centreline.X=[];
    centreline.Y=[];
    centreline.X=centreline1(1:end/2,1)';
    centreline.Y=centreline1((end/2)+1:end,1)';
    
    %   Read terminus data
    if terminus_format==1
        filenames=dir(strcat(terminus_directory,'*.shp'));
        termini=[];
        m=0;
        for n=1:length(filenames(:,1))
%             if ~isempty()
%             ind=strfind(filenames(),',');
%             end
            if ~isempty(str2num(filenames(n,1).name(2:4)))  %2:4 as avoids error if centreline in same directory AND is called 'flowline....shp' (will erroneously call matlab's 'flow' function)
                m=m+1;
                filenames(n,1).name
                termini{m,1}=shaperead(strcat(terminus_directory,filenames(n,1).name));
                termini_date(m,1)=str2num(filenames(n,1).name(1:4)); %year
                termini_date(m,2)=str2num(filenames(n,1).name(6:7)); %month
                termini_date(m,3)=str2num(filenames(n,1).name(9:10)); %day
                termini_date(m,4)=datenum(termini_date(m,1),termini_date(m,2),termini_date(m,3)); %serial date
                
                %gets rid of NaNs (if they exist)
                termini{m,1}.X(isnan(termini{m,1}.X))=[];
                termini{m,1}.Y(isnan(termini{m,1}.Y))=[];
                
                if termini{m,1}.X(1)>-180 && termini{m,1}.X(1)<180 && termini{m,1}.Y(1)>-90 ...
                    && termini{m,1}.Y(1)<90
                    [termini{m,1}.X,termini{m,1}.Y]=deg2utm(termini{m,1}.Y,termini{m,1}.X);
                end
                                
                %gets rid of duplicate coordinates (if any)
                terminus1=[];
                holder=[termini{m,1}.X,termini{m,1}.Y];
                terminus1=unique(holder,'rows','stable')';
                
                termini{m,1}.X=[];
                termini{m,1}.Y=[];
                
                termini{m,1}.X=terminus1(1,:)';
                termini{m,1}.Y=terminus1(2,:)';
                
            end
        end
    elseif terminus_format==2
        
        termini_dummy=[];
        termini=[];
        termini_dummy=shaperead(strcat(terminus_directory,terminus_filename));  %dummy variable to initially hold data
        [~,index] = sortrows({termini_dummy.Date}.'); termini_dummy = termini_dummy(index); clear index %sorts structure by date
        
        h=waitbar(0,'Loading shapefile data...')
        for n=1:length(termini_dummy)
            if exist('termini_dummy(n).Description')
                if strcmp(termini_dummy(n).Description,'EarthEngine')
                    termini_dummy(n).X=termini_dummy(n).Lon;
                    termini_dummy(n).Y=termini_dummy(n).Lat;
                end
            end
        %for n=1:length(termini_dummy(:,1))
            termini{n,1}=termini_dummy(n,1);
            termini_date(n,1)=str2num(termini_dummy(n,1).Date(1:4)); %year
            termini_date(n,2)=str2num(termini_dummy(n,1).Date(6:7)); %month
            termini_date(n,3)=str2num(termini_dummy(n,1).Date(9:10)); %day
            termini_date(n,4)=datenum(termini_date(n,1),termini_date(n,2),termini_date(n,3)); %serial date
            
            %gets rid of NaNs (if they exist)
            termini{n,1}.X(isnan(termini_dummy(n).X))=[];
            termini{n,1}.Y(isnan(termini_dummy(n).Y))=[];
            
            if termini{n,1}.X(1)>-180 && termini{n,1}.X(1)<180 && termini{n,1}.Y(1)>-90 ...
                && termini{n,1}.Y(1)<90
                [termini{n,1}.X,termini{n,1}.Y]=deg2utm(termini{n,1}.Y,termini{n,1}.X);
            end
            
            %gets rid of duplicate coordinates (if any)
                terminus1=[];
                holder=[termini{n,1}.X,termini{n,1}.Y];
                terminus1=unique(holder,'rows','stable')';
                
                termini{n,1}.X=[];
                termini{n,1}.Y=[];
                
                termini{n,1}.X=terminus1(1,:)';
                termini{n,1}.Y=terminus1(2,:)';
                
        
      %  end
    %For Google EarthEngine shapefiles
    if exist('termini_dummy(n).Description')
        if strcmp(termini_dummy(n).Description,'EarthEngine')
            termini{n}.ImagePath=termini_dummy(n).ImagePath;
            termini{n}.Satellite=termini_dummy(n).Satellite;
            termini{n}.Unclear=termini_dummy(n).Unclear;
        end
    end
    waitbar(n/length(termini_dummy))
        end
        termini
        close(h)
end