% Merge GeoJSON files, convert to format readable by terminus toolbox and
% convert to shapefile

clearvars termini merged

[filenames,folder_path]=uigetfile({'*.GeoJSON;*.JSON','JavaScript objects (*.GeoJSON,*.JSON)';...
    '*.shp','Shapefile (*.shp)'},'Select file(s) to convert/merge','MultiSelect','on');

merged=[];
if length(string(filenames))==1
    filenames={filenames};    
end
[~,~,suffix]=fileparts(filenames{1});
if strcmp(suffix,'.geojson')||strcmp(suffix,'.json')
    for n=1:length(filenames)
        input_file=fileread(strcat(folder_path,filenames{n}));
        input=jsondecode(input_file);
        merged=[merged;input];
    end

    a=0;
    b=0;
    for n=1:length(merged)
        if isfield(merged(n).features(1).properties,'satellite')
            for m=1:length(merged(n).features)
                if strcmp(merged(n).features(m).geometry.type,'LineString')
                    a=a+1;
                    margin(a).Name=merged(n).features(m).properties.Name;
                    margin(a).Date=merged(n).features(m).properties.date;
                    margin(a).Unclear=merged(n).features(m).properties.unclear;
                    margin(a).Description='EarthEngine';
                    margin(a).Satellite=merged(n).features(m).properties.satellite;
                    margin(a).Asc_Desc=merged(n).features(m).properties.Asc_Desc; 
                    margin(a).ImagePath=merged(n).features(m).properties.image_path;
                    margin(a).Notes=merged(n).features(m).properties.notes;  
                    margin(a).X=merged(n).features(m).geometry.coordinates(:,1);
                    margin(a).Y=merged(n).features(m).geometry.coordinates(:,2);
                    margin(a).Geometry='Line';
                elseif strcmp(merged(n).features(m).geometry.type,'Polygon')
                    b=b+1
                    margin1(b).Name=merged(n).features(m).properties.Name;
                    margin1(b).Date=merged(n).features(m).properties.date;
                    margin1(b).Unclear=merged(n).features(m).properties.unclear;
                    margin1(b).Description='EarthEngine';
                    margin1(b).Satellite=merged(n).features(m).properties.satellite;
                    margin1(b).Asc_Desc=merged(n).features(m).properties.Asc_Desc; 
                    margin1(b).ImagePath=merged(n).features(m).properties.image_path;
                    margin1(b).Notes=merged(n).features(m).properties.notes;  
                    margin1(b).X=merged(n).features(m).geometry.coordinates(:,1);
                    margin1(b).Y=merged(n).features(m).geometry.coordinates(:,2);
                    margin1(b).Geometry='Polygon';
                end
            end
        else
            margin2.Name='Name';
            margin2.X=merged.features.geometry.coordinates(:,1);
            margin2.Y=merged.features.geometry.coordinates(:,2);
            margin2.Geometry='Line';
        end
    end
elseif strcmp(suffix,'.shp')
    for n=1:length(filenames)
        merged{n}=shaperead(strcat(folder_path,filenames{n}));
    end
    a=0;
    for n=1:length(merged)
        for m=1:length(merged{n})
            a=a+1;
            margin(a).Name='Name';
            margin(a).Date=filenames{n}(1:10);
            margin(a).Description='Merged shapefile';
            margin(a).X=merged{n}.X;
            margin(a).Y=merged{n}.Y;
            margin(a).Geometry='Line';
        end
    end
end

[save_filename,save_folder_path]=uiputfile({'*.shp','Shapefile (*.shp)'},'Save shapefile as:');
if exist('margin')==1
    shapewrite(margin,strcat(save_folder_path,'Line_',save_filename));
end
if exist('margin1')==1
    shapewrite(margin1,strcat(save_folder_path,'Polygon_',save_filename));
end
if exist('margin2')==1
    shapewrite(margin2,strcat(save_folder_path,save_filename));
end

msgbox('Files successfully converted/merged. The single shapefile generated can now be used in GlaQiT.')