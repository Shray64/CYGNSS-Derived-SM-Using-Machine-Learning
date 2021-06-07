% Need to make 2 changes in this file before running
% first specify month number - line 24
% Second - specify month name - line 261 (last line)

clear all; close all; clc

flag_level = 3; 

%text_out_version = 'TxSON version 6_2, compiled 09/19/19'; %header line info on .txt output files
text_out_version = 'TxSON version 6_3, compiled 06/24/20'; %header line info on .txt output files

% load data from import_data_TxSON.m
load('TxSONv6_3.mat')

% [n_loc, n_cells] = size(cell_ind); %n_cells = # of grid cells of n_loc >=3
n_loc = 40
% n_cells = 16

grid_corners_9km = readtable('cygnss/fishnet_9km.csv')
n_cells = size(grid_corners_9km, 1)

opt.flag = flag_level;
opt.header = text_out_version;
date = datetime(Date)
month_ind = (date.Month == 11) %Have to change month number here
year_ind = (date.Year == 19)
ind_date = (month_ind & year_ind)

num_hours = sum(ind_date)

addpath('miscScripts') % fold of various scripts from Matlab Central, required for processing 

%%   Set QA level from QA_flag [tS(:), site(40), depth(4), flag_level(4)]
for i = 1:n_loc % site loop
    for j = 1:4 % depth loop
        qa_flag = squeeze(QA_flag(:,i,j,1:flag_level));
        QA_ind(:,i,j) = any(qa_flag,2);
    end
end
    
SWC(QA_ind) = NaN; % replace any flagged data with NaN
T(QA_ind) = NaN;    

clearvars qa_flag QA* i j 

% Get SWC values for required time period
SWC = SWC(ind_date,:,:)

%% ++++ ARITHMETIC MEAN ++++++++++++++++++++++++++++++++++++++++++++++++++
%  - for each grid cell indices, applying QA flags across all 4 soil depth (5, 10, 20 and 50 cm)
SWCa = NaN(num_hours,n_cells, 4);
SWCa_std = SWCa;  SWC_n = SWCa;
Ta = SWCa; Ta_std = SWCa;
PPTa = NaN(num_hours,n_cells); PPTa_std = PPTa;
% 
% for i = 1:n_cells % grid cell loop
%     swc = NaN(length(tS),4); % initial matrices
%     swc_std = swc; swc_n = swc;
%     ta = swc; ta_std = swc;
%     
%     for j = 1:4 % soil depth loop
%         swc(:,j) = nanmean(SWC(:,cell_ind(:,i),j),2);
%         swc_std(:,j) = nanstd(SWC(:,cell_ind(:,i),j)')';
%         swc_n(:,j) = nansum(~isnan(SWC(:,cell_ind(:,i),j)),2);
%         
%         ta(:,j) = nanmean(T(:,cell_ind(:,i),j),2);
%         ta_std(:,j) = nanstd(T(:,cell_ind(:,i),j)')';       
%     end
%     PPTa(:,i) = nanmean(ppt(:,cell_ind(:,i)),2);
%     PPTa_std(:,i) = nanstd(ppt(:,cell_ind(:,i))')';
%     
%     %export_TxSON(tS, PPTa(:,i), PPTa_std(:,i), swc_n, swc, swc_std, cell_ID{i},'SWC_A', opt)
%     
%     SWCa(:,i,:) = swc;
%     SWCa_std(:,i,:) = swc_std;
%     SWC_n(:,i,:) = swc_n; % number/scale is not changing between weighing functions
%     Ta(:,i,:) = ta;
%     Ta_std(:,i,:) = ta_std;
% end
% parfor i = 1:n_cells
%     export_TxSON(tS, PPTa(:,i), PPTa_std(:,i), SWC_n(:,i,:), SWCa(:,i,:), SWCa_std(:,i,:), cell_ID{i},'SWC_A', opt)
% end
% save('upscale_TxSON.mat', 'SWCa', 'SWCa_std', 'SWC_n','Ta','Ta_std','PPTa', 'PPTa_std')

%% +++++ INVERSE DISTANCE WEIGHING ++++++++++++++++++++++++++++++++++++++++
grid_interval = 100; %IDW interpolated grid size in meters 
r1 = 'fr'; %'fr' = fixed radius ;  'ng' = neighbours
r2 = 1500; %radius if r1 == 'fr' / number of neighbours == 'ng'
pow = 2; %IDW power
utm_z = '14 N'; %UTM zone for reprojections 

% create polygons (WGS) for each grid cell in cell_ind
% bbox = table2array(EASE2(:,6:9)); % corners

bbox = table2array(grid_corners_9km);

for i=1:n_cells
    poly_lon{i}(:)=[bbox(i,3) bbox(i,3) bbox(i,4) bbox(i,4) bbox(i,3)];
    poly_lat{i}(:)=[bbox(i,2) bbox(i,1) bbox(i,1) bbox(i,2) bbox(i,2)];
end

% Build grids for rasterized IDW
lat_range = [EASE2.LAT_corner_min(1) EASE2.LAT_corner_max(1)]; % corners of 36 km EASE2 grid in WGS
lon_range = [EASE2.LON_corner_min(1) EASE2.LON_corner_max(1)];

[x_range, y_range, utm_z] = wgs2utm(lat_range, lon_range);
x = [min(x_range):grid_interval:max(x_range)];
y = [min(y_range):grid_interval:max(y_range)];

[x_grid, y_grid] = meshgrid(x,y);
[size_x, size_y] = size(x_grid);

x_grid2 = reshape(x_grid, [],1);
y_grid2 = reshape(y_grid, [],1);
utm_z2 = sprintf('%s N', num2str(utm_z(1)));
utm_z_g = repmat(utm_z2, [length(x_grid2), 1]);

[LATg, LONg] = utm2deg(x_grid2,y_grid2, utm_z_g);

LATg = reshape(LATg,[size_x size_y]);
LONg = reshape(LONg,[size_x size_y]);

[xc, yc] = wgs2utm(LAT,LON); % convert station coordinates to UTM 

%% IDW time loop
%initialize matrices
SWCi = NaN(num_hours,n_cells, 4);
SWCi_std = SWCi;
PPTi = NaN(num_hours,n_cells); PPTi_std = PPTi;

% tic
% h1 = parfor_progressbar(744,'Computing IDW');
% 
% parfor i = 1:744 % i = time loop
% %     disp(i)
%     h1.iterate(1); %only use with parfor 
%     swci2 = zeros(4,n_cells);
%     swci2_std = swci2;
%     for j = 1:4 % j = depth loop
%         data = squeeze(SWC(i,:,j)); % send all data to produce gridded 36 km IDW cell per hour/depth
%         ind = ~isnan(data); % nan's need to be removed prior to IDW or it'll preferentially create a circle of them.
%         
%         [swci] = IDW(xc(ind,1), yc(ind,1), data(ind)', x,y, pow,r1,r2);
%         
%         for k = 1:n_cells % k = grid cell loop; parsing 36 km IDW into cells
%             in = inpolygon(LONg, LATg, poly_lon{k}, poly_lat{k});
%             swci2(j,k) = nanmean(swci(in),'all');
%             swci2_std(j,k) = nanstd(swci(in),0,'all');
%         end
%     end
%     
%     SWCi(i,:,:) = swci2';
%     SWCi_std(i,:,:) = swci2_std';
%     
%     %Precipitation loop
%     ppti2 = zeros(1,n_cells);
%     ppti2_std = ppti2;
%     
%     data = squeeze(ppt(i,:));
%     ind = ~isnan(data);
%     
%     [ppti] = IDW(xc(ind,1), yc(ind,1), data(ind)', x,y, pow,r1,r2);
%     for k = 1:n_cells % k = cell loop  
%         in = inpolygon(LONg, LATg, poly_lon{k}, poly_lat{k});
%         ppti2(k) = nanmean(ppti(in), 'all');
%         ppti2_std(k) = nanstd(ppti(in), 0, 'all');
%     end   
%     PPTi(i,:) = ppti2;
%     PPTi_std(i,:) = ppti2_std;  
% end
% toc
% close(h1)
% save('TxSON_upscaled/txson_jan_idw_3km.mat')

% parfor i = 1:n_cells
%     export_TxSON(tS, PPTi(:,i), PPTi_std(:,i), SWC_n(:,i,:), SWCi(:,i,:), SWCi_std(:,i,:), cell_ID{i},'SWC_I', opt)
% end
    
% save('upscale_TxSON.mat', 'SWCi', 'SWCi_std', 'PPTi', 'PPTi_std', '-append')

%% +++++ VORONOI ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Epsilon = 1E-8; %parameter for voronoi function in WGS
% Create a bounding box around EASE2 center of TxSON grid
dist_bound = 40000; %distance of larger boundary area in m
grid_res = .001; %grid resolution in degrees
hemi='N'; %N for northern hemisphere

% Build a grid for rasterized VORONOI index 1 is the 36 km pixel
[x1,y1,utm_z]=wgs2utm(EASE2.LAT_center(1),EASE2.LON_center(1));
% find the max and min utm coordinates of a dist_bound grid
x1max = x1+(dist_bound/2);  y1max = y1+(dist_bound/2);
x1min = x1-(dist_bound/2);  y1min = y1-(dist_bound/2);
%convert back to wgs
[Latmax,Lonmax]=utm2deg(x1max,y1max,sprintf('%s', [num2str(utm_z(1)),' ',hemi]));
[Latmin,Lonmin]=utm2deg(x1min,y1min,sprintf('%s', [num2str(utm_z(1)),' ',hemi]));

%create a grid of the dist_bound with a resolution of grid_res
long=Lonmin:grid_res:Lonmax;
latg=Latmin:grid_res:Latmax;
[LONgv, LATgv]=meshgrid(long, latg);

%% VORONOI time loop
%initialize matrices

SWCv= NaN(num_hours,n_cells, 4);
SWCv_std = SWCv;  SWCv_n = SWCv;
PPTv = NaN(num_hours,n_cells); 
PPTv_std = PPTv;  PPTv_n = PPTv;

h1 = parfor_progressbar(num_hours,'Computing VORONOI');
tic
parfor i=1:num_hours % i = time loop
    h1.iterate(1); %only use with parfor
    swcv2 = zeros(4,n_cells);
    swcv2_std = swcv2; swcv2_n = swcv2;
    for j=1:4 % j = depth loop
        data = squeeze(SWC(i,:,j));% send all data to produce gridded Voronoi station index per hour/depth
        [swcv] = voronoi_func_v4(LAT,LON,data',LATgv,LONgv,Epsilon);
        for k = 1:n_cells % k = grid cell loop; parsing Voronoi into cells
            disp(k)
            in = inpolygon(LONgv, LATgv, poly_lon{k}, poly_lat{k});
            station_unique=unique(swcv(in));
            weight=[];
            for w=1:length(station_unique)
                % weights of stations in the subsetted area
                weight(w)=length(find((swcv(in))==station_unique(w)))/numel(swcv(in));
            end
            swcv2_n(j,k)=length(station_unique);
            swcv2(j,k)=sum(data(station_unique).*weight);
            swcv2_std(j,k)=sqrt(var(data(station_unique),weight));
        end
    end
    
    SWCv_n(i,:,:) = swcv2_n';
    SWCv(i,:,:) = swcv2';
    SWCv_std(i,:,:) = swcv2_std';

    %Precipitation loop    
%     pptv2 = zeros(1,n_cells);
%     pptv2_std = pptv2;  pptv2_n = pptv2; 
%     data = squeeze(ppt(i,:)); % send all data to produce gridded Voronoi station index per hour/depth
%     [pptv] = voronoi_func_v4(LAT,LON,data',LATgv,LONgv,Epsilon);
%     for k = 1:n_cells % k = grid cell loop; parsing Voronoi into cells
%         in = inpolygon(LONgv, LATgv, poly_lon{k}, poly_lat{k});
%         station_unique=unique(pptv(in));
%         weight=[];
%         for w=1:length(station_unique)
%             % weights of stations in the subsetted area
%             weight(w)=length(find((pptv(in))==station_unique(w)))/numel(pptv(in));
%         end
%         pptv2_n(k)=length(station_unique);
%         pptv2(k)=sum(data(station_unique).*weight);
%         pptv2_std(k)=sqrt(var(data(station_unique),weight));
%     end       
%     PPTv_n(i,:) = pptv2_n;
%     PPTv(i,:) = pptv2;
%     PPTv_std(i,:) = pptv2_std;  
end
toc
close(h1)
save('TxSON_upscaled/txson_nov_voronoi_9km.mat') % Have to change month name here
