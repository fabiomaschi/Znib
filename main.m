%%% Znib Project
%%% Fabio Maschi
%%% 3 October 2018
%%% Version 0.5

clear all; close all;

%% Parameters
K_FILENAME = 'Picture2';
REF_LINES = ['A','B','B','F','H'];
REF_COLUM = [ 1 , 2 , 5 , 1 , 2 ];
B_INTERATIVE = false; % SELECT REFERENCE HALOS
%
K_SEGMENTATION = 1.5; % between 1.3 and 2
K_ROUND_TOLERANCE = 0.2;
K_CLOSE_KERNEL_SIZE = 2; % N*2 + 1
B_WITH_CLOSING = true; % CLOSE MORPHOLOGIC OPERATION
B_WITH_SAVING = true; % EXPORT INTERMEDIATE IMGS
%
K_NUCLEUS_SIZE = 100; % in pixels
K_HIT_SIZE = 10; % in pixels

%% Open image
img = double(rgb2gray(imread(strcat(K_FILENAME,'','.png'))));
txt = fopen(strcat(K_FILENAME,'','-result.txt'),'w');
[w,h] = size(img);

%% Segmentation
% histogram
hstgrm = zeros(1,255);
for x=1:w
    for y=1:h
        hstgrm(img(x,y)+1) = hstgrm(img(x,y)+1)+1;
    end
end
hstgrm = hstgrm / (w*h);

% % omega(i) = \sum_1^i { hstgrm(k) }
omega_0 = cumsum(hstgrm);
omega_1 = 1 - omega_0;

% % mu(i) = \sum_1^i { hstgrm(k)*i }
mu = cumsum(hstgrm .* (1:255));

mu_0 = mu ./ omega_0;
mu_1 = (mu(255) - mu) ./ (omega_1);

var2b = omega_0.*omega_1.*(mu_1 - mu_0).^K_SEGMENTATION;

[~, K] = max(var2b);
img_segmented = logical(img > K);
img = img/255;

if B_WITH_SAVING == true
    imwrite(~img_segmented, 'step_a-segmentation.jpg');
end
clear K mu mu_0 mu_1 var2b omega_0 omega_1 hstgrm;

%% Dilate and Erode (Close operation)
if B_WITH_CLOSING == true
    img_connected = single(~img_segmented);
    kernel = strel('disk',K_CLOSE_KERNEL_SIZE,0);
    kernel = single(kernel.getnhood());
    img_connected = conv2(img_connected, kernel, 'same');
    img_connected = single(~img_connected);
    img_connected = conv2(img_connected, kernel, 'same');
    img_connected = logical(img_connected);
else
    img_connected = logical(img_segmented);
end

if B_WITH_SAVING == true
    imwrite(~img_connected, 'step_b-close.jpg');
end
clear kernel;

%% Fill nucleus holes
CC = bwconncomp(img_connected);
numPixels = cellfun(@numel, CC.PixelIdxList);

for i=1:size(numPixels,2)
    if numPixels(i) <= K_NUCLEUS_SIZE
        img_connected(CC.PixelIdxList{i}) = 0;
    end
end
img_connected = ~img_connected;

if B_WITH_SAVING == true
    imwrite(img_connected, 'step_c-fill.jpg');
end
clear i numPixels CC;


%% ROUND CELLS
CC = bwconncomp(img_connected);
stats = regionprops(CC, 'Area','Perimeter');
pixels = regionprops(CC, 'PixelIdxList');
area_perim = [stats.Area; stats.Perimeter];
C = (4*pi*area_perim(1,:))./(area_perim(2,:).^2);
cellule_ronde = false(w,h);

for i = 1:CC.NumObjects 
   if (C(i) > (1-K_ROUND_TOLERANCE) && C(i) < (1+K_ROUND_TOLERANCE))
       cellule_ronde(CC.PixelIdxList{i}) = true;
   end
end

if B_WITH_SAVING == true
    mask = single(cat(3,img+img_connected,img+cellule_ronde,img));
    imwrite(mask, 'step_d-round.jpg');
end
clear C area_perim mask;

%% Grid detection
st = regionprops(cellule_ronde, 'Centroid');
c = vertcat(st.Centroid);
[~, x] = kmeans(c(:,1), 12);
[~, y] = kmeans(c(:,2), 8);
x = sort(x);
y = sort(y);

%% Select gridded cells
grid_cells = single(zeros(8,12));
img_stdenz = false(w, h);
for i = 1:length(x)
    mx = ceil(x(i)-K_HIT_SIZE:x(i)+K_HIT_SIZE);
    mx = min(h, max(1,mx));
    for j = 1:length(y)
        my = ceil(y(j)-K_HIT_SIZE:y(j)+K_HIT_SIZE);
        my = min(w, max(1,my));
        xy = sub2ind([w,h], my, mx);
        for k = 1:CC.NumObjects
            if(any(ismember(xy,pixels(k).PixelIdxList)))
                img_stdenz(CC.PixelIdxList{k}) = 1;
                grid_cells(j,i) = k;
                break;
            end
        end
    end
end

if B_WITH_SAVING == true
    mask = single(cat(3,img+img_connected,img+img_stdenz,img));
    mask(:,ceil(x),:) = 0;
    mask(:,ceil(x),2) = 1;
    mask(ceil(y),:,:) = 0;
    mask(ceil(y),:,2) = 1;
    imwrite(mask, 'step_e-grid.jpg');
	%line(c(:,1), c(:,2), 'LineStyle','none', 'Marker','+', 'Color','b')
end
clear st i j k my mx xy x y cellule_ronde mask;

%% Select References
ref_area = zeros(1,length(REF_LINES));
aux = 0;
img_refenz = false(w, h);

if B_INTERATIVE == true
    %# INTERATIVE MODE
    formatSpec = 'Ref #%d = %4d pixels\n';
    while(aux ~= N_REFERENCES)
        mask = single(cat(3,img+img_refenz,img+img_stdenz,img+img_stdenz));
        imshow(mask);
        [x, y] = ginput(1);
        xy = sub2ind([w,h], ceil(y), ceil(x));
        for j = 1:CC.NumObjects
            if(ismember(xy,pixels(j).PixelIdxList))
                img_stdenz(CC.PixelIdxList{j}) = 0;
                img_refenz(CC.PixelIdxList{j}) = 1;
                aux = aux + 1;
                ref_area(aux) = stats(j).Area;
                fprintf(txt,formatSpec, aux, stats(j).Area);
                break;
            end
        end
    end
else
    %# AUTOMATED
    formatSpec = 'Ref #%d %c%2d = %4d pixels\n';
    for i = 1:length(REF_LINES)
        letter = single(REF_LINES(i))-64;
        j = grid_cells(letter,REF_COLUM(i));
        if (j ~= 0)
            img_stdenz(CC.PixelIdxList{j}) = 0;
            img_refenz(CC.PixelIdxList{j}) = 1;
            aux = aux + 1;
            ref_area(aux) = stats(j).Area;
            fprintf(txt,formatSpec, aux, REF_LINES(i), REF_COLUM(i), stats(j).Area);
        end
    end
end
ref_mean = mean(ref_area);
ref_stdd = std(ref_area);
if B_WITH_SAVING == true
    mask = single(cat(3,img+img_refenz,img+img_stdenz,img+img_stdenz));
    imwrite(mask, 'step_f-ref.jpg');
end
clear c i j aux mask letter pixels ref_area;


%% Compute area
fprintf(txt,'\nAverage area = %4.2f\n', ref_mean);
fprintf(txt,'Std. deviation = %4.2f\n\n', ref_stdd);
formatSpec = '#%c%2d = %4d (%5.2f%%)\n';
str_bigger = '\nBigger halos:';
for y = 1:size(grid_cells,1)
    letter = char(y+64);
    for x = 1:size(grid_cells,2)
        aux_area = stats(grid_cells(y,x)).Area;
        if(aux_area < ref_mean)
            img_stdenz(CC.PixelIdxList{grid_cells(y,x)}) = 0;
        else
            str_bigger = strcat(str_bigger,sprintf(' #%c%d',letter,x));
        end
        fprintf(txt,formatSpec,letter,x, aux_area, aux_area/ref_mean*100);
    end
end
fprintf(txt,str_bigger);
fclose(txt);

if B_WITH_SAVING == true
    mask = single(cat(3,img+img_refenz,img+img_stdenz,img+img_stdenz));
    imwrite(mask, 'step_g-result.jpg');
end
mask = single(cat(3,img,img+img_stdenz,img+img_stdenz));
imwrite(mask, strcat(K_FILENAME,'','-result.jpg'));

clear i x y formatSpec str_bigger mask letter stats CC aux_area;