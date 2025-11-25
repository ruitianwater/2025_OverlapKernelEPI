function [eigen_maps,eigen_vals] = ESPIRiT_map(lowRes_dat,nk,ny,nx,thres_k)

%% 
[ny_low, nx_low, nc] = size(lowRes_dat); 
ns_y = ny_low-nk;
ns_x = nx_low-nk;

ny_zp = ceil((ny-nk)/2);
nx_zp = ceil((nx-nk)/2);

y_zp_odd = ceil((ny-nk)/2) ~= (ny-nk)/2;
x_zp_odd = ceil((nx-nk)/2) ~= (nx-nk)/2;

cali_matx = zeros(ns_y*ns_x,nk*nk*nc);

%% Formalize the calibration matrix

for index_coil = 1:nc
    
    x_cali = (1:nk*nk) + (index_coil-1)*nk*nk;
    
    for index_shift_x = 1:ns_x
        
        x_lowRes = (1:nk) + index_shift_x-1;
        
        for index_shift_y = 1:ns_y
            
            y_cali = index_shift_x + (index_shift_y-1)*ns_x;
            
            y_lowRes = (1:nk) + index_shift_y-1;
            
            dat_temp = lowRes_dat(y_lowRes,x_lowRes,index_coil);
            cali_matx(y_cali,x_cali) = reshape(dat_temp,[1 nk*nk]);
            
        end
        
    end
    
end

%% SVD of the calibration matrix

[U_cali,S_cali,V_cali] = svd(cali_matx,'econ');

cali_val = diag(S_cali);

id_trun = max(find(cali_val>=cali_val(1)*thres_k));

cali_filters = reshape(V_cali,[nk nk nc nk*nk*nc]);
cali_filters = cali_filters(:,:,:,1:id_trun);

cali_filters = permute(cali_filters,[1,2,4,3]); 
cali_filters = reshape(cali_filters,nk*nk*id_trun,nc);

[U_filter,S_filter,V_filter] = svd(cali_filters,'econ');
cali_filters = cali_filters*V_filter;

cali_filters = reshape(cali_filters,[nk nk id_trun nc]);
cali_filters = permute(cali_filters,[1,2,4,3]);

%% image domain operation

img_kernel = zeros(ny,nx,nc,id_trun);

for index_eigen = 1:id_trun
    
    cali_filters_temp = cali_filters(:,:,:,index_eigen);
    
    cali_filters_temp = conj(flip(flip(cali_filters_temp,1),2));
    
    cali_filters_temp_zp = padarray(cali_filters_temp,[ny_zp nx_zp],0,'both');
    
    if y_zp_odd == 1    
        cali_filters_temp_zp = cali_filters_temp_zp(1:end-1,:,:);    
    end 
    if x_zp_odd == 1 
        cali_filters_temp_zp = cali_filters_temp_zp(:,1:end-1,:);
    end
    
    img_kernel(:,:,:,index_eigen) = SymIfft(cali_filters_temp_zp,[1 2]);
        
end



%% SVD in the image domain

eigen_maps = zeros(ny,nx,nc,min(id_trun,nc));
eigen_vals = zeros(ny,nx,min(id_trun,nc));

for id_y = 1:ny
    for id_x =  1:nx
        
        img_kernel_yx = squeeze(img_kernel(id_y,id_x,:,:));

        [U_img,S_img,V_img] = svd(img_kernel_yx,'econ');
        
        ph_ref_rad = angle(U_img(1,:));
        ph_comp = repmat(exp(-1i*ph_ref_rad),[nc 1]);

        U_img = V_filter*(U_img.*ph_comp);

        S_img = diag(S_img);
        
        eigen_maps(id_y,id_x,:,:) = U_img;
        eigen_vals(id_y,id_x,:) = S_img;
        
    end
end



end

