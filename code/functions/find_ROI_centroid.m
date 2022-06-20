function [centroids, scaled_centroids, FOV_size] = ...
    find_ROI_centroid(mask, info)


s = regionprops(mask, 'centroid');
centroids = cat(1, s.Centroid);

ind = info.config.magnification;
x_scale = info.calibration(ind).x;
y_scale = info.calibration(ind).y;

FOV_size = [x_scale * 796, y_scale * 512];

scaled_centroids(:, 1) = centroids(:, 1) * x_scale;
scaled_centroids(:, 2) = centroids(:, 2) * y_scale;

end