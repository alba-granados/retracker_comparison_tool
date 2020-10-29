function [out_data,indx]=outliers_by_percentiles(y,outlier_percentil_low,outlier_percentil_high,IQR_times)
IQR=iqr(y);
percentiles = prctile(y,[outlier_percentil_low outlier_percentil_high]);
indx=(y<(percentiles(1)-IQR_times*IQR) | y>(percentiles(2)+IQR_times*IQR));
%force the outliers as missing data or NaN
out_data=y;
out_data(indx)=NaN;

end