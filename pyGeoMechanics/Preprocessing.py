import pandas as pd
import numpy as np

from scipy.ndimage import generic_filter, gaussian_filter1d
from scipy.signal import savgol_filter


class PreProcessor():

    def __init__(self,ingestion_key = None):
         self.dic_df = ingestion_key
         self.merged_df = self.data_prep(dic_df=self.dic_df)

    def is_numeric_log(self,series):
            return pd.api.types.is_numeric_dtype(series)

    def merge_logs(self,logs = None):
            
            # Sorting logs based on min depth
            list_sort = []
            col = logs[0].columns[1]
            # print(col)
            sample_size = round(logs[0]['DEPT'].iloc[1] - logs[0]['DEPT'].iloc[0],1)
            for count,idf in enumerate(logs):
                list_sort.append([count,idf['DEPT'].min(),idf['DEPT'].max()])
            sorted_list = sorted(list_sort, key=lambda x: x[1]) 
            # print("sorted_list:",sorted_list)

            # Consecutive merging
            for i in range(len(sorted_list)-1):

                first_log = logs[i]
                second_log = logs[i+1]
                result_log = None

                first_min_dept,first_max_dept = sorted_list[i][1],sorted_list[i][2]
                second_min_dept,second_max_dept = sorted_list[i+1][1],sorted_list[i+1][2]

                if second_min_dept>=first_min_dept and second_max_dept<=first_max_dept:
                    if (first_min_dept == second_min_dept) and (first_max_dept == second_max_dept):
                        not_null_1 = first_log[col].notna().sum()
                        not_null_2 = second_log[col].notna().sum()
                        if not_null_1 >= not_null_2:
                            result_log = first_log
                        else:
                            result_log = second_log
                    else:
                        result_log = first_log
                
                elif second_min_dept>=first_min_dept and second_max_dept>=first_max_dept:
                    result_log = pd.concat([first_log,second_log[second_log['DEPT']>=first_max_dept]])

                logs[i+1] = result_log
                logs[i+1].reset_index(inplace = True,drop =True)


            ov_df = logs[sorted_list[-1][0]]
            ov_df.replace(r'^\s*$', np.nan, regex=True, inplace=True) 
            ov_df = ov_df.sort_values('DEPT')
            ov_df = ov_df.astype(float)
            # ov_df[col][ov_df[col] < 0] = np.nan

            return ov_df

    def data_prep(self,track_columns=None, dic_df=None):

            min_depth_val = float('inf')
            sorted_df_list = []
            track_columns = set()
            common_depth = set()
            for i in dic_df:
                df = dic_df[i]
                if 'DEPT' in df.columns:
                    min_depth_val = df['DEPT'].iloc[0]
                    max_depth_val = df['DEPT'].iloc[-1]
                    sorted_df_list.append([i, min_depth_val, max_depth_val])
                    common_depth.update(df['DEPT'].values)
                    track_columns.update(list(df.columns))
            # print(len(common_depth))
            track_columns = list(track_columns)
            track_columns.remove('DEPT')
            sorted_list = sorted(sorted_df_list, key=lambda x: x[1])
            final_df = pd.DataFrame(columns=['DEPT'])
            final_df['DEPT'] = list(common_depth)
            final_df = final_df.sort_values(by='DEPT')
            # print(track_columns)
            for column in track_columns:
                logs = [] 
                for i in sorted_list:
                    df = dic_df[i[0]]

                    if column not in df.columns or not self.is_numeric_log(df[column]):
                        temp = df[['DEPT']].copy()
                        temp[column] = np.nan
                    else:
                        temp = df[['DEPT', column]].copy()
                    # print(temp[column])
                    if not temp[column].isna().all() or not temp[column].isnull().all():
                        logs.append(temp)
                # print(len(logs))
                if len(logs)>0:
                    ov_df = self.merge_logs(logs = logs)
                    final_df = pd.merge(final_df, ov_df, how='left', on='DEPT')
            for column in track_columns:
                final_df.rename(columns={column:column+'_merged'},inplace=True)
            final_df['DEPT'] = final_df['DEPT'].round(2)
            final_df.drop_duplicates(subset='DEPT',inplace=True)
            print("Successfully Merged all logs present in different LAS files")
            return final_df
    
    def smooth_log(self, data_frame = None, smooth_column = None, method='median', sample_size=11, ema_span=10, sg_polyorder=2, gaussian_sigma=2,depth_col='DEPT'):
        """
        Smooth a log using the specified method.
        
        Args:
            df: DataFrame with DEPT and data columns
            despiked_column: Name of the data column to smooth
            method: Smoothing method ('median', 'mean', 'ema', 'savgol', 'gaussian')
            sample_size: Window size for median, mean, and savgol (must be odd)
            ema_span: Span for exponential moving average
            sg_polyorder: Polynomial order for savgol filter
            gaussian_sigma: Standard deviation for gaussian filter
            depth_col: Name of depth column
            
        Returns:
            DataFrame with smoothed data
        """
        def nanmedian_filter(data, size):
            return generic_filter(data, lambda x: np.nanmedian(x), size=size, mode='nearest')
        
        def moving_average(data, size):
            return pd.Series(data).rolling(window=size, center=True, min_periods=1).mean().to_numpy()
        
        def exponential_moving_average(data, span):
            return pd.Series(data).ewm(span=span, adjust=False).mean().to_numpy()
        
        def savgol(data, size, polyorder):
            size = size if size % 2 != 0 else size + 1
            return savgol_filter(data, window_length=size, polyorder=polyorder)
        
        def gaussian(data, sigma):
            return gaussian_filter1d(data, sigma=sigma, mode='nearest')
        
        method_map = {
            'median': lambda x: nanmedian_filter(x, sample_size),
            'mean': lambda x: moving_average(x, sample_size),
            'ema': lambda x: exponential_moving_average(x, ema_span),
            'savgol': lambda x: savgol(x, sample_size, sg_polyorder),
            'gaussian': lambda x: gaussian(x, gaussian_sigma),
        }
        
        try:
            if method in method_map:
                result_df = data_frame.copy()
                result_df[smooth_column+'_smoothed'] = method_map[method](data_frame[smooth_column].to_numpy())
                print("Smoothing successfull")
                return result_df[[depth_col,smooth_column+'_smoothed']]
            else:
                raise ValueError(f"Unknown smoothing method: {method}")
        except Exception as e:
            print(f"Error during smoothing: {e}")
            return None
        
    def despike_log(self, data_frame = None, despiked_column = None, sigma=3.0, window_size=20, overlap_size=5, reduction_mode='remove', depth_col='DEPT'):
        """
        Despike a log using moving window approach with various reduction modes.
        
        Args:
            data_frame: DataFrame with DEPT and data columns
            despiked_column: Name of the data column to despike
            sigma: Standard deviation threshold
            window_size: Moving window size in depth units
            overlap_size: Window overlap in depth units
            reduction_mode: Type of despiking algorithm
            depth_col: Name of depth column
        
        Returns:
            Despiked DataFrame
        """
        if data_frame is None or data_frame.empty:
            return None
        
        result_df = data_frame.copy()
        result_df[f'{despiked_column}_despiked'] = data_frame[despiked_column].copy()
        
        if depth_col not in result_df.columns:
            raise ValueError(f"{depth_col} not in DataFrame")
        
        depths = result_df[depth_col].values
        min_depth = depths.min()
        max_depth = depths.max()
        
        window_starts = np.arange(min_depth, max_depth, max(window_size - overlap_size, 0.1))
        
        for start in window_starts:
            end = start + window_size
            window_mask = (result_df[depth_col] >= start) & (result_df[depth_col] < end)
            
            if not window_mask.any():
                continue
            
            window_series = result_df.loc[window_mask, f'{despiked_column}_despiked']
            window_data = window_series.dropna()
            
            if len(window_data) < 5:
                continue
            
            mean_val = window_data.mean()
            std_val = max(window_data.std(ddof=1), 1e-8)
            median_val = window_data.median()
            
            lower = mean_val - sigma * std_val
            upper = mean_val + sigma * std_val
            
            values = window_series.values.copy()
            spike_mask = (values < lower) | (values > upper)
            
            indices = result_df.index[window_mask]
            
            # Apply reduction mode
            if reduction_mode == 'clip':
                values = np.clip(values, lower, upper)
            
            elif reduction_mode == 'damp_half':
                values[spike_mask] = mean_val + 0.5 * (values[spike_mask] - mean_val)
            
            elif reduction_mode == 'damp_full':
                values[spike_mask] = mean_val
            
            elif reduction_mode == 'median_replace':
                values[spike_mask] = median_val
            
            elif reduction_mode == 'interp':
                temp = pd.Series(values, index=indices)
                temp[spike_mask] = np.nan
                temp = temp.interpolate(method='linear', limit_direction='both')
                values = temp.values
            
            elif reduction_mode == 'hampel':
                mad = np.median(np.abs(window_data - median_val))
                mad = max(mad, 1e-8)
                hampel_thresh = sigma * 1.4826 * mad
                hampel_mask = np.abs(values - median_val) > hampel_thresh
                values[hampel_mask] = median_val
            
            elif reduction_mode == 'physics_clip':
                values = np.clip(values, lower, upper)
            
            elif reduction_mode == 'remove':
                temp = pd.Series(values, index=indices,dtype=np.float32)
                temp[spike_mask] = np.nan
                values = temp.values
            
            result_df.loc[window_mask, f'{despiked_column}_despiked'] = values
        
        # Return only DEPT and despiked column
        print("Despiking Successfull")
        return result_df[[depth_col, f'{despiked_column}_despiked']].copy()
    
    def remove_outliers(self, data_frame = None, outlier_column = None, min_threshold=-np.inf, max_threshold=np.inf, depth_col='DEPT'):
        """
        Remove outliers from a log by setting values outside thresholds to NaN.
        
        Args:
            data_frame: DataFrame with DEPT and data columns
            outlier_column: Name of the data column
            min_threshold: Minimum allowable value
            max_threshold: Maximum allowable value
            depth_col: Name of depth column
        
        Returns:
            Cleaned DataFrame
        """
        if data_frame is None or data_frame.empty:
            return None
        
        result_df = data_frame.copy()
        
        # Create a copy of the data column
        cleaned_col = f'{outlier_column}_outlier'
        result_df[cleaned_col] = result_df[outlier_column].copy()
        
        # Apply thresholds - set outliers to NaN
        mask = (result_df[cleaned_col] < min_threshold) | (result_df[cleaned_col] > max_threshold)
        result_df.loc[mask, cleaned_col] = np.nan
        
        # Return only DEPT and cleaned column
        print("Outlier removal successfull")
        return result_df[[depth_col, cleaned_col]].copy()
    
    def interpolate(self,cleaned_df = None, track_columns = [], method = ''):
        
        methods = [
            'linear', 'time', 'index', 'values', 'nearest', 'zero',
            'spline', 'polynomial', 'pchip', 'akima',
            'barycentric', 'krogh', 'from_derivatives'
        ]
        if method in methods:
            try:
                for i in track_columns:
                    cleaned_df[i+'_interpolated'] = cleaned_df[i].interpolate(method=method, order=3,limit_direction = 'both')
                return cleaned_df
            except Exception as e:
                print(method,e)
        else:
            print(f"{method} not in available methods list",methods)