import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyGeoMechanics.Utils as ut

from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from scipy.ndimage import generic_filter, gaussian_filter1d

from matplotlib.ticker import LogLocator, ScalarFormatter


class PorePressurePredictor:
    def __init__(self, df: pd.DataFrame, 
                 depth_col='DEPT',
                 overburden_df = None ,
                 stratigraphy_df = None, 
                 strat_method = '3',
                 nct_type = None,

            ):
        # self.use_tvd = use_tvd
        self.ovd_df = overburden_df
        self.start_df =  stratigraphy_df
        self.final_df = df
        self.depth_col = depth_col 
        self.nct_type = nct_type

        # For nct calculation
        if strat_method == '3':
            self.shale_df = self.start_df[self.start_df['strat_code_method3'] == 1].reset_index(drop=True)
            self.shale_df.rename(columns = {'strat_code_method3':'strat_code'},inplace=True)
        elif strat_method == '2':
            self.shale_df = self.start_df[self.start_df['strat_code_method2'] == 1].reset_index(drop=True)
            self.shale_df.rename(columns = {'strat_code_method2':'strat_code'},inplace=True)
        else:
            self.shale_df = self.start_df[self.start_df['strat_code_method1'] == 1].reset_index(drop=True)
            self.shale_df.rename(columns = {'strat_code_method1':'strat_code'},inplace=True)

        self.generate_nct()

        self.calculate_pore_pressure()

        self.calculate_fracture_gradient()

        self.plot_pore_pressure()
        
        
    def pp_issue_tracker(self,df:pd.DataFrame):
        if len(df[(df.pore_pressure != df.p_hydro) & (df.pore_pressure.isna() == False)]) == 0:
            print("At all depth Pore pressure is Equal to Hydrostatic pressure")
            print("REASON:")
            if df.nct_val.isna().all():
                print("All nct values are Nan")
            else:
                total_nct_count = len(df[(df.nct_val.isna() == False)])
                sv_psi_na_count = len(df[(df.nct_val.isna() == False) & (df.SV_psi.isna() == True)])
                rhob_na_count = len(df[(df.nct_val.isna() == False) & (df.SV_psi.isna() == True) & (df.RHOB_merged.isna() == True)])
                print("Total Valid NCT values:",total_nct_count)
                print("Overburden Pressure Nan values count on valid nct values:",sv_psi_na_count)
                print("rhob nan values where overburden prssure is Nan:",rhob_na_count)
        else:
            print("Pore Pressure Calculation looks Good")
            # print(df[df.pore_pressure != df.p_hydro].to_string())


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


    def generate_nct(self):

        def get_nct_semilog(shale_df:pd.DataFrame,depth_col = 'DEPT', dtco_col = 'DTCO_merged'):

            shale_df = shale_df.copy()
            shale_df['nct_val'] = np.nan
            def exp_model(x, a, b, c):
                return a + b * np.exp(-c * x)

            # Sample data 
            x = shale_df[depth_col] 
            y = shale_df[dtco_col]

            mask = np.isfinite(x) & np.isfinite(y) 
            x_clean = x[mask] 
            y_clean = y[mask]

            # Initial guesses for a, b, c (important for convergence)
            initial_guess = [np.mean(y_clean), 1.0, 0.001]

            params, covariance = curve_fit(
                exp_model,
                x_clean,
                y_clean,
                p0=initial_guess,
                maxfev=10000
            )

            a, b, c =  params

            print(f"a = {a:.4f}, b = {b:.4f}, c = {c}")

            # x_fit = np.linspace(x_clean.min(), x_clean.max(), 500)
            y_fit = exp_model(x_clean, a, b, c)

            plt.figure(figsize=(5, 8))

            # Scatter: DTCO vs Depth
            plt.semilogx(
                y_clean,
                x_clean,
                'o',
                markersize=3,
                alpha=0.6,
                label="Data"
            )

            # Fitted curve
            plt.semilogx(
                y_fit,
                x_clean,
                color="red",
                linewidth=2,
                label="Exponential fit"
            )

            plt.xlabel("DTCO")
            plt.ylabel("Depth")

            # Depth increases downward (log-style convention)
            plt.gca().invert_yaxis()


            plt.grid(True, which="both", linestyle="--", alpha=0.6)
            plt.grid(True, which="major", linewidth=0.8)
            plt.grid(True, which="minor", linewidth=0.3)
            plt.legend()
            plt.tight_layout()

            

            ax = plt.gca()

            # Major ticks (10, 100, etc.)
            ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=100))

            # Minor ticks (2–9 between powers of 10)
            ax.xaxis.set_minor_locator(
                LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100)
            )

            # Show normal numbers instead of 10^x
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.xaxis.set_minor_formatter(ScalarFormatter())

            ax.tick_params(axis='x', which='major', length=6)
            ax.tick_params(axis='x', which='minor', length=3)
            plt.savefig('Data/nct.png')
            plt.show()

            shale_df.loc[mask, 'nct_val'] = y_fit

            return shale_df[[depth_col,'nct_val']]


        def get_nct_poly(shale_df:pd.DataFrame,degree:int = 2, depth_col = 'DEPT', dtco_col = 'DTCO_merged'):

            shale_df = shale_df.copy()

            shale_df['nct_val'] = None
            # Sample data 
            x = shale_df[depth_col] 
            y = shale_df[dtco_col]

            mask = np.isfinite(x) & np.isfinite(y)

            x_clean = x[mask]
            y_clean = y[mask]

            # Fit polynomial
            coeffs = np.polyfit(x_clean, y_clean, degree)

            # Create polynomial function
            poly = np.poly1d(coeffs)

            # Generate smooth x values for plotting
            x_smooth = np.linspace(x_clean.min(), x_clean.max(), 500)

            plt.scatter(x_clean, y_clean, s=5, alpha=0.4, label="Data")
            # plt.plot(x_clean, y_clean, alpha=0.4, label="Data")
            plt.plot(x_clean, poly(x_clean), color='red', linewidth=2, label="Polynomial fit")
            plt.legend()
            plt.show()

            shale_df.loc[mask, 'nct_val'] = poly(x_clean)

            return shale_df[[depth_col,'nct_val']]


        def get_nct_log(shale_df:pd.DataFrame,depth_col = 'DEPT', dtco_col = 'DTCO_merged'):
            shale_df = shale_df.copy()
            # Define the model function
            def linear_log_model(y, a, b):
                return a * np.log(y) + b
            
            shale_df['nct_val'] = None
            # Sample data 
            x = shale_df[depth_col] 
            y = shale_df[dtco_col]
            mask = np.isfinite(x) & np.isfinite(y) 
            x_clean = x[mask] 
            y_clean = y[mask]

            # Fit the curve
            params, covariance = curve_fit(linear_log_model, y_clean, x_clean)
            a, b = params

            print(f"a = {a}")
            print(f"b = {b}")

            # Solve for y given x
            def solve_for_y(x, a, b):
                return np.exp((x - b) / a)

            # Solve for all original x values
            y_predicted = solve_for_y(x_clean, a, b)

            # Compare with actual y values
            comparison = np.column_stack((x_clean, y_clean, y_predicted))
            print("\nComparison (first 10 rows):")
            print("x_actual | y_actual | y_predicted")
            print(comparison[:10])

            # Plot 1: Original fit (x vs log(y))
            y_sorted = np.sort(y_clean)
            x_fitted = linear_log_model(y_sorted, a, b)

            # Plot 2: Inverse prediction (y_actual vs y_predicted)
            plt.scatter(y_clean, x_clean, alpha=0.5, label='Predictions', s=20)
            plt.plot(y_predicted,x_clean,'r-', linewidth=2)
            plt.xlabel('y_actual (DTCO_merged)', fontsize=11)
            plt.ylabel('DEPT', fontsize=11)
            plt.title('Inverse Solution: y = exp((x - b)/a)', fontsize=12)
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.gca().invert_yaxis()

            plt.tight_layout()
            plt.show()

            shale_df.loc[mask, 'nct_val'] = y_predicted

        
            return shale_df[[depth_col,'nct_val']]
        
        if self.nct_type == 'semilog':
            self.nct_df = get_nct_semilog(self.shale_df)

        elif self.nct_type == 'poly':
            self.nct_df = get_nct_poly(self.shale_df)

        else:
            self.nct_df = get_nct_log(self.shale_df)
            
    def smoothing(self,cleaned_df=None, track_columns=['RHOB', 'DTCO', 'GR'], method='median', sample_size=11, ema_span=10, sg_polyorder=2, gaussian_sigma=2):
        
        cleaned_df = cleaned_df.copy().astype(float)
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
        if method in method_map:
            for col in track_columns:
                cleaned_df[f"{col}_smoothed"] = method_map[method](cleaned_df[col].to_numpy())
        return cleaned_df

    def calculate_pore_pressure(self):

        def attach_pore_pressure(pp_df: pd.DataFrame, sonic_expo: float = 1.5, sonic_fact: float = 1) -> pd.DataFrame:
            pp_df = pp_df.copy()
            
            # Validate required columns
            required_cols = ['SV_psi', 'p_hydro', 'DTCO_merged', 'nct_val']
            if not all(col in pp_df.columns for col in required_cols):
                raise ValueError(f"Missing required columns. Need: {required_cols}")
            
            # Check for division by zero
            if (pp_df['nct_val'].astype(float) == 0).any():
                raise ValueError("nct_val contains zero values - cannot proceed")
            
            
            pp_df['pore_pressure'] = np.where(pp_df['nct_val'].isna() == False ,pp_df['SV_psi'] - (pp_df['SV_psi'] - pp_df['p_hydro']) * sonic_fact * ((pp_df['DTCO_merged'] / pp_df['nct_val'].astype(float)) ** sonic_expo),pp_df['p_hydro'])
            return pp_df

        new_df = pd.merge(left=self.start_df,right=self.nct_df,how='left',on='DEPT')
        pp_df = new_df[['DEPT','DTCO_merged','nct_val']].copy().drop_duplicates(subset='DEPT').reset_index(drop=True)
        nct_ov_df =pd.merge(left=pp_df,right=self.ovd_df,how='left',on='DEPT')
        nct_ov_df['p_hydro'] = nct_ov_df[self.depth_col]*1.42
        final_result_df = attach_pore_pressure(nct_ov_df)
        try:
            final_result_df = pd.merge(left = final_result_df, right= self.final_df[[self.depth_col,'DTSM_merged']], on = self.depth_col, how = 'left')
        except:
            pass
        final_result_df['pore_pressure'] = np.where(final_result_df['pore_pressure']<final_result_df['p_hydro'],final_result_df['p_hydro'],final_result_df['pore_pressure'])

        self.pore_pressure_df = final_result_df

        self.pp_issue_tracker(self.pore_pressure_df)

    def calculate_fracture_gradient(self):

        def eaton_fracture_grad(smoothed_df:pd.DataFrame, a:int = 1):

            vp = 304.8/smoothed_df['DTCO_merged']
            vs = 304.8/smoothed_df['DTSM_merged']

            vdyn = (vp**2 - 2*vs**2)/(2*(vp**2 - vs**2))
            vstat = a*vdyn

            smoothed_df['fracture_pressure_eaton'] = smoothed_df['pore_pressure'] + (smoothed_df['SV_psi'] - smoothed_df['pore_pressure'])*(vstat/(1-vstat))

        def hubbert_willis_fracture_grad(smoothed_df:pd.DataFrame):

            smoothed_df['fracture_pressure_hw'] = 1/3*(smoothed_df['SV_psi'] + 2*smoothed_df['pore_pressure'])

        try:
            eaton_fracture_grad(self.pore_pressure_df)
        except Exception as e:
            print("Error in Calculating Eaton fracture gradient:",str(e))

        try:
            hubbert_willis_fracture_grad(self.pore_pressure_df)
        except Exception as e:
            print("Error in Calculating hubbert_willis_fracture_gradient:",str(e))
            
        try:
            self.pore_pressure_df = self.smoothing(cleaned_df=self.pore_pressure_df,track_columns=['pore_pressure','fracture_pressure_eaton','fracture_pressure_hw'],sample_size=200,method='mean')
        except Exception as e:
            print(print(e))
            
        ut.save_csv(filename='generated_data',dataframe=self.pore_pressure_df)

    def plot_pore_pressure(self):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10))

        # First plot
        ax1.plot(self.pore_pressure_df['pore_pressure'], self.pore_pressure_df[self.depth_col], label='Pore Pressure', color='green', alpha=0.5)
        ax1.plot(self.pore_pressure_df['p_hydro'], self.pore_pressure_df[self.depth_col], label='HydroStatic Pressure', color='blue')
        try:
            ax1.plot(self.pore_pressure_df['fracture_pressure_eaton'], self.pore_pressure_df[self.depth_col], label='Fracture Pressure Eaton', color='black')
        except:
            pass
        ax1.plot(self.pore_pressure_df['fracture_pressure_hw'], self.pore_pressure_df[self.depth_col], label='Fracture Pressure HW', color='brown')
        ax1.plot(self.pore_pressure_df['SV_psi'], self.pore_pressure_df[self.depth_col], label='Overburden Pressure', color='red')
        ax1.set_ylim(self.pore_pressure_df[self.depth_col].min()+10, self.pore_pressure_df[self.depth_col].max()+50)
        ax1.invert_yaxis()
        ax1.set_xlabel('Value')
        ax1.set_ylabel('Depth (ft)')
        ax1.set_title('Pore Pressure vs Depth')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Second plot (smoothed)
        ax2.plot(self.pore_pressure_df['pore_pressure_smoothed'], self.pore_pressure_df[self.depth_col], label='Pore Pressure', color='green', alpha=0.5)
        ax2.plot(self.pore_pressure_df['p_hydro'], self.pore_pressure_df[self.depth_col], label='HydroStatic Pressure', color='blue')
        try:
            ax2.plot(self.pore_pressure_df['fracture_pressure_eaton_smoothed'], self.pore_pressure_df[self.depth_col], label='Fracture Pressure Eaton', color='black')
        except:
            pass
        ax2.plot(self.pore_pressure_df['fracture_pressure_hw_smoothed'], self.pore_pressure_df[self.depth_col], label='Fracture Pressure HW', color='brown')
        ax2.plot(self.pore_pressure_df['SV_psi'], self.pore_pressure_df[self.depth_col], label='Overburden Pressure', color='red')
        ax2.set_ylim(self.pore_pressure_df[self.depth_col].min()+10, self.pore_pressure_df[self.depth_col].max()+50)
        ax2.invert_yaxis()
        ax2.set_xlabel('Value')
        ax2.set_ylabel('Depth (ft)')
        ax2.set_title('Pore Pressure (Smoothed) vs Depth')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        ut.make_folder()
        plt.savefig('Data/pore_pressure.png')
        plt.tight_layout()
        plt.show()





        
        


    
        





        



    