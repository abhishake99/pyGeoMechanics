import pandas as pd
import numpy as np
import pyGeoMechanics.Utils as ut

from pyGeoMechanics.Ingestion import Ingestor
from pyGeoMechanics.Preprocessing import PreProcessor
from pyGeoMechanics.Stratigraphy import Stratigraphy
from pyGeoMechanics.Overburden import OverburdensStressCalculator
from pyGeoMechanics.Pore_pressure import PorePressurePredictor
from pyGeoMechanics.Trajectory import Trajectory_calculator


def run_geomechanics(
    las_folder_path,
    trajectory_file_path,
    tvd_csv_path='tvd_calculated.csv',
    rkb=0.0,
    tvd_sample_rate=0.01,
    geological_zones=None,
    area_type='onshore',
    gap_fill_method='exponential',
    overburden_units=['psi', 'mpa'],
    strat_method='3',
    nct_type='semilog',
    outlier_ranges=None
):
    # Data Loading
    data_model = Ingestor(folder_path=las_folder_path) # Creates a dictionary from LAS files
    processor = PreProcessor(ingestion_key=data_model.ingestion_key) # Merges Logs and also holds other data preparation functions like despiking, interpolation, smoothing, outlier removal
    merged_df = processor.merged_df # Merged DataFrame
    ut.make_folder() # Creating "Data" folder in the current directory

    # TVD Generation
    try:
        tvd_df = pd.read_csv(f'Data/{tvd_csv_path}')  # Reads the TVD csv file if present else it calculates the TVD using the md, azi, inc columns
    except FileNotFoundError:
        traj_df = pd.read_csv(trajectory_file_path)
        trajectory_object = Trajectory_calculator(trajectory_df=traj_df)
        tvd_df = trajectory_object.calculate_tvd(depth_smaple_rate=tvd_sample_rate, merged_df=merged_df, rkb=rkb)
        tvd_df['DEPT'] = tvd_df['MD'].round(2)
        
        tvd_df.to_csv(f'Data/{tvd_csv_path}', index=False)
    merged_df = pd.merge(left=merged_df, right=tvd_df[['DEPT', 'TVDRT']], how='left', on='DEPT') 

    # Data Processing
    working_df = merged_df[['DEPT', 'TVDRT', 'DTCO_merged', 'RHOB_merged', 'GR_merged', 'DTSM_merged']].copy()
    processing_cols = ['DTCO_merged', 'RHOB_merged', 'GR_merged', 'DTSM_merged']
    # Default outlier ranges
    default_ranges = {
        'DTCO_merged': (40, 150),
        'DTSM_merged': (60, 360),
        'RHOB_merged': (1.2, 3),
        'GR_merged': (0, 250)
    }
    if outlier_ranges is None:
        outlier_ranges = default_ranges
    else:
        # Fill missing keys with defaults
        for col in processing_cols:
            if col not in outlier_ranges:
                outlier_ranges[col] = default_ranges[col]
    for col in processing_cols:
        try:
            despiked_df = processor.despike_log(data_frame=working_df[['DEPT', col]], despiked_column=col, window_size=20, overlap_size=4)
            working_df = pd.merge(left=working_df, right=despiked_df, how='left', on='DEPT')
        except Exception as e:
            print(str(e))
        min_thresh, max_thresh = outlier_ranges[col]
        try:
            outlier_df = processor.remove_outliers(
                data_frame=working_df[['DEPT', f"{col}_despiked"]],
                outlier_column=f"{col}_despiked",
                min_threshold=min_thresh,
                max_threshold=max_thresh
            )
            working_df = pd.merge(left=working_df, right=outlier_df, how='left', on='DEPT')
        except Exception as e:
            print(str(e))
    working_df = processor.interpolate(cleaned_df=working_df, track_columns=['GR_merged_despiked_outlier', 'RHOB_merged_despiked_outlier'], method='linear')
    for col in processing_cols:
        try:
            working_df.drop(columns=[col], inplace=True)
        except:
            pass
    working_df.rename(columns={
        'GR_merged_despiked_outlier_interpolated': 'GR_merged',
        'RHOB_merged_despiked_outlier_interpolated': 'RHOB_merged',
        'DTCO_merged_despiked_outlier': 'DTCO_merged',
        'DTSM_merged_despiked_outlier': 'DTSM_merged'
    }, inplace=True)
    working_df['flag'] = np.where(working_df['DTCO_merged'].isna(), 1, 0)
    working_df = processor.interpolate(cleaned_df=working_df, track_columns=['DTCO_merged'], method='linear')
    working_df.drop(columns=['DTCO_merged'], inplace=True)
    working_df.rename(columns={'DTCO_merged_interpolated': 'DTCO_merged'}, inplace=True)
    working_df['DTSM_merged'] = np.where(working_df.flag == 1, np.nan, working_df['DTSM_merged'])
    working_df['DTSM_merged'] = np.where(
        working_df.DTSM_merged.isna(),
        (304 * working_df.DTCO_merged) / (262.78 - (1.172 * working_df.DTCO_merged)),
        working_df['DTSM_merged']
    )
    working_df = working_df[['DEPT', 'TVDRT', 'DTCO_merged', 'RHOB_merged', 'GR_merged', 'DTSM_merged']].copy()

    # Stratigraphy Calculation
    if geological_zones is None:
        geological_zones = {"ZONE_1": {"dept": (0, 5000), "gr": (25, 120)}}
    strat = Stratigraphy(data_frame=working_df, gamma_ray_col='GR_merged', dtco_column='DTCO_merged', depth_col='DEPT')
    strat.method_3_vclay_based(geological_zones=geological_zones)
    strat.generate_band(method=strat_method)

    # Overburden Calculation
    rhob_df = working_df[['TVDRT', 'RHOB_merged', 'DEPT']].copy()
    calc_onshore = OverburdensStressCalculator(
        rhob_df,
        depth_col='TVDRT',
        density_col='RHOB_merged',
        area_type=area_type
    )
    result_onshore = calc_onshore.calculate(gap_fill_method=gap_fill_method, units=overburden_units)
    calc_onshore.plot_welllog_style(result_onshore)

    # Pore Pressure Calculation
    Pore_pressure_obj = PorePressurePredictor(
        df=working_df,
        overburden_df=result_onshore,
        stratigraphy_df=strat.data_frame,
        strat_method=strat_method,
        nct_type=nct_type,
        depth_col='TVDRT'
    )
    return {
        "working_df": working_df,
        "stratigraphy": strat.data_frame,
        "overburden": result_onshore,
        "pore_pressure_obj": Pore_pressure_obj
    }


if __name__ == "__main__":
    run_geomechanics()
