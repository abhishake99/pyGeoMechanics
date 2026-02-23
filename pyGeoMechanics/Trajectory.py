import wellpathpy as wp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
np.set_printoptions(suppress=True)


class Trajectory_calculator():
    def __init__(self,trajectory_df = None):
        self.trajectory_df = trajectory_df
        if self.trajectory_df is None:
            print("Please Provide Trajectory DataFrame Containing 'md', 'azi', 'inc' as columns")
            return
        
    def calculate_tvd(self, depth_smaple_rate = 0.01, merged_df = None, rkb = 0.0):
        print("Generating Trajectory")

        # header = wp.read_header_json(header_json)
        sample_rate = depth_smaple_rate

        # Get starting md from original dataframe
        start_md = self.trajectory_df['md'].iloc[0]

        # Generate new md values from 0 up to start_md (exclusive)
        new_md = np.arange(0, start_md, sample_rate)

        # Create new dataframe (fill other columns with NaN)
        new_df = pd.DataFrame({
            'md': [0.0],
            'inc':0.0,
            'azi':0.0
        })
        max_md =  max(merged_df['DEPT'].max(),self.trajectory_df['md'].max())

        new_df2 = pd.DataFrame({
            'md': [max_md],
            'inc':self.trajectory_df.iloc[-1].inc,
            'azi':self.trajectory_df.iloc[-1].azi
        })

        # Attach on top
        self.trajectory_df = pd.concat([new_df, self.trajectory_df], ignore_index=True)
        self.trajectory_df = pd.concat([self.trajectory_df,new_df2], ignore_index=True)
        self.trajectory_df = self.trajectory_df.astype(float)
        self.trajectory_df = self.trajectory_df.apply(pd.to_numeric, errors='coerce')
        self.trajectory_df.to_csv('BRK012_wellpath.csv',index=False)

        md, inc, azi = wp.read_csv('BRK012_wellpath.csv')

        os.remove('BRK012_wellpath.csv')

        dev = wp.deviation(md, inc, azi)

        depth_step = depth_smaple_rate
        depths = np.arange(0, dev.md[-1] + depth_step, depth_step)
        pos = dev.minimum_curvature().resample(depths = depths)
        resampled_dev = pos.deviation()
        pos_tvdss = pos.to_tvdss(datum_elevation=rkb)

        data = {
                'MD':resampled_dev.md, 
                'AZI':resampled_dev.azi,
                'INC':resampled_dev.inc,
                'TVDRT':pos.depth,
                'TVDSS':pos_tvdss.depth,
                'YLOC':pos.northing,
                'XLOC':pos.easting,
                # 'NORTHING': pos_wellhead.northing,
                # 'EASTING': pos_wellhead.easting
                }
        calc_df = pd.DataFrame(data)
        
        return calc_df