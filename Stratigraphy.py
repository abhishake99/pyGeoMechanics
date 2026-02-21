import pandas as pd
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# from matplotlib.collections import BrokenBarHCollection

warnings.filterwarnings('ignore')

class Stratigraphy():
    def __init__(self, data_frame=None, depth_col='DEPT', gamma_ray_col=None, dtco_column=None):
        self.gamma_ray_col = gamma_ray_col
        self.dtco_col = dtco_column
        self.depth_col = depth_col
        # Create a copy to avoid SettingWithCopy warnings on the original DF later
        self.data_frame = data_frame[[depth_col, self.gamma_ray_col, self.dtco_col,'TVDRT']].copy()

    def method_1_input_based(self, intervals=None):
        def assign_strat_code(depth, intervals):
            for (start, end), code in intervals.items():
                if start <= depth <= end:
                    return int(code)
            return None
        self.data_frame["strat_code_method1"] = self.data_frame[self.depth_col].apply(assign_strat_code, intervals=intervals)

    def method_2_gr_based(self, gamma_ray_cutoff=None):
        self.data_frame['strat_code_method2'] = np.where(self.data_frame[self.gamma_ray_col] > gamma_ray_cutoff, 1, 0)

    def method_3_vclay_based(self, geological_zones: dict = None):
        # (Your existing method_3 logic remains mostly the same, ensuring efficient concatenation)
        def vshale_gr_avail(gr_df, GR, DEPT, column, v_clay_cutoff=0.4):
            if gr_df is None or gr_df.empty:
                return gr_df
            
            gr_df = gr_df.copy()
            gr_min, gr_max = GR
            dept_min, dept_max = DEPT
            
            zone_df = gr_df[gr_df[self.depth_col].between(dept_min, dept_max)].copy()
            if zone_df.empty:
                return zone_df

            zone_df["v_clay"] = (zone_df[column] - gr_min) / (gr_max - gr_min)
            zone_df["v_clay"] = zone_df["v_clay"].clip(0, 1)
            zone_df["strat_code_method3"] = np.where(zone_df["v_clay"] >= v_clay_cutoff, 1, 0)
            return zone_df
        
        new_df_list = []
        for zone, params in geological_zones.items():
            temp_df = vshale_gr_avail(
                gr_df=self.data_frame,
                GR=params["gr"],
                DEPT=params["dept"],
                column=self.gamma_ray_col,
                v_clay_cutoff=0.4
            )
            new_df_list.append(temp_df)

        if new_df_list:
            new_df = pd.concat(new_df_list, ignore_index=True)
            self.data_frame = new_df.sort_values(self.depth_col).reset_index(drop=True)

    def generate_band(self, method=None):
        col_name = f'strat_code_method{method}'
        
        # 1. Prepare Data
        df = self.data_frame.dropna(subset=[col_name]).copy()
        df[col_name] = df[col_name].astype(int)
        df = df.sort_values(self.depth_col)

        if df.empty:
            print("No data to plot.")
            return

        # 2. Vectorized Grouping (Speed Up)
        # Identify contiguous blocks of the same lithology
        df['block_id'] = df[col_name].ne(df[col_name].shift()).cumsum()

        # Aggregate: Get Top, Base, and Code for each block
        blocks = df.groupby('block_id').agg(
            top=(self.depth_col, 'min'),
            base=(self.depth_col, 'max'),
            code=(col_name, 'first')
        ).reset_index(drop=True)

        # 3. Setup Definitions
        lithology_map = {
            0: "Non-Shale", 1: "Shale", 2: "Sand", 
            3: "Limestone", 4: "Dolomite", 5: "Coal", 6: "Halide"
        }
        color_map = {
            "Non-Shale": "lightblue", "Shale": "brown", "Sand": "gold",
            "Limestone": "lightgray", "Dolomite": "pink", "Coal": "black", "Halide": "purple"
        }

        blocks['lithology'] = blocks['code'].map(lithology_map)
        
        # 4. Plotting using ax.bar (Vectorized and Standard)
        fig, ax = plt.subplots(figsize=(3, 8))

        # Group by lithology so we can plot all "Shale" blocks in one command
        # This automatically handles the legend correctly.
        for lith_name, group in blocks.groupby('lithology'):
            if pd.isna(lith_name): continue
            
            color = color_map.get(lith_name, 'gray')
            
            # ax.bar parameters:
            # x      = X center position (0.5 to center it between 0 and 1)
            # height = Vertical thickness of the layer (base - top)
            # width  = Width of the bar (1.0 to fill the track)
            # bottom = The starting Y depth (top)
            ax.bar(
                x=[0.5] * len(group),          # Center all bars at x=0.5
                height=group['base'] - group['top'], 
                bottom=group['top'],
                width=1.0,
                color=color,
                label=lith_name,
                align='center',
                edgecolor='none' # Optional: removes thin lines between adjacent same-type blocks
            )

        # 5. Formatting
        ax.set_ylim(df[self.depth_col].max(), df[self.depth_col].min()) # Invert Y axis for depth
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_ylabel("Depth")
        ax.set_title("Stratigraphy vs Depth")

        ax.legend(loc='upper right', fontsize=8)

        plt.savefig('Data/Statigraphy.png')
        plt.tight_layout()
        plt.show()