import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

class OverburdensStressCalculator:
    def __init__(self, df, depth_col='DEPT', density_col='RHOB_merged', 
                 area_type='onshore', water_depth=0, surface_density=None):
        """
        Initialize the calculator.
        
        Parameters:
        - df: DataFrame with depth and density columns
        - depth_col: Column name for depth (TVD in ft)
        - density_col: Column name for bulk density (g/cc)
        - area_type: 'onshore' or 'offshore'
        - water_depth: Water depth in ft (for offshore) 
        - surface_density: Custom surface density (g/cc). If None, uses default. 
        """
        # Sort by depth and reset index to ensure proper ordering
        self.df = df.copy().sort_values(depth_col).reset_index(drop=True)
        # self.df = self.df[[depth_col,density_col]]
        self.depth_col = depth_col
        self.density_col = density_col
        self.area_type = area_type.lower()
        self.water_depth = water_depth
        
        # Set default surface densities
        if surface_density is None:
            if self.area_type == 'offshore':
                self.surface_density = 1.03  # Seawater density (g/cc)
                self.mudline_density = 1.65  # Mud line density (g/cc)
            else:  # onshore
                self.surface_density = 1.85  # Onshore sediment density (g/cc)
                self.mudline_density = self.surface_density
        else:
            self.surface_density = surface_density
            self.mudline_density = surface_density
    
    def fill_density_gap(self, method='exponential', depth_step=0.15, k=None):
        """
        Fill missing density from surface/mudline to first available log depth.

        Parameters
        ----------
        method : 'exponential' or 'linear'
        depth_step : Depth increment (ft)
        k : Compaction coefficient for exponential (auto if None)
        The compaction coefficient (k) in exponential density extrapolation is a geomechanical parameter that controls how fast bulk density increases 
        with depth due to sediment compaction.

        Returns
        -------
        DataFrame with filled density column
        """

        df = self.df.copy()

        # First valid density index
        first_idx = df[self.density_col].first_valid_index()
        if first_idx is None:
            raise ValueError("No valid density values found in dataframe.")

        first_depth = df.loc[first_idx, self.depth_col]
        first_density = df.loc[first_idx, self.density_col]

        # Depths to fill
        fill_depths = np.arange(0, first_depth, depth_step)

        density_filled = []

        # ---------------- OFFSHORE CASE ----------------
        if self.area_type == 'offshore':

            for z in fill_depths:
                if z <= self.water_depth:
                    density_filled.append(self.surface_density)
                else:
                    z_eff = z - self.water_depth
                    z_max = first_depth - self.water_depth

                    if method == 'exponential':
                        if k is None:
                            k = -np.log(0.05) / z_max

                        rho = self.mudline_density + (
                            (first_density - self.mudline_density)
                            * (1 - np.exp(-k * z_eff))
                        )
                    else:
                        rho = self.mudline_density + (
                            (first_density - self.mudline_density)
                            * (z_eff / z_max)
                        )

                    density_filled.append(rho)

        # ---------------- ONSHORE CASE ----------------
        else:
            if method == 'exponential':
                if k is None:
                    k = -np.log(0.05) / first_depth

                density_filled = (
                    self.surface_density
                    + (first_density - self.surface_density)
                    * (1 - np.exp(-k * fill_depths))
                )
            else:
                density_filled = (
                    self.surface_density
                    + (first_density - self.surface_density)
                    * (fill_depths / first_depth)
                )

        # Build extrapolated DataFrame
        extrap_df = pd.DataFrame({
            self.depth_col: fill_depths,
            self.density_col: density_filled
        })

        # Merge and return
        final_df = pd.concat([extrap_df, df], ignore_index=True)
        final_df = final_df.sort_values(self.depth_col).reset_index(drop=True)

        return final_df
    
    def calculate_sv_psi(self, filled_df):
        """
        Calculate overburden stress in PSI.
        Sv (psi) = 0.433 * Σ(ρ * Δh)
        where ρ is in g/cc and Δh is in ft
        """
        filled_df = filled_df.dropna(subset=[self.depth_col, self.density_col])
        depths = filled_df[self.depth_col].values
        densities = filled_df[self.density_col].values
        
        # Calculate depth increments
        depth_increments = np.diff(depths)
        
        # Use average density between points
        avg_densities = (densities[:-1] + densities[1:]) / 2
        
        # Cumulative overburden (0.433 conversion factor)
        overburden = np.cumsum(0.433 * 3.281 * avg_densities * depth_increments)
        
        # Create result dataframe
        result_df = filled_df.copy()
        result_df['SV_psi'] = np.insert(overburden, 0, 0)
        
        return result_df
    
    def calculate_sv_mpa(self, filled_df):
        """
        Calculate overburden stress in MPa.
        Sv (MPa) = (ρ * g * TVD) / 10^6
        where ρ is in kg/m³, g = 9.80665 m/s², TVD in meters
        """
        filled_df = filled_df.dropna(subset=[self.depth_col, self.density_col])
        depths_m = filled_df[self.depth_col].values * 0.3048  # Convert ft to m
        densities_kg_m3 = filled_df[self.density_col].values * 1000  # Convert g/cc to kg/m³
        
        g = 9.80665  # m/s²
        
        # Calculate depth increments
        depth_increments_m = np.diff(depths_m)
        
        # Use average density between points
        avg_densities = (densities_kg_m3[:-1] + densities_kg_m3[1:]) / 2
        
        # Cumulative overburden in Pa, convert to MPa
        overburden_pa = np.cumsum(avg_densities * g * depth_increments_m)
        overburden_mpa = overburden_pa / 1e6
        
        # Create result dataframe
        result_df = filled_df.copy()
        result_df['SV_mpa'] = np.insert(overburden_mpa, 0, 0)
        
        return result_df
    
    def calculate(self, gap_fill_method='exponential', units=['psi', 'mpa']):
        """
        Full calculation pipeline.
        
        Parameters:
        - gap_fill_method: 'exponential' or 'linear'
        - units: List of desired output units ['psi', 'mpa', or both]
        
        Returns:
        - DataFrame with SV calculations
        """
        # Fill density gap
        filled_df = self.fill_density_gap(method=gap_fill_method)
        
        # Calculate in requested units
        if 'psi' in units:
            filled_df = self.calculate_sv_psi(filled_df)
        
        if 'mpa' in units:
            filled_df = self.calculate_sv_mpa(filled_df)
        
        return filled_df
    
    def plot_welllog_style(self, result_df, figsize=(10, 14)):
        """
        Plot professional well log style with density on top and Sv curves below.
        """
        fig, axes = plt.subplots(1, 3, figsize=figsize, sharey=True)
        
        original_data = self.df.copy()
        top_original_depth = original_data[self.depth_col].min()
        gap_filled = result_df[result_df[self.depth_col] < top_original_depth]
        original_section = result_df[result_df[self.depth_col] >= top_original_depth]
        
        ax1 = axes[0]
        if len(gap_filled) > 0:
            ax1.plot(gap_filled[self.density_col], gap_filled[self.depth_col], 
                    'k--', linewidth=1.5, alpha=0.7, label='Synthetic')
        if len(original_section) > 0:
            ax1.plot(original_section[self.density_col], original_section[self.depth_col], 
                    'k-', linewidth=1.2)
            ax1.fill_betweenx(original_section[self.depth_col], 
                             original_section[self.density_col].min(),
                             original_section[self.density_col], 
                             alpha=0.3, color='black', label='Measured')
        ax1.invert_yaxis()
        ax1.set_xlabel('Density (g/cm³)', fontsize=11, fontweight='bold')
        ax1.set_ylabel('Depth (ft TVD)', fontsize=11, fontweight='bold')
        ax1.set_title('Bulk Density', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.4, linestyle=':', linewidth=0.8)
        ax1.legend(loc='upper right', fontsize=9)
        
        ax2 = axes[1]
        if 'SV_psi' in result_df.columns:
            ax2.plot(result_df['SV_psi'], result_df[self.depth_col], 'k-', linewidth=1.2)
            ax2.fill_betweenx(result_df[self.depth_col], 0, result_df['SV_psi'], 
                             alpha=0.3, color='green')
            ax2.invert_yaxis()
            ax2.set_xlabel('Sv (psi)', fontsize=11, fontweight='bold')
            ax2.set_title('Overburden Stress (PSI)', fontsize=12, fontweight='bold')
            ax2.grid(True, alpha=0.4, linestyle=':', linewidth=0.8)
        
        ax3 = axes[2]
        if 'SV_mpa' in result_df.columns:
            ax3.plot(result_df['SV_mpa'], result_df[self.depth_col], 'k-', linewidth=1.2)
            ax3.fill_betweenx(result_df[self.depth_col], 0, result_df['SV_mpa'], 
                             alpha=0.3, color='red')
            ax3.invert_yaxis()
            ax3.set_xlabel('Sv (MPa)', fontsize=11, fontweight='bold')
            ax3.set_title('Overburden Stress (MPa)', fontsize=12, fontweight='bold')
            ax3.grid(True, alpha=0.4, linestyle=':', linewidth=0.8)
        
        plt.suptitle('Well Log - Overburden Stress Analysis', fontsize=13, fontweight='bold', y=0.995)
        plt.savefig('Data/overburden_stress.png')
        plt.tight_layout()
        plt.show()
    
    def plot_results(self, result_df, show_original=True, figsize=(14, 8)):
        """
        Plot overburden stress results with density log.
        """
        fig = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(1, 3, figure=fig, hspace=0.3, wspace=0.3)
        
        original_data = self.df.copy()
        top_original_depth = original_data[self.depth_col].min()
        gap_filled = result_df[result_df[self.depth_col] < top_original_depth]
        original_section = result_df[result_df[self.depth_col] >= top_original_depth]
        
        ax1 = fig.add_subplot(gs[0, 0])
        if len(gap_filled) > 0:
            ax1.plot(gap_filled[self.density_col], gap_filled[self.depth_col], 
                    'b--', linewidth=2.5, label='Gap-filled (synthetic)', alpha=0.9)
            ax1.fill_betweenx(gap_filled[self.depth_col], gap_filled[self.density_col], 
                             alpha=0.2, color='blue', label='Gap-filled area')
        if len(original_section) > 0:
            orig_step = max(1, len(original_section) // 100)
            orig_plot = original_section.iloc[::orig_step]
            ax1.plot(orig_plot[self.density_col], orig_plot[self.depth_col], 
                    'b-', linewidth=2.5, label='Original log data')
            ax1.fill_betweenx(orig_plot[self.depth_col], orig_plot[self.density_col], 
                             alpha=0.25, color='blue')
        ax1.invert_yaxis()
        ax1.set_xlabel('Bulk Density (g/cc)', fontsize=11, fontweight='bold')
        ax1.set_ylabel('Depth (ft TVD)', fontsize=11, fontweight='bold')
        ax1.set_title('Density Log with Gap Filling', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3, which='both')
        ax1.legend(loc='best', fontsize=9)
        
        ax2 = fig.add_subplot(gs[0, 1])
        if 'SV_psi' in result_df.columns:
            step_size = max(1, len(result_df) // 200)
            plot_df = result_df.iloc[::step_size]
            ax2.plot(plot_df['SV_psi'], plot_df[self.depth_col], 'g-', linewidth=2.5)
            ax2.fill_betweenx(plot_df[self.depth_col], plot_df['SV_psi'], 
                             alpha=0.2, color='green')
            ax2.invert_yaxis()
            ax2.set_xlabel('Overburden Stress (psi)', fontsize=11, fontweight='bold')
            ax2.set_ylabel('Depth (ft TVD)', fontsize=11, fontweight='bold')
            ax2.set_title('Overburden Stress (PSI)', fontsize=12, fontweight='bold')
            ax2.grid(True, alpha=0.3, which='both')
        
        ax3 = fig.add_subplot(gs[0, 2])
        if 'SV_mpa' in result_df.columns:
            step_size = max(1, len(result_df) // 200)
            plot_df = result_df.iloc[::step_size]
            ax3.plot(plot_df['SV_mpa'], plot_df[self.depth_col], 'r-', linewidth=2.5)
            ax3.fill_betweenx(plot_df[self.depth_col], plot_df['SV_mpa'], 
                             alpha=0.2, color='red')
            ax3.invert_yaxis()
            ax3.set_xlabel('Overburden Stress (MPa)', fontsize=11, fontweight='bold')
            ax3.set_ylabel('Depth (ft TVD)', fontsize=11, fontweight='bold')
            ax3.set_title('Overburden Stress (MPa)', fontsize=12, fontweight='bold')
            ax3.grid(True, alpha=0.3, which='both')
        
        plt.suptitle('Overburden Stress (Sv) Analysis', fontsize=14, fontweight='bold', y=0.98)
        plt.tight_layout()
        plt.show()
    
    def plot_detailed(self, result_df, depth_range=None, figsize=(16, 10)):
        """
        Plot detailed view with multiple subplots for analysis.
        """
        top_original_depth = self.df[self.depth_col].min()
        gap_filled = result_df[result_df[self.depth_col] < top_original_depth]
        original_section = result_df[result_df[self.depth_col] >= top_original_depth]
        
        if depth_range:
            plot_df = result_df[(result_df[self.depth_col] >= depth_range[0]) & 
                               (result_df[self.depth_col] <= depth_range[1])]
            gap_filled = gap_filled[(gap_filled[self.depth_col] >= depth_range[0]) & 
                                   (gap_filled[self.depth_col] <= depth_range[1])]
        else:
            plot_df = result_df
        
        fig = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)
        
        ax1 = fig.add_subplot(gs[0, 0])
        if len(gap_filled) > 0:
            ax1.plot(gap_filled[self.density_col], gap_filled[self.depth_col], 
                    'b--', linewidth=2.5, label='Gap-filled (synthetic)', alpha=0.9)
            ax1.fill_betweenx(gap_filled[self.depth_col], gap_filled[self.density_col], 
                             alpha=0.2, color='blue')
        if len(plot_df) > 0:
            orig_step = max(1, len(plot_df) // 100)
            orig_plot = plot_df.iloc[::orig_step]
            ax1.plot(orig_plot[self.density_col], orig_plot[self.depth_col], 
                    'b-', linewidth=2.5, label='Original log data')
            ax1.fill_betweenx(orig_plot[self.depth_col], orig_plot[self.density_col], 
                             alpha=0.25, color='blue')
        ax1.invert_yaxis()
        ax1.set_xlabel('Density (g/cc)', fontsize=10, fontweight='bold')
        ax1.set_ylabel('Depth (ft TVD)', fontsize=10, fontweight='bold')
        ax1.set_title('Bulk Density Log', fontsize=11, fontweight='bold')
        ax1.grid(True, alpha=0.3, which='both')
        ax1.legend(loc='best', fontsize=9)
        
        ax2 = fig.add_subplot(gs[0, 1])
        if 'SV_psi' in plot_df.columns:
            step_size = max(1, len(plot_df) // 100)
            plot_psi = plot_df.iloc[::step_size]
            ax2.plot(plot_psi['SV_psi'], plot_psi[self.depth_col], 
                    'g-', linewidth=2.5, label='Sv (psi)')
            ax2.fill_betweenx(plot_psi[self.depth_col], plot_psi['SV_psi'], 
                             alpha=0.2, color='green')
            ax2.invert_yaxis()
            ax2.set_xlabel('Overburden Stress (psi)', fontsize=10, fontweight='bold')
            ax2.set_ylabel('Depth (ft TVD)', fontsize=10, fontweight='bold')
            ax2.set_title('Overburden Stress (PSI)', fontsize=11, fontweight='bold')
            ax2.grid(True, alpha=0.3, which='both')
        
        ax3 = fig.add_subplot(gs[1, 0])
        if 'SV_mpa' in plot_df.columns:
            step_size = max(1, len(plot_df) // 100)
            plot_mpa = plot_df.iloc[::step_size]
            ax3.plot(plot_mpa['SV_mpa'], plot_mpa[self.depth_col], 
                    'r-', linewidth=2.5, label='Sv (MPa)')
            ax3.fill_betweenx(plot_mpa[self.depth_col], plot_mpa['SV_mpa'], 
                             alpha=0.2, color='red')
            ax3.invert_yaxis()
            ax3.set_xlabel('Overburden Stress (MPa)', fontsize=10, fontweight='bold')
            ax3.set_ylabel('Depth (ft TVD)', fontsize=10, fontweight='bold')
            ax3.set_title('Overburden Stress (MPa)', fontsize=11, fontweight='bold')
            ax3.grid(True, alpha=0.3, which='both')
        
        ax4 = fig.add_subplot(gs[1, 1])
        if len(plot_df) > 1:
            step_size = max(1, len(plot_df) // 100)
            plot_grad = plot_df.iloc[::step_size]
            density_gradient = np.gradient(plot_grad[self.density_col].values, 
                                         plot_grad[self.depth_col].values)
            ax4.plot(density_gradient, plot_grad[self.depth_col], 
                    'purple', linewidth=2.5, label='Density Gradient')
            ax4.invert_yaxis()
            ax4.set_xlabel('Density Gradient (g/cc/ft)', fontsize=10, fontweight='bold')
            ax4.set_ylabel('Depth (ft TVD)', fontsize=10, fontweight='bold')
            ax4.set_title('Density Gradient', fontsize=11, fontweight='bold')
            ax4.grid(True, alpha=0.3, which='both')
            ax4.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
        
        title = 'Detailed Overburden Stress (Sv) Analysis'
        if depth_range:
            title += f' (Depth: {depth_range[0]:.0f} - {depth_range[1]:.0f} ft)'
        plt.suptitle(title, fontsize=13, fontweight='bold', y=0.995)
        
        
        plt.tight_layout()
        plt.show()
        
