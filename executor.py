from runner import run_geomechanics


outlier_ranges = {
    'DTCO_merged': (45, 140),
    'RHOB_merged': (1.5, 2.8),
    'GR_merged': (10, 200),
    'DTSM_merged': (70, 350)
}

run_geomechanics( las_folder_path=r"C:\Users\HP\Desktop\Enquest-petro-sols\Geomechanics-software\sample_data\2",
        trajectory_file_path=r'C:\Users\HP\Desktop\Enquest-petro-sols\Geomechanics-software\sample_data\1\1. Well Data\Well Trajectory\BRK012_Survey - Copy.csv',outlier_ranges=outlier_ranges)