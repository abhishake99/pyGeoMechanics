import glob
import lasio

class Ingestor:
    def __init__(self,folder_path = None):
        self.folder_path = folder_path
        self.ingestion_key = self.load_directory()

    def load_directory(self): 
        if self.folder_path is None:
             print("Provide a folder path please")
             return None
        else:
            glob_path = self.folder_path + '\\*.las'
            files = glob.glob(glob_path)
            dic_df = {}
            for file in files:
                las = lasio.read(file)
                df = las.df()
                sheet_name = file.split('\\')[-1].replace('.las', '')[:31] 
                check_df = df.reset_index()
                if 'DEPT' in check_df.columns:
                    df_sorted = check_df.sort_values(by='DEPT')
                elif 'DEPT:2' in check_df.columns:
                    df_sorted = check_df.sort_values(by='DEPT:2')
                elif 'DEPTH' in check_df.columns:
                        df_sorted = check_df.sort_values(by='DEPTH')
                else:
                    df_sorted = check_df 
                df_sorted.rename(columns={'DEPT:2':'DEPT','DT':'DTCO','DTC':'DTCO','DEPTH':'DEPT'}, inplace=True)
                dic_df[sheet_name] = df_sorted
            print("Returned a dictionary containing dataframes as values of las files names as keys")
            return dic_df
    


