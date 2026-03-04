from pathlib import Path
import pandas as pd

def make_folder(folder_name:str = "Data"):
    folder_path = Path(folder_name)
    try:
        folder_path.mkdir(parents=True, exist_ok=True)
        print(f"Folder '{folder_path}' created or already exists.")
    except Exception as e:
        print(f"An error occurred: {e}")
        
def save_csv(filename:str,dataframe:pd.DataFrame,folder_name:str = "Data"):
    
    dataframe.to_csv(f'{folder_name}/{filename}.csv',index=False)
    