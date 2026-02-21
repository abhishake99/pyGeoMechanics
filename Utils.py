from pathlib import Path
import pandas as pd

def make_folder():
    folder_path = Path("Data")
    try:
        folder_path.mkdir(parents=True, exist_ok=True)
        print(f"Folder '{folder_path}' created or already exists.")
    except Exception as e:
        print(f"An error occurred: {e}")
        
def save_csv(filename:str,dataframe:pd.DataFrame):
    
    dataframe.to_csv(f'Data/{filename}.csv',index=False)
    