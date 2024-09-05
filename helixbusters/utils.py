import pandas as pd
import os
import re
from Bio import SeqIO

def read_excel_column(file_path):
    """
    Reads a specified column from an Excel or CSV file and performs checks.

    Args:
        file_path (str): Path to the file (Excel or CSV).

    Returns:
        pandas.DataFrame: DataFrame containing the file data if successful, otherwise an empty DataFrame.
    """
    required_columns = ['Sample', 'Replicate', 'Group', 'PathReadF', 'SampleBarcode1']
    check_columns = ['Sample','PathReadF', 'SampleBarcode1']

    try:
        # Determine file type based on extension
        _, file_extension = os.path.splitext(file_path)

        if file_extension.lower() == '.xlsx' or file_extension.lower() == '.xls':
            # Read Excel file
            df = pd.read_excel(file_path, header=0)
        elif file_extension.lower() == '.csv':
            # Read CSV file
            df = pd.read_csv(file_path, header=0)
        else:
            raise ValueError(f"Unsupported file type: {file_extension}. Only .xlsx, .xls, and .csv are supported.")

        # Check if all required columns are present
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"The following required columns are missing: {', '.join(missing_columns)}")

        # Check if the DataFrame has exactly 5 columns
        if len(df.columns) != 5:
            raise ValueError(f"The file does not have exactly 5 columns. Found {len(df.columns)} columns.")

        # Check if required columns are unique
        for col in check_columns:
            if not df[col].is_unique:
                raise ValueError(f"The values in column '{col}' are not unique.")

        return df

    except Exception as e:
        print(f"An error occurred: {e}")
        return pd.DataFrame()