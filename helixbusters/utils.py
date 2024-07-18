import pandas as pd


def read_excel_column(file_path, sheet_name, column_name):
    """
    Reads a specified column from an Excel file.

    Args:
        file_path (str): Path to the Excel file.
        sheet_name (str): Name of the sheet to read from.
        column_name (str): Name of the column to extract.

    Returns:
        list: Values from the specified column.
    """
    try:
        # Read the Excel file
        df = pd.read_excel(file_path, sheet_name=sheet_name)

        # Check if the column exists
        if column_name not in df.columns:
            raise ValueError(f"Column '{column_name}' not found in the sheet '{sheet_name}'.")

        # Extract the column
        column_data = df[column_name].tolist()

        return column_data

    except Exception as e:
        print(f"An error occurred: {e}")
        return []

# Additional utility functions can be added here
