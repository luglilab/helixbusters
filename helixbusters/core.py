from helixbusters.utils import read_excel_column


class Helixbusters:
    def __init__(self, samplesheet, species, mismatch):
        self.samplesheet = samplesheet
        self.species = species

        # Validate species attribute
        if species.lower() not in ['mouse', 'human']:
            raise ValueError("Species must be 'mouse' or 'human'.")

        # Validate mismatch attribute
        if not 0 <= mismatch <= 3:
            raise ValueError("Mismatch must be between 0 and 3.")

        self.mismatch = mismatch

    def read_column_from_excel(self):
        """
        Reads a specified column from an Excel file and stores it in the instance.

        Args:
            sheet_name (str): Name of the sheet to read from.
            column_name (str): Name of the column to extract.

        Returns:
            list: Values from the specified column.
        """
        column_data = read_excel_column(self.samplesheet)
        return column_data
