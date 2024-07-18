import unittest
from helixbusters.utils import read_excel_column


class TestUtils(unittest.TestCase):

    def test_read_excel_column(self):
        # Use a sample Excel file for testing
        file_path = 'path/to/test/file.xlsx'

        # Call the function
        result = read_excel_column(file_path)

        # Check the result (update expected_data with what you expect from the test file)
        expected_data = ['value1', 'value2', 'value3']
        self.assertEqual(result, expected_data)


if __name__ == '__main__':
    unittest.main()