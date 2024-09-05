import unittest
from helixbusters.core import Helixbusters

class TestHelixbusters(unittest.TestCase):

    def test_read_column_from_excel(self):
        hb = Helixbusters(samplesheet='test')
        file_path = 'path/to/test/file.xlsx'

        # Call the method
        result = hb.read_column_from_excel(file_path)

        # Check the result (update expected_data with what you expect from the test file)
        expected_data = ['value1', 'value2', 'value3']
        self.assertEqual(result, expected_data)

if __name__ == '__main__':
    unittest.main()