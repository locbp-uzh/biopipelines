"""
Unit tests for the refactoring mixins.

Tests InputHandlerMixin, DatasheetNavigatorMixin, and FilePathDescriptor.
"""

import unittest
import os
import tempfile
import shutil
from unittest.mock import Mock, MagicMock

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from mixins import InputHandlerMixin, DatasheetNavigatorMixin, FilePathDescriptor, ResolvedInput
from base_config import DatasheetInfo, DatasheetContainer


class TestInputHandlerMixin(unittest.TestCase):
    """Test InputHandlerMixin functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.mixin = InputHandlerMixin()

    def test_resolve_input_with_list(self):
        """Test resolving list of file paths."""
        input_files = ['/path/to/file1.pdb', '/path/to/file2.pdb']
        resolved = self.mixin.resolve_input(input_files, 'structures')

        self.assertEqual(resolved.files, input_files)
        self.assertEqual(resolved.ids, ['file1', 'file2'])
        self.assertIsNone(resolved.datasheets)
        self.assertEqual(resolved.source, input_files)
        self.assertFalse(resolved.is_tool_output)

    def test_resolve_input_with_string_file(self):
        """Test resolving single file path string."""
        input_file = '/path/to/structure.pdb'
        resolved = self.mixin.resolve_input(input_file, 'structures')

        self.assertEqual(resolved.files, [input_file])
        self.assertEqual(resolved.ids, ['structure'])
        self.assertIsNone(resolved.datasheets)

    def test_resolve_input_with_dict(self):
        """Test resolving dictionary input."""
        input_dict = {
            'structures': ['/path/to/file1.pdb', '/path/to/file2.pdb'],
            'structure_ids': ['id1', 'id2'],
            'datasheets': {'main': '/path/to/datasheet.csv'}
        }
        resolved = self.mixin.resolve_input(input_dict, 'structures')

        self.assertEqual(resolved.files, input_dict['structures'])
        self.assertEqual(resolved.ids, input_dict['structure_ids'])
        self.assertEqual(resolved.datasheets, input_dict['datasheets'])

    def test_resolve_input_extracts_ids_from_filenames(self):
        """Test that IDs are extracted from filenames when not provided."""
        input_dict = {
            'structures': ['/path/to/struct1.pdb', '/path/to/struct2.pdb']
        }
        resolved = self.mixin.resolve_input(input_dict, 'structures')

        self.assertEqual(resolved.ids, ['struct1', 'struct2'])

    def test_resolve_multiple_inputs(self):
        """Test resolving multiple inputs at once."""
        input1 = ['/path/to/file1.pdb']
        input2 = ['/path/to/file2.pdb']

        results = self.mixin.resolve_multiple_inputs(input1, input2, expected_type='structures')

        self.assertEqual(len(results), 2)
        self.assertIsInstance(results[0], ResolvedInput)
        self.assertIsInstance(results[1], ResolvedInput)


class TestDatasheetNavigatorMixin(unittest.TestCase):
    """Test DatasheetNavigatorMixin functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.mixin = DatasheetNavigatorMixin()

    def test_get_datasheet_from_dict(self):
        """Test getting datasheet from dictionary."""
        source = {
            'structures': {'path': '/path/to/structures.csv', 'columns': ['id', 'file']},
            'sequences': {'path': '/path/to/sequences.csv', 'columns': ['id', 'seq']}
        }

        ds = self.mixin.get_datasheet(source, 'structures')
        self.assertEqual(ds['path'], '/path/to/structures.csv')

    def test_get_datasheet_with_fallbacks(self):
        """Test fallback behavior when primary name not found."""
        source = {
            'main': {'path': '/path/to/main.csv'}
        }

        # Request 'structures' but fall back to 'main'
        ds = self.mixin.get_datasheet(source, 'structures', fallback_names=['main'])
        self.assertEqual(ds['path'], '/path/to/main.csv')

    def test_get_datasheet_from_datasheet_container(self):
        """Test getting datasheet from DatasheetContainer."""
        ds_info = DatasheetInfo(
            name='structures',
            path='/path/to/structures.csv',
            columns=['id', 'file'],
            description='Test datasheet'
        )
        container = DatasheetContainer({'structures': ds_info})

        ds = self.mixin.get_datasheet(container, 'structures')
        self.assertEqual(ds.path, '/path/to/structures.csv')
        self.assertEqual(ds.name, 'structures')

    def test_get_datasheet_from_list(self):
        """Test getting datasheet from legacy list format."""
        source = ['/path/to/main.csv', '/path/to/other.csv']

        ds = self.mixin.get_datasheet(source, 'main')
        self.assertEqual(ds, '/path/to/main.csv')

    def test_get_datasheet_path(self):
        """Test extracting path from various datasheet formats."""
        # Test with DatasheetInfo
        ds_info = DatasheetInfo(
            name='test',
            path='/path/to/test.csv',
            columns=['id']
        )
        path = self.mixin.get_datasheet_path({'test': ds_info}, 'test')
        self.assertEqual(path, '/path/to/test.csv')

        # Test with dict format
        source = {'test': {'path': '/path/to/test2.csv'}}
        path = self.mixin.get_datasheet_path(source, 'test')
        self.assertEqual(path, '/path/to/test2.csv')

        # Test with string
        source = {'test': '/path/to/test3.csv'}
        path = self.mixin.get_datasheet_path(source, 'test')
        self.assertEqual(path, '/path/to/test3.csv')

    def test_get_all_datasheets_from_container(self):
        """Test getting all datasheets from container."""
        ds1 = DatasheetInfo(name='ds1', path='/path1.csv', columns=['id'])
        ds2 = DatasheetInfo(name='ds2', path='/path2.csv', columns=['id'])
        container = DatasheetContainer({'ds1': ds1, 'ds2': ds2})

        all_ds = self.mixin.get_all_datasheets(container)
        self.assertEqual(len(all_ds), 2)
        self.assertIn('ds1', all_ds)
        self.assertIn('ds2', all_ds)

    def test_datasheet_exists(self):
        """Test checking datasheet existence."""
        source = {
            'structures': {'path': '/path/to/structures.csv'}
        }

        self.assertTrue(self.mixin.datasheet_exists(source, 'structures'))
        self.assertFalse(self.mixin.datasheet_exists(source, 'nonexistent'))

    def test_get_datasheet_raises_on_missing(self):
        """Test that missing datasheet raises ValueError."""
        source = {'other': '/path/to/other.csv'}

        with self.assertRaises(ValueError) as context:
            self.mixin.get_datasheet(source, 'structures', fallback_names=[])

        self.assertIn('not found', str(context.exception))


class TestFilePathDescriptor(unittest.TestCase):
    """Test FilePathDescriptor functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def test_simple_filename(self):
        """Test descriptor with simple filename."""
        class TestTool:
            TOOL_NAME = "TestTool"
            output_folder = self.temp_dir
            test_file = FilePathDescriptor("test.csv")

        tool = TestTool()
        expected = os.path.join(self.temp_dir, "test.csv")
        self.assertEqual(tool.test_file, expected)

    def test_template_with_tool_name(self):
        """Test descriptor with template placeholder."""
        class TestTool:
            TOOL_NAME = "MyTool"
            output_folder = self.temp_dir
            test_file = FilePathDescriptor("{tool_name}_results.csv")

        tool = TestTool()
        expected = os.path.join(self.temp_dir, "MyTool_results.csv")
        self.assertEqual(tool.test_file, expected)

    def test_callable_template(self):
        """Test descriptor with callable template."""
        class TestTool:
            TOOL_NAME = "TestTool"
            output_folder = self.temp_dir
            name = "job123"
            test_file = FilePathDescriptor(lambda self: f"{self.name}_output.csv")

        tool = TestTool()
        expected = os.path.join(self.temp_dir, "job123_output.csv")
        self.assertEqual(tool.test_file, expected)

    def test_custom_folder_key(self):
        """Test descriptor with custom folder key."""
        class TestTool:
            TOOL_NAME = "TestTool"
            output_folder = self.temp_dir
            folders = {"HelpScripts": "/path/to/scripts"}
            test_script = FilePathDescriptor("script.py", folder_key="HelpScripts")

        tool = TestTool()
        expected = os.path.join("/path/to/scripts", "script.py")
        self.assertEqual(tool.test_script, expected)

    def test_caching(self):
        """Test that paths are cached after first access."""
        class TestTool:
            TOOL_NAME = "TestTool"
            output_folder = self.temp_dir
            test_file = FilePathDescriptor("test.csv")

        tool = TestTool()
        path1 = tool.test_file

        # Change output_folder
        tool.output_folder = "/different/path"

        # Path should still be cached to original value
        path2 = tool.test_file
        self.assertEqual(path1, path2)

    def test_manual_override(self):
        """Test that paths can be manually overridden."""
        class TestTool:
            TOOL_NAME = "TestTool"
            output_folder = self.temp_dir
            test_file = FilePathDescriptor("test.csv")

        tool = TestTool()
        original = tool.test_file

        # Manually override
        tool.test_file = "/custom/path/file.csv"
        self.assertEqual(tool.test_file, "/custom/path/file.csv")
        self.assertNotEqual(tool.test_file, original)


class TestResolvedInput(unittest.TestCase):
    """Test ResolvedInput container."""

    def test_resolved_input_creation(self):
        """Test creating ResolvedInput object."""
        files = ['/path/to/file1.pdb', '/path/to/file2.pdb']
        ids = ['file1', 'file2']
        datasheets = {'main': '/path/to/datasheet.csv'}
        source = Mock()

        resolved = ResolvedInput(files, ids, datasheets, source)

        self.assertEqual(resolved.files, files)
        self.assertEqual(resolved.ids, ids)
        self.assertEqual(resolved.datasheets, datasheets)
        self.assertEqual(resolved.source, source)

    def test_resolved_input_repr(self):
        """Test ResolvedInput string representation."""
        resolved = ResolvedInput(
            files=['file1.pdb', 'file2.pdb'],
            ids=['id1', 'id2'],
            datasheets=None,
            source="test"
        )

        repr_str = repr(resolved)
        self.assertIn('2', repr_str)  # Should show file count
        self.assertIn('ResolvedInput', repr_str)


if __name__ == '__main__':
    unittest.main()
