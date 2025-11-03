"""
Mixins for BioPipelines tool refactoring.

Provides reusable components for input handling, table navigation,
and other common patterns across tools.
"""

from .input_handler import InputHandlerMixin, ResolvedInput
from .table_navigator import TableNavigatorMixin
from .file_path_manager import FilePathDescriptor

__all__ = [
    'InputHandlerMixin',
    'ResolvedInput',
    'TableNavigatorMixin',
    'FilePathDescriptor',
]
