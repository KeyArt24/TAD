"""
Custom toolbar module for file operations.
"""

from PyQt6.QtWidgets import QToolBar, QFileDialog
from PyQt6.QtGui import QAction


class TBar(QToolBar):
    """Custom toolbar for file open and save operations."""

    def __init__(self):
        super().__init__()

        self.button_open = QAction('ОТКРЫТЬ ФАЙЛ', self)
        self.button_save = QAction('СОХРАНИТЬ ДАННЫЕ', self)
        self.button_information = QAction('ИНФОРМАЦИЯ О ПРОГРАММЕ', self)
        self.addAction(self.button_open)
        self.addAction(self.button_save)
        self.addAction(self.button_information)

        self.setStyleSheet('background-color: lightgray; border-radius: 4px;')

    def save_data(self, data=''):
        """Save data to a text file."""
        path = QFileDialog().getSaveFileName(
            None, "Сохранить файл", "", "Text Files (*.txt);;All Files (*)"
        )
        if len(path[0]) > 0:
            with open(path[0], 'w', encoding='utf-8') as file:
                file.write(data)

    def open_file(self):
        """Open and read a text file."""
        window = QFileDialog().getOpenFileNames()
        if len(window[0]) > 0 and window[0][0].endswith('.txt'):
            with open(window[0][0], 'r', encoding='utf-8', errors='ignore') as file:
                data = file.read()
                return data
        return False