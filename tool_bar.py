from PyQt6.QtWidgets import QToolBar, QFileDialog
from PyQt6.QtGui import QAction
import os


class TBar(QToolBar):
    def __init__(self):
        super().__init__()

        self.button_open = QAction('ОТКРЫТЬ ФАЙЛ', self)
        self.button_save = QAction('СОХРАНИТЬ ДАННЫЕ', self)
        self.addAction(self.button_open)
        self.addAction(self.button_save)

        self.setStyleSheet('background-color: lightgray; border-radius: 4px;')

    def save_data(self, data=''):
        path = QFileDialog().getSaveFileName()
        if len(path[0]) > 0:
            with open(path[0], 'w', encoding='utf-8') as file:
                file.write(data)

    def open_file(self):
        window = QFileDialog().getOpenFileNames()
        print(window)
