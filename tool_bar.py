from PyQt6.QtWidgets import QToolBar, QFileDialog
from PyQt6.QtGui import QAction

class TBar(QToolBar):
    def __init__(self):
        super().__init__()

        self.button_open = QAction('ОТКРЫТЬ ФАЙЛ', self)
        self.button_save = QAction('СОХРАНИТЬ ДАННЫЕ', self)
        self.addAction(self.button_open)
        self.addAction(self.button_save)

        self.button_save.triggered.connect(self.save_data)

        self.setStyleSheet('background-color: lightgray; border-radius: 4px;')



    def save_data(self):
        window = QFileDialog().getSaveFileName()
        print(window)