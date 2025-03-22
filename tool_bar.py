from PyQt6.QtWidgets import QToolBar




class TBar(QToolBar):
    def __init__(self):
        super().__init__()

        self.addAction('Новый файл')