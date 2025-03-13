import sys
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from PyQt6.QtWidgets import (QApplication, QWidget, QMainWindow, QPushButton, QHBoxLayout, QVBoxLayout,
                             QLabel, QSizePolicy, QTextEdit, QLineEdit, QScrollArea, QCheckBox, QSlider, QTabWidget)
from PyQt6.QtCore import Qt, QThread
from PyQt6.QtGui import QFont
from tm_oligo import deltaS_DNA, deltaH_DNA, temp_DNA_melt, DnaFraction, dimers_analyze
import numpy
import matplotlib
from threading import Thread
matplotlib.use('QtAgg')

list_samples = {}
seq_samples = {}
temp = numpy.array(range(0, 100)) + 273.15
temp_plot = numpy.array(range(0, 100))
concentration_DNA = 0.4
concentration_K = 50
concentration_Mg = 3


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi,
                     facecolor='lightgray', tight_layout=True)
        self.axes = fig.add_subplot()
        self.axes.plot(temp_plot, [1]*100)
        self.axes.grid(which='both')
        self.axes.set_xlabel('Температура в °C')
        self.axes.tick_params(axis='x', rotation=30)
        super(MplCanvas, self).__init__(fig)


class Sample(QLabel):
    def __init__(self):
        super().__init__()

        # виджеты
        self.box = QCheckBox()
        self.label = QLabel('SAMPLE')
        self.input = QLineEdit()
        self.label_features = QLabel(
            'Свойства: длина в нк. = 0, deltaS = 0, deltaH = 0, Tm = 0 °C, Ta ~ 0 °C')
        self.button_del = QPushButton(' УДАЛИТЬ ')
        self.button_homo_dimer = QPushButton(' Гомодимеры ')
        self.button_hetero_dimer = QPushButton(' Гетеродимеры ')
        self.lay = QHBoxLayout()
        self.lay_sample = QVBoxLayout()
        self.in_lay_sample = QHBoxLayout()

        # свойства виджетов
        self.button_del.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        self.button_homo_dimer.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        self.button_hetero_dimer.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        self.input.setStyleSheet(
            "border-style: solid; border-width: 1px; border-color: white; border-radius: 5px; background-color: white")
        self.button_del.setStyleSheet("background-color: white; color: black")
        self.button_homo_dimer.setStyleSheet(
            "background-color: white; color: black")
        self.button_hetero_dimer.setStyleSheet(
            "background-color: white; color: black")
        self.lay.setAlignment(Qt.AlignmentFlag.AlignTop)

        # подключение функций
        self.button_del.clicked.connect(self.sample_del)
        self.input.textChanged.connect(self.is_sample_selected)
        self.input.textChanged.connect(self.sample_features)
        self.button_homo_dimer.clicked.connect(self.homodimer_analyze)
        self.button_hetero_dimer.clicked.connect(self.heterodimer_analyze)

        # добавление виджетов
        self.in_lay_sample.addWidget(self.input)
        self.lay_sample.addLayout(self.in_lay_sample)
        self.lay_sample.addWidget(self.label_features)
        self.lay.addWidget(self.box)
        self.lay.addWidget(self.label)
        self.lay.addLayout(self.lay_sample)
        self.lay.addWidget(self.button_homo_dimer)
        self.lay.addWidget(self.button_hetero_dimer)
        self.lay.addWidget(self.button_del)

        # настройка главного окна
        self.setLayout(self.lay)
        self.setStyleSheet(
            "background-color: rgb(220,220,220); border-radius: 5px")
        self.setFixedHeight(60)

    def sample_del(self):
        self.box.setCheckState(Qt.CheckState.Unchecked)
        self.close()
        try:
            del list_samples[self.label.text()]
            del seq_samples[self.label.text()]
        except KeyError:
            pass

    def sample_features(self):
        concentration_DNA = w.data_conc_DNA.text()
        concentration_K = w.data_conc_K.text()
        concentration_Mg = w.data_conc_Mg.text()
        try:
            S = round(deltaS_DNA(self.input.text()))
            H = round(deltaH_DNA(self.input.text()))*1000
            Tm = temp_DNA_melt(self.input.text(), float(concentration_DNA), float(
                concentration_K), float(concentration_Mg))
            self.label_features.setText(f'Свойства: длина в нк. = {len(
                self.input.text())}, deltaS = {S}, deltaH = {H}, Tm = {Tm} °C, Ta ~ 0 °C')
        except KeyError:
            self.label_features.setText(
                'Неправильные символы или наличие пробелов в последовательности')
            self.setStyleSheet("background-color: rgb(227,38,54)")
        except ValueError:
            self.label_features.setText(
                'Значение концентрации ионов и ДНК должно быть больше 0')

    def sample_seq(self):
        concentration_K = w.data_conc_K.text()
        concentration_Mg = w.data_conc_Mg.text()
        try:
            melt_data = DnaFraction(concentration_DNA/100000, temp, deltaS_DNA(self.input.text(
            )), deltaH_DNA(self.input.text())*1000, float(concentration_K), float(concentration_Mg))
            list_samples.update({self.label.text(): melt_data})
            seq_samples.update({self.label.text(): self.input.text()})
        except KeyError:
            self.label_features.setText(
                'Неправильные символы или наличие пробелов в последовательности')
        except ValueError:
            pass

    def is_sample_selected(self):
        if self.box.isChecked() == True:
            self.setStyleSheet(
                'background-color: rgb(144, 238, 144); border-radius: 5px')
            self.sample_seq()
        else:
            self.setStyleSheet(
                "background-color: lightgray; border-radius: 5px")
            try:
                del list_samples[self.label.text()]
                del seq_samples[self.label.text()]
            except KeyError:
                pass

    def homodimer_analyze(self):
        w.widget_03.setText(f"{self.label.text()}\n")
        result = sorted(dimers_analyze(self.input.text(), self.input.text(
        )), key=lambda x: x[1].count('I'), reverse=True)
        for string in result:
            for str in string:
                w.widget_03.append(str)
            w.widget_03.append('---------------------------')
            w.widget_03.append('')

    def heterodimer_analyze(self):
        w.widget_03.setText(f"{self.label.text()}\n")
        result = []
        if len(seq_samples) > 0:
            for key, seq in seq_samples.items():
                if self.label.text() != key:
                    result.append([f'\n{self.label.text()} AND {key}\n'])
                    result += sorted(dimers_analyze(self.input.text(), seq),
                                     key=lambda x: x[1].count('I'), reverse=True)
            for string in result:
                for str in string:
                    w.widget_03.append(str)
                w.widget_03.append('---------------------------')


class Thread(QThread):
    def __init__(self, mainwindow):
        super().__init__()
        self.sc = mainwindow

    def run(self):
        self.sc.axes.cla()
        if len(list_samples) > 0:
            for key, values in list_samples.items():
                self.sc.axes.plot(temp_plot, values, label=key)
                self.sc.axes.legend()
        else:
            self.sc.axes.plot(temp_plot, [0 for i in range(100)])

        self.sc.axes.set_xlabel('Температура в °C')
        self.sc.axes.set_ylabel('Фракция двухцепочечной ДНК')
        self.sc.axes.grid(which='both', linestyle='--', drawstyle='steps')
        self.sc.axes.set_xticks(numpy.arange(0, 101, step=5))
        self.sc.axes.axhline(y=0.5, color='r')
        self.sc.draw()


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.resize(800, 600)

        # список виджетов
        self.label_conc_DNA = QLabel('Концентрация ДНК (мкМ)')
        self.data_conc_DNA = QLineEdit(str(concentration_DNA))
        self.label_conc_K = QLabel('Концентрация моновалентных инов (мМ)')
        self.data_conc_K = QLineEdit(str(concentration_K))
        self.label_conc_Mg = QLabel('Концентрация Мg (мМ)')
        self.data_conc_Mg = QLineEdit(str(concentration_Mg))

        self.slider_K = QSlider(Qt.Orientation.Horizontal)
        self.slider_K.setValue(concentration_K)
        self.slider_K.setMinimum(1)
        self.slider_K.setMaximum(100)
        self.slider_K.setSingleStep(1)
        self.slider_K.setTickPosition(QSlider.TickPosition.TicksBothSides)

        self.slider_Mg = QSlider(Qt.Orientation.Horizontal)
        self.slider_Mg.setValue(concentration_Mg*10)
        self.slider_Mg.setMinimum(10)
        self.slider_Mg.setMaximum(100)
        self.slider_Mg.setSingleStep(1)
        self.slider_Mg.setTickPosition(QSlider.TickPosition.TicksBothSides)

        lay_conditions = QVBoxLayout()
        widget_00 = QPushButton()
        widget_01 = QVBoxLayout()
        self.widget_03 = QTextEdit()
        widget_04 = QPushButton()
        button_add = QPushButton('Добавить образец')
        main_lay = QVBoxLayout()
        lay_UP = QHBoxLayout()
        self.scrollArea = QScrollArea()
        self.scrollWidget = QWidget()
        self.lay_DW = QVBoxLayout()
        self.lay_samples = QVBoxLayout()
        self.tab_lay = QTabWidget()

        # plot tm

        self.sc = MplCanvas(self, width=5, height=4, dpi=75)
        self.toolbar = NavigationToolbar2QT(self.sc)
        self.sc_thread = Thread(mainwindow=self.sc)

        # найтрока расположения и стиль виджетов
        self.data_conc_Mg.setStyleSheet(
            "background-color: white; border-style: solid; border-radius: 5px")
        self.data_conc_DNA.setStyleSheet(
            "background-color: white; border-style: solid; border-radius: 5px")
        self.data_conc_K.setStyleSheet(
            "background-color: white; border-style: solid; border-radius: 5px")
        lay_conditions.setAlignment(Qt.AlignmentFlag.AlignTop)
        widget_00.setSizePolicy(
            QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.widget_03.setSizePolicy(
            QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.widget_03.setStyleSheet(
            "background: white; background-color: white; border-style: solid; border-radius: 5px")
        self.widget_03.setFont(QFont('Courier new', 8))
        widget_04.setSizePolicy(
            QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.lay_DW.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.lay_samples.setAlignment(Qt.AlignmentFlag.AlignTop)
        button_add.setFixedHeight(50)
        button_add.setStyleSheet("background-color: lightgray")
        widget_00.setStyleSheet("background-color: lightgray")
        widget_04.setStyleSheet("background-color: lightgray")
        self.scrollArea.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setStyleSheet("""background-color: gray""")
        self.scrollArea.verticalScrollBar().setStyleSheet("""
    QScrollBar:vertical
    {
        background-color: gray;
        width: 15px;
        margin: 15px 3px 15px 3px;
        border: 1px transparent #2A2929;
        border-radius: 4px;
    }

    QScrollBar::handle:vertical
    {
        background-color: lightgray;         /* #605F5F; */
        min-height: 5px;
        border-radius: 4px;
    }

    QScrollBar::sub-line:vertical
    {
        margin: 3px 0px 3px 0px;
        border-image: url(:/qss_icons/rc/up_arrow_disabled.png);
        height: 10px;
        width: 10px;
        subcontrol-position: top;

    }

    QScrollBar::add-line:vertical
    {
        margin: 3px 0px 3px 0px;
        border-image: url(:/qss_icons/rc/down_arrow_disabled.png);
        height: 10px;
        width: 10px;
        subcontrol-position: bottom;
        subcontrol-origin: margin;
    }

    QScrollBar::up-arrow:vertical, QScrollBar::down-arrow:vertical
    {
        background: none;
    }

    QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical
    {
        background: none;
    }
""")

        # подключение функций
        button_add.clicked.connect(self.sample_add)
        self.slider_K.valueChanged.connect(self.changed_K)
        self.data_conc_K.textChanged.connect(self.changed_slide_K)
        self.slider_Mg.valueChanged.connect(self.changed_Mg)
        self.data_conc_Mg.textChanged.connect(self.changed_slide_Mg)

        self.count = 0

        self.widget_05 = QWidget()

        # подключение виджетов
        lay_conditions.addWidget(self.label_conc_DNA)
        lay_conditions.addWidget(self.data_conc_DNA)
        lay_conditions.addWidget(self.label_conc_K)
        lay_conditions.addWidget(self.data_conc_K)
        lay_conditions.addWidget(self.slider_K)
        lay_conditions.addWidget(self.label_conc_Mg)
        lay_conditions.addWidget(self.data_conc_Mg)
        lay_conditions.addWidget(self.slider_Mg)
        widget_04.setLayout(lay_conditions)
        widget_00.setLayout(widget_01)
        widget_01.addWidget(self.sc)

        self.tab_lay.setMovable(True)
        self.tab_lay.setStyleSheet(
            """
            QTabBar::tab-bar {
            alignment: center}
            QTabBar::tab  {background-color:rgb(220, 220, 220); border-radius: 4px; height: 20px; width: 150px; margin: 2px;}
            QTabWidget::pane {position; absolute}""")
        self.tab_lay.addTab(widget_00, 'График температуры')
        self.tab_lay.addTab(self.widget_03, 'Димеры')
        lay_UP.addWidget(self.tab_lay, 4)
        lay_UP.addWidget(widget_04, 1)
        self.scrollWidget.setLayout(self.lay_samples)
        self.scrollArea.setWidget(self.scrollWidget)
        self.lay_DW.addWidget(button_add)
        self.lay_DW.addWidget(self.scrollArea)

        main_lay.addLayout(lay_UP, 2)
        main_lay.addLayout(self.lay_DW, 2)

        # главное окно
        container = QWidget()
        container.setLayout(main_lay)
        self.setStyleSheet("background-color: rgb(110, 110, 110)")
        self.setCentralWidget(container)
        self.setWindowTitle('TAD')

        # потоки

    # функции

    def sample_add(self):
        self.count += 1
        self.sample = Sample()
        self.sample.label.setText(f'Олигонуклеотид {self.count}')
        self.lay_samples.addWidget(self.sample)
        self.sample.box.clicked.connect(self.thread_update_plot)
        self.sample.box.checkStateChanged.connect(self.thread_update_plot)
        self.sample.box.checkStateChanged.connect(
            self.sample.is_sample_selected)
        self.sample.input.textChanged.connect(self.thread_update_plot)
        self.sample.button_del.clicked.connect(self.thread_update_plot)
        self.data_conc_DNA.textChanged.connect(self.sample.sample_features)

        self.data_conc_K.textChanged.connect(self.sample.sample_features)
        self.data_conc_K.textChanged.connect(self.sample.is_sample_selected)
        self.data_conc_K.textChanged.connect(self.thread_update_plot)

        self.data_conc_Mg.textChanged.connect(self.sample.sample_features)
        self.data_conc_Mg.textChanged.connect(self.sample.is_sample_selected)
        self.data_conc_Mg.textChanged.connect(self.thread_update_plot)

    def thread_update_plot(self):
        self.sc_thread.start()

    def dimer_analyze(self):
        self.widget_03.setText(self.sample.input.text)

    def changed_K(self):
        concentration_K = self.slider_K.value()
        self.data_conc_K.setText(str(concentration_K))

    def changed_slide_K(self):
        value = self.data_conc_K.text()
        if value != '':
            self.slider_K.setValue(int(value))
        else:
            self.slider_K.setValue(1)

    def changed_Mg(self):
        concentration_Mg = self.slider_Mg.value()/10
        self.data_conc_Mg.setText(str(concentration_Mg))

    def changed_slide_Mg(self):
        value = self.data_conc_Mg.text()
        if value != '':
            value = float(value)*10
            self.slider_Mg.setValue(round(value))
        else:
            self.slider_Mg.setValue(10)


app = QApplication(sys.argv)
w = MainWindow()
w.show()
app.exec()
