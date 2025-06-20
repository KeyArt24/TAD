import sys
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from PyQt6.QtWidgets import (QApplication, QWidget, QMainWindow, QPushButton, QHBoxLayout, QVBoxLayout,
                             QLabel, QSizePolicy, QTextEdit, QLineEdit, QScrollArea, QCheckBox, QSlider, QTabWidget, QTextBrowser)
from PyQt6.QtCore import Qt, QThread
from PyQt6.QtGui import QFont, QIcon
from tm_oligo import deltaS_DNA, deltaH_DNA, deltaG_DNA, temp_DNA_melt, GC_features, DnaFraction, dimers_analyze, middles
from tool_bar import TBar
import numpy
import matplotlib
from threading import Thread
matplotlib.use('QtAgg')


list_samples = {}
dif_samples = {}
seq_samples = {}
count_samples = 0
temp = numpy.array([i * 0.1 for i in range(0, 1000)]) + 278.15
temp_plot = numpy.array([i * 0.1 for i in range(0, 1000)])
temp_dif = numpy.array(middles(temp)) - 278.15
concentration_DNA = 0.25
concentration_K = 50
concentration_Mg = 3


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.temp = []
        self.data = []
        self.name_x = ''
        self.name_y = ''

        fig = Figure(figsize=(width, height), dpi=dpi,
                     facecolor='lightgray', tight_layout=True)
        self.axes = fig.add_subplot()
        self.axes.plot(self.temp, self.data)
        self.axes.grid(which='both')
        self.axes.set_xlabel(self.name_x)
        self.axes.tick_params(axis='x', rotation=30)
        super(MplCanvas, self).__init__(fig)


class Sample(QLabel):
    def __init__(self):
        super().__init__()

        # виджеты
        self.box = QCheckBox()
        self.index = QLabel('0')
        self.label = QLineEdit(f'SAMPLE {w.index}')
        self.input = QLineEdit()
        self.label_features = QLineEdit(
            'Свойства: длина в нк. = 0, GC = 0%, deltaS = 0, deltaH = 0, deltaG = 0, Tm = 0 °C, Ta ~ 0 °C')
        self.button_del = QPushButton(' УДАЛИТЬ ')
        self.button_homo_dimer = QPushButton(' Гомодимеры ')
        self.button_hetero_dimer = QPushButton(' Гетеродимеры ')
        self.lay = QHBoxLayout()
        self.lay_sample = QVBoxLayout()
        self.in_lay_sample = QHBoxLayout()

        # свойства виджетов
        self.label.setFixedWidth(100)
        self.label.setFixedHeight(30)
        self.label.setStyleSheet("background-color: white")
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
        print(w.index)
        try:
            del list_samples[self.index.text()]
            del seq_samples[self.index.text()]
        except KeyError:
            pass

    def sample_features(self):
        concentration_DNA = w.data_conc_DNA.text()
        concentration_K = w.data_conc_K.text()
        concentration_Mg = w.data_conc_Mg.text()
        try:
            GC = GC_features(self.input.text())
            S = round(deltaS_DNA(self.input.text()))
            G = round(deltaG_DNA(self.input.text()))
            H = round(deltaH_DNA(self.input.text()))*1000
            Tm = temp_DNA_melt(self.input.text(), float(concentration_DNA), float(
                concentration_K), float(concentration_Mg))
            self.label_features.setText(
                f'Свойства: длина в нк. = {len(self.input.text())}, GC = {GC}%, deltaS = {S}, deltaH = {H}, deltaG = {G}, Tm = {Tm} °C, Ta ~ 0 °C')
        except KeyError:
            self.label_features.setText(
                'Неправильные символы или наличие пробелов в последовательности')
            self.setStyleSheet("background-color: rgb(227,38,54)")
        except ValueError:
            self.label_features.setText(
                'Значение концентрации ионов и ДНК должно быть больше 0')
        except ZeroDivisionError:
            pass

    def sample_seq(self):
        concentration_DNA = w.data_conc_DNA.text()
        concentration_K = w.data_conc_K.text()
        concentration_Mg = w.data_conc_Mg.text()
        try:
            melt_data = DnaFraction(float(concentration_DNA)/1E6, temp, deltaS_DNA(self.input.text(
            )), deltaH_DNA(self.input.text())*1E3, float(concentration_K), float(concentration_Mg))
            dif_data = numpy.diff(melt_data)
            list_samples.update(
                {self.index.text(): [self.label.text(), melt_data]})
            dif_samples.update(
                {self.index.text(): [self.label.text(), numpy.abs(dif_data)]})
            seq_samples.update(
                {self.index.text(): [self.label.text(), self.input.text()]})
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
                del list_samples[self.index.text()]
                del seq_samples[self.index.text()]
                del dif_samples[self.index.text()]
            except KeyError:
                pass

    def homodimer_analyze(self):
        w.tab_lay.setCurrentIndex(2)
        w.widget_03.setText(f"{self.label.text()}\n")
        result = sorted(dimers_analyze(self.input.text(), self.input.text(
        )), key=lambda x: x[1].count('I'), reverse=True)
        for line in result:
            for string in line:
                w.widget_03.append(str(string))
            w.widget_03.append('---------------------------')
            w.widget_03.append('')

    def heterodimer_analyze(self):
        w.tab_lay.setCurrentIndex(2)
        w.widget_03.setText(f"{self.label.text()}\n")
        result = []
        if len(seq_samples) > 0:
            for key, value in seq_samples.items():
                if self.label.text() != value[0]:
                    result.append([f'\n{self.label.text()} AND {value[0]}\n'])
                    result += sorted(dimers_analyze(self.input.text(), value[1]),
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
        try:
            self.sc.axes.cla()
            if len(self.sc.data) > 0:
                for key, values in self.sc.data.items():
                    self.sc.axes.plot(self.sc.temp, values[1], label=values[0])
                    self.sc.axes.legend()
            else:
                self.sc.axes.plot(
                    self.sc.temp, [0 for i in range(len(self.sc.temp))])

            self.sc.axes.set_xlabel(self.sc.name_x)
            self.sc.axes.set_ylabel(self.sc.name_y)
            self.sc.axes.grid(which='both', linestyle='--', drawstyle='steps')
            self.sc.axes.set_xticks(numpy.arange(0, 101, step=5))
            self.sc.draw()
        except RuntimeError:
            pass


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.resize(1024, 860)
        self.setWindowIcon(QIcon(r"D:\06_Python\TPD ver.1.03\HRM.png"))

        # список виджетов
        self.label_conc_DNA = QLabel('Концентрация ДНК (мкМ)')
        self.data_conc_DNA = QLineEdit(str(concentration_DNA))
        self.label_conc_K = QLabel('Концентрация моновалентных инов (мМ)')
        self.data_conc_K = QLineEdit(str(concentration_K))
        self.label_conc_Mg = QLabel('Концентрация Мg (мМ)')
        self.data_conc_Mg = QLineEdit(str(concentration_Mg))

        self.slider_conc = QSlider(Qt.Orientation.Horizontal)
        self.slider_conc.setValue(int(concentration_DNA*100))
        self.slider_conc.setMinimum(1)
        self.slider_conc.setMaximum(100)
        self.slider_conc.setSingleStep(1)
        self.slider_conc.setTickPosition(QSlider.TickPosition.TicksBothSides)

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
        widget_00 = QLabel()
        widget_01 = QVBoxLayout()
        self.widget_03 = QTextEdit()
        self.df_dx = QLabel()
        self.df_dx_sc2 = QVBoxLayout()
        widget_04 = QPushButton()
        button_add = QPushButton('Добавить образец')
        main_lay = QVBoxLayout()
        lay_UP = QHBoxLayout()
        self.scrollArea = QScrollArea()
        self.scrollWidget = QWidget()
        self.lay_DW = QVBoxLayout()
        self.lay_samples = QVBoxLayout()
        self.tab_lay = QTabWidget()

        self.index = 0

        # plot tm

        self.sc = MplCanvas(self, width=5, height=4, dpi=75)
        self.sc.temp = temp_plot
        self.sc.data = list_samples
        self.toolbar = NavigationToolbar2QT(self.sc)
        self.sc_thread = Thread(mainwindow=self.sc)

        # plot tm df/dx

        self.sc2 = MplCanvas(self, width=5, height=4, dpi=75)
        self.sc2.temp = temp_dif
        self.sc2.data = dif_samples
        self.toolbar2 = NavigationToolbar2QT(self.sc2)
        self.sc_thread2 = Thread(mainwindow=self.sc2)

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
        self.df_dx.setStyleSheet("background-color: lightgray")
        widget_04.setStyleSheet("background-color: lightgray")
        self.scrollArea.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setStyleSheet("""background-color: rgb(0, 32, 78)""")
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
        self.slider_conc.valueChanged.connect(self.change_conc)
        self.data_conc_DNA.textChanged.connect(self.changed_slide_conc)
        self.slider_K.valueChanged.connect(self.changed_K)
        self.data_conc_K.textChanged.connect(self.changed_slide_K)
        self.slider_Mg.valueChanged.connect(self.changed_Mg)
        self.data_conc_Mg.textChanged.connect(self.changed_slide_Mg)

        self.widget_05 = QWidget()

        # подключение виджетов
        lay_conditions.addWidget(self.label_conc_DNA)
        lay_conditions.addWidget(self.data_conc_DNA)
        lay_conditions.addWidget(self.slider_conc)
        lay_conditions.addWidget(self.label_conc_K)
        lay_conditions.addWidget(self.data_conc_K)
        lay_conditions.addWidget(self.slider_K)
        lay_conditions.addWidget(self.label_conc_Mg)
        lay_conditions.addWidget(self.data_conc_Mg)
        lay_conditions.addWidget(self.slider_Mg)
        widget_04.setLayout(lay_conditions)
        widget_00.setLayout(widget_01)
        widget_01.addWidget(self.sc)
        widget_01.addWidget(self.toolbar)
        self.df_dx_sc2.addWidget(self.sc2)
        self.df_dx_sc2.addWidget(self.toolbar2)
        self.df_dx.setLayout(self.df_dx_sc2)

        self.tab_lay.setMovable(True)
        self.tab_lay.setStyleSheet(
            """
            QTabBar::tab-bar {
            alignment: center}
            QTabBar::tab  {background-color:rgb(220, 220, 220); border-radius: 4px; height: 20px; width: 150px; margin: 2px;}
            QTabWidget::pane {position; absolute}""")

        self.tab_lay.addTab(widget_00, 'График температуры')
        self.tab_lay.addTab(self.df_dx, 'Функция dF/dx')
        self.tab_lay.addTab(self.widget_03, 'Димеры')

        lay_UP.addWidget(self.tab_lay, 4)
        lay_UP.addWidget(widget_04, 1)
        self.scrollWidget.setLayout(self.lay_samples)
        self.scrollArea.setWidget(self.scrollWidget)
        self.lay_DW.addWidget(button_add)
        self.lay_DW.addWidget(self.scrollArea)

        self.bar = TBar()
        self.bar.button_save.triggered.connect(self.save_data)
        self.bar.button_open.triggered.connect(self.open_data)
        self.bar.button_information.triggered.connect(self.information)

        self.lay_bar = QVBoxLayout()
        self.lay_bar.addWidget(self.bar)

        main_lay.addLayout(self.lay_bar)
        main_lay.addLayout(lay_UP, 2)
        main_lay.addLayout(self.lay_DW, 2)

        # главное окно
        container = QWidget()
        container.setLayout(main_lay)
        self.setStyleSheet("background-color: rgb(0, 32, 78)")
        self.setCentralWidget(container)
        self.setWindowTitle(
            'TAD: Termodynamic Analysis of DNA (KeyArt) v.1.1.3')

        # потоки

    # функции

    def save_data(self):
        self.bar.save_data(data=str(seq_samples) +
                           '\n' + self.widget_03.toPlainText())

    def open_data(self):
        data = self.bar.open_file()
        if data != False:
            for line in data.split('\n'):
                if len(line.split()) == 2:
                    self.sample_add()
                    self.sample.input.setText(line.split()[1])

    def information(self):
        self.new_window = QWidget()
        self.new_window.resize(600, 400)
        self.new_window.setWindowTitle('О программе')

        layout = QVBoxLayout()
        label = QLabel(
            "Здравствуйте!\nПрограмма предназначена для анализа ДНК последовательностей длиной до 100 нуклеотидов.\n В основе математических вычислений лежит модель ближайших соседей.\n"
            "Hatim T. Allawi and John SantaLucia, Jr. Thermodynamics and NMR of Internal G‚T Mismatches in DNA. Biochemistry 1997, 36, 10581-10594.\n\n"
            "Расчёт темпераутры плавления производиться по следующей формуле:\n"
            "Tm = (deltaH/(deltaS+1.987*log(CtDNA/1000000))) + (16.6*log10(salt/(1.0+0.7*salt))) - 273.15\n\n"
            "Расчёт терпмодинамического профиля по формуле:\n"
            "salt = (CtK/1000) + 4 * (CtMg/1000)**0.5\n"
            "CtKeq = Ct * numpy.exp(DeltaS/R - DeltaH/(R*T-16.6*log10(salt/(1.0+0.7*salt))))\n"
            "f = (1 + CtKeq - numpy.sqrt(1 + 2*CtKeq)) / CtKeq\n\n"
            "Версия 1.1.3 от 18 июня 2025")
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(label)
        self.new_window.setLayout(layout)
        self.new_window.show()

    def sample_add(self):
        self.index += 1
        self.sample = Sample()
        self.sample.index.setText(f'{self.index}')
        self.lay_samples.addWidget(self.sample)
        self.sample.box.clicked.connect(self.thread_update_plot)
        self.sample.box.checkStateChanged.connect(self.thread_update_plot)
        self.sample.box.checkStateChanged.connect(
            self.sample.is_sample_selected)

        self.sample.label.textChanged.connect(self.sample.sample_seq)
        self.sample.label.textChanged.connect(self.thread_update_plot)
        self.sample.input.textChanged.connect(self.thread_update_plot)
        self.sample.button_del.clicked.connect(self.thread_update_plot)

        self.data_conc_DNA.textChanged.connect(self.sample.sample_features)
        self.data_conc_DNA.textChanged.connect(self.sample.is_sample_selected)
        self.data_conc_DNA.textChanged.connect(self.thread_update_plot)

        self.data_conc_K.textChanged.connect(self.sample.sample_features)
        self.data_conc_K.textChanged.connect(self.sample.is_sample_selected)
        self.data_conc_K.textChanged.connect(self.thread_update_plot)

        self.data_conc_Mg.textChanged.connect(self.sample.sample_features)
        self.data_conc_Mg.textChanged.connect(self.sample.is_sample_selected)
        self.data_conc_Mg.textChanged.connect(self.thread_update_plot)

    def thread_update_plot(self):
        self.sc_thread.start()
        self.sc_thread2.start()

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

    def change_conc(self):
        concentration_DNA = self.slider_conc.value()/100
        self.data_conc_DNA.setText(str(concentration_DNA))

    def changed_slide_conc(self):
        value = self.data_conc_DNA.text()
        if value != '':
            value = float(value)*100
            self.slider_conc.setValue(round(value))
        else:
            self.slider_conc.setValue(10)

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
app.setStyle('Fusion')
app.exec()
