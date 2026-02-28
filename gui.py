import sys
from concurrent.futures import ProcessPoolExecutor

import matplotlib
import numpy
import seaborn as sns

from PyQt6.QtGui import QFont, QIcon, QAction, QPalette, QColor
from PyQt6.QtCore import Qt, QThread, pyqtSignal, pyqtSlot
from PyQt6.QtWidgets import (QApplication, QWidget,
                             QMainWindow, QPushButton,
                             QHBoxLayout, QVBoxLayout,
                             QLabel, QSizePolicy,
                             QTextEdit, QLineEdit,
                             QScrollArea, QCheckBox,
                             QSlider, QTabWidget, QMessageBox)

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib import pyplot
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure
from matplotlib import style

from tool_bar import TBar
from tm_oligo import deltaS_DNA, deltaH_DNA, deltaG_DNA, temp_DNA_melt, GC_features, DnaFraction, dimers_analyze, loops_analyze
from style import scroll_vertical, scroll_horizontal, tab_lay, css_buttons_sample


# Глобальные параметры
pyplot.rcParams.update({'font.size': 8})
pyplot.rcParams['agg.path.chunksize'] = 10000
style.use('fast')
matplotlib.use('Qt5Agg')

tm_samples = {}
dif_samples = {}
seq_samples = {}
count_samples = 0
tm_resolution = 500
temp = numpy.array([i * 0.2 for i in range(0, tm_resolution)]) + 276.15
temp_plot = [i * 0.2 for i in range(0, tm_resolution)]
temp_dif = [i * 0.2 for i in range(0, tm_resolution-1)]

concentration_DNA = 0.25
concentration_K = 50
concentration_Mg = 3

palette = sns.color_palette('deep', 50)
version = '1.1.8'


# Создание графика
class MplComplexCanvas(QWidget):
    cursor_moved = pyqtSignal(float, float)
    def __init__(self, parent=None):
        super(MplComplexCanvas, self).__init__()

        self.cursor_text = None
        self.data = {}
        self.temp = []
        self.legend_data = []


        self.name_x = 'Температура °C'
        self.name_y = ''
        self.on_text = False

        self.main_fig = Figure(dpi = 100, frameon = False)
        self.legend_fig = Figure(dpi = 100)

        self.main_fig.subplots_adjust(left=0.1, right=0.99)

        self.axes = self.main_fig.add_subplot()
        self.axes.grid(which='both')
        self.axes.set_xlabel(self.name_x)
        self.axes.set_ylabel(self.name_y)

        self.base_x_scale = self.axes.get_xlim()
        self.base_y_scale = self.axes.get_ylim()

        self.main_fig.text(0.18, 0.25, 'Гибридизация', ha='center', va='top', fontsize=8)
        self.main_fig.text(0.85, 0.25, 'Денатурация', ha='center', va='top', fontsize=8)

        self.main_fig.patch.set_alpha(0.0)
        self.axes.set_facecolor('#E8E9EB')
        self.axes.spines['top'].set_visible(False)
        self.axes.spines['right'].set_visible(False)
        self.axes.spines['bottom'].set_visible(False)
        self.axes.spines['left'].set_visible(False)
        self.legend_fig.patch.set_alpha(0.0)


        self.main_canvas = FigureCanvasQTAgg(self.main_fig)
        self.legend_canvas = FigureCanvasQTAgg(self.legend_fig)


    def legend_build(self):
        self.legend_fig.legend(handles = self.legend_data, loc = 'upper left')

    def setup_cursor_tracking(self):
        """Setup mouse event handlers"""
        self.cursor_position = None
        self.cursor_text = None


        # Connect mouse events
        self.main_canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.main_canvas.mpl_connect('button_press_event', self.on_mouse_click)
        self.main_canvas.mpl_connect('scroll_event', self.on_mouse_scroll)

    def on_mouse_move(self, event):
        """Handle mouse movement"""
        if event.inaxes:
            x, y = event.xdata, event.ydata

            # Remove previous cursor text
            if self.cursor_text:
                self.cursor_text.remove()

            # Display cursor position
            self.cursor_text = self.axes.text(
                0.02, 0.95,
                f'{self.name_x}={x:.2f}, y={y:.2f}',
                transform=self.axes.transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5),
            )

            # Draw vertical and horizontal lines
            self.draw_cursor_lines(x, y)

            # Emit signal with position
            self.cursor_moved.emit(x, y)

            self.main_canvas.draw_idle()  # Update canvas

    def on_mouse_click(self, event):
        """Handle mouse clicks"""
        if event.inaxes and event.button == 1:  # Left click
            x, y = event.xdata, event.ydata

            # Add a marker at click position
            self.axes.plot(x, y, 'go', markersize=5, alpha=0.5)

            # Add text annotation
            self.axes.annotate(f'({x:.2f}, {y:.2f})',
                             xy=(x, y),
                             xytext=(5, 5),
                             textcoords='offset points',
                             fontsize=8,
                             bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
            self.main_canvas.draw_idle()

    def on_mouse_scroll(self, event):
        """Handle mouse scroll for zoom"""
        if event.inaxes:
            base_scale = 1.1
            cur_xlim = self.axes.get_xlim()
            cur_ylim = self.axes.get_ylim()

            xdata = event.xdata
            ydata = event.ydata

            if event.button == 'up':
                # Zoom in
                scale_factor = 1 / base_scale
            elif event.button == 'down':
                # Zoom out
                scale_factor = base_scale
            else:
                scale_factor = 1

            # Set new limits
            self.axes.set_xlim([xdata - (xdata - cur_xlim[0]) * scale_factor, xdata + (cur_xlim[1] - xdata) * scale_factor])
            self.axes.set_ylim([ydata - (ydata - cur_ylim[0]) * scale_factor, ydata + (cur_ylim[1] - ydata) * scale_factor])

            self.main_canvas.draw_idle()

    def draw_cursor_lines(self, x, y):
        """Draw cursor crosshair lines"""
        # Remove old lines if they exist
        if hasattr(self, 'vline'):
            for line in self.axes.get_children():
                if line == self.vline:
                    line.remove()
        if hasattr(self, 'hline'):
            for line in self.axes.get_children():
                if line == self.hline:
                    line.remove()

        # Draw vertical line
        self.vline = self.axes.axvline(x=x, color='red', alpha=0.2, linestyle='-')

        # Draw horizontal line
        self.hline = self.axes.axhline(y=y, color='red', alpha=0.2, linestyle='-')

    def set_axes_text(self):
        self.axes.text(0.05, 0.5, "Гибридизация",
                              style ="italic",fontsize = 8, transform=self.axes.transAxes,)
        self.axes.text(0.85, 0.5, "Денатурация",
                              style="italic", fontsize=8, transform=self.axes.transAxes)

# Карточка образца ДНК
class Sample(QLabel):
    def __init__(self):
        super().__init__()

        # виджеты
        self.box = QCheckBox()
        self.box.setObjectName("label_checkbox")
        self.index = QLabel('0')
        self.label = QLineEdit(f'SAMPLE {w.index_sample}')
        self.label.setObjectName('label_name')
        self.input = QLineEdit()
        self.input.setObjectName('label_seq')
        self.label_features = QLineEdit('Свойства: длина в нк. = 0, '
                                        'GC = 0%, deltaS = 0, '
                                        'deltaH = 0, deltaG = 0, '
                                        'Tm = 0 °C')
        self.label_features.setObjectName("label_features")
        self.label_conc = QLabel('мкМ')
        self.field_conc = QLineEdit(w.data_conc_DNA.text())
        self.button_get_complementary = QPushButton('Комплементарная цепь')
        self.button_del = QPushButton(' УДАЛИТЬ ')
        self.button_loops = QPushButton('   Петли   ')
        self.button_homo_dimer = QPushButton(' Гомодимеры ')
        self.button_hetero_dimer = QPushButton(' Гетеродимеры ')

        self.check_field_conc = self.field_conc.text() == w.data_conc_DNA.text()

        # Основной слой содержащий все остальные области и размещающися кнопки в Sample
        self.lay = QHBoxLayout()
        # Слой содержащий виджет последовательности и свойства
        self.lay_sample = QVBoxLayout()
        # Cлой содержащий последовательность
        self.in_lay_sample = QHBoxLayout()

        # свойства виджетов
        self.label.setFixedWidth(100)
        self.label.setFixedHeight(30)
        self.input.setFixedHeight(20)
        self.field_conc.setFixedWidth(50)
        self.field_conc.setFixedHeight(20)
        self.label_conc.setFixedWidth(30)
        self.button_get_complementary.setFixedHeight(20)

        self.label.setStyleSheet("background-color: white")
        self.field_conc.setStyleSheet("background-color: white")
        self.button_get_complementary.setStyleSheet("background-color: white")
        self.button_get_complementary.setFixedWidth(150)
        self.button_del.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        self.button_loops.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        self.button_homo_dimer.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        self.button_hetero_dimer.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        self.input.setStyleSheet("border-style: solid; border-width: 1px; border-color: white; border-radius: 5px; background-color: white")
        self.button_del.setStyleSheet("background-color: white; color: black")
        self.button_loops.setStyleSheet(css_buttons_sample)
        self.button_homo_dimer.setStyleSheet(css_buttons_sample)
        self.button_hetero_dimer.setStyleSheet(css_buttons_sample)
        self.lay.setAlignment(Qt.AlignmentFlag.AlignLeft)

        # подключение функций
        self.button_del.clicked.connect(self.sample_del)
        self.input.editingFinished.connect(self.is_sample_selected)
        self.input.textChanged.connect(self.sample_features)
        self.field_conc.editingFinished.connect(self.func_field_conc)
        self.field_conc.editingFinished.connect(self.sample_features)
        self.field_conc.editingFinished.connect(self.is_sample_selected)
        self.button_loops.clicked.connect(self.loops_analyze)
        self.button_homo_dimer.clicked.connect(self.homodimer_analyze)
        self.button_hetero_dimer.clicked.connect(self.heterodimer_analyze)
        self.button_get_complementary.clicked.connect(self.get_complementary_sequence_5_3)


        # добавление виджетов
        self.in_lay_sample.addWidget(self.input)
        self.in_lay_sample.addWidget(self.field_conc)
        self.in_lay_sample.addWidget(self.label_conc)
        self.lay_sample.addLayout(self.in_lay_sample)
        self.lay_sample.addWidget(self.label_features)
        self.lay_sample.addWidget(self.button_get_complementary)

        self.lay.addWidget(self.box)
        self.lay.addWidget(self.label)
        self.lay.addLayout(self.lay_sample)
        self.lay.addWidget(self.button_loops)
        self.lay.addWidget(self.button_homo_dimer)
        self.lay.addWidget(self.button_hetero_dimer)
        self.lay.addWidget(self.button_del)

        # настройка главного окна
        self.setLayout(self.lay)
        self.setStyleSheet("background-color: rgb(220,220,220); border-radius: 5px")
        self.setFixedHeight(80)

    def func_field_conc(self):
        self.check_field_conc = self.field_conc.text() == w.data_conc_DNA.text()
        if set(self.field_conc.text()).issubset('0123456789.'):
            if float(self.field_conc.text()) == 0:
                self.field_conc.setText('0.1')
        else:
            self.field_conc.setText('0.1')

    def change_sample_conc(self):
        if self.check_field_conc:
            self.field_conc.setText(w.data_conc_DNA.text())

    def sample_del(self):
        self.box.setCheckState(Qt.CheckState.Unchecked)
        self.close()
        try:
            del tm_samples[self.index.text()]
            del seq_samples[self.index.text()]
        except KeyError:
            pass

    def sample_features(self):
        concentration_DNA = self.field_conc.text()
        concentration_K = w.data_conc_K.text()
        concentration_Mg = w.data_conc_Mg.text()
        try:
            GC = GC_features(self.input.text())
            S = round(deltaS_DNA(self.input.text()))
            G = round(deltaG_DNA(self.input.text()))
            H = round(deltaH_DNA(self.input.text()))*1000
            Tm = round(
                temp_DNA_melt(self.input.text(),
                float(concentration_DNA),
                float(concentration_K),
                float(concentration_Mg)), 2)

            self.label_features.setText(
                f'Свойства: длина в нк. = {len(self.input.text())}, GC = {GC}%, deltaS = {S}, deltaH = {H}, deltaG = {G}, Tm = {Tm} °C')
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
        concentration_DNA = self.field_conc.text()
        concentration_K = w.data_conc_K.text()
        concentration_Mg = w.data_conc_Mg.text()
        try:
            melt_data = DnaFraction(
                float(concentration_DNA)/1E6,
                temp,
                deltaS_DNA(self.input.text()),
                deltaH_DNA(self.input.text())*1E3,
                float(concentration_K),
                float(concentration_Mg))

            dif_data = numpy.diff(melt_data)
            tm_samples.update(
                {self.index.text(): [self.label.text(), melt_data]})
            dif_samples.update(
                {self.index.text(): [self.label.text(), numpy.abs(dif_data)]})
            seq_samples.update(
                {self.index.text(): [self.label.text(), self.input.text()]})
        except KeyError:
            self.label_features.setText(
                'Неправильные символы или наличие пробелов в последовательности')
        except ValueError:
            print('Ошибка значений')

    def is_sample_selected(self):
        if self.box.isChecked() == True:
            self.setStyleSheet('background-color: rgb(144, 238, 144); border-radius: 5px')
            self.sample_seq()
        else:
            self.setStyleSheet("background-color: lightgray; border-radius: 5px")
            try:
                del tm_samples[self.index.text()]
                del seq_samples[self.index.text()]
                del dif_samples[self.index.text()]
            except KeyError:
                pass

    def homodimer_analyze(self):

        w.tab_lay.setCurrentIndex(2)

        w.field_text_dimers.setText(f"{self.label.text()}"
                            f"\t {self.label_features.text()}"
                            f"\t {self.input.text()}\n")
        w.field_text_dimers.append(f"ГОМОДИМЕРЫ\n")

        result = sorted(dimers_analyze(self.input.text(), self.input.text()), key=lambda x: abs(float(x[3].split()[1])), reverse=True)
        for line in result:
            for string in line[0:4]:
                w.field_text_dimers.append(str(string))
            w.field_text_dimers.append('---------------------------')
            w.field_text_dimers.append('')

    def loops_analyze(self):

        w.tab_lay.setCurrentIndex(3)

        w.field_text_loops.setText(f"{self.label.text()}"
                            f"\t {self.label_features.text()}"
                            f"\t {self.input.text()}\n")
        w.field_text_loops.append(f"ПЕТЛИ\n")

        result = loops_analyze(self.input.text())
        try:
            for line in result:
                for string in line[0:3]:
                    w.field_text_loops.append(str(string))
                w.field_text_loops.append(f'{str(line[3])} kкал/моль')
                w.field_text_loops.append('---------------------------')
                w.field_text_loops.append('')
        except:
            print('ОШИБКА В LOOPS ANALYZE')

    def heterodimer_analyze(self):
        # seq_samples -> dict{index: [name, sequence]}
        w.tab_lay.setCurrentIndex(2)
        w.field_text_dimers.setText(f'{w.data_prep_for_save()}\n')
        w.field_text_dimers.append(f"ГЕТЕРОДИМЕРЫ\n")
        result = []
        if len(seq_samples) > 0:
            # seq_samples -> dict{index: [name, sequence]}
            for key, value in seq_samples.items():
                if self.label.text() != value[0] and value[1]:
                    result.append([f'\n{self.label.text()} AND {value[0]}\n'])
                    result += sorted(dimers_analyze(self.input.text(),
                                     value[1]), key=lambda x: abs(float(x[3].split()[1])), reverse=True)
                elif not value[1]:
                    result.append([f'\n{self.label.text()} AND {value[0]}\n', f'ОТСУТСТВУЕТ ПОСЛЕДОВАТЕЛЬНОСТЬ {value[0]}\n'])
            for line in result:
                for string in line[0:4]:
                    w.field_text_dimers.append(string)
                w.field_text_dimers.append('---------------------------')

    def get_complementary_sequence_5_3(self):
        result = ''
        nucleotides_complementary = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        for symbol in self.input.text()[::-1]:
            result += nucleotides_complementary[symbol]
        self.input.setText(result)

    def change_input_seq(self):
        input_text = self.input.text().strip().replace(' ', '').upper()
        if set(input_text).issubset('CATG'):
            self.input.setText(input_text)
        else:
            self.input.setText(''.join([symbol for symbol in input_text if symbol in 'ACGT']))

# Обновление данных графика в отдельном потоке
class Thread_Canvas(QThread):
    finished = pyqtSignal()
    def __init__(self, canvas: MplComplexCanvas):
        super().__init__()
        self.canvas = canvas

    @pyqtSlot()
    def run(self):
        self.canvas.legend_fig.clear()
        self.canvas.axes.cla()
        # self.canvas.data -> dict{index: [name, array]}
        if len(self.canvas.data) > 0:
            [self.canvas.axes.plot(self.canvas.temp, values[1], label = values[0], color = palette[int(key)%50]) for key, values in self.canvas.data.items()]

            self.canvas.legend_data = self.canvas.axes.lines
            self.canvas.legend_build()

        self.canvas.axes.set_xlabel(self.canvas.name_x)
        self.canvas.axes.set_ylabel(self.canvas.name_y)
        self.canvas.axes.grid(which='both', linestyle='--', drawstyle='steps')
        self.canvas.axes.set_xticks(numpy.arange(0, 101, step=5))
        self.canvas.main_canvas.draw_idle()
        self.canvas.legend_canvas.draw_idle()
        self.canvas.setup_cursor_tracking()

    def __del__(self):
        self.quit()
        self.wait(1000)

# Окно главного приложения
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.index_sample = 0
        self.setMinimumSize(1248, 860)
        self.showMaximized()
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

        palette = self.slider_conc.palette()
        palette.setColor(QPalette.ColorRole.Highlight, QColor(0, 50, 100))
        self.slider_conc.setPalette(palette)
        self.slider_K.setPalette(palette)
        self.slider_Mg.setPalette(palette)

        lay_conditions = QVBoxLayout()
        widget_00 = QLabel()
        widget_01 = QVBoxLayout()
        self.field_text_dimers = QTextEdit()
        self.field_text_loops = QTextEdit()
        self.df_dx = QLabel()
        self.df_dx_sc2 = QVBoxLayout()
        widget_04 = QWidget()

        button_add = QPushButton('Добавить образец')
        self.main_lay = QVBoxLayout()
        self.lay_UP = QHBoxLayout()
        self.scrollArea = QScrollArea()
        self.scrollWidget = QWidget()
        self.lay_DW = QVBoxLayout()
        self.lay_samples = QVBoxLayout()
        self.tab_lay = QTabWidget()

        self.legend_scroll = QScrollArea()
        self.legend_widget = QWidget()
        self.legend_layout = QVBoxLayout()

        self.legend_widget.setMaximumWidth(180)
        self.scroll_widget_height = 500
        self.scroll_widget_weight = 100

        self.index = 0

        # plot tm

        self.sc = MplComplexCanvas()
        self.sc.main_fig.suptitle('Профиль плавления ДНК')
        self.sc.on_text = True
        self.sc.temp = temp_plot
        self.sc.data = tm_samples
        self.sc.name_y = 'Фракция двухцепочечных молекул'
        self.sc_worker1 = Thread_Canvas(self.sc)


        # plot tm df/dx

        self.sc2 = MplComplexCanvas()
        self.sc2.main_fig.suptitle('Производная функция профиля плавления ДНК')
        self.sc2.temp = temp_dif
        self.sc2.data = dif_samples
        self.sc2.name_y = 'Скорость перехода молекул \n между состояниями'
        self.sc_worker2 = Thread_Canvas(self.sc2)


        # найтрока расположения и стиль виджетов
        self.data_conc_Mg.setStyleSheet(
            "background-color: white; border-style: solid; border-radius: 5px")
        self.data_conc_DNA.setStyleSheet(
            "background-color: white; border-style: solid; border-radius: 5px")
        self.data_conc_K.setStyleSheet(
            "background-color: white; border-style: solid; border-radius: 5px")
        lay_conditions.setAlignment(Qt.AlignmentFlag.AlignTop)
        widget_00.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        
        self.field_text_dimers.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.field_text_dimers.setStyleSheet("background: white; background-color: #E8E9EB; border-style: solid; border-radius: 5px")
        self.field_text_dimers.setFont(QFont('Courier new', 8))
        self.field_text_loops.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.field_text_loops.setStyleSheet("background: white; background-color: #E8E9EB; border-style: solid; border-radius: 5px")
        self.field_text_loops.setFont(QFont('Courier new', 8))
        
        
        widget_04.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        widget_00.setStyleSheet("background: white; background-color: white; border-style: solid; border-radius: 5px")
        self.lay_DW.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.lay_samples.setAlignment(Qt.AlignmentFlag.AlignTop)
        button_add.setFixedHeight(50)
        button_add.setStyleSheet("background-color: #E8E9EB")
        widget_00.setStyleSheet("background-color: #E8E9EB; border-style: solid; border-radius: 5px")
        self.df_dx.setStyleSheet("background-color: #E8E9EB; border-style: solid; border-radius: 5px")
        widget_04.setStyleSheet("background-color: #E8E9EB; border-style: solid; border-radius: 5px")

        self.scrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setStyleSheet("""background-color: rgb(0, 32, 78)""")
        self.legend_scroll.verticalScrollBar().setStyleSheet(scroll_vertical)
        self.legend_scroll.horizontalScrollBar().setStyleSheet(scroll_horizontal)
        self.scrollArea.verticalScrollBar().setStyleSheet(scroll_vertical)
        self.field_text_dimers.verticalScrollBar().setStyleSheet(scroll_vertical)
        self.field_text_loops.verticalScrollBar().setStyleSheet(scroll_vertical)

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
        widget_01.addWidget(self.sc.main_canvas)

        self.df_dx_sc2.addWidget(self.sc2.main_canvas)
        self.df_dx.setLayout(self.df_dx_sc2)

        self.tab_lay.setMovable(True)
        self.tab_lay.setStyleSheet(tab_lay)

        self.tab_lay.addTab(widget_00, 'График температуры')
        self.tab_lay.addTab(self.df_dx, 'Функция dF/dx')
        self.tab_lay.addTab(self.field_text_dimers, 'Димеры')
        self.tab_lay.addTab(self.field_text_loops, 'Петли')

        self.legend_layout.addWidget(self.sc.legend_canvas)
        self.legend_widget.setLayout(self.legend_layout)
        self.legend_scroll.setWidget(self.legend_widget)

        self.lay_UP.addWidget(widget_04, 3)
        self.lay_UP.addWidget(self.tab_lay, 9)
        self.lay_UP.addWidget(self.legend_scroll, 2)

        self.scrollWidget.setLayout(self.lay_samples)
        self.scrollArea.setWidget(self.scrollWidget)
        self.lay_DW.addWidget(button_add)
        self.lay_DW.addWidget(self.scrollArea)

        # ДОБАВЛЕНИЕ ОСНОВНЫХ ОБЛАСТЕЙ
        # ВЕРХНЯЯ ОБЛАСТЬ
        self.main_lay.addLayout(self.lay_UP, 2)
        # НИЖНЯЯ ОБЛАСТЬ
        self.main_lay.addLayout(self.lay_DW, 2)

        # главное окно
        container = QWidget()
        container.setLayout(self.main_lay)
        self.setWindowIcon(QIcon('TAD.ico'))
        self.bar = TBar()
        self.main_toolbar()
        self.setStyleSheet("background-color: rgb(0, 32, 78)")
        self.setCentralWidget(container)
        self.setWindowTitle(f'TAD: Termodynamic Analysis of DNA (KeyArt) v.{version}')


        # подключение функций
        button_add.clicked.connect(self.sample_add)
        self.slider_conc.sliderReleased.connect(self.change_conc)
        self.data_conc_DNA.editingFinished.connect(self.changed_slide_conc)
        self.slider_K.sliderReleased.connect(self.changed_K)
        self.data_conc_K.editingFinished.connect(self.changed_slide_K)
        self.slider_Mg.sliderReleased.connect(self.changed_Mg)
        self.data_conc_Mg.editingFinished.connect(self.changed_slide_Mg)

    def main_toolbar(self):
        main_toolbar = self.addToolBar('Основные инструменты')
        main_toolbar.setStyleSheet("background-color: #E8E9EB")

        icon_open_file = QAction(QIcon(), 'Загрузить данные', self)
        icon_open_file.triggered.connect(self.open_data)
        main_toolbar.addAction(icon_open_file)

        icon_save_file = QAction(QIcon(), 'Сохранить результаты', self)
        icon_save_file.triggered.connect(self.save_data)
        main_toolbar.addAction(icon_save_file)

        icon_information = QAction(QIcon(), 'Информация о программе', self)
        icon_information.triggered.connect(self.information)
        main_toolbar.addAction(icon_information)

    def main_toolbar_addition(self):
        main_toolbar = self.addToolBar('Дополнительные инструменты')
        main_toolbar.setStyleSheet("background-color: lightgray")

        icon_open_file = QAction(QIcon(), 'Сбросить масштаб', self)
        main_toolbar.addAction(icon_open_file)

    # функции
    def save_data(self):
        data_dimers = self.field_text_dimers.toPlainText()
        data_loops = self.field_text_loops.toPlainText()
        if len(data_dimers) == 0 or len(data_loops) == 0:
            msgbox = QMessageBox()
            msgbox.setWindowTitle('Сохранение данных')
            msgbox.setText('Нет данных. Для сохранения нужно провести анализ. Отметьте необходимые образцы, и нажмите кнопки петли, гомодимеры и гетеродимеры')
            msgbox.exec()
        else:
            if 'ГЕТЕРОДИМЕРЫ' in data_dimers:
                self.bar.save_data(data = data_dimers)
            elif 'ГОМОДИМЕРЫ' in data_dimers:
                self.bar.save_data(data=data_dimers + data_loops)

    def open_data(self):
        data = self.bar.open_file()

        if data:
            for line in data.split('\n'):
                name_and_seq = line.split()
                if name_and_seq and set(name_and_seq[1].upper()).issubset('ACGT'):
                    self.sample_add()
                    self.sample.label.setText(name_and_seq[0])
                    self.sample.input.setText(name_and_seq[1])

    def information(self):
        self.new_window = QWidget()
        self.new_window.resize(600, 400)
        self.new_window.setWindowTitle('О программе')

        layout = QVBoxLayout()
        label = QTextEdit()
        label.setMarkdown("Здравствуйте!\n\n Программа предназначена для анализа ДНК последовательностей длиной до 100 нуклеотидов.\n В основе математических вычислений лежит модель ближайших соседей.\n"
            "Hatim T. Allawi and John SantaLucia, Jr. Thermodynamics and NMR of Internal G‚T Mismatches in DNA. Biochemistry 1997, 36, 10581-10594.\n\n"
            "Расчёт температуры плавления производится по следующей формуле:\n\n"
            "Tm = (deltaH/(deltaS+1.987*log(ConcDNA/1000000))) + (16.6*log10(salt/(1.0+0.7*salt))) - 273.15\n\n"
            "Расчёт термодинамического профиля по формуле:\n\n"
            "salt = (ConcK/1000) + 4 * (ConcMg/1000)**0.5\n\n"
            "CtKeq = ConcDNA * (DeltaS/R - DeltaH/(R*T-16.6*log10(salt/(1.0+0.7*salt))))\n\n"
            "f = (1 + CtKeq - (1 + 2 * CtKeq)) / CtKeq\n\n"
                          "__________________________________________________\n\n"
            f"Версия {version} от 27 февраля 2026\n\n"
            "Разработчик Артём Махнёв\n\n"
            "Почта: artyom_alex@rambler.ru")
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        label.setStyleSheet("background-color: lightgray; border-radius: 5px")
        layout.addWidget(label)
        self.new_window.setLayout(layout)
        self.new_window.show()

    def sample_add(self):
        self.index_sample += 1
        self.sample = Sample()
        self.sample.index.setText(f'{self.index_sample}')
        #self.legend_widget.setMinimumWidth(self.scroll_widget_weight)
        self.lay_samples.addWidget(self.sample)

        self.sample.input.editingFinished.connect(self.sample.change_input_seq)
        self.data_conc_DNA.textChanged.connect(self.sample.change_sample_conc)

        self.sample.box.checkStateChanged.connect(self.sample.is_sample_selected)
        self.data_conc_DNA.textChanged.connect(self.sample.is_sample_selected)
        self.data_conc_Mg.textChanged.connect(self.sample.is_sample_selected)
        self.data_conc_K.textChanged.connect(self.sample.is_sample_selected)

        self.data_conc_DNA.textChanged.connect(self.sample.sample_features)
        self.data_conc_Mg.textChanged.connect(self.sample.sample_features)
        self.data_conc_K.textChanged.connect(self.sample.sample_features)

        self.sample.box.checkStateChanged.connect(self.change_size_widget_legend)

        self.sample.box.checkStateChanged.connect(self.thread_update_plot)
        self.sample.field_conc.editingFinished.connect(self.thread_update_plot)
        self.sample.button_del.clicked.connect(self.thread_update_plot)
        self.data_conc_K.textChanged.connect(self.thread_update_plot)
        self.data_conc_DNA.textChanged.connect(self.thread_update_plot)
        self.data_conc_Mg.textChanged.connect(self.thread_update_plot)
        self.sample.input.editingFinished.connect(self.thread_update_plot)

    def thread_update_plot(self):
        self.sc_worker1.start()
        self.sc_worker2.start()

        self.sc_worker1.quit()
        self.sc_worker2.quit()

    def dimer_analyze(self):
        self.field_text_dimers.setText(self.sample.input.text)

    def changed_K(self):
        concentration_K = self.slider_K.value()
        self.data_conc_K.setText(str(concentration_K))

    def changed_slide_K(self):
        value = self.data_conc_K.text()
        if value != '' and set(value).issubset('1234567890.'):
            self.slider_K.setValue(int(value))
        else:
            self.slider_K.setValue(1)

    def change_conc(self):
        concentration_DNA = self.slider_conc.value()/100
        self.data_conc_DNA.setText(str(concentration_DNA))

    def changed_slide_conc(self):
        value = self.data_conc_DNA.text()
        if value != '' and set(value).issubset('1234567890.'):
            value = float(value)*100
            self.slider_conc.setValue(round(value))
        else:
            self.slider_conc.setValue(10)

    def changed_Mg(self):
        concentration_Mg = self.slider_Mg.value()/10
        self.data_conc_Mg.setText(str(concentration_Mg))

    def changed_slide_Mg(self):
        value = self.data_conc_Mg.text()
        if value != '' and set(value).issubset('1234567890.'):
            value = float(value)*10
            self.slider_Mg.setValue(round(value))
        else:
            self.slider_Mg.setValue(10)

    def change_size_widget_legend(self):
        for a in self.sc.legend_fig.get_default_bbox_extra_artists():
            self.legend_widget.setMinimumHeight(int(a.get_window_extent().height)+50)
            self.legend_widget.setMinimumWidth(int(a.get_window_extent().width)+50)

    def data_prep_for_save(self):
        widgets = [self.lay_samples.itemAt(i).widget() for i in range(self.lay_samples.count())]
        result = [f"{widget.findChild(QLineEdit, "label_name").text()} "
               f"\t {widget.findChild(QLineEdit, "label_features").text()} "
               f"\t {widget.findChild(QLineEdit, "label_seq").text()}"
               for widget in widgets if widget.findChild(QLineEdit, "label_seq").text() and widget.findChild(QCheckBox, "label_checkbox").isChecked()]

        return '\n'.join(result)








# Запуск приложения
if __name__ == "__main__":
    print("ЗАПУСК ПРОГРАММЫ")
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon('TAD.ico'))
    w = MainWindow()
    w.show()
    app.setStyle('Fusion')
    app.exec()
