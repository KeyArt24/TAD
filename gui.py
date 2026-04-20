"""
GUI module for DNA thermodynamic analysis application.
"""

import sys
import logging
from typing import Dict, List, Optional, Any

import numpy
import seaborn as sns
from matplotlib import style, pyplot
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from PyQt6.QtGui import QFont, QIcon, QAction, QPalette, QColor
from PyQt6.QtCore import Qt, QThread, pyqtSignal, pyqtSlot
from PyQt6.QtWidgets import (
    QApplication, QWidget, QMainWindow, QPushButton,
    QHBoxLayout, QVBoxLayout, QLabel, QSizePolicy,
    QTextEdit, QLineEdit, QScrollArea, QCheckBox,
    QSlider, QTabWidget, QMessageBox, QStackedWidget
)

from tool_bar import TBar
from tm_oligo import (
    deltaS_DNA, deltaH_DNA, deltaG_DNA, temp_DNA_melt,
    GC_features, DnaFraction, dimers_analyze, loops_analyze
)
from style import scroll_vertical, scroll_horizontal, tab_lay, css_buttons_sample, css_lines_filter

# Global constants
pyplot.rcParams.update({'font.size': 8})
pyplot.rcParams['agg.path.chunksize'] = 10000
style.use('fast')

TM_SAMPLES: Dict[str, list] = {}
DIF_SAMPLES: Dict[str, list] = {}
SEQ_SAMPLES: Dict[str, list] = {}
COUNT_SAMPLES: int = 0
TM_RESOLUTION: int = 500
TEMP_K = numpy.array([i * 0.2 + 273.15 for i in range(TM_RESOLUTION)])
TEMP_C = [i * 0.2 for i in range(TM_RESOLUTION)]
TEMP_DIF = [i * 0.2 for i in range(TM_RESOLUTION - 1)]

DEFAULT_CONC_DNA: float = 0.25
DEFAULT_CONC_K: float = 50
DEFAULT_CONC_MG: float = 3
MAX_DG_DIMERS: int = 0
MIN_DG_DIMERS: int = 0
MAX_DG_LOOPS: int = 0
MIN_DG_LOOPS: int = 0

PALETTE = sns.color_palette('deep', 50)
VERSION: str = '1.3.4'

logging.basicConfig(
    format='%(filename)s[LINE:%(lineno)d]# %(levelname)-8s [%(asctime)s]  %(message)s',
    level=logging.WARNING
)


class MplComplexCanvas(QWidget):
    """Matplotlib canvas widget for DNA melt curves."""

    cursor_moved = pyqtSignal(float, float)

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)

        self.data: Dict[str, list] = {}
        self.temp: list = []
        self.legend_data: list = []

        self.name_x: str = 'Температура °C'
        self.name_y: str = ''

        self.main_fig = Figure(dpi=100, frameon=False)
        self.legend_fig = Figure(dpi=100)

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
        for spine in ['top', 'right', 'bottom', 'left']:
            self.axes.spines[spine].set_visible(False)
        self.legend_fig.patch.set_alpha(0.0)

        self.main_canvas = FigureCanvasQTAgg(self.main_fig)
        self.legend_canvas = FigureCanvasQTAgg(self.legend_fig)

        self.diff_loop: bool = False
        self.cursor_tracking: bool = False

        self.vline = None
        self.hline = None
        self.cursor_text = None

    def legend_build(self) -> None:
        """Build legend from current plot data."""
        self.legend_fig.legend(handles=self.legend_data, loc='upper left')

    def setup_cursor_tracking(self) -> None:
        """Setup mouse event handlers for cursor tracking."""
        self.main_canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.main_canvas.mpl_connect('button_press_event', self.on_mouse_click)
        self.main_canvas.mpl_connect('scroll_event', self.on_mouse_scroll)

    def on_mouse_move(self, event) -> None:
        """Handle mouse movement for cursor position display."""
        if event.inaxes and self.cursor_tracking:
            x, y = event.xdata, event.ydata

            if self.cursor_text:
                self.cursor_text.remove()

            self.cursor_text = self.axes.text(
                0.02, 0.95,
                f'{self.name_x}={x:.2f}, y={y:.2f}',
                transform=self.axes.transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5),
            )

            self.draw_cursor_lines(x, y)
            self.cursor_moved.emit(x, y)
            self.main_canvas.draw_idle()

    def on_mouse_click(self, event) -> None:
        """Handle mouse clicks for adding markers."""
        if event.inaxes and event.button == 1 and self.cursor_tracking:
            x, y = event.xdata, event.ydata

            self.axes.plot(x, y, 'go', markersize=5, alpha=0.5)
            self.axes.annotate(
                f'({x:.2f}, {y:.2f})',
                xy=(x, y),
                xytext=(5, 5),
                textcoords='offset points',
                fontsize=8,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7)
            )
            self.main_canvas.draw_idle()

    def on_mouse_scroll(self, event) -> None:
        """Handle mouse scroll for zoom functionality."""
        if event.inaxes:
            base_scale = 1.1
            cur_xlim = self.axes.get_xlim()
            cur_ylim = self.axes.get_ylim()

            xdata, ydata = event.xdata, event.ydata

            if event.button == 'up':
                scale_factor = 1 / base_scale
            elif event.button == 'down':
                scale_factor = base_scale
            else:
                scale_factor = 1

            self.axes.set_xlim([
                xdata - (xdata - cur_xlim[0]) * scale_factor,
                xdata + (cur_xlim[1] - xdata) * scale_factor
            ])
            self.axes.set_ylim([
                ydata - (ydata - cur_ylim[0]) * scale_factor,
                ydata + (cur_ylim[1] - ydata) * scale_factor
            ])

            self.main_canvas.draw_idle()

    def draw_cursor_lines(self, x: float, y: float) -> None:
        """Draw cursor crosshair lines."""
        if self.vline:
            try:
                self.vline.remove()
            except ValueError:
                pass

        if self.hline:
            try:
                self.hline.remove()
            except ValueError:
                pass

        self.vline = self.axes.axvline(x=x, color='red', alpha=0.2, linestyle='-')
        self.hline = self.axes.axhline(y=y, color='red', alpha=0.2, linestyle='-')

    def set_axes_text(self) -> None:
        """Set annotation text on axes."""
        self.axes.text(
            0.05, 0.5, "Гибридизация",
            style="italic", fontsize=8, transform=self.axes.transAxes
        )
        self.axes.text(
            0.85, 0.5, "Денатурация",
            style="italic", fontsize=8, transform=self.axes.transAxes
        )


class Sample(QLabel):
    """DNA sample card widget."""

    def __init__(self, main_window: 'MainWindow', index: int):
        super().__init__()
        self.main_window = main_window
        self.index_val = index

        self._init_ui()
        self._connect_signals()

    def _init_ui(self) -> None:
        """Initialize UI components."""
        # Widgets
        self.box = QCheckBox()
        self.box.setObjectName("label_checkbox")
        self.index_label = QLabel(str(self.index_val))
        self.name_edit = QLineEdit(f'SAMPLE {self.index_val}')
        self.name_edit.setObjectName('label_name')
        self.seq_edit = QLineEdit()
        self.seq_edit.setObjectName('label_seq')
        self.features_label = QLineEdit(
            'Свойства: длина в нк. = 0, GC = 0%, deltaS = 0, '
            'deltaH = 0, deltaG = 0, Tm = 0 °C'
        )
        self.features_label.setObjectName("label_features")
        self.conc_label = QLabel('мкМ')
        self.conc_edit = QLineEdit(self.main_window.conditions.conc_dna_edit.text())
        self.complementary_btn = QPushButton('Комплементарная цепь')
        self.del_btn = QPushButton(' УДАЛИТЬ ')
        self.loops_btn = QPushButton('   Петли   ')
        self.homodimer_btn = QPushButton(' Гомодимеры ')
        self.heterodimer_btn = QPushButton(' Гетеродимеры ')

        self.check_field_conc = (
            self.conc_edit.text() == self.main_window.conditions.conc_dna_edit.text()
        )

        # Layouts
        self.lay = QHBoxLayout()
        self.lay_sample = QVBoxLayout()
        self.in_lay_sample = QHBoxLayout()

        # Styling
        self.name_edit.setFixedWidth(100)
        self.name_edit.setFixedHeight(30)
        self.seq_edit.setFixedHeight(20)
        self.conc_edit.setFixedWidth(50)
        self.conc_edit.setFixedHeight(20)
        self.conc_label.setFixedWidth(30)
        self.complementary_btn.setFixedHeight(20)

        self.name_edit.setStyleSheet("background-color: white")
        self.conc_edit.setStyleSheet("background-color: white")
        self.complementary_btn.setStyleSheet("background-color: white")
        self.complementary_btn.setFixedWidth(150)
        self.del_btn.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        self.loops_btn.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        self.homodimer_btn.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        self.heterodimer_btn.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)

        self.seq_edit.setStyleSheet(
            "border-style: solid; border-width: 1px; border-color: white; "
            "border-radius: 5px; background-color: white"
        )
        self.del_btn.setStyleSheet("background-color: white; color: black")
        self.loops_btn.setStyleSheet(css_buttons_sample)
        self.homodimer_btn.setStyleSheet(css_buttons_sample)
        self.heterodimer_btn.setStyleSheet(css_buttons_sample)

        self.lay.setAlignment(Qt.AlignmentFlag.AlignLeft)

        # Build layout
        self.in_lay_sample.addWidget(self.seq_edit)
        self.in_lay_sample.addWidget(self.conc_edit)
        self.in_lay_sample.addWidget(self.conc_label)

        self.lay_sample.addLayout(self.in_lay_sample)
        self.lay_sample.addWidget(self.features_label)
        self.lay_sample.addWidget(self.complementary_btn)

        self.lay.addWidget(self.box)
        self.lay.addWidget(self.name_edit)
        self.lay.addLayout(self.lay_sample)
        self.lay.addWidget(self.loops_btn)
        self.lay.addWidget(self.homodimer_btn)
        self.lay.addWidget(self.heterodimer_btn)
        self.lay.addWidget(self.del_btn)

        self.setLayout(self.lay)
        self.setStyleSheet("background-color: rgb(220,220,220); border-radius: 5px")
        self.setFixedHeight(80)

    def _connect_signals(self) -> None:
        """Connect all signal handlers."""
        self.del_btn.clicked.connect(self.sample_del)
        self.seq_edit.editingFinished.connect(self.is_sample_selected)
        self.seq_edit.textChanged.connect(self.sample_features)
        self.conc_edit.editingFinished.connect(self.func_field_conc)
        self.conc_edit.editingFinished.connect(self.sample_features)
        self.conc_edit.editingFinished.connect(self.is_sample_selected)
        self.loops_btn.clicked.connect(self.loops_analyze)
        self.homodimer_btn.clicked.connect(self.homodimer_analyze)
        self.heterodimer_btn.clicked.connect(self.heterodimer_analyze)
        self.complementary_btn.clicked.connect(self.get_complementary_sequence_5_3)

    def func_field_conc(self) -> None:
        """Validate concentration input."""
        self.check_field_conc = (
            self.conc_edit.text() == self.main_window.conditions.conc_dna_edit.text()
        )
        if set(self.conc_edit.text()).issubset('0123456789.'):
            if float(self.conc_edit.text()) == 0:
                self.conc_edit.setText('0.1')
        else:
            self.conc_edit.setText('0.1')

    def change_sample_conc(self) -> None:
        """Update sample concentration from global settings."""
        if self.check_field_conc:
            self.conc_edit.setText(self.main_window.conditions.conc_dna_edit.text())

    def sample_del(self) -> None:
        """Delete sample from all data structures."""
        self.box.setCheckState(Qt.CheckState.Unchecked)
        self.close()
        TM_SAMPLES.pop(self.index_label.text(), None)
        SEQ_SAMPLES.pop(self.index_label.text(), None)

    def _get_concentration_values(self):
        """Get current concentration values."""
        return {
            'dna': self.conc_edit.text(),
            'k': self.main_window.conditions.conc_k_edit.text(),
            'mg': self.main_window.conditions.conc_mg_edit.text()
        }

    def sample_features(self) -> None:
        """Calculate and display DNA sample features."""
        conc = self._get_concentration_values()
        try:
            seq = self.seq_edit.text()
            length = len(seq)
            gc = GC_features(seq)
            s = deltaS_DNA(seq)
            g = deltaG_DNA(seq)
            h = deltaH_DNA(seq)
            tm = temp_DNA_melt(
                Conc_DNA=float(conc['dna']),
                Length_Seq=length,
                Conc_K=float(conc['k']),
                Conc_Mg=float(conc['mg']),
                dH=float(h),
                dS=float(s)
            )

            self.features_label.setText(
                f'Свойства: длина в нк. = {length}, GC = {gc}%, '
                f'deltaS = {round(s, 2)}, deltaH = {round(h, 2)}, '
                f'deltaG (37 °C) = {tm["dG"]}, Tm = {tm["Tm"]} °C'
            )
        except KeyError:
            self.features_label.setText(
                'Неправильные символы или наличие пробелов в последовательности'
            )
            self.setStyleSheet("background-color: rgb(227,38,54)")
        except ValueError:
            self.features_label.setText('Значение концентрации ионов и ДНК должно быть больше 0')
        except Exception as e:
            logging.error(f"Error in sample_features: {e}")

    def sample_seq(self) -> None:
        """Calculate melt curve data for sample."""
        conc = self._get_concentration_values()
        seq = self.seq_edit.text()
        length = len(seq)

        try:
            melt_data = DnaFraction(
                Conc_DNA=float(conc['dna']),
                T=TEMP_K,
                Length_Seq=length,
                DeltaS=deltaS_DNA(seq),
                DeltaH=deltaH_DNA(seq),
                Conc_K=float(conc['k']),
                Conc_Mg=float(conc['mg'])
            )

            loop = self.get_max_loop()
            dif_data_loop = []
            if loop is not None:
                melt_data_loop = DnaFraction(
                    Conc_DNA=float(conc['dna']),
                    T=TEMP_K,
                    Length_Seq=length,
                    DeltaS=loop[1],
                    DeltaH=loop[0],
                    Conc_K=float(conc['k']),
                    Conc_Mg=float(conc['mg']),
                    type='loop'
                )
                dif_data_loop = numpy.diff(melt_data_loop)

            dif_data = numpy.diff(melt_data)
            TM_SAMPLES[self.index_label.text()] = [self.name_edit.text(), melt_data]
            DIF_SAMPLES[self.index_label.text()] = [
                self.name_edit.text(),
                numpy.abs(dif_data),
                numpy.abs(dif_data_loop)
            ]
            SEQ_SAMPLES[self.index_label.text()] = [self.name_edit.text(), seq]

        except KeyError:
            self.features_label.setText(
                'Неправильные символы или наличие пробелов в последовательности'
            )
        except ValueError:
            logging.error('Ошибка значений в sample_seq')

    def is_sample_selected(self) -> None:
        """Handle sample selection/deselection."""
        if self.box.isChecked():
            self.setStyleSheet('background-color: rgb(144, 238, 144); border-radius: 5px')
            self.sample_seq()
        else:
            self.setStyleSheet("background-color: lightgray; border-radius: 5px")
            TM_SAMPLES.pop(self.index_label.text(), None)
            SEQ_SAMPLES.pop(self.index_label.text(), None)
            DIF_SAMPLES.pop(self.index_label.text(), None)

    def _update_analysis_text(self, text_widget: QTextEdit, title: str) -> None:
        """Update analysis text widget with sample info."""
        text_widget.setText(
            f"{self.name_edit.text()}\t {self.features_label.text()}\t {self.seq_edit.text()}\n"
        )
        text_widget.append(f"{title}\n")

    def homodimer_analyze(self) -> None:
        """Analyze homodimers for this sample."""
        conc = self._get_concentration_values()
        self.main_window.widget_options.setCurrentIndex(1)
        self.main_window.tab_lay.setCurrentIndex(2)
        self._update_analysis_text(self.main_window.field_text_dimers, "ГОМОДИМЕРЫ")

        for dimer in dimers_analyze(self.seq_edit.text(), self.seq_edit.text(),
                                    max_positive_dG=0, max_negatrive_dG=0):
            tm_dimer = temp_DNA_melt(
                Conc_DNA=float(conc['dna']),
                Conc_K=float(conc['k']),
                Conc_Mg=float(conc['mg']),
                dH=float(dimer["total_dH_dS"][0]),
                dS=float(dimer["total_dH_dS"][1]),
                Length_Seq=dimer["length_dimer"]
            )
            self._append_dimer_result(dimer, tm_dimer)

    def _append_dimer_result(self, dimer: dict, tm_dimer: dict) -> None:
        """Append dimer analysis result to text widget."""
        filter_min = float(self.main_window.filter.line_dimers_negative_max.text())
        filter_max = float(self.main_window.filter.line_dimers_positive_max.text())
        dg = tm_dimer['dG']

        if (filter_min <= dg <= filter_max) or (filter_min == 0 and filter_max == 0):
            self.main_window.field_text_dimers.append(dimer['UP-chain'])
            self.main_window.field_text_dimers.append(dimer['bounds'])
            self.main_window.field_text_dimers.append(dimer['DW-chain'])
            self.main_window.field_text_dimers.append('\n')
            self.main_window.field_text_dimers.append(
                f'ΔG {dg} ккал/моль, Tm: {tm_dimer["Tm"]} °C'
            )
            self.main_window.field_text_dimers.append('=' * 60)

    def loops_analyze(self, state=None) -> None:
        """Analyze hairpin loops for this sample."""
        conc = self._get_concentration_values()
        self.main_window.widget_options.setCurrentIndex(1)
        self.main_window.tab_lay.setCurrentIndex(3)
        self._update_analysis_text(self.main_window.field_text_loops, "ПЕТЛИ")

        try:
            for loop in loops_analyze(self.seq_edit.text(), max_positive_dG=0, max_negatrive_dG=0):
                self._append_loop_result(loop, conc)
        except Exception as e:
            logging.error(f"ОШИБКА В LOOPS ANALYZE: {e}")

    def _append_loop_result(self, loop: dict, conc: dict) -> None:
        """Append loop analysis result to text widget."""
        length = loop["length_loop"] + loop["length_stem"]
        tm_stem = temp_DNA_melt(
            Conc_DNA=float(conc['dna']),
            Conc_K=float(conc['k']),
            Conc_Mg=float(conc['mg']),
            dH=float(loop["stem_dH_dS"][0]),
            dS=float(loop["stem_dH_dS"][1]),
            Length_Seq=length,
            type='loop'
        )
        tm_stem_loop = temp_DNA_melt(
            Conc_DNA=float(conc['dna']),
            Conc_K=float(conc['k']),
            Conc_Mg=float(conc['mg']),
            dH=float(loop["total_dH_dS"][0]),
            dS=float(loop["total_dH_dS"][1]),
            Length_Seq=length,
            type='loop'
        )
        tm_loop = temp_DNA_melt(
            Conc_DNA=float(conc['dna']),
            Conc_K=float(conc['k']),
            Conc_Mg=float(conc['mg']),
            dH=float(loop["loop_dH_dS"][0]),
            dS=float(loop["loop_dH_dS"][1]),
            Length_Seq=length,
            type='loop'
        )

        filter_min = float(self.main_window.filter.line_loops_negative_max.text())
        filter_max = float(self.main_window.filter.line_loops_positive_max.text())

        if (filter_min <= tm_stem['dG'] <= filter_max) or (filter_min == 0 and filter_max == 0):
            self.main_window.field_text_loops.append(loop["5'-side"])
            self.main_window.field_text_loops.append(loop["bounds"])
            self.main_window.field_text_loops.append(loop["3'-side"])
            self.main_window.field_text_loops.append('\n')
            self.main_window.field_text_loops.append(
                f'Cтебель-петля: ΔG {tm_stem_loop["dG"]} kкал/моль, '
                f'Tm: {tm_stem_loop["Tm"]} °C'
            )
            self.main_window.field_text_loops.append(
                f'Стебель: ΔG {tm_stem["dG"]} kкал/моль, Tm: {tm_stem["Tm"]} °C'
            )
            self.main_window.field_text_loops.append(f'Петля: ΔG {tm_loop["dG"]} kкал/моль')
            self.main_window.field_text_loops.append('=' * 60)
            self.main_window.field_text_loops.append('')

    def get_max_loop(self) -> Optional[list]:
        """Get the loop with maximum melting temperature."""
        conc = self._get_concentration_values()
        max_tm = 0
        hairpin_with_max_tm_dh_ds = None

        try:
            for loop in loops_analyze(self.seq_edit.text(), max_positive_dG=0, max_negatrive_dG=0):
                length = loop["length_loop"] + loop["length_stem"]
                tm_stem_loop = temp_DNA_melt(
                    Conc_DNA=float(conc['dna']),
                    Conc_K=float(conc['k']),
                    Conc_Mg=float(conc['mg']),
                    dH=float(loop["total_dH_dS"][0]),
                    dS=float(loop["total_dH_dS"][1]),
                    Length_Seq=length,
                    type='loop'
                )
                if tm_stem_loop['Tm'] > max_tm:
                    max_tm = tm_stem_loop['Tm']
                    hairpin_with_max_tm_dh_ds = loop["total_dH_dS"]
        except Exception as e:
            logging.error(f"ОШИБКА В GET_MAX_LOOP: {e}")

        return hairpin_with_max_tm_dh_ds

    def heterodimer_analyze(self) -> None:
        """Analyze heterodimers between this sample and others."""
        conc = self._get_concentration_values()
        self.main_window.widget_options.setCurrentIndex(1)
        self.main_window.tab_lay.setCurrentIndex(2)

        self.main_window.field_text_dimers.setText(f'{self.main_window.data_prep_for_save()}\n')
        self.main_window.field_text_dimers.append("ГЕТЕРОДИМЕРЫ\n")

        if SEQ_SAMPLES:
            for key, value in SEQ_SAMPLES.items():
                if self.name_edit.text() != value[0] and len(value[1]) != 0:
                    self._analyze_heterodimer_with_sample(value, conc)
                elif len(value[1]) == 0:
                    self.main_window.field_text_dimers.append(
                        f'\n{self.name_edit.text()} AND {value[0]}\n'
                        f'\nОТСУТСТВУЕТ ПОСЛЕДОВАТЕЛЬНОСТЬ {value[0]}\n'
                    )

    def _analyze_heterodimer_with_sample(self, other_sample: list, conc: dict) -> None:
        """Analyze heterodimers with another sample."""
        self.main_window.field_text_dimers.append(f'\n{self.name_edit.text()} AND {other_sample[0]}\n')
        for dimer in dimers_analyze(self.seq_edit.text(), other_sample[1],
                                    max_positive_dG=0, max_negatrive_dG=0):
            tm_dimer = temp_DNA_melt(
                Conc_DNA=float(conc['dna']),
                Conc_K=float(conc['k']),
                Conc_Mg=float(conc['mg']),
                dH=float(dimer["total_dH_dS"][0]),
                dS=float(dimer["total_dH_dS"][1]),
                Length_Seq=dimer["length_dimer"]
            )
            self._append_dimer_result(dimer, tm_dimer)

    def get_complementary_sequence_5_3(self) -> None:
        """Generate complementary sequence (5'->3')."""
        nucleotides_complementary = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        result = ''.join(nucleotides_complementary[s] for s in self.seq_edit.text()[::-1])
        self.seq_edit.setText(result)

    def change_input_seq(self) -> None:
        """Validate and clean input sequence."""
        input_text = self.seq_edit.text().strip().replace(' ', '').upper()
        if set(input_text).issubset('CATG'):
            self.seq_edit.setText(input_text)
        else:
            self.seq_edit.setText(''.join(s for s in input_text if s in 'ACGT'))


class ThreadCanvas(QThread):
    """Thread for updating canvas plots."""

    finished = pyqtSignal()

    def __init__(self, canvas: MplComplexCanvas):
        super().__init__()
        self.canvas = canvas

    @pyqtSlot()
    def run(self) -> None:
        """Run canvas update in separate thread."""
        self.canvas.legend_fig.clear()
        self.canvas.axes.cla()

        if self.canvas.data:
            for key, values in self.canvas.data.items():
                color = PALETTE[int(key) % 50]
                self.canvas.axes.plot(
                    self.canvas.temp, values[1], label=values[0], color=color
                )
            if self.canvas.diff_loop:
                for key, values in self.canvas.data.items():
                    if len(values[2]) > 0:
                        color = PALETTE[int(key) % 50]
                        self.canvas.axes.plot(
                            self.canvas.temp, values[2], label=values[0],
                            linestyle='dashed', color=color
                        )
            self.canvas.legend_data = self.canvas.axes.lines
            self.canvas.legend_build()

        self.canvas.axes.set_xlabel(self.canvas.name_x)
        self.canvas.axes.set_ylabel(self.canvas.name_y)
        self.canvas.axes.grid(which='both', linestyle='--', drawstyle='steps')
        self.canvas.axes.set_xticks(numpy.arange(0, 101, step=5))
        self.canvas.main_canvas.draw_idle()
        self.canvas.legend_canvas.draw_idle()
        self.canvas.setup_cursor_tracking()

    def __del__(self) -> None:
        self.quit()
        self.wait(1000)


class FilterConditions(QWidget):
    """Filter widget for dimer and loop analysis."""

    def __init__(self):
        super().__init__()
        self._init_ui()
        self.config()

    def _init_ui(self) -> None:
        """Initialize UI components."""
        self.label_deltaG_dimers = QLabel('<= ∆G <=')
        self.label_deltaG_loops = QLabel('<= ∆G <=')

        self.label_dimers = QLabel('ДИМЕРЫ-DIMERS')
        self.label_dimers_negative_max = QLabel('Отрицательное значение ∆G')
        self.label_dimers_positive_max = QLabel('Положительное значение ∆G')
        self.line_dimers_negative_max = QLineEdit('0')
        self.line_dimers_positive_max = QLineEdit('0')
        self.lay_values_dimers = QHBoxLayout()

        self.label_loops = QLabel('ПЕТЛИ-HAIRPINS (∆G Стебля)')
        self.label_loops_negative_max = QLabel('Отрицательное значение ∆G')
        self.label_loops_positive_max = QLabel('Положительное значение ∆G')
        self.line_loops_negative_max = QLineEdit('0')
        self.line_loops_positive_max = QLineEdit('0')
        self.lay_values_loops = QHBoxLayout()

        self.description = QLabel(
            '\n\nФильтр результатов анализа петель и димеров.\n'
            'Работает по принципу выбора диапазона.'
            '\nЕсли в полях максимального положительного\nи отрицательного значения стоят нули (0),'
            '\nтогда буду выданы все значимые структуры.'
            '\nПри заданых значениях, например:'
            '\nПоложительное значение ∆G 2'
            '\nОтрицательное значение ∆G -3'
            '\nРезультат, все образования с ∆G от -3 до 2 ккал/моль'
            '\nДля работы фильтра, нужно задать значения '
            '\nи снова нажать кнопку анализа (петли или димеры).'
        )

        self.main_lay = QVBoxLayout()

        self.main_lay.addWidget(self.label_dimers)
        self.lay_values_dimers.addWidget(self.line_dimers_negative_max)
        self.lay_values_dimers.addWidget(self.label_deltaG_dimers)
        self.lay_values_dimers.addWidget(self.line_dimers_positive_max)
        self.main_lay.addLayout(self.lay_values_dimers)

        self.main_lay.addWidget(self.label_loops)
        self.lay_values_loops.addWidget(self.line_loops_negative_max)
        self.lay_values_loops.addWidget(self.label_deltaG_loops)
        self.lay_values_loops.addWidget(self.line_loops_positive_max)
        self.main_lay.addLayout(self.lay_values_loops)

        self.main_lay.addWidget(self.description)
        self.setLayout(self.main_lay)

    def config(self) -> None:
        """Configure filter widget styling and connections."""
        self.main_lay.setAlignment(Qt.AlignmentFlag.AlignTop)

        self.line_dimers_negative_max.setStyleSheet(css_lines_filter)
        self.line_dimers_positive_max.setStyleSheet(css_lines_filter)
        self.line_loops_negative_max.setStyleSheet(css_lines_filter)
        self.line_loops_positive_max.setStyleSheet(css_lines_filter)
        self.description.setStyleSheet("color: black;")

        for widget in [self.line_dimers_negative_max, self.line_dimers_positive_max,
                       self.line_loops_negative_max, self.line_loops_positive_max]:
            widget.textChanged.connect(self.change_input)

    def change_input(self, text: str) -> None:
        """Ensure filter input is not empty."""
        widget = self.sender()
        if not text:
            widget.setText('0')


class ConditionsWidget(QWidget):
    """Experimental conditions widget."""

    def __init__(self):
        super().__init__()
        self._init_ui()
        self.config()

    def _init_ui(self) -> None:
        """Initialize UI components."""
        self.conc_dna_label = QLabel('Концентрация ДНК (мкМ)')
        self.conc_dna_edit = QLineEdit(str(DEFAULT_CONC_DNA))
        self.conc_k_label = QLabel('Концентрация моновалентных инов (мМ)')
        self.conc_k_edit = QLineEdit(str(DEFAULT_CONC_K))
        self.conc_mg_label = QLabel('Концентрация Мg (мМ)')
        self.conc_mg_edit = QLineEdit(str(DEFAULT_CONC_MG))

        self.slider_conc = QSlider(Qt.Orientation.Horizontal)
        self.slider_conc.setValue(int(DEFAULT_CONC_DNA * 100))
        self.slider_conc.setMinimum(1)
        self.slider_conc.setMaximum(100)
        self.slider_conc.setSingleStep(1)
        self.slider_conc.setTickPosition(QSlider.TickPosition.TicksBothSides)

        self.slider_k = QSlider(Qt.Orientation.Horizontal)
        self.slider_k.setValue(int(DEFAULT_CONC_K))
        self.slider_k.setMinimum(1)
        self.slider_k.setMaximum(100)
        self.slider_k.setSingleStep(1)
        self.slider_k.setTickPosition(QSlider.TickPosition.TicksBothSides)

        self.slider_mg = QSlider(Qt.Orientation.Horizontal)
        self.slider_mg.setValue(int(DEFAULT_CONC_MG * 10))
        self.slider_mg.setMinimum(10)
        self.slider_mg.setMaximum(100)
        self.slider_mg.setSingleStep(1)
        self.slider_mg.setTickPosition(QSlider.TickPosition.TicksBothSides)

        self.main_lay = QVBoxLayout()
        self.main_lay.addWidget(self.conc_dna_label)
        self.main_lay.addWidget(self.conc_dna_edit)
        self.main_lay.addWidget(self.slider_conc)
        self.main_lay.addWidget(self.conc_k_label)
        self.main_lay.addWidget(self.conc_k_edit)
        self.main_lay.addWidget(self.slider_k)
        self.main_lay.addWidget(self.conc_mg_label)
        self.main_lay.addWidget(self.conc_mg_edit)
        self.main_lay.addWidget(self.slider_mg)

        self.setLayout(self.main_lay)

    def config(self) -> None:
        """Configure conditions widget styling."""
        palette = self.slider_conc.palette()
        palette.setColor(QPalette.ColorRole.Highlight, QColor(0, 50, 100))

        self.slider_conc.setPalette(palette)
        self.slider_k.setPalette(palette)
        self.slider_mg.setPalette(palette)

        style_sheet = "background-color: white; border-style: solid; border-radius: 5px"
        self.conc_mg_edit.setStyleSheet(style_sheet)
        self.conc_dna_edit.setStyleSheet(style_sheet)
        self.conc_k_edit.setStyleSheet(style_sheet)

        self.main_lay.setAlignment(Qt.AlignmentFlag.AlignTop)


class MainWindow(QMainWindow):
    """Main application window."""

    def __init__(self):
        super().__init__()
        self.index_sample = 0
        self.samples: List[Sample] = []

        self._init_ui()
        self._connect_signals()

    def _init_ui(self) -> None:
        """Initialize main window UI."""
        self.setMinimumSize(1248, 860)
        self.showMaximized()
        self.setWindowIcon(QIcon(r"D:\06_Python\TPD ver.1.03\HRM.png"))

        # Widgets
        self.filter = FilterConditions()
        self.conditions = ConditionsWidget()
        self.widget_options = QStackedWidget()
        self.widget_options.addWidget(self.conditions)
        self.widget_options.addWidget(self.filter)

        widget_00 = QLabel()
        widget_01 = QVBoxLayout()

        self.field_text_dimers = QTextEdit()
        self.field_text_loops = QTextEdit()
        self.df_dx = QLabel()
        self.df_dx_sc2 = QVBoxLayout()

        button_add = QPushButton('Добавить образец')
        self.main_lay = QVBoxLayout()
        self.lay_up = QHBoxLayout()
        self.scroll_area = QScrollArea()
        self.scroll_widget = QWidget()
        self.lay_down = QVBoxLayout()
        self.lay_samples = QVBoxLayout()
        self.tab_lay = QTabWidget()

        self.legend_scroll = QScrollArea()
        self.legend_widget = QWidget()
        self.legend_layout = QVBoxLayout()

        self.legend_widget.setMaximumWidth(180)

        # Plot widgets
        self.sc = MplComplexCanvas()
        self.sc.main_fig.suptitle('Профиль плавления ДНК')
        self.sc.temp = TEMP_C
        self.sc.data = TM_SAMPLES
        self.sc.name_y = 'Фракция двухцепочечных молекул'
        self.sc_worker1 = ThreadCanvas(self.sc)

        self.sc2 = MplComplexCanvas()
        self.sc2.main_fig.suptitle('Производная функция профиля плавления ДНК')
        self.sc2.diff_loop = False
        self.sc2.temp = TEMP_DIF
        self.sc2.data = DIF_SAMPLES
        self.sc2.name_y = 'Скорость перехода молекул \n между состояниями'
        self.sc_worker2 = ThreadCanvas(self.sc2)

        # Styling
        widget_00.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)

        style_base = "background: white; background-color: #E8E9EB; border-style: solid; border-radius: 5px"
        self.field_text_dimers.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.field_text_dimers.setStyleSheet(style_base)
        self.field_text_dimers.setFont(QFont('Courier new', 8))

        self.field_text_loops.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.field_text_loops.setStyleSheet(style_base)
        self.field_text_loops.setFont(QFont('Courier new', 8))

        self.widget_options.setStyleSheet(style_base)
        widget_00.setStyleSheet("background: white; background-color: white; border-style: solid; border-radius: 5px")

        self.lay_down.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.lay_samples.setAlignment(Qt.AlignmentFlag.AlignTop)

        button_add.setFixedHeight(50)
        button_add.setStyleSheet("background-color: #E8E9EB")
        widget_00.setStyleSheet("background-color: #E8E9EB; border-style: solid; border-radius: 5px")
        self.df_dx.setStyleSheet("background-color: #E8E9EB; border-style: solid; border-radius: 5px")

        self.scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setStyleSheet("background-color: rgb(0, 32, 78)")

        for scroll in [self.legend_scroll.verticalScrollBar(),
                       self.scroll_area.verticalScrollBar(),
                       self.field_text_dimers.verticalScrollBar(),
                       self.field_text_loops.verticalScrollBar()]:
            scroll.setStyleSheet(scroll_vertical)

        self.legend_scroll.horizontalScrollBar().setStyleSheet(scroll_horizontal)

        # Layout assembly
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

        self.lay_up.addWidget(self.widget_options, 2)
        self.lay_up.addWidget(self.tab_lay, 9)
        self.lay_up.addWidget(self.legend_scroll, 2)

        self.scroll_widget.setLayout(self.lay_samples)
        self.scroll_area.setWidget(self.scroll_widget)
        self.lay_down.addWidget(button_add)
        self.lay_down.addWidget(self.scroll_area)

        self.main_lay.addLayout(self.lay_up, 2)
        self.main_lay.addLayout(self.lay_down, 2)

        container = QWidget()
        container.setLayout(self.main_lay)
        self.setWindowIcon(QIcon(r'C:\Users\Stupnikova\PycharmProjects\PythonProject\TAD.ico'))
        self.setStyleSheet("background-color: rgb(0, 32, 78)")
        self.bar = TBar()
        self.setCentralWidget(container)
        self.setWindowTitle(f'TAD: Termodynamic Analysis of DNA (KeyArt) v.{VERSION}')

        # Connect button
        button_add.clicked.connect(self.sample_add)

    def _connect_signals(self) -> None:
        """Connect all signal handlers."""
        self.tab_lay.currentChanged.connect(self.on_tab_changed)

        self.conditions.slider_conc.actionTriggered.connect(self.change_conc)
        self.conditions.conc_dna_edit.editingFinished.connect(self.changed_slide_conc)
        self.conditions.slider_k.sliderReleased.connect(self.changed_k)
        self.conditions.conc_k_edit.editingFinished.connect(self.changed_slide_k)
        self.conditions.slider_mg.sliderReleased.connect(self.changed_mg)
        self.conditions.conc_mg_edit.editingFinished.connect(self.changed_slide_mg)

        self._setup_toolbar()

    def _setup_toolbar(self) -> None:
        """Setup main toolbar."""
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

        # Additional toolbar
        additional_toolbar = self.addToolBar('Дополнительные инструменты')
        additional_toolbar.setStyleSheet("background-color: #E8E9EB")

        self.icon_tracking_cursor = QAction(QIcon(), 'Трэкинг курсора - Выкл', self)
        self.icon_tracking_cursor.triggered.connect(self.switch_tracking_cursor)
        additional_toolbar.addAction(self.icon_tracking_cursor)

        self.icon_projection_loop = QAction(QIcon(), 'Проекция петли - Выкл', self)
        self.icon_projection_loop.triggered.connect(self.projection_loop)
        additional_toolbar.addAction(self.icon_projection_loop)

    def save_data(self) -> None:
        """Save analysis results to file."""
        data_dimers = self.field_text_dimers.toPlainText()
        data_loops = self.field_text_loops.toPlainText()

        if not data_dimers or not data_loops:
            msgbox = QMessageBox()
            msgbox.setWindowTitle('Сохранение данных')
            msgbox.setText(
                'Нет данных. Для сохранения нужно провести анализ. '
                'Отметьте необходимые образцы, и нажмите кнопки петли, гомодимеры и гетеродимеры'
            )
            msgbox.exec()
        else:
            if 'ГЕТЕРОДИМЕРЫ' in data_dimers:
                self.bar.save_data(data=data_dimers)
            elif 'ГОМОДИМЕРЫ' in data_dimers:
                self.bar.save_data(data=data_dimers + data_loops)

    def open_data(self) -> None:
        """Open and load data from file."""
        data = self.bar.open_file()

        if data:
            for line in data.split('\n'):
                parts = line.split()
                if parts and set(parts[1].upper()).issubset('ACGT'):
                    self.sample_add()
                    self.samples[-1].name_edit.setText(parts[0])
                    self.samples[-1].seq_edit.setText(parts[1])
        else:
            msgbox = QMessageBox()
            msgbox.setWindowTitle('ВНИМАНИЕ!!!')
            msgbox.setText(
                ' Неверный формат данных.\n Файл данных должен быть формата .txt.\n '
                'В котором в виде списка внесены:\n Имя Последовательность\n'
                ' Forward_F1 AGCTAGCGAGCGAGCGAGCAC\n'
                ' Revers_R1 AGCTAGCAGGAGAGCGAGGCGAGC'
            )
            msgbox.exec()

    def information(self) -> None:
        """Show information dialog."""
        self.new_window = QWidget()
        self.new_window.resize(600, 400)
        self.new_window.setWindowTitle('О программе')

        layout = QVBoxLayout()
        label = QTextEdit()
        label.setMarkdown(
            "Здравствуйте!\n\n Программа предназначена для анализа ДНК последовательностей "
            "длиной до 100 нуклеотидов.\n В основе математических вычислений лежит модель ближайших соседей.\n"
            "Hatim T. Allawi and John SantaLucia, Jr. Thermodynamics and NMR of Internal G‚T Mismatches in DNA. "
            "Biochemistry 1997, 36, 10581-10594.\n\n"
            "Расчёт температуры плавления производится по следующей формуле:\n\n"
            "Температура плавления ДНК и димеров: \n\n"
            "Tm = dH / (dS+1.987*numpy.log(Conc_DNA/4)) + Salt_Correction - 273.15\n\n"
            "Температура плавления петли: \n\n"
            "Tm = dH / dS + Salt_Correction - 273.15\n\n"
            "Расчёт термодинамического профиля по формуле:\n\n"
            "Conc_Salt = (Conc_K / 1000) + 4 * (Conc_Mg / 1000) ** 0.5\n\n"
            "Effective_Conc_Salt = Conc_Salt/(1.0+0.7*Conc_Salt)\n\n"
            "Salt_Correction = 0.368 * (Length_Seq - 1) * log(Effective_Conc_Salt)\n\n"
            "Keq = numpy.exp(DeltaS / R - DeltaH / (R * (T - Salt_Correction)))\n\n"
            "ДНК и димеры:\n\n"
            "CtKeq = (Conc_DNA / 4) * Keq\n\n"
            "f = (4 * numpy.round(CtKeq, 5) + 1 - numpy.sqrt(8 * numpy.round(CtKeq, 5) + 1)) / (4 * CtKeq)\n\n"
            "Петли:\n\n"
            "CtKeq = Keq\n\n"
            "f = CtKeq / (1 + CtKeq)\n\n"
            "=====================================================\n\n"
            f"Версия {VERSION} от 01 апреля 2026\n\n"
            "Разработчик Артём Махнёв\n\n"
            "Почта: artyom_alex@rambler.ru"
        )
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        label.setStyleSheet("background-color: lightgray; border-radius: 5px")
        layout.addWidget(label)
        self.new_window.setLayout(layout)
        self.new_window.show()

    def sample_add(self) -> None:
        """Add new sample to the interface."""
        self.index_sample += 1
        sample = Sample(self, self.index_sample)
        self.samples.append(sample)
        self.lay_samples.addWidget(sample)

        sample.seq_edit.editingFinished.connect(sample.change_input_seq)
        self.conditions.conc_dna_edit.textChanged.connect(sample.change_sample_conc)

        sample.box.checkStateChanged.connect(sample.is_sample_selected)
        for widget in [self.conditions.conc_dna_edit, self.conditions.conc_mg_edit,
                       self.conditions.conc_k_edit]:
            widget.textChanged.connect(sample.is_sample_selected)
            widget.textChanged.connect(sample.sample_features)

        sample.box.checkStateChanged.connect(self.change_size_widget_legend)
        sample.box.checkStateChanged.connect(self.thread_update_plot)
        sample.conc_edit.editingFinished.connect(self.thread_update_plot)
        sample.del_btn.clicked.connect(self.thread_update_plot)

        for widget in [self.conditions.conc_k_edit, self.conditions.conc_dna_edit,
                       self.conditions.conc_mg_edit]:
            widget.textChanged.connect(self.thread_update_plot)
        sample.seq_edit.editingFinished.connect(self.thread_update_plot)

    def thread_update_plot(self) -> None:
        """Update plots in separate threads."""
        self.sc_worker1.start()
        self.sc_worker2.start()
        self.sc_worker1.quit()
        self.sc_worker2.quit()

    def changed_k(self) -> None:
        """Handle K+ concentration slider change."""
        concentration_k = self.conditions.slider_k.value()
        self.conditions.conc_k_edit.setText(str(concentration_k))

    def changed_slide_k(self) -> None:
        """Handle K+ concentration text change."""
        value = self.conditions.conc_k_edit.text()
        if value and set(value).issubset('1234567890.'):
            self.conditions.slider_k.setValue(int(float(value)))
        else:
            self.conditions.slider_k.setValue(1)

    def change_conc(self) -> None:
        """Handle DNA concentration slider change."""
        concentration_dna = self.conditions.slider_conc.value() / 100
        self.conditions.conc_dna_edit.setText(str(concentration_dna))

    def changed_slide_conc(self) -> None:
        """Handle DNA concentration text change."""
        value = self.conditions.conc_dna_edit.text()
        if value and set(value).issubset('1234567890.'):
            value = float(value) * 100
            self.conditions.slider_conc.setValue(round(value))
        else:
            self.conditions.slider_conc.setValue(10)

    def changed_mg(self) -> None:
        """Handle Mg2+ concentration slider change."""
        concentration_mg = self.conditions.slider_mg.value() / 10
        self.conditions.conc_mg_edit.setText(str(concentration_mg))

    def changed_slide_mg(self) -> None:
        """Handle Mg2+ concentration text change."""
        value = self.conditions.conc_mg_edit.text()
        if value and set(value).issubset('1234567890.'):
            value = float(value) * 10
            self.conditions.slider_mg.setValue(round(value))
        else:
            self.conditions.slider_mg.setValue(10)

    def change_size_widget_legend(self) -> None:
        """Adjust legend widget size."""
        for artist in self.sc.legend_fig.get_default_bbox_extra_artists():
            self.legend_widget.setMinimumHeight(int(artist.get_window_extent().height) + 50)
            self.legend_widget.setMinimumWidth(int(artist.get_window_extent().width) + 50)

    def data_prep_for_save(self) -> str:
        """Prepare data string for saving."""
        result = []
        for sample in self.samples:
            if sample.seq_edit.text() and sample.box.isChecked():
                result.append(
                    f"{sample.name_edit.text()}\t {sample.features_label.text()}\t {sample.seq_edit.text()}"
                )
        return '\n'.join(result)

    def on_tab_changed(self, index: int) -> None:
        """Handle tab change event."""
        if index in [2, 3]:
            self.widget_options.setCurrentIndex(1)
        elif index in [0, 1]:
            self.widget_options.setCurrentIndex(0)

    def projection_loop(self) -> None:
        """Toggle loop projection display."""
        if self.sc2.diff_loop:
            self.sc2.diff_loop = False
            self.icon_projection_loop.setText('Проекция петли - Выкл')
        else:
            self.sc2.diff_loop = True
            self.icon_projection_loop.setText('Проекция петли - Вкл')
        self.thread_update_plot()

    def switch_tracking_cursor(self) -> None:
        """Toggle cursor tracking."""
        if self.sc.cursor_tracking and self.sc2.cursor_tracking:
            self.sc.cursor_tracking = False
            self.sc2.cursor_tracking = False
            self.icon_tracking_cursor.setText('Трэкинг курсора - Выкл')
        else:
            self.sc.cursor_tracking = True
            self.sc2.cursor_tracking = True
            self.icon_tracking_cursor.setText('Трэкинг курсора - Вкл')
        self.thread_update_plot()


if __name__ == "__main__":
    print('ЗАПУСК ПРОГРАММЫ')
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon(r'C:\Users\Stupnikova\PycharmProjects\PythonProject\TAD.ico'))
    window = MainWindow()
    window.show()
    app.setStyle('Fusion')
    app.exec()