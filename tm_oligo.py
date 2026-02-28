from math import *
import numpy

NN_list_G = {'AA': -1, 'TT': -1, 'CC': -1.84, 'GG': -1.84, 'AT': -0.88, 'TA': -0.58, 'AC': -1.44, 'CA':
             -1.45, 'AG': -1.28, 'GA': -1.30, 'CG': -2.17, 'GC': -2.24, 'TC': -1.30, 'CT': -1.28, 'TG': -1.45, 'GT': -1.44}

NN_list_S = {'AA': -22.2, 'TT': -22.2, 'CC': -19.9, 'GG': -19.9, 'AT': -20.4, 'TA': -21.3, 'AC': -22.4, 'CA': -22.7,
             'AG': -21.0, 'GA': -22.2, 'CG': -27.2, 'GC': -24.4, 'TC': -22.2, 'CT': -21.0, 'TG': -22.7, 'GT': -22.4}

NN_list_H = {'AA': -7.9, 'TT': -7.9, 'CC': -8.0, 'GG': -8.0, 'AT': -7.2, 'TA': -7.2, 'AC': -8.4, 'CA':
             -8.5, 'AG': -7.8, 'GA': -8.2, 'CG': -10.6, 'GC': -9.8, 'TC': -8.2, 'CT': -7.8, 'TG': -8.5, 'GT': -8.4}


def DnaFraction(Ct, T, DeltaS, DeltaH, CtK=50, CtMg=3, R=1.987):
    """
    Returns the fraction of dsDNA in a solution containing equal concentrations
    of two complementary ssDNA oligos, as a function of total [DNA],
    temperature, entropy of annealing, and enthalpy of annealing.

    Default units are mole, calorie, kelvin. For other unit systems, supply
    the appropriate value for the optional gas constant parameter.

    T can be a single value or a numpy.array of values.
    """
    # Compute Ct * Keq
    salt = (CtK/1000) + 4 * (CtMg/1000)**0.5
    CtKeq = Ct * numpy.exp(DeltaS/R - DeltaH /(R*T-16.6*log10(salt/(1.0+0.7*salt))))
    CtKeq1 = numpy.round(CtKeq, 5)
    # Compute f
    f = (1 + CtKeq1 - numpy.sqrt(1 + 2*CtKeq1)) / CtKeq

    return f

def middles(arr):
    """
    Returns a new numpy array, 1 element shorter than the input array,
    whose values come from averaging each two adjacent values in the input.
    """
    result = []
    for i in range(0, len(arr) - 1):
        element = (arr[i] + arr[i+1]) / 2.0
        result.append(element)

    return numpy.array(result)

def deltaS_DNA(sequence):
    sum_S = 0
    if len(sequence) > 0:
        if str(sequence[0][-1]).upper() == 'A' or str(sequence[0][-1]).upper() == 'T':
            sum_S += 0
        elif str(sequence[0][-1]).upper() == 'C' or str(sequence[0][-1]).upper() == 'G':
            sum_S += 0
        for i in range(len(sequence)-1):
            sum_S += NN_list_S[str(sequence[i:i+2]).upper()]
        sum_S += -5.7
    return sum_S

def deltaH_DNA(sequence):
    sum_H = 0
    if len(sequence) > 0:
        if str(sequence[0][-1]).upper() == 'A' or str(sequence[0][-1]).upper() == 'T':
            sum_H += 0
        elif str(sequence[0][-1]).upper() == 'C' or str(sequence[0][-1]).upper() == 'G':
            sum_H += 0
        for i in range(len(sequence)-1):
            sum_H += NN_list_H[str(sequence[i:i+2]).upper()]
        sum_H += 0.2
    return sum_H

def deltaG_DNA(sequence):
    sum_G = 0
    if len(sequence) > 0:
        if str(sequence[0][-1]).upper() == 'A' or str(sequence[0][-1]).upper() == 'T':
            sum_G += 0
        elif str(sequence[0][-1]).upper() == 'C' or str(sequence[0][-1]).upper() == 'G':
            sum_G += 0
        for i in range(len(sequence)-1):
            sum_G += NN_list_G[str(sequence[i:i+2]).upper()]
        sum_G += 0
    return sum_G

def GC_features(sequence: str):
    countG = sequence.upper().count('G')
    countC = sequence.upper().count('C')
    result = (countG+countC)/len(sequence)
    return round(result*100)

def dimers_analyze(seq1: str, seq2: str)-> list:
    result = []
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    rows = len(seq1)
    cols = len(seq2)
    seq2 = seq2[::-1]
    for i in range(1, cols + rows):
        sub_result = []
        bounds = ''
        deltaG = ''
        deltaH = ''
        deltaS = ''
        start_column = max(0, i - rows)
        end_column = min(i, cols)
        for j in range(start_column, end_column):
            n = rows-i-1+j+1
            m = j
            if any([(seq1[n] == 'A' and seq2[m] == 'T') or (seq1[n] == 'T' and seq2[m] == 'A'),
                    (seq1[n] == 'C' and seq2[m] == 'G') or (seq1[n] == 'G' and seq2[m] == 'C')]):
                bounds += 'I'
                deltaG += str(n) + " "
                deltaH += str(n) + " "
                deltaS += str(n) + " "
            else:
                bounds += '-'
                deltaG += '-'

        bounds = bounds.replace('-I-', '-*-').replace('-I-', '-*-')
        if bounds[0:2] == 'I-':
            bounds = '*' + bounds[1::]
        if bounds[-2::] == '-I':
            bounds = bounds[0:-1] + '*'
        positions_dimers = [position.split() for position in deltaG.split('-') if len(position.strip()) > 2]
        if len(deltaG) != 0:
            deltaG = sum([sum([NN_list_G[seq1.replace('+', '')[int(numbers):int(numbers)+2]] for numbers in position
                               if len(seq1[int(numbers):int(numbers)+2]) > 1]) for position in positions_dimers])
            deltaH = sum([sum([NN_list_H[seq1.replace('+', '')[int(numbers):int(numbers)+2]] for numbers in position
                               if len(seq1[int(numbers):int(numbers)+2]) > 1]) for position in positions_dimers])
            deltaS = sum([sum([NN_list_S[seq1.replace('+', '')[int(numbers):int(numbers)+2]] for numbers in position
                               if len(seq1[int(numbers):int(numbers)+2]) > 1]) for position in positions_dimers])

            deltaG = round(deltaG, 2)
            deltaH = round(deltaH, 2)
            deltaS = round(deltaS, 2)
        else:
            deltaG = 0

        if bounds.count('I') > 2 and deltaG < 0:
            sub_result.append(f"{' ' * (m) + "5'-" + seq1.replace('+', '')}-3'")
            sub_result.append(f"   {' ' * max(n, m) + bounds}")
            sub_result.append(f"{' ' * (n) + "3'-" + seq2.replace('+', '')}-5'")
            sub_result.append(f"\n deltG {deltaG} ккал/моль")
            sub_result.append(f"\n deltH {deltaH} ккал/моль")
            sub_result.append(f"\n deltS {deltaS} ккал/моль")
            result.append(sub_result)

    return result

def temp_DNA_melt(sequence, CtDNA, CtK, CtMg):
    Tm = 0
    if len(sequence) > 0:
        salt = (CtK/1000) + 4 * (CtMg/1000)**0.5
        deltaH = deltaH_DNA(sequence)*1000
        deltaS = deltaS_DNA(sequence)
        Tm = (deltaH/(deltaS+1.987*log(CtDNA/1000000))) + \
            (16.6*log10(salt/(1.0+0.7*salt))) - 273.15
    return Tm

def loops_analyze(seq):
    result = list()
    string = list(seq)
    loop = list()

    nn_params = {
        'AT': -1.00, 'TA': -0.88, 'GC': -2.17, 'CG': -1.84,
        'AA': -0.93, 'TT': -1.11, 'AG': -1.10, 'GA': -1.43,
        'CT': -1.44, 'TC': -1.28, 'GT': -1.92, 'TG': -1.45,
        'CA': -1.28, 'AC': -1.44, 'GG': -2.08, 'CC': -1.91,
    }

    # Энергия инициализации и петель
    init_terminal = 0.98
    loop_energies = {3: 5.70, 4: 5.30, 5: 5.60, 6: 5.70, 7: 5.80, 8: 6.00}


    while string:
        overhang = '│'
        #loop.append(string.pop(0))
        for shift in range(2):
            bounds = ''
            positions_stem = []
            positions_loop = []
            for i in range(min([len(string), len(loop)])):
                if any([string[i] == 'A' and loop[::-1][i] == 'T', string[i] == 'T' and loop[::-1][i] == 'A',
                        string[i] == 'C' and loop[::-1][i] == 'G', string[i] == 'G' and loop[::-1][i] == 'C']):
                    bounds += 'I'
                    positions_stem.append(i)
                    positions_loop.append('n')
                else:
                    bounds += '-'
                    positions_stem.append('n')
                    positions_loop.append(i)

            bounds = bounds.replace('-I-', '-*-').replace('-I-', '-*-')
            if bounds[0:2] == 'I-':
                bounds = '*' + bounds[1::]
            if bounds[-2::] == '-I':
                bounds = bounds[0:-1] + '*'


            if bounds.count('II') >= 1:

                bounds_positions = {'loop': positions_loop, 'stem': positions_stem}

                # Расчёт длины петли
                loop_symbols = ([loop[::-1][value] for value in bounds_positions['loop'] if value != 'n'] +
                                [string[value] for value in bounds_positions['loop'] if value != 'n'])
                loop_GC = loop_symbols.count('G') + loop_symbols.count('C')
                loop_len = len(loop_symbols)
                if overhang in 'ACTG':
                    loop_len += 1
                if overhang in 'CG':
                    loop_GC += 1

                # Расчёт энергии стебля (сумма динуклеотидов)
                stem_energy = init_terminal
                choose_string = {len(string): string, len(loop): loop}

                for value in bounds_positions['stem']:
                    if value != 'n':
                        pair = ''.join(choose_string[min(choose_string)][value: value+2])
                        if pair in nn_params:
                            stem_energy += nn_params[pair]
                        elif pair in ['A', 'T']:
                            stem_energy += -1.0
                        elif pair in ['C', 'G']:
                            stem_energy += -2.0

                # Энергия петли
                loop_energy = 0
                if loop_len in loop_energies:
                    loop_energy = loop_energies[loop_len]
                elif loop_len != 0 or loop_len not in loop_energies:
                    # Аппроксимация для больших петель
                    loop_energy = 6.0 + 0.3 * (loop_len - 8) if loop_len > 8 else 6.0

                # Коррекция на GC-состав петли
                if loop_len != 0:
                    gc_content = loop_GC/ loop_len
                    loop_energy -= gc_content * 0.3  # GC-пары стабилизируют

                # Общая энергия
                total_delta_g = stem_energy + loop_energy
                total_delta_g = round(total_delta_g, 2)

                if total_delta_g < 2:

                    result.append([f"┌{''.join(loop[::-1])}-5'",
                                   f'{overhang}{bounds}',
                                   f"└{''.join(string)}-3'",
                                   total_delta_g])

            if shift == 0 and string:
                overhang = string.pop(0)

        loop.append(overhang)

    result = sorted(result, key=lambda x: x[3])
    return result




