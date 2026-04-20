"""
Thermodynamic calculations for DNA oligonucleotides.
Contains functions for DNA melting temperature, dimers, and loops analysis.
"""

from math import log, log10
import numpy

# Словари термодинамических параметров (NN model)
NN_LIST_G = {
    'AA': -1, 'TT': -1, 'CC': -1.84, 'GG': -1.84,
    'AT': -0.88, 'TA': -0.58, 'AC': -1.44, 'CA': -1.45,
    'AG': -1.28, 'GA': -1.30, 'CG': -2.17, 'GC': -2.24,
    'TC': -1.30, 'CT': -1.28, 'TG': -1.45, 'GT': -1.44
}

NN_LIST_S = {
    'AA': -22.2, 'TT': -22.2, 'CC': -19.9, 'GG': -19.9,
    'AT': -20.4, 'TA': -21.3, 'AC': -22.4, 'CA': -22.7,
    'AG': -21.0, 'GA': -22.2, 'CG': -27.2, 'GC': -27.4,
    'TC': -22.2, 'CT': -21.0, 'TG': -22.7, 'GT': -22.4
}

NN_LIST_H = {
    'AA': -7.9, 'TT': -7.9, 'CC': -8.0, 'GG': -8.0,
    'AT': -7.2, 'TA': -7.2, 'AC': -8.4, 'CA': -8.5,
    'AG': -7.8, 'GA': -8.2, 'CG': -10.6, 'GC': -10.6,
    'TC': -8.2, 'CT': -7.8, 'TG': -8.5, 'GT': -8.4
}


def DnaFraction(Conc_DNA, T, DeltaS, DeltaH, Length_Seq=0, Conc_K=50,
                Conc_Tris=0, Conc_Mg=3, R=1.987, type='hybridization'):
    """
    Calculate DNA melting curve fraction.

    Parameters
    ----------
    Conc_DNA : float
        DNA concentration in uM
    T : numpy.array
        Temperature array in Kelvin
    DeltaS : float
        Entropy change in cal/(mol*K)
    DeltaH : float
        Enthalpy change in kcal/mol
    Length_Seq : int
        Sequence length
    Conc_K : float
        Potassium concentration in mM
    Conc_Tris : float
        Tris concentration in mM (not used)
    Conc_Mg : float
        Magnesium concentration in mM
    R : float
        Gas constant in cal/(mol*K)
    type : str
        'hybridization' or 'loop'

    Returns
    -------
    numpy.array
        Fraction of double-stranded molecules
    """
    # Перевод DeltaH в кал/моль
    DeltaH = DeltaH * 1000
    # Перевод концентрации ДНК в моль/л
    Conc_DNA = Conc_DNA * 1e-6 / 4
    # Учёт концентрации солей
    if Conc_Mg < 0:
        Conc_Mg = 0
    Conc_Salt = (Conc_K / 1000) + 4 * (Conc_Mg / 1000) ** 0.5
    Effective_Conc_Salt = Conc_Salt / (1.0 + 0.7 * Conc_Salt)
    Salt_Correction = 0.368 * (Length_Seq - 1) * log(Effective_Conc_Salt)

    if type == 'loop':
        Salt_Correction = 16.6 * log10(Effective_Conc_Salt)
        Keq = numpy.exp(DeltaS / R - DeltaH / (R * (T - Salt_Correction)))
        CtKeq = Keq
        f = CtKeq / (1 + CtKeq)
    elif type == 'hybridization':
        Keq = numpy.exp((DeltaS + Salt_Correction) / R - DeltaH / (R * T))
        CtKeq = Conc_DNA * Keq
        f = (4 * numpy.round(CtKeq, 5) + 1 -
             numpy.sqrt(8 * numpy.round(CtKeq, 5) + 1)) / (4 * CtKeq)

    return f


def middles(arr):
    """
    Returns a new numpy array, 1 element shorter than the input array,
    whose values come from averaging each two adjacent values in the input.
    """
    result = []
    for i in range(0, len(arr) - 1):
        element = (arr[i] + arr[i + 1]) / 2.0
        result.append(element)
    return numpy.array(result)


def deltaS_DNA(sequence):
    """Calculate entropy change (Delta S) for DNA sequence."""
    sum_S = 0
    if len(sequence) > 0:
        if sequence[0].upper() in ['A', 'T'] or sequence[-1].upper() in ['A', 'T']:
            sum_S += -1.4
        elif str(sequence[0][-1]).upper() == 'C' or str(sequence[0][-1]).upper() == 'G':
            sum_S += -2.8
        for i in range(len(sequence) - 1):
            sum_S += NN_LIST_S[str(sequence[i:i + 2]).upper()]
    return sum_S


def deltaH_DNA(sequence):
    """Calculate enthalpy change (Delta H) for DNA sequence."""
    sum_H = 0
    if len(sequence) > 0:
        if sequence[0].upper() in ['A', 'T'] or sequence[-1].upper() in ['A', 'T']:
            sum_H += 0
        elif str(sequence[0][-1]).upper() == 'C' or str(sequence[0][-1]).upper() == 'G':
            sum_H += -0.1
        for i in range(len(sequence) - 1):
            sum_H += NN_LIST_H[str(sequence[i:i + 2]).upper()]
    return sum_H


def deltaG_DNA(sequence):
    """Calculate Gibbs free energy change (Delta G) for DNA sequence."""
    sum_G = 0.98
    if len(sequence) > 0:
        if str(sequence[0][-1]).upper() == 'A' or str(sequence[0][-1]).upper() == 'T':
            sum_G += 0
        elif str(sequence[0][-1]).upper() == 'C' or str(sequence[0][-1]).upper() == 'G':
            sum_G += 0
        for i in range(len(sequence) - 1):
            sum_G += NN_LIST_G[str(sequence[i:i + 2]).upper()]
    return sum_G


def GC_features(sequence: str):
    """Calculate GC content percentage for DNA sequence."""
    countG = sequence.upper().count('G')
    countC = sequence.upper().count('C')
    result = (countG + countC) / len(sequence)
    return round(result * 100)


def dimers_analyze(seq1: str, seq2: str, max_negatrive_dG=0, max_positive_dG=0):
    """
    Analyze dimer formations between two DNA sequences.

    Parameters
    ----------
    seq1 : str
        First DNA sequence
    seq2 : str
        Second DNA sequence
    max_negatrive_dG : float
        Minimum Delta G filter value
    max_positive_dG : float
        Maximum Delta G filter value

    Returns
    -------
    list
        List of dimer dictionaries
    """
    result = []
    n_dimers = 0
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    rows = len(seq1)
    cols = len(seq2)
    seq2_rev = seq2[::-1]

    for i in range(1, cols + rows):
        bounds = ''
        positions_bounds_seq1 = ''
        positions_bounds_seq2 = ''
        start_column = max(0, i - rows)
        end_column = min(i, cols)

        for j in range(start_column, end_column):
            n = rows - i - 1 + j + 1
            m = j
            if any([
                (seq1[n] == 'A' and seq2_rev[m] == 'T'),
                (seq1[n] == 'T' and seq2_rev[m] == 'A'),
                (seq1[n] == 'C' and seq2_rev[m] == 'G'),
                (seq1[n] == 'G' and seq2_rev[m] == 'C')
            ]):
                bounds += 'I'
                positions_bounds_seq1 += str(n) + " "
                positions_bounds_seq2 += str(m) + " "
            else:
                bounds += '-'
                positions_bounds_seq1 += "-"
                positions_bounds_seq2 += "-"

        bounds = bounds.replace('-I-', '-*-').replace('-I-', '-*-')
        if bounds[0:2] == 'I-':
            bounds = '*' + bounds[1::]
        if bounds[-2::] == '-I':
            bounds = bounds[0:-1] + '*'

        positions_bounds_seq1 = [
            pos.split() for pos in positions_bounds_seq1.split('-')
            if len(pos.strip()) > 2
        ]
        positions_bounds_seq2 = [
            pos.split() for pos in positions_bounds_seq2.split('-')
            if len(pos.strip()) > 2
        ]

        deltaG = 0.98
        deltaH = 0
        deltaS = -1.4

        if bounds.count('II') > 1:
            for chain_points in range(len(positions_bounds_seq1)):
                for num in range(len(positions_bounds_seq1[chain_points])):
                    nums = positions_bounds_seq1[chain_points][num:num + 2]
                    if len(nums) > 1:
                        pair = seq1[int(nums[0])] + seq1[int(nums[1])]
                    else:
                        continue

                    deltaG += NN_LIST_G[pair]
                    deltaH += NN_LIST_H[pair]
                    deltaS += NN_LIST_S[pair]

            n_dimers += 1
            lenght_dimer = bounds.count('I')
            up_chain = " " * m + "5'-" + seq1 + "-3'"
            dw_chain = " " * n + "3'-" + seq2_rev + "-5'"
            dimer = {
                'name': f"dimer {n_dimers}",
                "UP-chain": up_chain,
                "bounds": "   " + " " * max(n, m) + bounds,
                "DW-chain": dw_chain,
                "total_dG": round(deltaG, 2),
                "total_dH_dS": [round(deltaH, 2), round(deltaS, 2)],
                "length_dimer": lenght_dimer
            }
            result.append(dimer)

    result = sorted(result, key=lambda x: x['total_dG'])
    if max_positive_dG == 0 and max_negatrive_dG == 0:
        return result

    return list(filter(
        lambda x: max_negatrive_dG <= x['total_dG'] <= max_positive_dG,
        result
    ))


def temp_DNA_melt(Conc_DNA, Length_Seq=0, Conc_K=50, Conc_Mg=3,
                  dH=0, dS=0, type='hybridization'):
    """
    Calculate DNA melting temperature and Delta G.

    Parameters
    ----------
    Conc_DNA : float
        DNA concentration in uM
    Length_Seq : int
        Sequence length
    Conc_K : float
        Potassium concentration in mM
    Conc_Mg : float
        Magnesium concentration in mM
    dH : float
        Enthalpy change in kcal/mol
    dS : float
        Entropy change in cal/(mol*K)
    type : str
        'hybridization' or 'loop'

    Returns
    -------
    dict
        Dictionary with 'Tm' and 'dG' values
    """
    Tm = 0
    dH_cal = dH * 1000
    Conc_DNA_mol = Conc_DNA * 1e-6 / 4

    if Conc_Mg < 0:
        Conc_Mg = 0

    Conc_Salt = (Conc_K / 1000) + 4 * (Conc_Mg / 1000) ** 0.5
    Effective_Conc_Salt = Conc_Salt / (1.0 + 0.7 * Conc_Salt)
    Salt_Correction = 0.368 * (Length_Seq - 1) * log(Effective_Conc_Salt)
    dS_Conc_DNA_Correction = dS + 1.987 * numpy.log(Conc_DNA_mol)

    if any([dH_cal != 0, dS != 0]):
        if type == 'hybridization':
            Tm = dH_cal / (dS_Conc_DNA_Correction + Salt_Correction) - 273.15
            dG_val = (dH_cal - (37 + 273.15) * (dS_Conc_DNA_Correction + Salt_Correction)) / 1000
        elif type == 'loop':
            Salt_Correction = 16.6 * log10(Effective_Conc_Salt)
            Tm = dH_cal / dS + Salt_Correction - 273.15
            dG_val = (dH_cal - (37 + 273.15) * dS + Salt_Correction) / 1000

    return {'Tm': round(Tm, 2), 'dG': round(dG_val, 2)}


def loops_analyze(seq, max_negatrive_dG=0, max_positive_dG=0):
    """
    Analyze hairpin loop formations in DNA sequence.

    Parameters
    ----------
    seq : str
        DNA sequence
    max_negatrive_dG : float
        Minimum Delta G filter value
    max_positive_dG : float
        Maximum Delta G filter value

    Returns
    -------
    list
        List of loop dictionaries
    """
    result = []
    right_side_seq = list(seq)
    left_side_seq = []
    n_loops = 0

    NN_PARAMS = {
        'AT': -1.00, 'TA': -0.88, 'GC': -2.17, 'CG': -1.84,
        'AA': -0.93, 'TT': -1.11, 'AG': -1.10, 'GA': -1.43,
        'CT': -1.44, 'TC': -1.28, 'GT': -1.92, 'TG': -1.45,
        'CA': -1.28, 'AC': -1.44, 'GG': -2.08, 'CC': -1.91,
    }

    # Энергия инициализации и петель
    init_terminal = 0.98
    loop_delta_G = {3: 3.0, 4: 3.5, 5: 4.0, 6: 4.5, 7: 5.0, 8: 5.4, 9: 5.8, 10: 6.1}

    while right_side_seq:
        overhang = '│'
        for shift in range(2):
            bounds = ''
            positions_stem = []
            positions_loop = []
            for i in range(min([len(right_side_seq), len(left_side_seq)])):
                if any([
                    right_side_seq[i] == 'A' and left_side_seq[::-1][i] == 'T',
                    right_side_seq[i] == 'T' and left_side_seq[::-1][i] == 'A',
                    right_side_seq[i] == 'C' and left_side_seq[::-1][i] == 'G',
                    right_side_seq[i] == 'G' and left_side_seq[::-1][i] == 'C'
                ]):
                    bounds += 'I'
                    positions_stem.append(i)
                    positions_loop.append('n')
                else:
                    bounds += '-'
                    positions_stem.append('n')
                    positions_loop.append(i)

            bounds = '-' + bounds[1::]
            if positions_stem and positions_loop:
                positions_stem[0] = 'n'
                positions_loop[0] = 0
            bounds = bounds.replace('-I-', '-*-').replace('-I-', '-*-')

            if bounds[0:2] == 'I-':
                bounds = '*' + bounds[1::]
            if bounds[-2::] == '-I':
                bounds = bounds[0:-1] + '*'

            if bounds.count('II') >= 1:
                n_loops += 1

                bounds_positions = {
                    'loop': positions_loop,
                    'stem': positions_stem,
                    'overhang': overhang
                }
                end_position_stem = max([
                    s for s in bounds_positions['stem'] if s not in ['n', '*']
                ])

                # Параметры петли
                loop_params = {
                    3: {'dh': -2, 'ds': -10},
                    4: {'dh': -1.5, 'ds': -15},
                    5: {'dh': -1.0, 'ds': -18},
                    6: {'dh': -1.0, 'ds': -21},
                    7: {'dh': -1.0, 'ds': -24},
                    8: {'dh': -1.0, 'ds': -27},
                    9: {'dh': -1.0, 'ds': -30},
                    10: {'dh': -1.0, 'ds': -33},
                }

                # Параметры инициализации
                init_dh = 0.2
                init_ds = -5.7

                # Расчёт длины петли
                loop_symbols = (
                    [left_side_seq[::-1][value] for value in bounds_positions['loop']
                     if value != 'n' and value < end_position_stem][::-1] +
                    [bounds_positions['overhang'] if bounds_positions['overhang'] in 'ACTG' else ''] +
                    [right_side_seq[value] for value in bounds_positions['loop']
                     if value != 'n' and value < end_position_stem]
                )
                loop_GC = loop_symbols.count('G') + loop_symbols.count('C')
                loop_len = len(''.join(loop_symbols))

                if loop_len >= 3 or (loop_len == 1 and overhang in 'ACTG' and bounds.count('II') >= 2):
                    # Энергия петли
                    loop_energy = 0
                    if loop_len in loop_delta_G:
                        loop_energy = loop_delta_G[loop_len]
                    elif loop_len > 1 and loop_len not in loop_delta_G:
                        loop_energy = 6.0 + 0.3 * (loop_len - 8) if loop_len > 8 else 6.0
                    elif loop_len == 1:
                        loop_energy = loop_delta_G[3]
                        bounds = '-' + bounds[1::]

                    total_dh_loop = init_dh
                    total_ds_loop = init_ds

                    if loop_len in loop_params:
                        total_dh_loop += loop_params[loop_len]['dh']
                        total_ds_loop += loop_params[loop_len]['ds']
                    elif loop_len > 1 and loop_len not in loop_params:
                        base_dh = 5.0 + 0.2 * (loop_len - 8)
                        base_ds = 14.0 + 0.4 * (loop_len - 8)
                        total_dh_loop += base_dh
                        total_ds_loop += base_ds
                    elif loop_len == 1:
                        total_dh_loop += loop_params[3]['dh']
                        total_ds_loop += loop_params[3]['ds']
                        loop_len = 3

                    # Коррекция на GC-состав петли
                    gc_content = loop_GC / loop_len
                    loop_energy -= gc_content * 0.3
                    total_dh_loop -= gc_content * 0.5
                    total_ds_loop -= gc_content * 1.2

                    # АНАЛИЗ СТЕБЛЯ
                    stem_energy = init_terminal
                    total_dh_stem = 0.2
                    total_ds_stem = -5.7
                    choose_chain = {
                        len(right_side_seq): right_side_seq,
                        len(left_side_seq): left_side_seq[::-1]
                    }
                    stem_len = 0
                    for value in bounds_positions['stem']:
                        if value != 'n':
                            pair = ''.join(choose_chain[min(choose_chain)][value:value + 2])
                            stem_len += 1
                            if pair in NN_PARAMS:
                                stem_energy += NN_LIST_G[pair]
                                total_dh_stem += NN_LIST_H[pair]
                                total_ds_stem += NN_LIST_S[pair]
                            elif pair in ['A', 'T']:
                                stem_energy += -1.0
                            elif pair in ['C', 'G']:
                                stem_energy += -2.0

                    # Общая энергия
                    loop_energy = round(loop_energy, 2)
                    stem_energy = round(stem_energy, 2)
                    total_dh_stem = round(total_dh_stem, 2)
                    total_ds_stem = round(total_ds_stem, 2)
                    total_dh_loop = round(total_dh_loop, 2)
                    total_ds_loop = round(total_ds_loop, 2)
                    total_dh = total_dh_loop + total_dh_stem
                    total_ds = total_ds_loop + total_ds_stem
                    total_delta_g = total_dh / (total_ds / 1000) - 273.15
                    total_delta_g = round(total_delta_g, 2)

                    result.append({
                        'name': f"loop {n_loops}",
                        "5'-side": f"┌{''.join(left_side_seq[::-1])}-5'",
                        "bounds": f'{overhang}{bounds}',
                        "3'-side": f"└{''.join(right_side_seq)}-3'",
                        "loop_seq": ''.join(loop_symbols),
                        "total_dG": total_delta_g,
                        "stem_dG": stem_energy,
                        "loop_dG": loop_energy,
                        "stem_dH_dS": [total_dh_stem, total_ds_stem],
                        "loop_dH_dS": [total_dh_loop, total_ds_loop],
                        "total_dH_dS": [total_dh, total_ds],
                        "length_loop": loop_len,
                        "length_stem": stem_len
                    })

            if shift == 0 and right_side_seq:
                overhang = right_side_seq.pop(0)

        left_side_seq.append(overhang)

    result = sorted(result, key=lambda x: x['total_dG'], reverse=True)
    if max_positive_dG == 0 and max_negatrive_dG == 0:
        return result

    return list(filter(
        lambda x: max_negatrive_dG <= x['total_dG'] <= max_positive_dG,
        result
    ))


def test():
    """Test function for loops_analyze."""
    for loop in loops_analyze('TGAAGTACACCGGAATTGCCAGGAGACTTCA'):
        print(loop["5'-side"])
        print(loop['bounds'])
        print(loop["3'-side"])
        print(loop["loop_seq"])