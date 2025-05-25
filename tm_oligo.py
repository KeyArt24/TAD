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
    CtKeq = Ct * numpy.exp(DeltaS/R - DeltaH /
                           (R*T-16.6*log10(salt/(1.0+0.7*salt))))
    # Compute f
    f = (1 + CtKeq - numpy.sqrt(1 + 2*CtKeq)) / CtKeq

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
        sum_S += 0
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
        sum_H += 0
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


# def temp_DNA_melt(sequence, CtDNA, CtK, CtMg):
    if len(sequence) > 0:
        salt = (CtK/1000) + 4 * (CtMg/1000)**0.5
        deltaH = deltaH_DNA(sequence)*1000
        deltaG = deltaG_DNA(sequence)*1000
        Tm = ((298.2*deltaH)/(deltaH-deltaG+(1.99*298.2)*log(CtDNA/1000000))
              ) + (16.6*log10(salt/(1.0+0.7*salt))) - 269.3
        return Tm


def GC_features(sequence: str):
    countG = sequence.upper().count('G')
    countC = sequence.upper().count('C')
    result = (countG+countC)/len(sequence)
    return round(result*100)


def dimers_analyze(seq1: str, seq2: str):
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
        start_column = max(0, i - rows)
        end_column = min(i, cols)
        for j in range(start_column, end_column):
            n = rows-i-1+j+1
            m = j
            if (seq1[n] == 'A' and seq2[m] == 'T') or (seq1[n] == 'T' and seq2[m] == 'A'):
                bounds += 'I'
                deltaG += str(n) + " "
            elif (seq1[n] == 'C' and seq2[m] == 'G') or (seq1[n] == 'G' and seq2[m] == 'C'):
                bounds += 'I'
                deltaG += str(n) + " "
            else:
                bounds += '-'
                deltaG += '-'

        bounds = bounds.replace('-I-', '-*-')
        deltaG = [position.split() for position in deltaG.split('-')
                  if len(position.strip()) > 2]
        if len(deltaG) != 0:
            deltaG = sum([sum([NN_list_G[seq1.replace('+', '')[int(numbers):int(numbers)+2]]
                         for numbers in position if len(seq1[int(numbers):int(numbers)+2]) > 1]) for position in deltaG])

        if bounds.count('I') > 2:
            sub_result.append(
                f"{' ' * (m) + "5'-" + seq1.replace('+', '')}-3'")
            sub_result.append(f"   {' ' * max(n, m) + bounds}")
            sub_result.append(
                f"{' ' * (n) + "3'-" + seq2.replace('+', '')}-5'")
            sub_result.append(
                f"\n deltG {deltaG} ккал/моль")
            result.append(sub_result)

    return result


def temp_DNA_melt(sequence, CtDNA, CtK, CtMg):
    if len(sequence) > 0:
        salt = (CtK/1000) + 4 * (CtMg/1000)**0.5
        deltaH = deltaH_DNA(sequence)*1000
        deltaG = deltaG_DNA(sequence)*1000
        deltaS = deltaS_DNA(sequence)
        Tm = (deltaH/(deltaS+1.987*log(CtDNA/1000000))) + \
            (16.6*log10(salt/(1.0+0.7*salt))) - 273.15
        return Tm
