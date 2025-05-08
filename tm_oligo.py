from math import *
import numpy

NN_list_G = {'AA': -0.55, 'TT': -0.55, 'CC': -1.25, 'GG': -1.25, 'AT': -0.28, 'TA': -0.16, 'AC': -0.89, 'CA':
             -1.00, 'AG': -0.91, 'GA': -0.87, 'CG': -1.25, 'GC': -1.31, 'TC': -0.87, 'CT': -0.91, 'TG': -1.00, 'GT': -0.89}

NN_list_S = {'AA': -19.2, 'TT': -19.2, 'CC': -8.9, 'GG': -8.9, 'AT': -29.4, 'TA': -13.3, 'AC': -26.8, 'CA': -38.8,
             'AG': -7.9, 'GA': -13.0, 'CG': -16.1, 'GC': -9.3, 'TC': -13.0, 'CT': -7.9, 'TG': -38.8, 'GT': -26.8}

NN_list_H = {'AA': -6.5, 'TT': -6.5, 'CC': -4.0, 'GG': -4.0, 'AT': -9.4, 'TA': -4.3, 'AC': -13.1, 'CA':
             -13.1, 'AG': -3.4, 'GA': -4.9, 'CG': -6.4, 'GC': -4.2, 'TC': -4.9, 'CT': -3.4, 'TG': -13.1, 'GT': -9.2}


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


def temp_DNA_melt(sequence, CtDNA, CtK, CtMg):
    if len(sequence) > 0:
        salt = (CtK/1000) + 4 * (CtMg/1000)**0.5
        deltaH = deltaH_DNA(sequence)*1000
        deltaG = deltaG_DNA(sequence)*1000
        Tm = ((298.2*deltaH)/(deltaH-deltaG+(1.99*298.2)*log(CtDNA/1000000))
              ) + (16.6*log10(salt/(1.0+0.7*salt))) - 269.3
        return Tm


def GC_features(sequence):
    countG = sequence.count('G')
    countC = sequence.count('C')
    result = (countG+countC)/len(sequence)
    return result


def dimers_analyze(seq1, seq2):
    result = []
    rows = len(seq1)
    cols = len(seq2)
    seq2 = seq2[::-1]
    for i in range(1, cols + rows):
        sub_result = []
        bounds = ''
        start_column = max(0, i - rows)
        end_column = min(i, cols)
        for j in range(start_column, end_column):
            n = rows-i-1+j+1
            m = j
            if (seq1[n] == 'A' and seq2[m] == 'T') or (seq1[n] == 'T' and seq2[m] == 'A'):
                bounds += 'I'
            elif (seq1[n] == 'C' and seq2[m] == 'G') or (seq1[n] == 'G' and seq2[m] == 'C'):
                bounds += 'I'
            else:
                bounds += '-'

        bounds = bounds.replace('-I-', '-*-')

        if bounds.count('I') > 2:
            sub_result.append(f"{' ' * (m) + "5'-" + seq1}-3'")
            sub_result.append(f"   {' ' * max(n, m) + bounds}")
            sub_result.append(f"{' ' * (n) + "3'-" + seq2}-5'")
            result.append(sub_result)

    return result
