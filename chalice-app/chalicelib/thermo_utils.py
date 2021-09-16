import random as rnd

# Quick, approximate Tm calculation for a DNA sequence without using nupack
def TM(sequence):
    if len(sequence) < 1:
        return 1.0
    GC_count = float(sequence.count('G') + sequence.count('C'))
    AT_count = float(sequence.count('A') + sequence.count('T'))
    if len(sequence) > 13:
        # Wallace Formula - Wallace RB et al. (1979) Nucleic Acids Res 6:3543-3557, PMID 158748
        return 64.9 +41.0*(GC_count-16.4)/(GC_count+AT_count)
    else:
        # Marmur Formula
        return float(4*GC_count + 2*AT_count)

# Returns the reverse complement of the input DNA sequence
def reverse_complement(sequence):
    sequence = sequence[::-1].upper()
    sequence = sequence.replace('A','1').replace('T','2').replace('G','3').replace('C','4')
    sequence = sequence.replace('1','T').replace('2','A').replace('3','C').replace('4','G')
    return sequence

# Truncate a DNA sequence to a specific length
def truncate(sequence, end, trunc_size=1):
    length = len(sequence)
    if end == 3:
        sequence = sequence[:length-trunc_size]
    elif end == 5:
        sequence = sequence[trunc_size:]
    else:
        print("Improper args given to truncate()!")
    return sequence

# Randomly truncate a DNA sequence
def generate_truncations(seqlength,minlength):
    # Probe truncations
    trunc_length = rnd.randint(0,seqlength-minlength)
    five_trunc = rnd.randint(0,trunc_length)
    output = [five_trunc, trunc_length - five_trunc]
    # Probe complement truncation
    if rnd.randint(0,1) == 0:
        try:
            output += [output[1], rnd.randint(output[0],seqlength-minlength-trunc_length), 3]
        except(ValueError):
            output += [output[1], output[0], 3]
    else:
        try:
            output += [rnd.randint(output[1],seqlength-minlength-trunc_length), output[0], 5]
        except(ValueError):
            output += [output[1], output[0], 5]
    #Sink truncations
    trunc_length = rnd.randint(0,seqlength-minlength)
    five_trunc = rnd.randint(0,trunc_length)
    output += [five_trunc, trunc_length - five_trunc]
    # Sink complement truncation
    trunc_length = rnd.randint(trunc_length,seqlength-minlength)
    five_trunc = rnd.randint(output[6],trunc_length)
    output += [five_trunc, trunc_length - five_trunc]
    return output