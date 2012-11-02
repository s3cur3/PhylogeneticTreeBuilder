def getConsensusLetter(char0, char1):
    '''
    @return The FASTA consensus letter of 2 mismatched characters
    '''
    if len(char0) != 1 or len(char1) != 1:
        print("Error in characters!")
        exit()

    if char0 == char1:
        return char0

    purines = "AG"
    if char0 in purines and char1 in purines:
        return 'R'

    pyrimidines = "CTU"
    if char0 in pyrimidines and char1 in pyrimidines:
        return 'Y'

    withAminoGroups = "AC"
    if char0 in withAminoGroups and char1 in withAminoGroups:
        return 'M'

    strongInteraction = "CG"
    if char0 in strongInteraction and char1 in strongInteraction:
        return 'S'

    weakInteraction = "ATU"
    if char0 in weakInteraction and char1 in weakInteraction:
        return 'W'

    ketones = "GTU"
    if char0 in ketones and char1 in ketones:
        return 'K'

    # If we get here, one or more characters were not A, C, G, T
    actu = "YM" + pyrimidines + withAminoGroups
    if char0 in actu and char1 in actu:
        return 'H'

    acg = "MS" + strongInteraction + withAminoGroups
    if char0 in acg and char1 in acg:
        return 'V'

    agtu = "RW" + weakInteraction + purines
    if char0 in agtu and char1 in agtu:
        return 'D'

    cgtu = "SK" + strongInteraction + ketones
    if char0 in cgtu and char1 in cgtu:
        return 'B'

    indeterminate = "-X"
    if char0 == 'N' or char1 == 'N' \
        and char0 not in indeterminate and char1 not in indeterminate:
        return 'N'

    # At this point, one or more characters is either X, or -, so it's indeterminate
    return 'X'