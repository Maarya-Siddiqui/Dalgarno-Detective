import pandas as pd
# Maarya Siddiqui
# CPS501 - F2024

seq = {"thrA":"GTAACGAGGTAACAACCatg",
        "yaaJ":"CAATAAGAGGGATATGCatg",
        "satP":"ATGATTTTTGAGGAATTatg",
        "amyA":"TTGAATGGAGTCGAAGAatg",
        "bioA":"TGGTTTACAAGTCGATTatg",
        "ding":"GTATTTCAGGTTTTCTCatg",
        "folA":"TTTATCGGGAAATCTCAatg",
        "murP":"TAGACAAGGAATAACCCatg",
        "rhlB":"TGATATTCTACCACACTatg",
        "proP":"TGATATTCTACCACACTatg",
        "yjdN":"TAACTGAGGTCACCATCatg",
        "loiP":"TGTGCGGAAGAGAAAACatg",
        "mscS":"TTGAAAAGGAATATTGAatg",
        "insH11":"ATAGAACAGGGTTCATCatg",
        "gadE":"CGGCTAAGGAGCAAGTTatg",
        "iraD":"CGTCAGGAGTGCGCAAAatg",
        "hypT":"TGACAGAGTAAAACGTAatg",
        "rapA":"TATAGAGCCGAACACATatg",
        "tsr":"CACAGGAAAGAGAAACCatg",
        "sucB":"TAAATAAAGGATACACAatg"}

TARGET = "AGGAGG"

def shineDalgarno(seq):
    '''Looks for 'AGGAGG' starting at the first letter onwards, counting 
    the mismatches each time, while keeping track of the position where the best match occurs. 

    Args: Dictionary of gene name and upstream sequence to locate AGGAGG in

    Creates DataFrame using Pandas containing the following info: 
        - Gene Name
        - Shine-Dalgarno Sequence
        - Number of mismatches in best case
        - Seperation between Shine and Start
        - 17 Upstream Bases

    Returns: List containing the Shine sequence for each gene
    '''

    data = []
    shines = []
    for gene, sequence in seq.items():

        final_mm=0
        sd_pos = -1

        for i in range(len(sequence)-len(TARGET)+1):
            count=0
            for j in range(len(TARGET)):                            
                if TARGET[j] != sequence[i+j]:
                    count += 1
            # first comparison --> note number of mismatches for baseline comparison
            if i == 0:                                              
                final_mm=count
            else:
                if count < final_mm:
                    final_mm=count
                    sd_pos=i

        # Calculating length of seperation Shine to Start codon
        start = sd_pos+len(TARGET)
        stop = sequence.find("atg")
        seperation = stop-start if stop > start else None

        # Determining shine sequence
        shine = sequence[sd_pos:(sd_pos+6)]
        shines.append(shine)
        up_b = sequence[:-3]

        # Storing info for each sequence in a list, then converting it to a Pandas Dataframe
        info = {'Gene':gene, 'Shine':shine, 'Mismatches':final_mm, 'Seperation': seperation, '17 Upstream Bases': up_b}
        data.append(info)

        table=pd.DataFrame(data)
        table.index = range(1, len(table)+1)

    print("\n")
    print(table)
    print("\n")

    # Calculating the mean and standard deviation of the numerical columns (Mismatches and Seperation)
    print("STATISTICS")
    stats = table.describe().loc[['mean', 'std']]
    print(stats)
    print("\n")

    return shines

def get_consensus(shines):
    '''This function calculates the percentage of each letter at each position in all shine sequences obtained, and returns the
    consensus Shine-Dalgarno sequence'''
    
    # Preparing a matrix transpose to compare all the letters in same position
    matrix = []
    for seq in shines:
        matrix.append([letter for letter in seq])
    transpose = list(map(list, zip(*matrix)))

    consensus = ""
    total = len(shines)
    for i in transpose:
        a_count, t_count, c_count, g_count = 0, 0, 0, 0
        for letter in i:
            if letter == "A":
                    a_count += 1
            elif letter == "T":
                t_count += 1
            elif letter == "C":
                c_count += 1
            elif letter == "G":
                g_count += 1
        # Determine average for each letter in the given position, i and convert them to percentages
        a_avg = round(a_count/total,2)*100
        t_avg = round(t_count/total,2)*100
        c_avg = round(c_count/total,2)*100
        g_avg = round(g_count/total,2)*100

        # Based on the maximum value in the dictionary, obtain the corresponding key and concatenate it to consensus string for final Shine Delgarno sequence
        data = {"A":a_avg, "T": t_avg, "C": c_avg, "G":g_avg}
        print(str(data)+"\t"+max(data, key=data.get))
        final_letter = max(data, key=data.get)
        consensus+=final_letter

    return consensus


def main():
    '''Mainline for program'''
    shines = shineDalgarno(seq)
    final_shine = get_consensus(shines)
    print("\nCONSENSUS SHINE DELGARNO SEQUENCE: "+final_shine)

main()



