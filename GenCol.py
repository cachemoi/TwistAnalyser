import random

def GenerateStrandGXY(length):

    strand_length = range(length)

    aa_probs = {'L': 3.35615362884089, 'I': 1.589963710163489, 'T': 2.087807704977877, 'H': 0.6760663957883509,
                'F': 1.1763417532513833, 'S': 3.3691109171343636, 'R': 4.7552093171840255, 'C': 0.022846209524532153,
                'P': 21.08678244043536, 'Z': 0.0038505968425105886, 'Y': 0.27857212746094984, 'V': 2.2430820247461045,
                'M': 0.9600225957747107, 'X': 0.0, 'N': 1.5045087984115857, 'D': 3.384144345152309,
                'W': 0.01593087589089932, 'E': 4.851703587206036, 'Q': 3.381052927323254, 'K': 4.25069111103824,
                'G': 33.86050181266035, 'A': 7.145657120192756}


    g_strand = ""

    count = 1

    for i in strand_length:
        print(count)
        if count == 3:
            count = 1
            g_strand += "G"
            continue
        else:
            count += 1

        print("aa added")

        possible_aa = []
        generating = True

        while generating:

            temp = random.random() * 100
            print(temp)

            for aa in aa_probs:
                if temp <= aa_probs[aa]:
                    possible_aa.append(aa)

            if possible_aa:
                print("succeeded")
                generating =False
                g_strand += random.choice(possible_aa)
            else:
                generating =True

    print(g_strand)

def GenerateStrandNoGXY(length):

    strand_length = range(length)

    aa_probs = {'K': 4.5691937022889535, 'W': 0.459212404800139, 'I': 2.586021574823164, 'E': 5.264748403057716,
                'Y': 1.2336124200299565, 'D': 4.313070095440298, 'M': 1.1710329204189691, 'L': 4.7317372968599845,
                'X': 0.052521008403361345, 'G': 25.275004243192615, 'N': 2.4167359608785555, 'F': 1.9745792281975412,
                'Z': 0.0028129790855004996, 'Q': 3.804279081972282, 'P': 16.659719318954426, 'S': 4.66288442917996,
                'V': 3.3419053385685875, 'C': 0.990903930140905, 'T': 3.567952664801641, 'A': 7.13612636494623,
                'H': 1.07123068268337, 'R': 4.714715951275851}

    g_strand = ""

    for i in strand_length:

        print("aa added")

        possible_aa = []
        generating = True

        while generating:

            temp = random.random() * 100
            print(temp)

            for aa in aa_probs:
                if temp <= aa_probs[aa]:
                    possible_aa.append(aa)

            if possible_aa:
                print("succeeded")
                generating = False
                g_strand += random.choice(possible_aa)
            else:
                generating = True

    print(g_strand)

if __name__ == "__main__":
    GenerateStrandNoGXY(4318)