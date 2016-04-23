import re
import os
from matplotlib import pyplot as plt
from matplotlib import style
from collections import defaultdict
import random
from itertools import permutations
import  numpy

class Strand:
    def __init__(self,aa):
        self.aa = aa
        self.P_values = []

    def reset(self):
        self.aa = ''
        self.P_values = []

class TotalStrand:
    def __init__(self):
        self.P_values = []
        self.P_number = 0
        self.angles = []
        self.translation = []
        self.average_translation = 0
        self.average_angle = 0
        self.step_total = 0
        self.all_count = 0
        self.aa_counts = {}
        self.step_distrib = {'G22': 0,
                             'G21': 0,
                             'G20': 0,
                             'G11': 0,
                             'G10': 0,
                             'G00': 0}
        self.step_distrib_percentage = {'G22': 0,
                                        'G21': 0,
                                        'G20': 0,
                                        'G11': 0,
                                        'G10': 0,
                                        'G00': 0}


    def reset(self):
        self.P_values = []
        self.P_number = 0
        self.angles = []
        self.translation = []
        self.average_translation = 0
        self.average_angle = 0
        self.step_total = 0
        self.all_count = 0
        self.aa_counts = {}
        self.step_distrib = {'G22': 0,
                             'G21': 0,
                             'G20': 0,
                             'G11': 0,
                             'G10': 0,
                             'G00': 0}
        self.step_distrib_percentage = {'G22': 0,
                                        'G21': 0,
                                        'G20': 0,
                                        'G11': 0,
                                        'G10': 0,
                                        'G00': 0}

class MetaData:
    def __init__(self):
        self.all_angle_diff = []
        self.all_t_diff = []

    def reset(self):
        self.all_angle_diff = []
        self.all_t_diff = []

strand1 = Strand("")
strand2 = Strand("")
strand3 = Strand("")

total = TotalStrand()
meta_data = MetaData()

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

        possible_aa = []
        generating = True

        while generating:

            temp = random.random() * 100

            for aa in aa_probs:
                if temp <= aa_probs[aa]:
                    possible_aa.append(aa)

            if possible_aa:
                generating =False
                g_strand += random.choice(possible_aa)
            else:
                generating =True

    return g_strand

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

        possible_aa = []
        generating = True

        while generating:

            temp = random.random() * 100

            for aa in aa_probs:
                if temp <= aa_probs[aa]:
                    possible_aa.append(aa)

            if possible_aa:
                generating = False
                g_strand += random.choice(possible_aa)
            else:
                generating = True

    return g_strand

def IncrementStrands():

    strand1.aa = strand1.aa[2:]
    strand2.aa = strand2.aa[1:-1]
    strand3.aa = strand3.aa[:-2]

def FindDomain(Gvalue):

    Gasked = "G.{" + str(Gvalue) + "}"

    strand1.aa = "".join(re.findall(3 * Gasked, strand1.aa))
    strand2.aa = "".join(re.findall(3 * Gasked, strand2.aa))
    strand3.aa = "".join(re.findall(3 * Gasked, strand3.aa))

def InitiateStrandsHomotrimer(sequence):
    strand1.aa = sequence
    strand2.aa = sequence
    strand3.aa = sequence

def InitiateStrandHeterotrimer(sequence1, sequence2, sequence3):
    strand1.aa = sequence1
    strand2.aa = sequence2
    strand3.aa = sequence3

class FileInfo:
    def __init__(self):
        self.folder = os.path.join('C:\\','Users','user','PycharmProjects','TwistAnalyser')
        self.file_name =''
        self.file_path= os.path.join('C:\\','Users','user','PycharmProjects','TwistAnalyser')

    def SetFile(self,file,filename):
        self.file_name = filename
        self.file_path = os.path.join('C:\\','Users','user','PycharmProjects','TwistAnalyser',file, filename)

    def reset(self):
        self.folder = os.path.join('C:\\','Users','user','PycharmProjects','TwistAnalyser')
        self.file_name =''
        self.file_path= os.path.join('C:\\','Users','user','PycharmProjects','TwistAnalyser')

def CountP():

    for c in strand1.aa:
        if c == "P":
            strand1.P_values.append(1)
        else:
            strand1.P_values.append(0)

    for c in strand2.aa:
        if c == "P":
            strand2.P_values.append(1)
        else:
            strand2.P_values.append(0)

    for c in strand3.aa:
        if c == "P":
            strand3.P_values.append(1)
        else:
            strand3.P_values.append(0)

    for v1, v2, v3 in zip(strand1.P_values, strand2.P_values, strand3.P_values):
        sum_P = v1 + v2 + v3
        total.P_values.append(sum_P)

def StepIteration():

    def G22():

        total.angles.append(-102.6)
        total.translation.append(2.84)
        total.step_distrib['G22'] += 1

    def G21():
        total.angles.append(-104.5)
        total.translation.append(2.86)
        total.step_distrib['G21'] += 1


    def G20():
        # no values for type 2???
        total.angles.append(-105.0)
        total.translation.append(2.86)
        total.step_distrib['G20'] += 1


    def G11():
        total.angles.append(-105.0)
        total.translation.append(2.86)
        total.step_distrib['G11'] += 1

    def G10():
        total.angles.append(-105.9)
        total.translation.append(2.89)
        total.step_distrib['G10'] += 1

    def G00():
        total.angles.append(-108.1)
        total.translation.append(2.90)
        total.step_distrib['G00'] += 1

    step_options = {'22': G22,
                    '21': G21,
                    '12': G21,
                    '20': G20,
                    '02': G20,
                    '11': G11,
                    '10': G10,
                    '01': G10,
                    '00': G00
                    }

    for index, value in enumerate(total.P_values):
        if index < len(total.P_values) - 1:
            temp = str(total.P_values[index]) + str(total.P_values[index + 1])
            if temp in step_options:
                step_options[temp]()
                total.P_number += value
            else:
                continue
        else:
            break

def CalcAverages():

    if len(total.angles) == 0:
        print("divided by zero")

    total.average_angle = numpy.mean(total.angles)
    total.average_translation = numpy.mean(total.translation)

    total.step_total = total.step_distrib['G22'] + total.step_distrib['G21'] + total.step_distrib['G20'] + \
                       total.step_distrib['G11'] + total.step_distrib['G10'] + total.step_distrib['G00']

def DrawGraph(run):

    position = list(range(len(total.angles)))

    plt.figure(1)

    style.use('ggplot')

    #local angle variations

    plt.subplot(221)
    plt.plot(position,total.angles,'b')
    plt.ylabel('Angles (K°)')
    plt.xlabel('Position')
    plt.title('Triple helical twist diagram')
    plt.axis([0,len(total.angles),-109,-102])

    #local translation variations

    plt.subplot(222)
    plt.plot(position,total.translation, 'b')
    plt.ylabel('Translation (t)')
    plt.xlabel('Position')
    plt.title('Translation diagram')
    plt.axis([0,len(total.translation),2.92,2.82])

    #aa composition

    barx1 = []
    bary1 = []

    for aa in total.aa_counts:
        barx1.append(aa)
        bary1.append(float("{0:.2f}".format(round(total.aa_counts[aa],2))))

    plt.subplot(224)
    plt.bar(range(len(bary1)),bary1,align='center')
    plt.xticks(range(len(bary1)),barx1)
    plt.ylabel('% Composition')
    plt.xlabel('Amino Acid')
    plt.title('Amino Acid composition')

    #step distibution bar chart

    barx2 = []
    bary2 = []

    for step in total.step_distrib:
        barx2.append(step)
        #turn it to percent
        total.step_distrib[step] = total.step_distrib[step]*100/total.step_total

        bary2.append(float("{0:.2f}".format(round(total.step_distrib[step],2))))

    plt.subplot(223)
    plt.bar(range(len(bary2)),bary2,align='center')
    plt.xticks(range(len(bary2)),barx2)
    plt.ylabel('% Frequency')
    plt.xlabel('Step')
    plt.title('Step Distribution')

    plt.tight_layout()

    file_open = FileInfo()
    file_open.SetFile("Graphs",run +".png")
    plt.savefig(file_open.file_path)

    plt.close()

def ClearMemory():

    strand1.reset()
    strand2.reset()
    strand3.reset()
    total.reset()

def main(run_number):

    angle_diff = 0
    t_diff = 0

    for i in range(run_number):

        InitiateStrandsHomotrimer(GenerateStrandNoGXY(4318))

        FindDomain(5)
        IncrementStrands()

        comb_angles = []
        comb_translation = []

        angles_max = 0
        angles_min = 0

        t_max = 0
        t_min = 0

        for combination in permutations((strand1.aa, strand2.aa, strand3.aa)):

            strand1.aa = combination[0]
            strand2.aa = combination[1]
            strand3.aa = combination[2]

            CountP()
            StepIteration()
            CalcAverages()

            comb_angles.append(total.average_angle)
            comb_translation.append(total.average_translation)
        else:

            angles_max = max(comb_angles)
            angles_min = min(comb_angles)

            t_max = max(comb_translation)
            t_min = min(comb_translation)

            angle_diff = angles_max - angles_min
            t_diff = t_max - t_min

            meta_data.all_angle_diff.append(angle_diff)
            meta_data.all_t_diff.append(t_diff)

        ClearMemory()

        print('Run' + str(i))

    angle_diff = numpy.mean(meta_data.all_angle_diff)
    t_diff = numpy.mean(meta_data.all_t_diff)

    print(angle_diff)
    print(t_diff)

main(500)