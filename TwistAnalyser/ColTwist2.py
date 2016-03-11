import os
import re
from matplotlib import pyplot as plt
from matplotlib import style
import numpy as np

class Strand:
    def __init__(self,aa):
        self.aa = aa
        self.P_values = []

    def reset(self):
        self.aa = ''
        self.P_values = []

class TotalStrand():
    def __init__(self):

        self.P_values = []
        self.P_number = 0
        self.angles = []
        self.translation = []
        self.average_translation =0
        self.average_angle =0
        self.step_total =0
        self.all_count =0
        self.aa_counts = {}
        self.step_distrib = {'G22' : 0,
        'G21' : 0,
        'G20' : 0,
        'G11' : 0,
        'G10' : 0,
        'G00' : 0}

    def reset(self):
        self.P_values = []
        self.P_number = 0
        self.angles = []
        self.translation = []
        self.average_translation =0
        self.average_angle =0
        self.step_total =0
        self.all_count =0
        self.aa_counts = {}
        self.step_distrib = {'G22' : 0,
        'G21' : 0,
        'G20' : 0,
        'G11' : 0,
        'G10' : 0,
        'G00' : 0}

class TwistData():
    def __init__(self):
        self.triplet_coord = np.zeros((3,1))
        self.axis_direction = np.zeros((3,1))


twist_calc = TwistData()
print(twist_calc.triplet_coord)

class FileInfo():
    def __init__(self):
        self.folder = os.path.join('C:\\','Users','user','PycharmProjects','TwistAnalyser','Data')
        self.file_name =''
        self.file_path= os.path.join('C:\\','Users','user','PycharmProjects','TwistAnalyser','Data')

    def SetFile(self,value):
        self.file_name = value
        self.file_path = os.path.join('C:\\','Users','user','PycharmProjects','TwistAnalyser','Data', value)

    def reset(self):
        self.folder = os.path.join('C:\\','Users','user','PycharmProjects','TwistAnalyser','Data')
        self.file_name =''
        self.file_path= os.path.join('C:\\','Users','user','PycharmProjects','TwistAnalyser','Data')

strand1 = Strand("")
strand2 = Strand("")
strand3 = Strand("")
total = TotalStrand()
file_open = FileInfo()

def OpenFile():

    with open(file_open.file_path, 'r', encoding='utf-8-sig') as file:
        count = 1

        def strand1line():
            strand1.aa += line

        def strand2line():
            strand2.aa += line

        def strand3line():
            strand3.aa += line

        file_options = {1 : strand1line,
               2 : strand2line,
               3 : strand3line,
        }

        for line in file:
            if not line.strip():
                count = 1
            else:
                line = line.replace(" ","")
                line = line.replace("\n", "")
                file_options[count]()
                count +=1

def PrepStrands():

    strand1.aa = "".join(re.findall(5 * "G.{2}", strand1.aa))
    strand2.aa = "".join(re.findall(5 * "G.{2}", strand2.aa))
    strand3.aa = "".join(re.findall(5 * "G.{2}", strand3.aa))

    strand1.aa = strand1.aa[2:]
    strand2.aa = strand2.aa[1:-1]
    strand3.aa = strand3.aa[:-2]

def CountP():

    for c in strand1.aa:
        if c == "P":
            strand1.P_values.append(1)
        else:
            strand1.P_values.append(0)

    for c in strand2.aa:
        if c=="P":
            strand2.P_values.append(1)
        else:
            strand2.P_values.append(0)

    for c in strand3.aa:
        if c=="P":
            strand3.P_values.append(1)
        else:
            strand3.P_values.append(0)

    for v1, v2, v3 in zip(strand1.P_values, strand2.P_values,strand3.P_values):

        sum_P = v1 + v2 + v3
        total.P_values.append(sum_P)

def StepIteration():

    def G22():

        total.angles.append(-102.6)
        total.translation.append(2.84)
        total.step_distrib['G22']+=1

    def G21():
        total.angles.append(-104.5)
        total.translation.append(2.86)
        total.step_distrib['G21']+=1

    def G20():
        #no values for type 2???
        total.angles.append(-105.0)
        total.translation.append(2.86)
        total.step_distrib['G20']+=1

    def G11():
        total.angles.append(-105.0)
        total.translation.append(2.86)
        total.step_distrib['G11']+=1

    def G10():
        total.angles.append(-105.9)
        total.translation.append(2.89)
        total.step_distrib['G10']+=1

    def G00():
        total.angles.append(-108.1)
        total.translation.append(2.90)
        total.step_distrib['G00']+=1

    step_options = {'22' : G22,
        '21' : G21,
        '12' : G21,
        '20' : G20,
        '02' : G20,
        '11' : G11,
        '10' : G10,
        '01' : G10,
        '00' : G00
    }

    for index, value in enumerate(total.P_values):
        if index<len(total.P_values)-1:
            temp = str(total.P_values[index]) + str(total.P_values[index+1])
            step_options[temp]()
            total.P_number += value
        else:
            break

def CalcAverages():
    sum_angles = 0
    sum_translations = 0

    for angle in total.angles:
        sum_angles += angle

    for distance in total.translation:
        sum_translations += distance

    total.average_angle=sum_angles/len(total.angles)
    total.average_translation = sum_translations/len(total.translation)

    total.step_total = total.step_distrib['G22']+total.step_distrib['G21']+total.step_distrib['G20']+total.step_distrib['G11']+total.step_distrib['G10']+total.step_distrib['G00']

def CalcComposition():

    singleletter_dict = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
                         'Trp': 'W', 'Thr': 'T', 'Asn': 'N', 'Pro': 'P', 'Phe': 'F',
                         'Ala': 'A', 'Gly': 'G', 'Ile': 'I', 'Leu': 'L', 'His': 'H',
                         'Arg': 'R', 'Met': 'M', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y',}

    total.all_count = len(strand1.aa)+len(strand2.aa)+len(strand3.aa)

    for aa in singleletter_dict.values():
        total.aa_counts[aa] = 0

    for v1, v2, v3 in zip(strand1.aa, strand2.aa,strand3.aa):
        total.aa_counts[v1] += 1
        total.aa_counts[v2] += 1
        total.aa_counts[v3] += 1

    for aa in total.aa_counts:
        total.aa_counts[aa]= total.aa_counts[aa]*100/total.all_count

def DrawGraph():
    position = list(range(len(total.angles)))

    plt.figure(1)

    style.use('ggplot')

    plt.subplot(221)
    plt.plot(position,total.angles,'b')
    plt.ylabel('Angles (K°)')
    plt.xlabel('Position')
    plt.title('Triple helical twist diagram')
    plt.axis([0,len(total.angles),-109,-102])

    plt.subplot(222)
    plt.plot(position,total.translation, 'b')
    plt.ylabel('Translation (t)')
    plt.xlabel('Position')
    plt.title('Translation diagram')
    plt.axis([0,len(total.translation),2.92,2.82])

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

    plt.show()

def SaveResults():
    file_open.file_name = re.sub(r"\.(.*)", "", file_open.file_name)

    file_open.SetFile(file_open.file_name + " results.txt")

    print(file_open.file_path)

    with open(file_open.file_path, 'w', encoding='utf-8-sig') as file:
        file.truncate()
        file.write("Step\tFrequency (%)" + "\n")
        for key in total.step_distrib:
            file.write(key + " \t"+ str(float("{0:.2f}".format(round(total.step_distrib[key],2))))+ "\n")

        file.write("Total " + "\t" + str(total.step_total) + "\n"+ "\n"
        "Average Angle " + "\t" + str(float("{0:.2f}".format(round(total.average_angle,2)))) + "\n"
        "Average Translation " + "\t" +  str(float("{0:.2f}".format(round(total.average_translation,2)))) + "\n"+ "\n"
        "Strand 1 length " + "\t" + str(len(strand1.aa)) + "\n"
        "Strand 2 length " + "\t" + str(len(strand2.aa)) + "\n"
        "Strand 3 length " + "\t" + str(len(strand3.aa)) + "\n"+ "\n"
        "Number of P in Strand " + "\t" + str(total.P_number) + "\n"+ "\n"
        "aa" + "\t" + "%composition" +"\n")

        for aa in total.aa_counts:
            file.write(aa + "\t" + "{0:.2f}".format(round(total.aa_counts[aa],2))+"\n")

def ClearMemory():
    strand1.reset()
    strand2.reset()
    strand3.reset()
    total.reset()

def Main():

    for fn in os.listdir(file_open.folder):

        if re.search(r"\bresults\b", fn):
            continue
        else:

            file_open.SetFile(fn)

            OpenFile()
            PrepStrands()
            CountP()
            StepIteration()
            CalcAverages()
            CalcComposition()
            DrawGraph()
            SaveResults()

            print(file_open.file_name)
            print(total.step_distrib)
            print(strand1.aa)
            print(strand2.aa)
            print(strand3.aa)
            print(strand1.P_values)
            print(strand2.P_values)
            print(strand3.P_values)
            print(total.P_values)
            print(total.P_number)
            print(total.angles)
            print(total.average_angle)
            print(total.average_translation)
            print("\n")
        ClearMemory()

Main()