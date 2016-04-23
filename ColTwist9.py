#!/home/max/anaconda3/bin/python3

import cgitb
cgitb.enable()
import sys, json
import os
import re
from matplotlib import pyplot as plt
from matplotlib import style
from collections import defaultdict

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
        self.step_distrib_percentage = {'G22' : 0,
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
        self.step_distrib_percentage = {'G22' : 0,
        'G21' : 0,
        'G20' : 0,
        'G11' : 0,
        'G10' : 0,
        'G00' : 0}

class MetaData:
    def __init__(self):
        self.all_angles_averages = []
        self.all_average_angles = 0
        self.all_translation_averages = []
        self.all_average_translation = 0
        self.all_aa_count = {}
        self.all_aa_count_percentage = {}
        self.all_aa_total = 0
        self.all_steps = {'G22' : 0,
        'G21' : 0,
        'G20' : 0,
        'G11' : 0,
        'G10' : 0,
        'G00' : 0}
        self.all_steps_total =0
        self.all_steps_percentage = {'G22' : 0,
        'G21' : 0,
        'G20' : 0,
        'G11' : 0,
        'G10' : 0,
        'G00' : 0}

        self.species_angle =defaultdict(list)
        self.species_translation =defaultdict(list)
        self.species_average_angle = defaultdict(list)
        self.species_average_translation = defaultdict(list)

        #aa_count  is a percentage
        self.species_aa_count = defaultdict(list)
        self.species_aa_total = {}
        self.species_aa_average = {}

        self.species_steps = {}
        self.species_steps_total ={}
        self.species_steps_percentage = defaultdict(dict)

    def reset(self):
        self.all_angles_averages = []
        self.all_average_angles = 0
        self.all_translation_averages = []
        self.all_average_translation = 0

        self.all_aa_count = {}
        self.all_aa_count_percentage = {}
        self.all_aa_total = 0

        self.all_steps = {'G22' : 0,
        'G21' : 0,
        'G20' : 0,
        'G11' : 0,
        'G10' : 0,
        'G00' : 0}
        self.all_steps_total =0
        self.all_steps_percentage = {'G22' : 0,
        'G21' : 0,
        'G20' : 0,
        'G11' : 0,
        'G10' : 0,
        'G00' : 0}

        self.species_angle =defaultdict(list)
        self.species_translation =defaultdict(list)
        self.species_average_angle = defaultdict(list)
        self.species_average_translation = defaultdict(list)

        #aa count is a percentage
        self.species_aa_count = defaultdict(list)
        self.species_aa_total = {}
        self.species_aa_average = {}

        self.species_steps = {}
        self.species_steps_total ={}
        self.species_steps_percentage = defaultdict(dict)

    def add_species_step(self,species):
        self.species_steps[species] = {'G22' : 0,
        'G21' : 0,
        'G20' : 0,
        'G11' : 0,
        'G10' : 0,
        'G00' : 0}

        self.species_steps_percentage[species] = {'G22' : 0,
        'G21' : 0,
        'G20' : 0,
        'G11' : 0,
        'G10' : 0,
        'G00' : 0}

        self.species_steps_total[species] = 0

class FileInfo:
    def __init__(self):
        self.folder = os.path.join('var','www','html')
        self.file_name =''
        self.file_path= os.path.join('var','www','html')

    def SetFile(self,file,filename):
        self.file_name = filename
        self.file_path = os.path.join('var','www','html',file, filename)

    def reset(self):
        self.folder = os.path.join('var','www','html')
        self.file_name = ''
        self.file_path = os.path.join('var','www','html')

class Sequences:
    def __init__(self):
        self.sequence = {}
        self.species = {}

class CommandOUT:
    def __init__(self):
        self.MD = False
        self.SPMD = False
        self.SP = []
        self.P = []
        self.OAA = 0
        self.OAT = 0

    def reset(self):
        self.MD = False
        self.SPMD = False
        self.SP = []
        self.P = []
        self.OAA = 0
        self.OAT = 0

    def to_JSON(self):
        return json.dumps(self, default=lambda o: o.__dict__,sort_keys=True, indent=4)

strand1 = Strand("")
strand2 = Strand("")
strand3 = Strand("")
total = TotalStrand()
meta_data = MetaData()
sequences = Sequences()

def PullPDB(codes=None, DOM = False):

    if codes is None:
        asked = []
    else:
        asked = codes

    if DOM is True:
        domains = True
    else:
        domains =False

    import json
    import sys

    PY3 = sys.version > '3'

    if PY3:
        import urllib.request as urllib2
    else:
        import urllib2

    SERVER_URL = "https://www.ebi.ac.uk/pdbe/api"

    FASTA_CALL = "pdb/entry/files/"
    SPECIE_CALL = "pdb/entry/experiment"

    class SequenceList:
        def __init__(self):
            self.url = {}

    sequence_list = SequenceList()

    def make_request(url, data):

        request = urllib2.Request(url)


        url_file = urllib2.urlopen(request, data)


        return url_file.read().decode()

    def get_request(url, arg, pretty=False):
        full_url = "%s/%s/%s?pretty=%s" % (SERVER_URL, url, arg, str(pretty).lower())

        return make_request(full_url, None)

    def post_request(url, data, pretty=False):

        full_url = "%s/%s/?pretty=%s" % (SERVER_URL, url, str(pretty).lower())

        if isinstance(data, (list, tuple)):
            data = ",".join(data)

        return make_request(full_url, data.encode())

    def PrepareSequence():
        for c in sequences.sequence:
            temp = sequences.sequence[c].splitlines()
            temp_seq = ""

            for index, value in enumerate(temp):
                if index%2 == 1:
                    temp_seq += "|"
                    temp_seq += value

            sequences.sequence[c] = temp_seq

    sequence_list = SequenceList()

    def main(codes, DOM):

        if DOM is False:
            asking = codes
        else:
            asking= ['1CAG','1CGD','1A3I','1A3J','1G9W','1ITT','1K6F','2CUO','3AH9','1BKV','1QSU',
                      '1NAY','1Q7D','1V7H','4OY5','1V4F','1V6Q','3B0S','1X1K','1WZB','1YM8','1EI8',
                      '2G66','2DRT','2DRX','3DMW','2D3F','2D3H','3A0A','2KLW','3A08','3A19','3IPN',
                      '3ADM','3A1H','3A0M','3POD','3PON','3P46','3T4F','3U29','3B2C','4AXY','4DMT',
                      '4GYX','3ABN','3WN8']

        response_fasta = post_request(FASTA_CALL, asking, True)

        if response_fasta:
            data = json.loads(response_fasta)

            collected_codes = []

            for value in data:
                collected_codes.append(value)

            for code in collected_codes:
                for value in data[code]['molecule']['downloads']:
                    sequence_list.url[code] = value['url']

                    FASTA_req = urllib2.urlopen(value['url'])

                    FASTA = FASTA_req.read().decode('utf-8')

                    sequences.sequence[code] = FASTA

        PrepareSequence()

        response_specie = post_request(SPECIE_CALL, asking, True)


        if response_specie:

            data = json.loads(response_specie)
            collected_codes = []

            for value in data:
                collected_codes.append(value)

            for code in collected_codes:
                for prot_dat in data[code]:
                    if not prot_dat["expression_host_scientific_name"]:
                        sequences.species[code] = "Source Unknown"
                    else:
                        sequences.species[code] = prot_dat["expression_host_scientific_name"][0]["scientific_name"]

    main(asked, domains)

def PullUP(domains =False, input=None):

    if input is None:
        input = []

    import urllib.error
    import urllib.request
    import urllib.parse
    import re

    url = 'http://www.uniprot.org/uniprot/'

    class Codes :
        def __init__(self):
            self.codes =[]

    codes = Codes()

    def GetDomain():
        params = {
        'from' :'ACC',
        'to':'P_REFSEQ_AC',
        'format':'list',
        'query': 'domain:fibrillar collagen'
        }

        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        while True:
            try:
                request = urllib.request.urlopen(url, data)
            except urllib.error.HTTPError as err:
                  if err.code == 503:
                    continue
            else:
                domain_codes = request.read().decode('utf-8')
                newlines = re.compile(r"\n")

                query_id = re.split(newlines, domain_codes)
                del query_id[-1]

                for code in query_id:
                    codes.codes.append(code)
                break

    def PrepareSequence(input,code):
        newlines = re.compile(r"\n")

        split_FASTA = re.split(newlines, input)

        specie = re.findall(r"OS=(\w+\s\w+)", split_FASTA[0])
        specie = str(specie[0])

        del split_FASTA[0]

        split_FASTA = "".join(split_FASTA)

        sequences.sequence[code]=split_FASTA

        sequences.species[code] =specie

    def GetFASTA():
        for code in codes.codes:
            params = {
            'from':'ACC',
            'to':'P_REFSEQ_AC',
            'format':'fasta',
            'query':code
            }

            data = urllib.parse.urlencode(params)
            data = data.encode('utf-8')
            while True:
                try:
                    request = urllib.request.urlopen(url, data)

                except urllib.error.HTTPError as err:
                      if err.code == 503:
                        continue
                else:
                    FASTA = request.read().decode('utf-8')
                    PrepareSequence(FASTA,code)
                    break

    def GetByID(code):
        params = {
        'from':'ACC',
        'to':'P_REFSEQ_AC',
        'format':'fasta',
        'query':code
        }

        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        while True:
            try:
                request = urllib.request.urlopen(url, data)

            except urllib.error.HTTPError as err:
                  if err.code == 503:
                    continue
            else:
                FASTA = request.read().decode('utf-8')
                PrepareSequence(FASTA,code)
                break

    if domains == True:
        GetDomain()
        GetFASTA()
    else:
        for code in input:
            GetByID(code)

def ProcessFASTA(FASTA):

    strand1.aa = FASTA
    strand2.aa = FASTA
    strand3.aa = FASTA

def ProcessFASTAPDB(FASTA):

    temp = FASTA.split("|")
    temp = list(filter(None, temp))

    if len(temp) == 3:
        strand1.aa = temp[0]
        strand2.aa = temp[1]
        strand3.aa = temp[2]
    elif len(temp) == 2:
        strand1.aa = temp[0]
        strand2.aa = temp[0]
        strand3.aa = temp[1]
    elif len(temp) == 1:
        strand1.aa = temp[0]
        strand2.aa = temp[0]
        strand3.aa = temp[0]

def ProcessFASTAUI(FASTA):

    temp = FASTA.split("|")
    temp = list(filter(None, temp))

    if len(temp) == 3:
        strand1.aa = temp[0]
        strand2.aa = temp[1]
        strand3.aa = temp[2]
    elif len(temp) == 2:
        strand1.aa = temp[0]
        strand2.aa = temp[0]
        strand3.aa = temp[1]
    elif len(temp) == 1:
        strand1.aa = temp[0]
        strand2.aa = temp[0]
        strand3.aa = temp[0]

def IncrementStrands():

    strand1.aa = strand1.aa[2:]
    strand2.aa = strand2.aa[1:-1]
    strand3.aa = strand3.aa[:-2]

def FindDomain(Gvalue):

    Gasked = "G.{" + str(Gvalue) + "}"

    strand1.aa = "".join(re.findall(3 * Gasked, strand1.aa))
    strand2.aa = "".join(re.findall(3 * Gasked, strand2.aa))
    strand3.aa = "".join(re.findall(3 * Gasked, strand3.aa))

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

def StepIteration(code):

    if not sequences.species[code] in meta_data.species_steps:
        meta_data.add_species_step(sequences.species[code])

    def G22():

        total.angles.append(-102.6)
        total.translation.append(2.84)
        total.step_distrib['G22']+=1

        meta_data.all_steps['G22']+=1
        meta_data.species_steps[sequences.species[code]]['G22']+= 1

    def G21():
        total.angles.append(-104.5)
        total.translation.append(2.86)
        total.step_distrib['G21']+=1

        meta_data.all_steps['G21']+=1
        meta_data.species_steps[sequences.species[code]]['G21']+= 1

    def G20():
        #no values for type 2???
        total.angles.append(-105.0)
        total.translation.append(2.86)
        total.step_distrib['G20']+=1

        meta_data.all_steps['G20']+=1
        meta_data.species_steps[sequences.species[code]]['G20']+= 1

    def G11():
        total.angles.append(-105.0)
        total.translation.append(2.86)
        total.step_distrib['G11']+=1

        meta_data.all_steps['G11']+=1
        meta_data.species_steps[sequences.species[code]]['G11']+= 1

    def G10():
        total.angles.append(-105.9)
        total.translation.append(2.89)
        total.step_distrib['G10']+=1

        meta_data.all_steps['G10']+=1
        meta_data.species_steps[sequences.species[code]]['G10']+= 1

    def G00():
        total.angles.append(-108.1)
        total.translation.append(2.90)
        total.step_distrib['G00']+=1

        meta_data.all_steps['G00']+=1
        meta_data.species_steps[sequences.species[code]]['G00']+= 1

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
            if temp in step_options:
                step_options[temp]()
                total.P_number += value
            else:
                continue
        else:
            break

def CalcAverages():
    sum_angles = 0
    sum_translations = 0

    for angle in total.angles:
        sum_angles += angle

    for distance in total.translation:
        sum_translations += distance

    else:
        total.average_angle=sum_angles/len(total.angles)
        total.average_translation = sum_translations/len(total.translation)

    total.step_total = total.step_distrib['G22']+total.step_distrib['G21']+total.step_distrib['G20']+\
                       total.step_distrib['G11']+total.step_distrib['G10']+total.step_distrib['G00']

def CalcMetaData():

    sum_angles = 0
    sum_translation = 0

    for angle, translation in zip(meta_data.all_angles_averages, meta_data.all_translation_averages):
        sum_angles += angle
        sum_translation += translation

    meta_data.all_average_angles = sum_angles / len(meta_data.all_angles_averages)
    meta_data.all_average_translation = sum_translation / len(meta_data.all_translation_averages)

    for key in total.aa_counts:
        meta_data.all_aa_count[key] += total.aa_counts[key]

    for step in meta_data.all_steps:
        meta_data.all_steps_total += meta_data.all_steps[step]

    for step in meta_data.all_steps:
        meta_data.all_steps_percentage[step] += (100 * meta_data.all_steps[step]) / meta_data.all_steps_total

    for aa in meta_data.all_aa_count:
        meta_data.all_aa_total += meta_data.all_aa_count[aa]

    for aa in meta_data.all_aa_count:
        meta_data.all_aa_count_percentage[aa] = (100 * meta_data.all_aa_count[aa]) / meta_data.all_aa_total

def CalcMetaDataSP():

    for code in sequences.species:

        specie = sequences.species[code]

        sum_angles_specie = 0
        sum_translation_specie = 0

        for translation_temp, angle_temp in zip(meta_data.species_average_translation[specie],
                                                meta_data.species_average_angle[specie]):
            sum_angles_specie += angle_temp
            sum_translation_specie += translation_temp

        meta_data.species_average_angle[specie] = [sum_angles_specie / len(meta_data.species_average_angle[specie])]
        meta_data.species_average_translation[specie] = [sum_translation_specie / len(meta_data.species_average_translation[specie])]

        #duplicate

        meta_data.species_steps_total[specie] = meta_data.species_steps[specie]['G22'] + meta_data.species_steps[specie]['G21'] + \
                                                meta_data.species_steps[specie]['G20'] + meta_data.species_steps[specie]['G11'] + \
                                                meta_data.species_steps[specie]['G10'] + meta_data.species_steps[specie]['G00']

        for step in meta_data.species_steps[specie]:
            meta_data.species_steps_percentage[specie][step] = (100 * meta_data.species_steps[specie][step])/meta_data.species_steps_total[specie]

def AddToMeta():
    #overall

    for aa in total.aa_counts:
        meta_data.all_aa_count[aa] += total.aa_counts[aa]

    meta_data.all_angles_averages.append(total.average_angle)
    meta_data.all_translation_averages.append(total.average_translation)

def AddToMetaSP(code):
    #species data

    specie = sequences.species[code]

    for angle, translation in  zip(total.angles , total.translation):

        meta_data.species_angle[specie].append(angle)
        meta_data.species_translation[specie].append(translation)

    meta_data.species_average_angle[specie].append(total.average_angle)
    meta_data.species_average_translation[specie].append(total.average_translation)

    if specie not in meta_data.species_aa_count:
        meta_data.species_aa_count[specie] = total.aa_counts
        meta_data.species_aa_total[specie] = total.all_count
    else:
        for aa in total.aa_counts:
            meta_data.species_aa_count[specie][aa] += total.aa_counts[aa]
            meta_data.species_aa_count[specie][aa] /= 2
        meta_data.species_aa_total[specie] += total.all_count

    #specie name is already included in the dict for some reason? must initialise before
    #duplicate

    if all(step_value == 0 for step_value in meta_data.species_steps_percentage[specie].values()):
        for step in total.step_distrib_percentage:
            meta_data.species_steps_percentage[specie][step] += total.step_distrib_percentage[step]
    else:
        for step in total.step_distrib_percentage:
            meta_data.species_steps_percentage[specie][step] += total.step_distrib_percentage[step]
            meta_data.species_steps_percentage[specie][step] /= 2

    if specie not in meta_data.species_aa_average:
        meta_data.species_aa_average[specie] = total.all_count
    else:
        meta_data.species_aa_average[specie] = (total.all_count + meta_data.species_aa_average[specie])/2

def CalcComposition(code):

    singleletter_dict = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
                         'Trp': 'W', 'Thr': 'T', 'Asn': 'N', 'Pro': 'P', 'Phe': 'F',
                         'Ala': 'A', 'Gly': 'G', 'Ile': 'I', 'Leu': 'L', 'His': 'H',
                         'Arg': 'R', 'Met': 'M', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y',
                         'Glx' : 'Z', 'Unk': 'X'}

    total.all_count = len(strand1.aa)+len(strand2.aa)+len(strand3.aa)

    for aa in singleletter_dict.values():
        total.aa_counts[aa] = 0
        if aa not in meta_data.all_aa_count:
            meta_data.all_aa_count[aa] = 0
            meta_data.all_aa_count_percentage[aa] = 0

    for v1, v2, v3 in zip(strand1.aa, strand2.aa,strand3.aa):
        total.aa_counts[v1] += 1
        total.aa_counts[v2] += 1
        total.aa_counts[v3] += 1

    else:
        for aa in total.aa_counts:
            total.aa_counts[aa]= total.aa_counts[aa]*100/total.all_count

    # do step distribution

    for step in total.step_distrib:

        #turn it to percent
        total.step_distrib_percentage[step] = total.step_distrib[step]*100/total.step_total

def DrawGraph(code):
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
    file_open.SetFile("Graphs",code +".png")
    plt.savefig(file_open.file_path)

    plt.close()

def PlotMetaData():
    file_open = FileInfo()

    plt.figure(1)

    #step percentage

    bary3 = []
    barx3 =[]

    for steps in meta_data.all_steps_percentage:
        step_value = meta_data.all_steps_percentage[steps]
        bary3.append(float("{0:.2f}".format(round(step_value, 2))))
        barx3.append(steps)

    plt.subplot(223)
    plt.bar(range(len(bary3)),bary3, align='center')
    plt.xticks(range(len(bary3)), barx3)
    plt.ylim(0, 100)
    plt.ylabel('% Frequency')
    plt.title('Step Distribution')

    #aa comp

    barx4 = []
    bary4 = []

    for aa in meta_data.all_aa_count_percentage:
        barx4.append(aa)
        bary4.append(float("{0:.2f}".format(round(meta_data.all_aa_count_percentage[aa], 2))))

    plt.subplot(224)
    plt.bar(range(len(bary4)), bary4, align='center')
    plt.xticks(range(len(bary4)), barx4)
    plt.ylabel('% Composition')
    plt.xlabel('Amino Acid')
    plt.title('Amino Acid composition')

    #angle distrib

    as_array = []

    for angle in meta_data.all_angles_averages:
        as_array.append(float("{0:.2f}".format(round(angle, 2))))

    bary9 = as_array

    plt.subplot(222)
    plt.boxplot(bary9, showmeans=True, widths =0.2)
    plt.ylim(-110,-100)
    plt.tick_params(axis='x', which='both', bottom='off', labelbottom='off')
    plt.xlabel("n= " +str(len(bary9)))
    plt.ylabel('Angle (K)')
    plt.title('Average angle distributions')

    # translation distrib

    formatted = []

    for angle in meta_data.all_translation_averages:
        formatted.append(float("{0:.2f}".format(round(angle,2))))

    bary11 = formatted

    plt.subplot(221)
    plt.boxplot(bary11, showmeans=True, widths=0.2)
    plt.tick_params(axis='x', which='both', bottom='off', labelbottom='off')
    plt.ylim(2.70, 3)
    plt.xlabel("n= " + str(len(bary11)))
    plt.ylabel('Translation (t)')
    plt.title('Average translation distribution')

    plt.tight_layout()

    file_open.SetFile("Graphs", "Metadata.png")
    plt.savefig(file_open.file_path)

    plt.close()

def PlotMetaDataWithSpecies():
    file_open = FileInfo()

    plt.figure(1)

    style.use('ggplot')

    #boxplot with the species step percentage

    barx3 = []
    bary3 = []

    for specie in meta_data.species_steps_percentage:

        as_array = []

        for steps in meta_data.species_steps_percentage[specie]:

            step_value = meta_data.species_steps_percentage[specie][steps]
            as_array.append(float("{0:.2f}".format(round(step_value,2))))

        bary3.append(as_array)

        shortened_specie_name = "".join(item[0].upper() for item in re.findall("\w+", specie))

        barx3.append(shortened_specie_name)

    barx3_lab = list(range(len(barx3)))

    for value in barx3_lab:
        barx3_lab[value] = value + 1

    plt.subplot(224)
    plt.boxplot(bary3, labels=barx3, showmeans=True, widths =0.2)
    plt.xticks(barx3_lab,barx3, rotation = 'vertical')
    plt.ylim(0,60)
    plt.ylabel('% Frequency')
    plt.xlabel('Specie')
    plt.title('Step Distribution')

    #boxplot for the species angle values

    barx9 = []
    bary9 = []

    for specie in meta_data.species_angle:

        as_array = []

        for angle in meta_data.species_angle[specie]:
            as_array.append(float("{0:.2f}".format(round(angle,2))))

        bary9.append(as_array)

        shortened_specie_name = "".join(item[0].upper() for item in re.findall("\w+", specie))

        barx9.append(shortened_specie_name)

    barx9_lab = list(range(len(barx9)))

    for value in barx9_lab:
        barx9_lab[value] = value + 1

    plt.subplot(222)
    plt.boxplot(bary9, labels=barx9, showmeans=True, widths =0.2)
    plt.xticks(barx9_lab,barx9, rotation = 'vertical')
    plt.ylabel('Angle (K)')
    plt.xlabel('Specie')
    plt.title('Angle distributions')

    #bar chart for the length of the collagen domain in different specie

    barx10 = []
    bary10 = []

    for specie in meta_data.species_aa_average:

        length_value = meta_data.species_aa_average[specie]
        bary10.append(float("{0:.2f}".format(round(length_value,2))))

        shortened_specie_name = "".join(item[0].upper() for item in re.findall("\w+", specie))
        barx10.append(shortened_specie_name)

    plt.subplot(223)
    plt.bar(range(len(bary10)),bary10, align='center')
    plt.xticks(range(len(bary10)), barx10, rotation = 'vertical')
    plt.ylabel('Domain length (aa)')
    plt.xlabel('Specie')
    #actually number of aa total, not average length
    plt.title('Average domain length')

    # boxplot for the translation observed in species

    barx11= []
    bary11 = []

    for specie in meta_data.species_translation:

        as_array = []

        for angle in meta_data.species_translation[specie]:
            as_array.append(float("{0:.2f}".format(round(angle,2))))

        bary11.append(as_array)

        shortened_specie_name = "".join(item[0].upper() for item in re.findall("\w+", specie))

        barx11.append(shortened_specie_name)

    barx11_lab = list(range(len(barx9)))

    for value in barx11_lab:
        barx11_lab[value] = value + 1

    plt.subplot(221)
    plt.boxplot(bary11, labels=barx11, showmeans=True, widths =0.2)
    plt.xticks(barx11_lab,barx11, rotation = 'vertical')
    plt.ylabel('Translation (t)')
    plt.xlabel('Specie')
    plt.title('Translation distribution')

    plt.tight_layout()

    file_open.SetFile("Graphs", "Species Metadata.png")
    plt.savefig(file_open.file_path)

    plt.close()

def PlotSpecieData(code):
    specie = sequences.species[code]

    file_open = FileInfo()

    plt.figure(1)

    style.use('ggplot')

    #translation averages distribution

    formatted = []

    for t in meta_data.species_average_translation[specie]:
        formatted.append(float("{0:.2f}".format(round(t, 2))))

    bary5 = formatted

    plt.subplot(221)
    plt.boxplot(bary5, showmeans=True, widths =0.2)
    plt.tick_params(axis='x', which='both', bottom='off', labelbottom='off')
    plt.ylim(2.7, 3)
    plt.xlabel("n= " + str(len(bary5)))
    plt.ylabel('t (Å)')
    plt.title('Average translation distribution')

    # angles averages distribution

    as_array = []

    for angle in meta_data.species_average_angle[specie]:
        as_array.append(float("{0:.2f}".format(round(angle, 2))))

    bary12 = as_array

    plt.subplot(222)
    plt.boxplot(bary12, showmeans=True, widths=0.2)
    plt.tick_params(axis='x', which='both', bottom='off', labelbottom='off')
    plt.ylim(-100,-110)
    plt.xlabel("n= " + str(len(bary12)))
    plt.ylabel('K°')
    plt.title('Average angle distribution')

    #aa count bar chart

    barx6 = []
    bary6 = []

    for aa in meta_data.species_aa_count[specie]:
        barx6.append(aa)
        bary6.append(float("{0:.2f}".format(round(meta_data.species_aa_count[specie][aa],2))))

    plt.subplot(224)
    plt.bar(range(len(bary6)),bary6 , align='center')
    plt.xticks(range(len(bary6)),barx6)
    plt.ylabel('% Composition')
    plt.xlabel('Amino Acid')
    plt.title('Amino Acid composition')

    # step percentage for the specie

    barx7 = []
    bary7 = []

    for step in meta_data.species_steps_percentage[specie]:
        barx7.append(step)
        bary7.append(float("{0:.2f}".format(round(meta_data.species_steps_percentage[specie][step],2))))

    plt.subplot(223)
    plt.bar(range(len(bary7)),bary7,align='center')
    plt.xticks(range(len(bary7)),barx7)
    plt.ylabel('% Frequency')
    plt.xlabel('Step')
    plt.title('Step Distribution')

    plt.tight_layout()

    file_open.SetFile("Graphs", specie + " Metadata.png")
    plt.savefig(file_open.file_path)

    plt.close()

def SaveResults(code):

    file_open = FileInfo()
    file_open.SetFile("Data", code +(".txt"))

    with open(file_open.file_path, 'w', encoding='utf-8-sig') as file:
        file.truncate()

        file.write(str(sequences.species[code])+"\n")
        file.write("\nStep\tFrequency (%)" + "\n")
        for key in total.step_distrib_percentage:
            file.write(key + " \t"+ str(float("{0:.2f}".format(round(total.step_distrib_percentage[key],2))))+ "\n")

        file.write("\nTotal number of steps" + "\t" + str(total.step_total) + "\n"+ "\n"
        "Average Angle " + "\t" + str(float("{0:.2f}".format(round(total.average_angle,2)))) + "\n"
        "Average Translation " + "\t" +  str(float("{0:.2f}".format(round(total.average_translation,2)))) + "\n"+ "\n"
        "Strand 1 length " + "\t" + str(len(strand1.aa)) + "\n"
        "Strand 2 length " + "\t" + str(len(strand2.aa)) + "\n"
        "Strand 3 length " + "\t" + str(len(strand3.aa)) + "\n"+ "\n"
        "Number of P in Strand " + "\t" + str(total.P_number) + "\n"+ "\n"
        "aa" + "\t" + "%composition" +"\n")

        for aa in total.aa_counts:
            file.write(aa + "\t" + "{0:.2f}".format(round(total.aa_counts[aa],2))+"\n")

def SaveMetaData():

    #save for overall

    file_open = FileInfo()
    file_open.SetFile("Metadata","Overall.txt")

    with open(file_open.file_path, 'w', encoding='utf-8-sig') as file:

        file.truncate()

        file.write("\nn= " + str(len(meta_data.all_angles_averages)) + "\n")
        file.write("\n" + "Average angle: " + "\t" + str(meta_data.all_average_angles) +
                   "\n" + "Average Translation: " + "\t" + str(meta_data.all_average_translation) + "\n")

        file.write("\nStep\tFrequency" + "\n")
        for steps in meta_data.all_steps:
            file.write(steps + " \t" + str(meta_data.all_steps[steps]) +"\n")

        file.write("Total steps" + " \t" + str(meta_data.all_steps_total) + "\n\n")

        file.write("\nStep\tFrequency (%)" + "\n")
        for steps in meta_data.all_steps_percentage:
            step_value = meta_data.all_steps_percentage[steps]
            file.write(steps + " \t" + str(float("{0:.2f}".format(round(step_value, 2))))+"\n")

        file.write("\naa\t%composition\n")
        for aa in meta_data.all_aa_count_percentage:
            file.write(aa + "\t" + "{0:.2f}".format(round(meta_data.all_aa_count_percentage[aa], 2)) + "\n")

    file_open.reset()

    #save for species

    for specie in meta_data.species_steps_percentage:

        file_specie = FileInfo()
        file_specie.SetFile("Metadata",specie +(".txt"))

        with open(file_specie.file_path, 'w', encoding='utf-8-sig') as file:

            file.truncate()
            file.write("\n" + specie)

            file.write("\n\n n= " + str(len(meta_data.species_average_translation[specie]))+"\n")

            file.write("\n\n" + "Average angle: " + "\t" + str(meta_data.species_average_angle[specie][0]) +
                       "\n" + "Average Translation: " + "\t" + str(meta_data.species_average_translation[specie][0]) + "\n")

            file.write("\nStep\tFrequency" + "\n")
            for steps in meta_data.species_steps[specie]:
                file.write(steps + " \t" + str(meta_data.species_steps[specie][steps]) +"\n")

            file.write("Total steps" + " \t" + str(meta_data.species_steps_total[specie]) + "\n\n")

            file.write("\nStep\tFrequency (%)" + "\n")
            for steps in meta_data.species_steps_percentage[specie]:
                step_value = meta_data.species_steps_percentage[specie][steps]
                file.write(steps + " \t" + str(float("{0:.2f}".format(round(step_value, 2)))) +"\n")

            file.write("\naa" + "\t" + "%composition" + "\n")

            for aa in meta_data.species_aa_count[specie]:
                file.write(aa + "\t" + "{0:.2f}".format(round(meta_data.species_aa_count[specie][aa], 2)) + "\n")

def ClearMemory():

    strand1.reset()
    strand2.reset()
    strand3.reset()
    total.reset()

def main(UP=False, PDB=False, UI="", GXY_num=None, SP=False,OV=False, inc=False, ND=None, DOM=False):

    if ND is None:
        ND = []

    global strand1
    global strand2
    global strand3
    global total
    global meta_data
    global sequences

    if UP is True:

        if DOM is True:
            PullUP(domains=True)

        #does not check if ND is also empty (flaw)

        else:
            PullUP(input=ND)

        for code in sequences.sequence:

            commandOUT.SP.append(sequences.species[code])
            commandOUT.P.append(code)

            ProcessFASTA(sequences.sequence[code])

            if GXY_num:
                FindDomain(GXY_num)

            if inc is True:
                IncrementStrands()

            CountP()
            StepIteration(code)
            CalcAverages()
            CalcComposition(code)

            AddToMeta()
            AddToMetaSP(code)

            DrawGraph(code)

            if SP is True:
                PlotSpecieData(code)

            SaveResults(code)

            ClearMemory()

        CalcMetaData()
        CalcMetaDataSP()

        if OV is True:
            PlotMetaData()
        if SP is True:
            PlotMetaDataWithSpecies()

        if OV or SP:
            SaveMetaData()

        commandOUT.OAA = meta_data.all_average_angles
        commandOUT.OAT = meta_data.all_average_translation

    if PDB is True:

        if DOM is True:
            PullPDB(DOM=True)
        else:
            PullPDB(codes=ND)

        for code in sequences.sequence:

            commandOUT.SP.append(sequences.species[code])
            commandOUT.P.append(code)

            ProcessFASTAPDB(sequences.sequence[code])

            if GXY_num:
                FindDomain(GXY_num)

            if inc is True:
                IncrementStrands()

            CountP()
            StepIteration(code)
            CalcAverages()
            CalcComposition(code)
            AddToMeta()
            AddToMetaSP(code)
            DrawGraph(code)

            #Does not always work due to eclectic PDB classification

            if SP is True:
                PlotSpecieData(code)

            SaveResults(code)

            ClearMemory()

        CalcMetaData()
        CalcMetaDataSP()

        if SP is True:
            PlotMetaDataWithSpecies()
        if OV is True:
            PlotMetaData()

        if OV or SP:
            SaveMetaData()

        commandOUT.OAA = meta_data.all_average_angles
        commandOUT.OAT = meta_data.all_average_translation

    if UI:

        temp = UI.split("?")
        list(filter(None, temp))

        code = temp[0]
        specie = temp[1]
        FASTA = temp[2]

        sequences.sequence[code]=FASTA
        sequences.species[code] = specie

        for code in sequences.sequence:

            commandOUT.SP.append(sequences.species[code])
            commandOUT.P.append(code)

            ProcessFASTAUI(sequences.sequence[code])

            if GXY_num:
                FindDomain(GXY_num)

            if inc is True:
                IncrementStrands()

            CountP()
            StepIteration(code)
            CalcAverages()
            CalcComposition(code)
            AddToMeta()
            AddToMetaSP(code)
            DrawGraph(code)

            # Does not always work due to eclectic PDB classification

            if SP is True:
                PlotSpecieData(code)

            SaveResults(code)

            ClearMemory()

        CalcMetaData()
        CalcMetaDataSP()

        if SP is True:
            PlotMetaDataWithSpecies()
        if OV is True:
            PlotMetaData()

        if OV or SP:
            SaveMetaData()

        commandOUT.OAA = meta_data.all_average_angles
        commandOUT.OAT = meta_data.all_average_translation

def RunOptions(options):

    if "UI" in options:
        UI= options["UI"]
    else:
        UI = ""

    if "UP" in options:
        UP = True
    else:
        UP = False

    if "PDB" in options:
        PDB = True
    else:
        PDB = False

    if "GXY_num" in options:
        GXY_num = options["GXY_num"]
    else:
        GXY_num = None

    if "SP" in options:
        SP = True
        commandOUT.SPMD = True
    else:
        SP = False

    if "OV" in options:
        OV = True
        commandOUT.MD = True
    else:
        OV = False

    if "inc" in options:
        inc = True
    else:
        inc = False

    if "ND" in options:
        ND = options["ND"]
    else:
        ND = None

    if "DOM" in options:
        DOM = True
    else:
        DOM = False

    main(UP=UP, PDB=PDB, UI=UI, GXY_num=GXY_num, SP=SP,
         OV=OV, inc=inc, ND=ND, DOM=DOM)

commandIN =sys.stdin.read()
commandIN = json.loads(commandIN)

commandOUT = CommandOUT()

RunOptions(commandIN)

output = commandOUT

print('Content-Type: application/json')
print()
print(type(input))
print(type(output))
print(output.to_JSON())

#Some runs which worked

#main(UP=True, DOM= True, GXY_num=5, inc=True, SP=True, OV=True)

#main(UP=True, DOM=True, GXY_num=5, inc=True, SP=True, OV=True)

#main(UP=True, ND=["P02745","P20909"], GXY_num=5, inc=True, SP=True, OV=True)

#main(UP=True, DOM=True, GXY_num=5, inc=True, SP=True, OV=True)

#main(PDB=True, DOM = True, inc=True, OV=True)

#main(PDB=True, DOM = True, inc=True, OV=True)

#main(PDB=True, inc=True, OV=True, ND=["3pod","4dmt"])

#main(UI="CODE?SPECIE?GPPGPPG|PGPPGPP|PPGPPGP", inc=True, SP=True, OV=True)
