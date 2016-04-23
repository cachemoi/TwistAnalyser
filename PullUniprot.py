import urllib.error
import urllib.request
import urllib.parse
import re

url = 'http://www.uniprot.org/uniprot/'

class Codes :
    def __init__(self):
        self.codes =[]

class Sequences :
    def __init__(self):
        self.sequence = {}
        self.species = {}


sequences = Sequences()
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
                print("503 error")
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
                    print("503 error")
                    continue
            else:
                FASTA = request.read().decode('utf-8')
                PrepareSequence(FASTA,code)
                break

        print("Downloading FASTA")

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
                print("503 error")
                continue
        else:
            FASTA = request.read().decode('utf-8')
            print(FASTA)
            PrepareSequence(FASTA,code)
            break

    print("Downloading")

#GetDomain()
#print("codes acquired")
#GetFASTA()

GetByID("P02745")
GetByID("P20909")
print("Fasta acquired")