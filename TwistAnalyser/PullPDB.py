import json
import sys
import timeit

start = timeit.default_timer()

PY3 = sys.version > '3'

if PY3:
    import urllib.request as urllib2
else:
    import urllib2

class SequenceList:
    def __init__(self):
        self.url = {}
class Sequences :
    def __init__(self):
        self.sequence = {}

SERVER_URL = "https://www.ebi.ac.uk/pdbe/api"

SUMMARY = "pdb/entry/files/"

def make_request(url, data):
    request = urllib2.Request(url)

    try:
        url_file = urllib2.urlopen(request, data)
    except urllib2.HTTPError as e:
        if e.code == 404:
            print("[NOTFOUND %d] %s" % (e.code, url))
        else:
            print("[ERROR %d] %s" % (e.code, url))

        return None

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
sequences = Sequences()

def Main():
    input= ['1CAG','1CGD','1A3I','1A3J','1G9W','1ITT','1K6F','2CUO','3AH9','1BKV','1QSU',
              '1NAY','1Q7D','1V7H','4OY5','1V4F','1V6Q','3B0S','1X1K','1WZB','1YM8','1EI8',
              '2G66','2DRT','2DRX','3DMW','2D3F','2D3H','3A0A','2KLW','3A08','3A19','3IPN',
              '3ADM','3A1H','3A0M','3POD','3PON','3P46','3T4F','3U29','3B2C','4AXY','4DMT',
              '4GYX','3ABN','3WN8']

    response = post_request(SUMMARY, input, True)

    if response:
        data = json.loads(response)

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


Main()

stop = timeit.default_timer()

print ("\n" +  "Time for request: " + str(stop - start) + " seconds" + '\n')
