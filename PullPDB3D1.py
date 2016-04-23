import json
import re
import sys
import timeit
from collections import defaultdict

import numpy as np
import io

start = timeit.default_timer()

PY3 = sys.version > '3'

if PY3:
    import urllib.request as urllib2
else:
    import urllib2

class SequenceList:
    def __init__(self):
        self.url = {}
class Proteins :
    def __init__(self):
        self.sequence = {}
        self.coord = {}
        self.species = {}
        self.atoms = defaultdict(list)

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
    for c in proteins.sequence:
        temp = proteins.sequence[c].splitlines()
        temp_seq = ""

        for index, value in enumerate(temp):
            if index%2 == 1:
                temp_seq += "|"
                temp_seq += value

        proteins.sequence[c] = temp_seq

sequence_list = SequenceList()
proteins = Proteins()

def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)

def get_coordinates_pdb(code, PDB, ignore_hydrogens):
    """
    Get coordinates from the first chain in a pdb file
    and return a vectorset with all the coordinates.
    """
    # PDB files tend to be a bit of a mess. The x, y and z coordinates
    # are supposed to be in column 31-38, 39-46 and 47-54, but this is not always the case.
    # Because of this the three first columns containing a decimal is used.
    # Since the format doesn't require a space between columns, we use the above
    # column indices as a fallback.
    x_column = None
    V = []
    # Same with atoms and atom naming. The most robust way to do this is probably
    # to assume that the atomtype is given in column 3.
    atoms = []

    buf = io.StringIO(PDB)

    for line in buf:
            if line.startswith("TER") or line.startswith("END"):
                break
            if line.startswith("ATOM"):
                tokens = line.split()

                # Try to get the atom type
                try:
                    atom = tokens[2][0]
                    if ignore_hydrogens and atom == "H":
                        continue
                    elif atom in ["H", "C", "N", "O", "S", "P"]:
                        atoms.append(atom)
                except:
                        exit("Error parsing atomtype for the following line: \n%s" % line)

                if x_column == None:
                    try:
                        # look for x column
                        for i, x in enumerate(tokens):
                            if "." in x and "." in tokens[i+1] and "." in tokens[i+2]:
                                x_column = i
                                break
                    except IndexError:
                        exit("Error parsing coordinates for the following line: \n%s" % line)
                # Try to read the coordinates
                try:
                    V.append(np.asarray(tokens[x_column:x_column+3],dtype=float))
                except:
                    # If that doesn't work, use hardcoded indices
                    try:
                        x = line[30:38]
                        y = line[38:46]
                        z = line[46:54]
                        V.append(np.asarray([x,y,z],dtype=float))
                    except:
                        exit("Error parsing input for the following line: \n%s" % line)


    V = np.asarray(V)
    """
    np.set_printoptions(precision=4,
                       threshold=10,
                       linewidth=150)
    """
    proteins.coord[code] = V
    proteins.atoms[code] = atoms

def Main():
    input= ['1CAG', '1CGD']

    response = post_request(SUMMARY, input, True)

    if response:
        data = json.loads(response)

        collected_codes = []

        for value in data:
            collected_codes.append(value)

        for code in collected_codes:

            for value in data[code]['PDB']['downloads']:
                if value['label'] == 'PDB file':

                    sequence_list.url[code] = value['url']

                    PDB_req = urllib2.urlopen(value['url'])

                    PDB = PDB_req.read().decode('utf-8')

                    get_coordinates_pdb(code, PDB, True)

                    proteins.sequence[code] = PDB

    PrepareSequence()


Main()

for code in proteins.coord:
    print(rmsd(proteins.coord['1cag'],proteins.coord['1cgd']))

for code in proteins.atoms:
    print(proteins.atoms[code])

stop = timeit.default_timer()

print ("\n" +  "Time for request: " + str(stop - start) + " seconds" + '\n')
