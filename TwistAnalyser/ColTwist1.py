import sys
import requests
import urllib.request as urllib2
import re

#use beautifulsoup to make url requests
#djago cherrypy and flask are popular web framework
#matplotlib to create graphs
#numpy for matrices
#scipy is better numpy


def open_file():
    with open('Collagen-2_mouse.txt', 'r', encoding='utf-8-sig') as file:
        count = 1
        def strand1line():
            global strand1
            strand1 += line

        def strand2line():
            global strand2
            strand2 += line

        def strand3line():
            global strand3
            strand3 += line

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

def init_strands():

    global strand1
    global strand2
    global strand3
    global P_strand1
    global P_strand2
    global P_strand3
    global P_total

    strand1 = "".join(re.findall(5 * "G.{2}", strand1))
    strand2 = "".join(re.findall(5 * "G.{2}", strand2))
    strand3 = "".join(re.findall(5 * "G.{2}", strand3))

    strand1 = strand1[2:]
    strand2 = strand2[1:-1]
    strand3 = strand3[:-2]

    P_strand1 = []
    P_strand2 = []
    P_strand3 = []

    P_total = []


def count_P():
    for c in strand1:
        if c == "P":
            P_strand1.append(1)
        else:
            P_strand1.append(0)

    for c in strand2:
        if c=="P":
            P_strand2.append(1)
        else:
            P_strand2.append(0)

    for c in strand3:
        if c=="P":
            P_strand3.append(1)
        else:
            P_strand3.append(0)

    for v1, v2, v3 in zip(P_strand1, P_strand2, P_strand3):

        total = v1 + v2 + v3
        P_total.append(total)

def step_iteration():
    global angles
    global translation
    global G22_number
    global G21_number
    global G20_number
    global G11_number
    global G10_number
    global G00_number
    global P_number

    P_number = 0

    angles = []
    translation = []

    G22_number = 0
    G21_number = 0
    G20_number = 0
    G11_number = 0
    G10_number = 0
    G00_number = 0


    def G22():
        global angles
        angles.append(-102.6)
        translation.append(2.84)
        global G22_number
        G22_number+=1

    def G21():
        global angles
        angles.append(-104.5)
        translation.append(2.86)
        global G21_number
        G21_number+=1

    def G20():
        #no values for type 2???

        global angles
        angles.append(-105.0)
        translation.append(2.86)
        global G20_number
        G20_number+=1

    def G11():
        global angles
        angles.append(-105.0)
        translation.append(2.86)
        global G11_number
        G11_number+=1

    def G10():
        global angles
        angles.append(-105.9)
        translation.append(2.89)
        global G10_number
        G10_number+=1

    def G00():
        global angles
        angles.append(-108.1)
        translation.append(2.90)
        global G00_number
        G00_number+=1

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

    for index, value in enumerate(P_total):
        if index<len(P_total)-1:
            temp = str(P_total[index]) + str(P_total[index+1])
            step_options[temp]()
            P_number += value
        else:
            break

def calc_averages():
    total_angle = 0
    total_translation = 0

    for angle in angles:
        total_angle += angle

    for distance in translation:
        total_translation += distance

    global average_translation
    global average_angle
    global step_total

    average_angle=total_angle/len(angles)
    average_translation = total_translation/len(translation)

    step_total = G22_number+G21_number+G20_number+G11_number+G10_number+G00_number

def save_results():
    with open('ColTwist1.txt', 'w', encoding='utf-8-sig') as file:
        file.truncate()
        file.write("Step\tFrequency" + "\n"
        "G22 " + "\t"+ str(G22_number)+ "\n"
        "G21 " + "\t"+ str(G21_number)+ "\n"
        "G20 " + "\t"+ str(G20_number)+ "\n"
        "G11 " + "\t"+ str(G11_number)+ "\n"
        "G10 " + "\t"+ str(G10_number)+ "\n"
        "G00 " + "\t"+ str(G00_number)+ "\n"
        "Total " + "\t" + str(step_total) + "\n"+ "\n"
        "Average Angle " + "\t" + str(average_angle) + "\n"
        "Average Translation " + "\t" + str(average_translation) + "\n"+ "\n"
        "Strand 1 length " + "\t" + str(len(strand1)) + "\n"
        "Strand 2 length " + "\t" + str(len(strand2)) + "\n"
        "Strand 3 length " + "\t" + str(len(strand3)) + "\n"+ "\n"
        "Number of P in Strand " + "\t" + str(P_number) + "\n")

strand1 = ""
strand2 = ""
strand3 = ""

def main():

    open_file()
    init_strands()

    print(strand1)
    print(strand2)
    print(strand3)

    count_P()

    print(P_strand1)
    print(P_strand2)
    print(P_strand3)
    print(P_total)

    step_iteration()

    print (angles)
    print (translation)

    calc_averages()
    save_results()

main()