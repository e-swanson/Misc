#!/usr/bin/python3

import sys, os

length = {}
taxon = {}
read = {}

transcrips = []

g_length = "0.VIRGO2.geneLength.txt"
tax = "1.VIRGO2.taxon.txt"
reads = "Mapping_stats.csv"

dir = os.listdir(sys.argv[1])
os.mkdir(sys.argv[2])
file_count = 0
for i in dir:
    file_count += 1

with open(g_length, "r") as file:
    for line in file:
        l = line.split("\t")
        length[l[0]] = l[1].strip("\n")
        transcrips.append(l[0])
print("Made gene length dictionary")

with open(tax, "r") as file:
    for line in file:
        l = line.split("\t")
        taxon[l[1]] = [l[2].strip("\n")]
print("Made taxon dictionary")

with open(reads, "r") as file:
    for line in file:
        l = line.split(",")
        pid = l[0].split("_")
        pid = pid[0]
        read[pid] = {"Frag": l[5].replace("bp", ""), "aligned": l[7]}
print("Made mapping stats dictionary")

progress = 0
for file in dir:
    temp = {}
    pid = file.split("_")
    pid = pid[0]
    out = open(str(sys.argv[2]) + "/" + str(pid) + ".tsv", "w")
    out.write("target_id\tlength\teff_length\test_counts\ttpm\n")
    with open(str(sys.argv[1]) + "/" + str(file), "r") as f:
        for line in f:
            l = line.split("\t")
            temp[l[0]] = l[1].strip("\n")
    for item in transcrips:
        try:
            count = temp[item]
            l = length[item]
            f = read[pid]["Frag"]
            ef_l = float(l) - (float(f))
            if ef_l < 0:
                ef_l = 1
            tpm = float(count)/float(read[pid]["aligned"])
            out.write(str(item) + "\t" + str(l) + "\t" + str(ef_l) + "\t" + str(count) + "\t" + str(tpm) + "\n")
        except:
            out.write(str(item) + "\tNA\tNA\tNA\tNA\n")
    progress += 1
    print(str(progress) + " out of " + str(file_count) + " files finished...")

print("Done")
