__author__ = 'euphrasieservant'
# population de départ : N séquences d'ARNt de longueur l =
# chaques séquences est constituées de l nucléotides choisis au hasard
import random

# constitution de la population de départ
tirage = 0
nuc = ""
atgc = ["A","T","U", "C"]
sequences = []

l = 75  # longueur de la séquence
for i in range(10000):
    seq = []
    for y in range(l):
        tirage = random.randint(4)
        nuc = atgc[tirage]
        seq.append(nuc)
    sequences.append(seq)

for i in range(10):
    print(sequences[i])

# phénotype

# appariements de la séquence cible
