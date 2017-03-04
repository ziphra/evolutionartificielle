__author__ = 'euphrasieservant'
# population de départ : N séquences d'ARNt de longueur l =
# chaques séquences est constituées de l nucléotides choisis au hasard
import random


# constitution de la population de départ
def PopDeDepart(l):

    atgc = ["A","U","G", "C"]
    sequences = []

    for i in range(10000):
        seq = []
        for y in range(l):
            tirage = random.randint(0, 3)
            nuc = atgc[tirage]
            seq.append(nuc)
        sequences.append(seq)
    return sequences


# phénotype
# appariements de la séquence cible, sous forme d'une liste de paires de positions
# fonctions qui génère la liste des appariements
def ListeAppariements(repli):
    listeapp = []
    for i in range(0, 4):
        pos = []
        a = y = repli[i][0]
        b = repli[i][2]
        z = repli[i][1]
        pos.append(a)
        pos.append(b)
        listeapp.append(pos)
        for j in range(y, z):
            pos = []
            a += 1
            b -= 1
            pos.append(a)
            pos.append(b)
            listeapp.append(pos)
    return listeapp


#permet de verifier si 2 nucléotides à deux positions précises de la séquences sont complémentaires ou non
def complementaire(sequences,pos1, pos2):
    if (sequences[pos1] == "A" and sequences[pos2] == "U") or (sequences[pos2] == "A" and sequences[pos1] == "U"):
        complement = True
    elif (sequences[pos1] == "G" and sequences[pos2] == "C") or (sequences[pos2] == "G" and sequences[pos1] == "C"):
        complement = True
    else:
        complement = False
    return complement


# score
def score(listeapp, seq):
    score = 0
    for i in range(len(listeapp)):
        pos1 = listeapp[i][0]
        pos2 = listeapp[i][1]
        if complementaire(seq, pos1, pos2) == True:
            score += 1
            print(pos1, pos2, seq[pos1], seq[pos2])
    return score




#P R O G R A M M E P R I N C I P A L E

#population de depart
l = 75  # longueur de la séquence

sequences = PopDeDepart(l)

#appariements
repli1 = [1, 7, 72]
repli2 = [10, 13, 25]
repli3 = [27, 31, 43]
repli4 = [49, 53, 65]
repli = [repli1] + [repli2] + [repli3] + [repli4]

listeapp = (ListeAppariements(repli))


#score
# test : donne le score pour les séquences 1 à 10
for i in range(10):
    seq = sequences[i]
    print(score(listeapp, seq))

