__authors__ = 'euphrasieservant', 'agathewallet'
# -*- coding:utf-8 *-*
# population de départ : N séquences d'ARNt de longueur l =
# chaque séquence est constitué de l nucléotides choisis au hasard
import random
import matplotlib.pyplot as plt


# constitution de la population de départ
def PopDeDepart(l,N):

    augc = ["A","U","G", "C"]
    sequences = []

    for i in range(N):
        seq = []
        for y in range(l):
            tirage = random.randint(0, 3)
            nuc = augc[tirage]
            seq.append(nuc)
        sequences.append(seq)
    return sequences


# PHÉNOTYPE
# fonction qui génère la liste des appariements
# sous forme d'une liste de paires de positions
def ListeAppariements(repli):
    listeapp = []
    for i in range(0, len(repli)):
        a = repli[i][0]  #premier élément du repli dans le sens 5'-3'   (1 dans le cas du repli1)
        b = repli[i][2]  #premier élement du repli dans le sens 3'-5'   (72 dans le cas du repli1)
        c = repli[i][1]  #dernier element du repli dans le sens 5'-3'   (7 dans le cas du repli1)
        for j in range(a, c+1):
            pos = []
            pos.append(a)
            pos.append(b)
            listeapp.append(pos)
            a += 1
            b -= 1
    return listeapp


#permet de verifier si 2 nucléotides à deux positions précises de la séquences sont complémentaires ou non
def complementaire(sequences,pos1, pos2):
    if (sequences[pos1] == "A" and sequences[pos2] == "U") or (sequences[pos2] == "A" and sequences[pos1] == "U"):
        complement = True
    elif (sequences[pos1] == "G" and sequences[pos2] == "C") or (sequences[pos2] == "G" and sequences[pos1] == "C"):
        complement = True
    elif (sequences[pos1] == "G" and sequences[pos2] == "U") or (sequences[pos2] == "G" and sequences[pos1] == "U"):
        complement = True
    else:
        complement = False
    return complement


# score
# range les scores dans une liste
def score(listeapp, sequences, consensus):
    ListeScore = []
    for k in range(0,len(sequences)): # len sequences = n
        score = 0
        for i in range(0,(len(listeapp))):
            pos1 = listeapp[i][0]
            pos2 = listeapp[i][1]
            if complementaire(sequences[k], pos1, pos2) == True:
                score += 1
        for cle in consensus:  # pour A, T , G, C, U qui sont les clés de mon dic
            value = consensus.get(cle)  # on recupere les valeurs : donc les positions ex : (33, 8) pour U
            for o in range(len(value)):  # pour chacune de ces valeurs
                n = value[o]  # on recupere dans n une des valeurs : ex : value[O] = 33
                if sequences[k][n] == cle:  # si le nucleotide de la sequence k, à la position n = cle (ex : si nucleotide de la sequence k = U)
                    score += 1  # alors le score augmente

        ListeScore.append(score)
    return ListeScore

#reproduction
#le score est proportionel à la probabilité de se reproduire
#score de 5 = 5/21 chances de se reproduire (car 21 = nombre d'appariements de la cible)

def reproductionbis(sequences, N):  # prends plus de temps que l'autre fonction reproduction
    enfants = []
    ListeScore = score(listeapp, sequences, consensus)
    for i in range(0, (len(ListeScore))):
        scorei = ListeScore[i]
        a = random.randint(0, 42)  # 21 appariements + 21 consensus
        if a <= scorei:
            NewSeq = mutation(sequences[i], probamutation)
            enfants.append(NewSeq)
    while len(enfants) < N:        # tant que la [enfants] ne contient pas autant d'elements que la liste de départ, repeter
        b = random.randint(0, (len(enfants)-1))  #tire un nombre entre 0 et (len de enfants)-1
        ListeScore = score(listeapp, enfants, consensus)
        for j in range((ListeScore[b]//2)):
            enfants.append(enfants[b])
    return [enfants, ListeScore]  # on retourne le score et les enfants dans une liste pour pouvoir les retourner et les utiliser séparement.




def reproduction(sequences, N):
    enfants = []
    ListeScore = score(listeapp, sequences, consensus)
    while len(enfants) < (N+1):
        j = random.randint(0,(N-1))
        scorej = ListeScore[j]
        for i in range(scorej//2):
            NewSeq = mutation(sequences[j], probamutation)
            enfants.append(NewSeq)

    return [enfants, ListeScore]







# induit une mutation avec une probabilité de 1/100
def mutation(seq,probamutation):  # seq = un seul ARNt
    tirage = random.randint(1,100)
    augc = ["A","U","G", "C"]
    if tirage <= probamutation:  # 1 chance sur 100

        pos_mutation = random.randint(0, (len(seq)-1))  # une position choisie au hasard dans la séquence de longueur l
        nuc = augc[random.randint(0, 3)]  # un nucléotide de remplacement choisi au hasard
        while nuc == seq[pos_mutation]:  # au cas ou le nucléotide de remplacement est le même que l'original
            nuc = augc[random.randint(0, 3)]
        seq[pos_mutation] = nuc

    return seq

def MoyenneScore(ListeScore):
    moy = 0
    for i in range(0, (len(ListeScore))):
        moy += ListeScore[i]
    if len(ListeScore) == 0:
        moy = 0
    else:
        moy = moy/(len(ListeScore))
    return moy

#P R O G R A M M E  P R I N C I P A L

#population de depart
l = 75 # longueur de la séquence
N = 100 #nombre de chaines
enfants = PopDeDepart(l,N)

#sequences consensus

consensus = {"U":(33, 8), "C":(11, 25, 32, 48, 60, 62), "T":(11, 25, 32, 48, 60, 62), "A":(14, 21, 58, 9, 10, 15, 24, 37, 52, 57), "G":(19, 18, 53,9, 10, 15, 24, 37, 52, 57)}

#appariements
repli1 = [1, 7, 72]
repli2 = [10, 13, 25]
repli3 = [27, 31, 43]
repli4 = [49, 53, 65]
repli = [repli1] + [repli2] + [repli3] + [repli4]

listeapp = (ListeAppariements(repli)) #liste de liste contenant une paire de position : la position des appariements


# mutation
probamutation = 1

#reproduction

nbr_generation = 350
generations_enfants = []
ListeMoy = []


for i in range(nbr_generation):
    EnfantsScore = reproduction(enfants, N)
    ListeScore = EnfantsScore[1]
    enfants = EnfantsScore[0]

    moy = MoyenneScore(ListeScore)
    ListeMoy.append(moy)



# G R A P H I Q U E
abs = []
for i in range(nbr_generation):
    abs.append(i+1)

plt.plot(abs,ListeMoy)
plt.ylabel('Score Moyen')
plt.xlabel('Nombre De Generation')
plt.show()
