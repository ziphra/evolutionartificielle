__author__ = 'euphrasieservant'
# -*- coding:utf-8 *-*
# population de départ : N séquences d'ARNt de longueur l =
# chaque séquence est constitué de l nucléotides choisis au hasard
import random


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
    elif (sequences[pos1] == "G" and sequences[pos2] == "U") or (sequences[pos2] =="G" and sequences[pos1] == "U"):
        complement = True
    else:
        complement = False
    return complement


# score
# range les scores dans une liste
def score(listeapp, sequences):
    ListeScore = []
    for k in range(0,len(sequences)): #len sequences = n
        score = 0
        for i in range(0,(len(listeapp))):
            pos1 = listeapp[i][0]
            pos2 = listeapp[i][1]
            if complementaire(sequences[k], pos1, pos2) == True:
                score += 1

        ListeScore.append(score)
    return ListeScore

#reproduction
#le score est proportionel à la probabilité de se reproduire
#score de 5 = 5/21 chances de se reproduire (car 21 = nombre d'appariements de la cible)

def reproduction(ListeScore, enfants):
    c = 0  #un compteur qui nous indique combien de nouvelles séquences ont été crées
    for i in range (0,(len(ListeScore))):
        scorei = ListeScore[i]
        n = random.randint(0,20)
        if n <= scorei:
            NewSeq = sequences[i]
            mutation(NewSeq,probamutation)
            enfants.append(NewSeq)
            c+=1
    #print(c)

    return enfants




# induit une mutation avec une probabilité de 1/100
def mutation(seq,probamutation):  #seq = un seul ARNt
    tirage = random.randint(1,100)
    augc = ["A","U","G", "C"]
    if tirage == probamutation: #  1 chance sur 100
        pos_mutation = random.randint(0,(len(seq)-1))  # une position choisie au hasard dans la séquence de longueur l
        nuc = augc[random.randint(0,3)]  # un nucléotide de remplacement choisi au hasard
        while nuc == seq[pos_mutation]:  # au cas ou le nucléotide de remplacement est le même que l'original
            nuc = augc[random.randint(0,3)]
        seq[pos_mutation] = nuc
        print(seq)
    return seq



#P R O G R A M M E  P R I N C I P A L

#population de depart
l = 10 # longueur de la séquence
N = 100  #nombre de chaines
sequences = PopDeDepart(l,N)

#appariements
repli1 = [1, 7, 72]
repli2 = [10, 13, 25]
repli3 = [27, 31, 43]
repli4 = [49, 53, 65]
repli = [repli1] + [repli2] + [repli3] + [repli4]

listeapp = (ListeAppariements(repli))
print(listeapp)



# mutation
probamutation = 100

#reproduction
#Après avoir compté le nombre de nouvelles séquences créées après la première genération, et l'avoir comparé avec les scores,
# le resultat de cette fonction semble cohérent.
# => Les scores oscillent autour de 5. Si en moyenne, une séquence a un score de 5 (sur un maximum de 21)
# la probabilité de se reproduire est de 6/21
# (6/21) * 10 000 = 2 857 => à peu près égale aux nombres de nouvelles séquences créées

nbr_generation = 4
generations_enfants = []

for i in range(nbr_generation) :
    ListeScore = score(listeapp, sequences)
    enfants = reproduction(ListeScore, sequences)
    generations_enfants.append(enfants)
    sequences = enfants
    print (enfants)





#print(generations_enfants)


