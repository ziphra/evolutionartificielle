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
def complementaire(sequences, pos1, pos2):
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
    for k in range(0, len(sequences)): # len sequences = n
        score = 0
        for i in range(0, (len(listeapp))):
                pos1 = listeapp[i][0]
                pos2 = listeapp[i][1]

                if (len(sequences[k])) <= pos1 and (len(sequences[k])) <= pos2: # conditions pour éviter les erreurs de types
                                                                                # liste index out of range
                    if complementaire(sequences[k], pos1, pos2) == True:
                        score += 1
        for cle in consensus:  # pour A, G, C, U qui sont les clés de mon dic
            value = consensus.get(cle)  # on recupere les valeurs : donc les positions ex : (33, 8) pour U
            for o in range(len(value)):  # pour chacune de ces valeurs
                n = value[o]  # on recupere dans n une des valeurs : ex : value[O] = 33
                if (len(sequences[k])-1) >= n and sequences[k][n] == cle :  # si le nucleotide de la sequence k, à la position n = cle
                                            # (ex : si nucleotide de la sequence k = U)
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
            NewSeq = mutation(sequences[i], substitution, deletion, insertion)
            enfants.append(NewSeq)
    while len(enfants) < N:        # tant que la [enfants] ne contient pas autant d'elements que la liste de départ,
                                    #  repeter
        b = random.randint(0, (len(enfants)-1))  #tire un nombre entre 0 et (len de enfants)-1
        ListeScore = score(listeapp, enfants, consensus)
        for j in range((ListeScore[b]//2)):
            enfants.append(enfants[b])
    return [enfants, ListeScore]  # on retourne le score et les enfants dans une liste pour pouvoir les retourner et les
                                  #  utiliser séparement.


# avec l'ajout des consensus, le score max devient 42 (21 appariements + 21 consensus)
# et de la recombinaison
def reproduction(sequences, N, moy):
    enfants = []

    ListeScore = score(listeapp, sequences, consensus)
    while len(enfants) < (N):

        a = random.randint(0, 42)  # 21 appariements + 21 consensus
        j = random.randint(0, (N-1))  # un nombre au hasard entre 0 N-1 soit pour chaque indice de la liste
                                    # ex : une liste de N = 75 sequences indicées de 0 à N-1 = 74)
        scorej = ListeScore[j]      # on tire la séquence associé au score grace à l'indice du score
        if (a <= scorej) or (scorej >= moy):             # scorj chance sur 42 que a soit <= scorej
                                                           #  on rajoute une sélection : scorej >= à la moyenne
                                                            # des scores précédents
            NewSeq = mutation(sequences[j], substitution, deletion, insertion)  # mutation ...

            enfants.append(NewSeq) # stockage des séquences dans une liste
    ListeScore = score(listeapp, enfants, consensus)

    TirageRecomb = random.randint(0,100)
    if TirageRecomb == 16:             # 1 chance sur 100
        c = random.randint(0, (len(enfants)-1))
        d = random.randint(0, (len(enfants)-1))
        while c == d:
            c = random.randint(0, (len(enfants)-1))
            d = random.randint(0, (len(enfants)-1))
        sequence1 = enfants[c]
        sequence2 = enfants[d]
        recombinés = recombinaison(sequence1, sequence2)
        enfants[c] = recombinés[0]
        enfants[d] = recombinés[1]

    return [enfants, ListeScore]  # on retourne une liste avec les enfants et les scores pour pouvoir retourner
                                    # ces deux elements en même temps et les utiliser séparément



# induit une mutation avec une probabilité de 1/100
def mutation(seq, substitution, deletion, insertion):  # seq = un seul ARNt
    tirage = random.randint(1, 300) # il ya une chance sur 100 d'avoir une mutation de 3 types différents
                                    # une chance sur 300 pour chaque
    augc = ["A","U","G", "C"]
    pos_mutation = random.randint(0, (len(seq)-1)) # une position choisie au hasard dans la séquence de longueur l

    if tirage == substitution:  # 3 chance sur 300 = 1/100
        nuc = augc[random.randint(0, 3)]  # un nucléotide de remplacement choisi au hasard
        while nuc == seq[pos_mutation]:  # au cas ou le nucléotide de remplacement est le même que l'original
            nuc = augc[random.randint(0, 3)]
        seq[pos_mutation] = nuc

    elif tirage == deletion:
        del seq[pos_mutation]

    elif tirage == insertion:
        nuc = augc[random.randint(0, 3)]  # un nucléotide choisi au hasard
        seq.insert(pos_mutation, nuc)

    return seq

# moyenne
def MoyenneScore(ListeScore):
    moy = 0
    for i in range(0, (len(ListeScore))):
        moy += ListeScore[i]
    if len(ListeScore) == 0:
        moy = 0
    else:
        moy = moy/(len(ListeScore))
    return moy


# recombinaison
def recombinaison(sequence1, sequence2):
    pos1 = random.randint(0, (len(sequence1)-1))  # on choisit 2 position au hasard dans 2 séquences qui seront
                                                    # les points d'inititaions de la recombinaison
    pos2 = random.randint(0, (len(sequence2)-1))

    a = random.randint(10, 35)  # on choisit la longueur des morceaux d'ADN qui vont se recombiner
    b = random.randint(10, 35)
    if pos1 < a:    # car si pos1 > a, alors a sera en dehors de la liste
        seq1 = sequence1[pos1:(pos1+a)]  # seq1 = morceau de séquence qui va être recombiné
    else:
        seq1 = sequence1[pos1:(pos1+a)]
    if pos2 < b:
        seq2 = sequence2[pos1:(pos1+b)]
    else:
        seq2 = sequence2[pos1:(pos1+b)]

    sequence1bis = sequence1[0:(pos1-1)] + seq2 + sequence1[a:(len(sequence1)-1)]
    sequence2bis = sequence2[0:(pos2-1)] + seq1 + sequence2[a:(len(sequence2)-1)]
    return [sequence1bis, sequence2bis]



#P R O G R A M M E  P R I N C I P A L

#population de depart
l = 75 # longueur de la séquence
N = 100 #nombre de chaines
enfants = PopDeDepart(l,N)

#sequences consensus

consensus = {"C": (11, 25, 32, 48, 60, 62, 74, 75), "U": (8, 11, 25, 32, 33, 48, 60, 62), "A": (14, 21, 58, 9, 10, 15, 24, 37, 52, 57, 75), "G":(19, 18, 53,9, 10, 15, 24, 37, 52, 57)}

#appariements
repli1 = [1, 7, 72]
repli2 = [10, 13, 25]
repli3 = [27, 31, 43]
repli4 = [49, 53, 65]
repli = [repli1] + [repli2] + [repli3] + [repli4]

listeapp = (ListeAppariements(repli)) #liste de liste contenant une paire de position : la position des appariements
                                    # de la sequence cible


# mutation
substitution = 1
deletion = 2
insertion = 3

#reproduction
nbr_generation = 100
generations_enfants = []
ListeMoy = []
moy = 0

for i in range(nbr_generation):
    print(moy)
    EnfantsScore = reproduction(enfants, N, moy)  # enfants ici = PopDeDépart
    ListeScore = EnfantsScore[1]        # en postion 1 dans la liste on a retourné Liste score
                                        # return [enfants, ListeScore]
    moy = MoyenneScore(ListeScore)       # on calcul la moyenne pour tracer le graphique
    ListeMoy.append(moy)                # on stocke les moyennes dans une liste
    enfants = EnfantsScore[0]           # en position 0, la liste des nv sequences





# G R A P H I Q U E
abs = []
for i in range(nbr_generation):
    abs.append(i+1)     # numero de la génération

plt.plot(abs, ListeMoy)          # (numero de la génération en x, moyenne en y)
plt.ylabel('Score Moyen')           # titre axe des y
plt.xlabel('Nombre De Generation')      # titre axe des x
plt.show()      # affichage
