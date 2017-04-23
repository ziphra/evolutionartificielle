__author__ = 'euphrasieservant'
__authors__ = 'euphrasieservant', 'agathewallet'
# -*- coding:utf-8 *-*
import random
import copy
import matplotlib.pyplot as plt


# constitution de la population de départ
# param l : nombre de nucléotides dans une séquence
# param N : nombre de séquences
# return : une liste de N séquences de l nucléotides
# (une séquence est une liste de nucléotides
# une nucléotide est représentée par une chaine de 1 caratère: A, U, G ou C)
def PopDeDepart(l, N):

    augc = ["A", "U", "G", "C"]
    sequences = []
    # On crée N séquences en passant N fois dans cette boucle
    for i in range(N):
        # seq = une liste initialement vide qui contiendra UNE séquence
        seq = []
        # La séquence contiendra l nucléotide en passant l fois dans la boucle
        for y in range(l):
            nuc = augc[random.randint(0, 3)]
            seq.append(nuc)
        sequences.append(seq)
    return sequences


# PHÉNOTYPE
# Génération de la liste des appariements
# sous forme d'une liste de paires de positions
# param repli : liste contenant 3 éléments : le premier  et le dernier
# nucléotide dans le sens 5'-3' du repli, et le premier dans le sens 3'-5'
# return : une liste de paires de positions (chaque paire
# représente les positions des nucléotides qui s'apparient entre eux)
# pour la formation des boucles
def ListeAppariements(repli):
    listeapp = []
    # Pour tous les repli (repli1, repli2, repli3, repli4)
    for i in range(0, len(repli)):
        # Premier élément du repli dans le sens 5'-3'   (1 dans le cas du repli1)
        a = repli[i][0]
        # dernier élement du repli dans le sens 3'-5'   (72 dans le cas du repli1)
        b = repli[i][2]
        # 2e element du repli dans le sens 5'-3'   (7 dans le cas du repli1)
        c = repli[i][1]
        # Du premier élément (7), au dernier INCLU (d'ou c+1) dans le sens 5'-3'
        for j in range(a, c+1):
            pos = []
            # a = 1, puis 2, 3...7
            pos.append(a)
            # b = 72,puis 71, 70....
            pos.append(b)
            # [1,72], [2,71]
            listeapp.append(pos)
            # On incrémente a et b de 1 pour passer aux nucléotides suivants
            a += 1
            b -= 1
    return listeapp


# Vérifier si 2 nucléotides à deux positions précises de la séquences
# sont complémentaires ou non
# param sequence : une liste de nucléotides
# param pos1 : position (entier) du 1er nucléotide
# param pos2 : position (entier) du 2eme nucléotide
# return true si les 2 nucléotides sont complémentaires
def complementaire(sequence, pos1, pos2):
    if (len(sequence)-1 < pos1):
        # pos1 en dehors de la séquence
        return False

    if (len(sequence)-1 < pos2):
        # pos2 en dehors de la séquence
        return False

    complement = False
    if (sequence[pos1] == "A"):
        # A s'apparie avec U
        complement = (sequence[pos2] == "U")
    elif (sequence[pos1] == "U"):
        # U s'apparie avec A ou G
        complement = ((sequence[pos2] == "A") or (sequence[pos2] == "G"))
    elif (sequence[pos1] == "G"):
       # U s'apparie avec C ou U
        complement = ((sequence[pos2] == "C") or (sequence[pos2] == "U"))
    elif (sequence[pos1] == "C"):
        # C s'apparie avec G
        complement = (sequence[pos2] == "G")

    return complement


# SCORE
# calcule les scores des séquences
# param listeapp : listes des paires de nucléotides s'appariant dans la séquence cible
# param sequences : une population de séquences (liste de séquences)
# param consensus : un dictionnaire dont les clefs sont les nucléotides consensus,
# et les valeurs sont les positions où on les rencontre dans la séquence
# return : scores des séquences (liste de même taille que la liste de séquences: un score par séquence)

#score V1 :
def scoreV1(listeapp, sequences, consensus):
     ListeScore = []                             # Une liste initialement vide qui contiendra tous les scores
     for k in range(0, len(sequences)):          # Pour chaque nucléotides de la séquence
         score = 0                               # Le score est initialisé à 0
         for i in range(0, (len(listeapp))):     # Pour chaque appariement
                 pos1 = listeapp[i][0]           # La position des nucléotides s'appariants dans la formation des boucles
                 pos2 = listeapp[i][1]
                 if (len(sequences[k])) >= pos1 and (len(sequences[k])) >= pos2:     # conditions pour éviter les erreurs
                     print(len(sequences[k]))                                                                  # de types liste index out of range
                                                                                     # si la séquence se raccourcie
                     if complementaire(sequences[k], pos1, pos2) == True:            # Si les nucs peuvent s'apparier
                         score += 1                                                  # On incrémente le score de 1
         for cle in consensus:                   # Pour A, G, C, U qui sont les clés de mon dic (voir prog principal)
             value = consensus.get(cle)          # On recupere les valeurs : donc les positions ex : (33, 8) pour U
             for o in range(len(value)):         # Pour chacune de ces valeurs (len(value))
                 n = value[o]                    # On recupere dans n une des valeurs : ex : value[O] = 33

                 if (len(sequences[k])-1) >= n and sequences[k][n] == cle :          # Conditions pour éviter les erreurs
                                                                                     # de types liste index out of range
                                                                                     # si la séquence se raccourcie
                                                                                     # Si le nucleotide de la sequence k,
                                                                                     # à la position n = cle
                                                                                     # Ex : sequences[k][33] == U
                     score += 1                  # alors le score augmente + 1
         ListeScore.append(score)                # On ajoute le score obtenu à ListeScore
     return ListeScore




# scoreV2 : permet de tester le phénotype sur des séquences plus longues que 75 nucleotides
# run trop long
def scoreV2(listeapp, sequences, consensus):
    ListeScore = []
    # Pour chaque nucléotides de la séquence
    for k in range(0, len(sequences)):
        score = 0
        #si la séquence a 75 ou - de nucleotides
        if len(sequences[k]) <= 76:
            # Pour chaque appariement
            for i in range(0, (len(listeapp))):
                     # La position des nucléotides s'appariants dans la formation des boucles
                    pos1 = listeapp[i][0]
                    pos2 = listeapp[i][1]

                    # conditions pour éviter les erreurs de types liste index out of range
                    if (len(sequences[k])) >= pos1 and (len(sequences[k])) >= pos2:

                        if complementaire(sequences[k], pos1, pos2) == True:
                            score += 1
            # Pour A, G, C, U qui sont les clés de mon dic
            for cle in consensus:
                #ex : (33, 8) pour U
                value = consensus.get(cle)
                for o in range(len(value)):
                    #ex : value[O] = 33
                    n = value[o]
                    # Conditions pour éviter les erreurs de types liste index out of range si la séquence se raccourcie
                    # si le nucleotide de la sequence k, à la position n = cle
                    if (len(sequences[k])-1) >= n and sequences[k][n] == cle:
                        score += 1
            ListeScore.append(score)
        else:
            # pour les séquences plus longues que 75 nucs => analyse des scores selon diff cadres de lecture
            scorexl = phenotypexl(listeapp, sequences[k], consensus)
            ListeScore.append(scorexl)
    print(ListeScore)
    return ListeScore


# calcul le meilleur score d'une séquence > 75nucs selon diff cadres de lecture
# param listeapp : listes des paires de nucléotides s'appariant dans la séquence cible
# param sequence : une séquence = une liste de nucs
# param consensus : un dictionnaire dont les clefs sont les nucléotides consensus,
# et les valeurs sont les positions où on les rencontre dans la séquence
# return : le score maximum de cette fonction selon le cadre de lecture
def phenotypexl(listeapp, seq, consensus):
    scoresinter = []
    for k in range(0, len(seq)):
        # cadre = nombre de cadlres de lecture possibles
        cadre = len(seq)-76
        for j in range(cadre):
            score = 0
            for cle in consensus:
                value = consensus.get(cle)
                for o in range(len(value)):
                    # on demarre le cadre de lecture un nuc plus tard
                    n = value[o] + j
                    if seq[n] == cle:
                        score += 1
            for k in range(0, (len(listeapp))):
                    # on demarre le cadre de lecture un nuc plus tard
                    pos1 = listeapp[k][0] + j
                    pos2 = listeapp[k][1] + j
                    if complementaire(seq, pos1, pos2) == True:
                        score += 1
            # les scores de chaque cadres de lecture sont stockés dans cette liste
            scoresinter.append(score)
    # on selectionne le meilleur score
    score = Fctmax(scoresinter)
    return score

#REPRODUCTION
# génère une nouvelle population de séquences, en favorisant les séquences de plus grand score
# en induisant ou pas des mutations, et des recombinaisons
# param sequences : liste des séquences générées précédemment
# param N : taille de la population = nombre de séquences
# return : une nouvelle population, et la liste des scores associé à chaque séquence
def reproduction(sequences, N):
    enfants = []
    # On fait la liste des scores des séquences de la pop en question
    ListeScore = scoreV2(listeapp, sequences, consensus)
    #ListeScore = scoreV1(listeapp, sequences, consensus)
    while len(enfants) < (N):
        a = random.randint(0, 45)
        j = random.randint(0, (N-1))
        # On tire la séquence associée au score grace à l'indice du score
        scorej = ListeScore[j]
        # scorj chance sur 45 que a soit <= scorej
        if scorej > a:
            NewSeq = mutation(sequences[j])
            enfants.append(NewSeq)
    TirageRecomb = random.randint(0, 100)
    if TirageRecomb == 16:
        c = random.randint(0, (len(enfants)-1))
        d = random.randint(0, (len(enfants)-1))
        # Au cas ou c == d et donc qu'on tire deux fois la même séquence,
        while c == d:
            c = random.randint(0, (len(enfants)-1))
            d = random.randint(0, (len(enfants)-1))
        sequence1 = enfants[c]
        sequence2 = enfants[d]
        recombines = recombinaison(sequence1, sequence2)
        enfants[c] = recombines[0]
        enfants[d] = recombines[1]
    ListeScore = scoreV2(listeapp, sequences, consensus)
    # on retourne une liste avec les enfants et les scores pour pouvoir
    # retourner ces deux elements en même temps et les utiliser séparément
    return [enfants, ListeScore]


# MUTATION
# induit une mutation ou pas
# param seq : une séquence
# return : la même séquence mutée ou pas
def mutation(seq):   # seq = un seul ARNt
    # 1 chance sur 100 d'avoir une mutation, 1/300 chances pour chaque type
    tirage = random.randint(1, 100)
    seqbis = copy.deepcopy(seq)
    augc = ["A", "U", "G", "C"]
    pos_mutation = random.randint(0, (len(seq)-1))
    #substitution
    if tirage == 1:
        print("s")
        nuc = augc[random.randint(0, 3)]
        # Au cas ou le nuc de remplacement est le même que l'original
        while nuc == seqbis[pos_mutation]:
            nuc = augc[random.randint(0, 3)]
        seqbis[pos_mutation] = nuc
    #deletion
    elif tirage == 2:
        print("d")
        del seqbis[pos_mutation]
    #insertion
    elif tirage == 3:
        print("i")
        nuc = augc[random.randint(0, 3)]                # un nucléotide choisi au hasard
        seqbis.insert(pos_mutation, nuc)
    else:
        print("pas de mutation")
    return seqbis


# MOYENNE
# calcul la moyenne des score pour une population
# param ListeScore : la liste des scores pour chaque séquence de la population
# return : la moyenne des scores de la ListeScore
def MoyenneScore(ListeScore):
    moy = 0
    for i in range(0, (len(ListeScore))):
        moy += ListeScore[i]
    else:
        moy = moy/(len(ListeScore))
    return moy

# RECOMBINAISON
# effectue ou non une recombinaison entre deux séquences
# param sequence1 : une séquence
# param sequence2 : une autre séquence
# return : ces 2 séquences recombinées ou pas
def recombinaison(sequence1, sequence2):
    # On tire soit un 1 soit un 0 au hasard, ce qui va definir si la recombinaison
    # aura lieu sur l'extrémité 5' ou 3'
    ext = random.randint(0, 1)
    # Extrémité 5'
    if ext == 0:
        # On tire 2 nombres au hasard qui vont définir le site de cassure dans la première moitié de la sequence
        pos1 = random.randint(0, (len(sequence1)-1)//2)
        pos2 = random.randint(0, (len(sequence2)-1)//2)
        seq1 = sequence1[0:pos1]
        seq2 = sequence2[0:pos2]
        sequence1bis = seq2 + sequence1[pos1:]
        sequence2bis = seq1 + sequence2[pos2:]

    # Extrémité 3'
    if ext == 1:
        pos1 = random.randint((len(sequence1)-1)//2, (len(sequence1)-1))
        pos2 = random.randint((len(sequence2)-1)//2, (len(sequence2)-1))
        seq1 = sequence1[pos1:]
        seq2 = sequence2[pos2:]
        sequence1bis = sequence1[0: pos1] + seq2
        sequence2bis = sequence2[0: pos2] + seq1

    return [sequence1bis, sequence2bis]

def Fctmax(ListeScore):
    c = ListeScore[0]
    for X in ListeScore:
        if X > c:
            c = X
    return c


#P R O G R A M M E  P R I N C I P A L

#population de depart
# l = longueur de la séquence
l = 75
# N = nombre de chaines
N = 1000
enfants = PopDeDepart(l, N)

#sequences consensus
consensus = {"C": (11, 25, 32, 48, 60, 62, 74, 75), "U": (8, 11, 25, 32, 33, 48, 60, 62), "A": (14, 21, 58, 9, 10, 15, 24, 37, 52, 57, 75), "G":(19, 18, 53,9, 10, 15, 24, 37, 52, 57)}

#appariements : repli pour ListeAppariements
repli1 = [1, 7, 72]
repli2 = [10, 13, 25]
repli3 = [27, 31, 43]
repli4 = [49, 53, 65]
repli = [repli1] + [repli2] + [repli3] + [repli4]

#liste de liste contenant 2 positions : la position des appariements de la sequence cible
listeapp = (ListeAppariements(repli))


#reproduction
nbr_generation = 300
# liste qui contiendra les moyennes des scores au fil des générations
ListeMoy = []
moy = 0
#liste des scores max
LSMax = []

for i in range(nbr_generation):
    # enfants ici = PopDeDépart
    EnfantsScore = reproduction(enfants, N)
    # en postion 1 dans la liste on a retourné Liste score: (return [enfants, ListeScore])
    ListeScore = EnfantsScore[1]
    moy = MoyenneScore(ListeScore)
    ListeMoy.append(moy)
    #ajout du score maximal de la génération à la liste de scores maximaux de toutes les générations
    LSMax.append(Fctmax(ListeScore))
    # en position 0, on avait retourné la liste des nv sequences
    enfants = EnfantsScore[0]





# G R A P H I Q U E
abs = []
for i in range(nbr_generation):
     # numero de la génération
    abs.append(i+1)

# (numero de la génération en x, moyenne en y)
plt.plot(abs, ListeMoy, label="Score moyen")
# titre axe des y
plt.ylabel('Score')
plt.plot(abs, LSMax, label="Score max")
# titre axe des x
plt.xlabel('Nombre De Generation')
# affichage
plt.show()
