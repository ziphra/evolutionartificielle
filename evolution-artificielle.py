__authors__ = 'euphrasieservant', 'agathewallet'
# -*- coding:utf-8 *-*
# population de départ : N séquences d'ARNt de longueur l =
# chaque séquence est constitué de l nucléotides choisis au hasard
import random
import matplotlib.pyplot as plt


# constitution de la population de départ
def PopDeDepart(l, N):

    augc = ["A", "U", "G", "C"]             # Nucléotides
    sequences = []                          # Une liste initialement vide qui contiendra LES séquences de la pop de départ

    for i in range(N):                      # On crée N séquences en passant N fois dans cette boucle
        seq = []                            # Une liste initialement vide qui contiendra UNE séquence
        for y in range(l):                  # La séquence contiendra l nucléotide en passant l fois dans la boucle
            tirage = random.randint(0, 3)   # On choisit un nombre au hasard compris entre 0 et 3
            nuc = augc[tirage]              # qui representera l'indexe d'un nucléotide dans la liste augc
            seq.append(nuc)                 # On ajoute ce nucléotide à la sequence
        sequences.append(seq)               # À l'issu de cette boucle on ajoute la séquence à la population
    return sequences


# PHÉNOTYPE
# fonction qui génère la liste des appariements
# sous forme d'une liste de paires de positions
def ListeAppariements(repli):               # Repli est une liste contenant 3 éléments : le premier  et le dernier
                                            # nucléotide dans le sens 5'-3' du repli, et le premier dans le sens 3'-5'
    listeapp = []                           # Une liste initialement vide qui contiendra les paires de positions des
                                            # nucléotides appariés dans la formation des boucles
    for i in range(0, len(repli)):          # Pour tous les repli (repli1, repli2, repli3, repli4)
        a = repli[i][0]                     # Premier élément du repli dans le sens 5'-3'   (1 dans le cas du repli1)
        b = repli[i][2]                     # Premier élement du repli dans le sens 3'-5'   (72 dans le cas du repli1)
        c = repli[i][1]                     # Dernier element du repli dans le sens 5'-3'   (7 dans le cas du repli1)
        for j in range(a, c+1):             # Du premier élément (7), au dernier INCLU (d'ou c+1) dans le sens 5'-3'
            pos = []                        # Une liste initialement vide qui contiendra les positions des 2 nucléotides
                                            # qui s'apparient
            pos.append(a)
            pos.append(b)
            listeapp.append(pos)            # On ajoute cette liste dans listeapp
            a += 1                          # On incrémente a et b de 1 pour passer aux nucléotides suivants
            b -= 1
    return listeapp


# permet de verifier si 2 nucléotides à deux positions précises de la séquences sont complémentaires ou non
def complementaire(sequence, pos1, pos2):   # pos1 et pos2 sont les positions des nucleotides dont on cherche
                                            # à savoir la complémentarité
    # A et U ou U et A s'apparient
    if ((len(sequence)-1) >= (pos1) and (len(sequence)-1) >= (pos2)) and ((sequence[pos1] == "A" and sequence[pos2] == "U") or (sequence[pos2] == "A" and sequence[pos1] == "U")):
        complement = True

    # G et C ou C et G s'apparient
    elif ((len(sequence)-1) >= (pos1) and (len(sequence)-1) >= (pos2)) and ((sequence[pos1] == "G" and sequence[pos2] == "C") or (sequence[pos2] == "G" and sequence[pos1] == "C")):
        complement = True

    # G et U ou U et G s'apparient
    elif ((len(sequence)-1) >= (pos1) and (len(sequence)-1) >= (pos2)) and ((sequence[pos1] == "G" and sequence[pos2] == "U") or (sequence[pos2] == "G" and sequence[pos1] == "U")):
        complement = True

    # si à pos1 et à pos2 on trouve 2 nucs ne s'appariant pas on retourne False
    else:
        complement = False
    return complement                           # On retourne True si les nucs s'apparient, False sinon


# SCORE
# range les scores dans une liste
def score(listeapp, sequences, consensus):
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

#REPRODUCTION
def reproduction(sequences, N):
    enfants = []                                # Une liste vide qui contiendra les séquences qui se sont dupliquées
    ListeScore = score(listeapp, sequences, consensus)  # On fait la liste des scores des séquences de la pop en question
    while len(enfants) < (N):                   # Tant que la population ne dépasse pas N

        a = random.randint(0, 45)               # 21 appariements + 24 consensus
        j = random.randint(0, (N-1))            # Un nombre au hasard entre 0, N-1 soit pour chaque indice de la liste
        scorej = ListeScore[j]                  # On tire la séquence associée au score grace à l'indice du score
        if scorej >= a :                        # scorj chance sur 42 que a soit <= scorej
            NewSeq = mutation(sequences[j])  # mutation ou pas...
            enfants.append(NewSeq)              # stockage des séquences dans une liste
    ListeScore = score(listeapp, enfants, consensus)  # calcul des nouveaux score

    TirageRecomb = random.randint(0, 100)       # On tire un nombre au hasard
    if TirageRecomb == 16:                      # 1 chance sur 100 de tomber sur 16
        c = random.randint(0, (len(enfants)-1)) # On choisit 2 séquences au hasard en tirant un nombre aléatoire
        d = random.randint(0, (len(enfants)-1))
        while c == d:                           # Au cas ou c == d et donc qu'on tire deux fois la même séquence,
                                                # On retire aléatoirement 2 nombre jusqu'a ce qu'ils soient différents
            c = random.randint(0, (len(enfants)-1))
            d = random.randint(0, (len(enfants)-1))
        sequence1 = enfants[c]                  # ce nombre représente l'index de la séquence dans la liste
        sequence2 = enfants[d]
        recombinés = recombinaison(sequence1, sequence2)
        enfants[c] = recombinés[0]
        enfants[d] = recombinés[1]

    return [enfants, ListeScore]                # on retourne une liste avec les enfants et les scores pour pouvoir
                                                # retourner ces deux elements en même temps et les utiliser séparément





# induit une mutation avec une probabilité de 1/100
def mutation(seq):   # seq = un seul ARNt
    substitution = 1
    deletion = 2
    insertion = 3
    tirage = random.randint(1, 300)                     # il ya 1 chance /100 d'avoir une mutation de 3 types différents
                                                        # 1 chance /300 pour chaque type
                                                        # On tire aléatoirement un nombre entre 1, 300
    augc = ["A","U","G", "C"]
    pos_mutation = random.randint(0, (len(seq)-1))      # Une position choisie au hasard dans la séquence de longueur l

    if tirage == substitution:                          # 3 chance /300 = 1/100
        nuc = augc[random.randint(0, 3)]                # Un nucléotide de remplacement choisi au hasard
        while nuc == seq[pos_mutation]:                 # Au cas ou le nuc de remplacement est le même que l'original
            nuc = augc[random.randint(0, 3)]            # on retire un nuc jusqu'a ce qu'il soit différent de l'original
        seq[pos_mutation] = nuc                         # On remplace

    elif tirage == deletion:                            # 3 chance /300 = 1/100
        del seq[pos_mutation]

    elif tirage == insertion:                           # 3 chance /300 = 1/100
        nuc = augc[random.randint(0, 3)]                # un nucléotide choisi au hasard
        seq.insert(pos_mutation, nuc)                   # On l'insert

    return seq

# MOYENNE
def MoyenneScore(ListeScore):
    moy = 0                                             # On initialise la moyenne = 0
    for i in range(0, (len(ListeScore))):               # Pour chaque score dans ListeScore
        moy += ListeScore[i]                            # On ajoute à la moyenne le ième score
    else:                                               # Quand on a terminé de parcourir la liste,
        moy = moy/(len(ListeScore))                     # On termine le calcul en divisant par le nombre de sequences
                                                        # definit par le nombre d'éléments dans la liste
    return moy


def recombinaison(sequence1, sequence2):
    ext = random.randint(0, 1)              # On tire soit un 1 soit un 0 au hasard, ce qui va definir si la recombinaison
                                            # aura lieu sur l'extrémité 5' ou 3'

    if ext == 0:                                            # Extrémité 5'
        pos1 = random.randint(0, (len(sequence1)-1)//2)     # On tire 2 nombres au hasard qui vont définir le site de
        pos2 = random.randint(0, (len(sequence2)-1)//2)     # cassure pour la recombinaison
        seq1 = sequence1[0:pos1]                            # On coupe les séquences
        seq2 = sequence2[0:pos2]                            # seq = le morceau de séquence qui va se recombiner
        sequence1bis = seq2 + sequence1[pos1:]              # On assemble les nouvelles séquences
        sequence2bis = seq1 + sequence2[pos2:]

    if ext == 1:                                            # Extrémité 3'
        pos1 = random.randint((len(sequence1)-1)//2, (len(sequence1)-1))    # de la même manière
        pos2 = random.randint((len(sequence2)-1)//2, (len(sequence2)-1))
        seq1 = sequence1[pos1:]
        seq2 = sequence2[pos2:]
        sequence1bis = sequence1[0: pos1] + seq2
        sequence2bis = sequence2[0: pos2] + seq1

    return [sequence1bis, sequence2bis]



#P R O G R A M M E  P R I N C I P A L

#population de depart
l = 75          # longueur de la séquence
N = 100         # nombre de chaines
enfants = PopDeDepart(l,N)

#sequences consensus
consensus = {"C": (11, 25, 32, 48, 60, 62, 74, 75), "U": (8, 11, 25, 32, 33, 48, 60, 62), "A": (14, 21, 58, 9, 10, 15, 24, 37, 52, 57, 75), "G":(19, 18, 53,9, 10, 15, 24, 37, 52, 57)}

#appariements
repli1 = [1, 7, 72]
repli2 = [10, 13, 25]
repli3 = [27, 31, 43]
repli4 = [49, 53, 65]
repli = [repli1] + [repli2] + [repli3] + [repli4]

listeapp = (ListeAppariements(repli))  #liste de liste contenant 2 positions : la position des appariements de la sequence cible


# mutation
# on définit le tirage qui induira une mutation dans la fonction mutation
substitution = 1
deletion = 2
insertion = 3

#reproduction
nbr_generation = 100
ListeMoy = []
moy = 0

for i in range(nbr_generation):                 # On répète autant de fois que de générations désirées le processus programmé
    EnfantsScore = reproduction(enfants, N)     # enfants ici = PopDeDépart
    ListeScore = EnfantsScore[1]                # en postion 1 dans la liste on a retourné Liste score
                                                # return [enfants, ListeScore]
    moy = MoyenneScore(ListeScore)              # on calcul la moyenne pour tracer le graphique
    ListeMoy.append(moy)                        # on stocke les moyennes dans une liste
    enfants = EnfantsScore[0]                   # en position 0, la liste des nv sequences





# G R A P H I Q U E
abs = []
for i in range(nbr_generation):
    abs.append(i+1)                         # numero de la génération

plt.plot(abs, ListeMoy)                     # (numero de la génération en x, moyenne en y)
plt.ylabel('Score Moyen')                   # titre axe des y
plt.xlabel('Nombre De Generation')          # titre axe des x
plt.show()                                  # affichage
