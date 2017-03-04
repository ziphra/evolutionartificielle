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
        tirage = random.randint(0, 3)
        nuc = atgc[tirage]
        seq.append(nuc)
    sequences.append(seq)

# afficher les 10 premières séquences
for i in range(10):
    print(sequences[i])

# phénotype
# appariements de la séquence cible, sous forme d'une liste de paires de positions
# fonctions qui génère la liste des appariements
def ListeAppariements(a,b, pos1, pos2):
    app = []  # là ou sera rangé les paires de position correspondant aux appariements
    for i in range (a,b):
        pos = [] #là ou sera rangé la paire de positions d'un appariement à chaque tour de boucle. Doit être vidé à chaque boucle pour pouvoir stocker la nouvelle position, qui est enregistré dans la liste app.
        pos1 += 1
        pos2 += 1
        pos.append(pos1)
        pos.append(pos2)
        app.append(pos)
    return app


# premier repliement de la séquence cible : 
a = pos1 = 1
pos2 = 66
b = 7
appariement1 = ListeAppariements(a,b, pos1, pos2)
print (appariement1)