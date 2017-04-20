# Évolution Artificielle

## 2V483 project 

Dans le cadre d'un projet pour l'unité d'enseignement complémentaire *"programmation avancée en python pour la bioinformatique"*, nous avons tenté de simuler le processus de l'évolution biologique. 
Pour ce faire, nous avons pris comme modèle une molécule d'ARN de transfert (ARNt). 
L'ARNt est une petite molécule d'ADN qui permet de traduire un codon en son acide aminé correspondant lors de la traduction de l'ADN en protéines. 
Grâce au langage de programmation python, nous avons codé différents processus à l'origine de l'évolution biologique 
intervenant lors de la reproduction tel que : 
* les mutations (délétions, insertions, ou substitutions)
* la recombinaison entre séquences 
* et la sélection des séquences les mieux adaptés, en les comparant à une "séquence cible", c'est à dire dont la configuration est optimale.

Ainsi nous avons simulé la reproduction d'une population d'ARNt. Nous avons ensuite observé à l'aide d'une courbe tracée grâce à la bibliothèque python Matplotlib, l'évolution de nos séquences en analysant leur score. Le score des séquences étant une mesure de la ressemblance de nos séquences générées à la "séquence cible".  

