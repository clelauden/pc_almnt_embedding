# Script d'alignement de protéine à partir de valeurs d'embedding
# Selon 3 méthodes : Smith-Waterman, Needleman-Wunsch et Semi-Global

# Chargement des Modules
import numpy as np

def AjoutColonneLigne(matrice):
# Ajout de la ligne vide et la colonne vide en tête à la matrice de dot product
     addrow = np.zeros((1, matrice.shape[1]))
     matrice_1 = np.concatenate((addrow,matrice),axis = 0)
     addcol = np.zeros((matrice_1.shape[0], 1))
     matrice_2 = np.concatenate((addcol,matrice_1),axis = 1)
     return matrice_2

def matrix_sw(seq1, seq2, siml, gap=0):
# Production la matrice de score avec la méthode de Smith-Waterman
    sim = AjoutColonneLigne(siml)
    seq1.insert(0, '')
    seq2.insert(0, '')
    H = np.zeros((len(seq1), len(seq2)))
    for i in range(1,len(seq1)):
        for j in range (1,len(seq2)):
            asso = H[i-1, j-1] + (sim[i,j] if seq1[i-1] == seq2[j-1] else - abs(siml[i,j]))
            dele = H[i-1,j] - gap
            inse = H[i, j-1] - gap
            H[i,j] = max(asso,dele,inse, 0)
    return H

def matrix_nw(seq1nw, seq2nw, simlnw, gapnw=0):
# Production la matrice de score avec la méthode de Needleman-Wunsch
    simnw = AjoutColonneLigne(simlnw)
    seq1nw.insert(0, '')
    seq2nw.insert(0, '')
    Hnw = np.zeros((len(seq1nw), len(seq2nw)))
    for i in range(1,len(seq1nw)):
        for j in range (1,len(seq2nw)):
            assonw = Hnw[i-1, j-1] + (simnw[i,j] if seq1nw[i-1] == seq2nw[j-1] else - abs(simnw[i,j]))
            delenw = Hnw[i-1,j] - gapnw
            insenw = Hnw[i, j-1] - gapnw
            Hnw[i,j] = max(assonw,delenw,insenw)
    return Hnw


def backtrak_sw(seq1mn, seq2mn, M):
# Réalisation du backtracking et de l'alignement avec la méthode de Smith-Waterman
    seq1mn.insert(0, '')
    seq2mn.insert(0, '')
    coord = np.max(M)
    c , d = np.where(M==coord)
    i = np.max(c)
    j = np.max(d)
    al1 = ''
    al2 = ''
    while M[i,j] > 0:
        if M[i-1,j-1] > M[i-1,j] and M[i-1,j-1] > M[i,j-1] or M[i-1,j-1] == M[i-1,j] or M[i-1,j-1] == M[i,j-1]:
            if i == 0 or j == 0:
                i = i-1
                j = j-1
            else :
                al1 = al1 + seq1mn[i]
                al2 = al2 + seq2mn[j]
                i = i-1
                j = j-1
        elif M[i-1,j] > M[i-1,j-1] and M[i-1,j] > M[i,j-1]:
            if i == 0 or j == 0:
                i = i-1
                j = j
            else:
                al1 = al1 + '-'
                al2 = al2 + seq2mn[j]
                i = i-1
                j = j
        elif M[i,j-1] > M[i-1,j-1] and M[i,j-1] > M[i-1,j]:
            if i == 0 or j == 0:
                i = i
                j = j-1
            else :
                al1 = al1 + seq1mn[i]
                al2 = al2 + '-'
                i = i
                j = j-1
    return al1, al2


def backtrak_nw(seq1mn, seq2mn, M):
# Réalisation du backtracking et de l'alignement avec la méthode de Needlman-Wunsch
    seq1mn.insert(0, '')
    seq2mn.insert(0, '')
    i = 0
    j = 0
    al1 = ''
    al2 = ''
    while i < len(seq1mn)-1 and j < len(seq2mn)-1:
        if M[i-1,j-1] > M[i-1,j] and M[i-1,j-1] > M[i,j-1] or M[i-1,j-1] == M[i-1,j] or M[i-1,j-1] == M[i,j-1]:
            if i == 0 or j == 0:
                i = i+1
                j = j+1
            else :
                al1 = al1 + seq1mn[i]
                al2 = al2 + seq2mn[j]
                i = i+1
                j = j+1
        elif M[i-1,j] > M[i-1,j-1] and M[i-1,j] > M[i,j-1]:
            if i == 0 or j == 0:
                i = i+1
                j = j
            else :
                al1 = al1 + '-'
                al2 = al2 + seq2mn[j]
                i = i+1
                j = j
        elif M[i,j-1] > M[i-1,j-1] and M[i,j-1] > M[i-1,j]:
            if i == 0 or j == 0:
                i = i
                j = j+1
            else :
                al1 = al1 + seq1mn[i]
                al2 = al2 + '-'
                i = i
                j = j+1
    return al1, al2

def backtrak_sg(seq1mn, seq2mn, M):
# Réalisation du backtracking et de l'alignement avec la méthode de Semi-global
    seq1mn.insert(0, '')
    seq2mn.insert(0, '')
    coord = np.max(M[:,len(seq2mn)-1])
    c , d = np.where(M==coord)
    i = np.max(c)
    j = len(seq2mn)-1
    al1 = ''
    al2 = ''
    while i > 0 or j > 0 :
        if M[i-1,j-1] > M[i-1,j] and M[i-1,j-1] > M[i,j-1] or M[i-1,j-1] == M[i-1,j] or M[i-1,j-1] == M[i,j-1]:
            if i == 0 or j == 0:
                i = i-1
                j = j-1
            else :
                al1 = al1 + seq1mn[i]
                al2 = al2 + seq2mn[j]
                i = i-1
                j = j-1
        elif M[i-1,j] > M[i-1,j-1] and M[i-1,j] > M[i,j-1]:
            if i == 0 or j == 0:
                i = i-1
                j = j
            else:
                al1 = al1 + '-'
                al2 = al2 + seq2mn[j]
                i = i-1
                j = j
        elif M[i,j-1] > M[i-1,j-1] and M[i,j-1] > M[i-1,j]:
            if i == 0 or j == 0:
                i = i
                j = j-1
            else :
                al1 = al1 + seq1mn[i]
                al2 = al2 + '-'
                i = i
                j = j-1
    return al1, al2

def alignement (prot1, prot2, dotp):
# Lacement de l'alignement choisi par l'input
    if input_alignement == 'SW':
        result = matrix_sw(prot1, prot2, dotp)
        bt = backtrak_sw(prot1, prot2, result)
    elif input_alignement == 'NW':
        result2 = matrix_nw(prot1, prot2, dotp)       
        bt = backtrak_nw(prot1, prot2, result2)
    elif input_alignement == 'SG':
        result3 = matrix_nw(prot1, prot2, dotp)       
        bt = backtrak_sg(prot1, prot2, result3)
    return bt

def result_file(prot1, prot2, dotp):
# Création du fichier de sortie avec les alignements
    almt1 , almt2 = alignement (prot1, prot2, dotp)
    almt1_rev = almt1[::-1]
    almt2_rev = almt2[::-1]
    new_file = open('alignement.txt', 'w')
    new_file.write(almt1_rev)
    new_file.write("\n")
    new_file.write(almt2_rev)
    new_file.close()
    print('File created')
    return


# Input des séquences protéiques et des embedding
input_emb1 = input ('Enter the name of the file with the 1st embedding:')
input_emb2 = input ('Enter the name of the file with the 2nd embedding:')
input_prot1 = input ('Enter the name of the file with the 1nd protein:')
input_prot2 = input ('Enter the name of the file with the 2nd protein:')

# Input du type d'alignment
input_alignement = input ('Alignment method (SW, NW or SG):')

#Lire des fichier d'embedding
emb1=np.loadtxt(input_emb1)
emb2=np.loadtxt(input_emb2)

# Lecture et mise en forme de séquences de protéines
with open(input_prot1) as f:
    pr1 = f.readlines()
    pro1 = pr1[1:]

with open(input_prot2) as f:
    pr2 = f.readlines()
    pro2 = pr2[1:]

# On transforme les séquences en des listes avec autant d'éléments qu'il y a d'AA dans la séquence
listAA1 = ''.join(pro1)
prot1 = list(listAA1)
listAA2 = ''.join(pro2)
prot2 = list(listAA2)

# Calcul du dot product entre les 2 matrices d'embedding
# On commence par retourner la 2ème matrice
emb2_t = emb2.T
# Puis on réalise le dot product
dotp = np.dot(emb1, emb2_t, out = None)

# Lancement de la fonction produisant l'alignement
result_file(prot1, prot2, dotp)


#6PF2K_1bif.t5emb
#SKI_1shka.t5emb
#6PF2K_1BIF.fasta
#SKI_1SHKA.fasta
