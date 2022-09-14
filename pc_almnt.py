# Script d'alignement de protéine à partir de valeurs d'embedding
# Selon 3 méthodes : Smith-Waterman, Needleman-Wunsch et Semi-Global

# Chargement des Modules
import numpy as np

def add_col_row(matrice):
# Ajout de la ligne vide et la colonne vide en tête à la matrice de dot product
# Input = matrice de dot product
     add_row = np.zeros((1, matrice.shape[1]))
     matrice_1 = np.concatenate((add_row,matrice),axis = 0)
     add_col = np.zeros((matrice_1.shape[0], 1))
     matrice_2 = np.concatenate((add_col,matrice_1),axis = 1)
     return matrice_2

def matrice_sw(seq1sw, seq2sw, simlsw, gap=0):
# Production la matrice de score avec la méthode de Smith-Waterman
# Input = sequences des protéine 1 et 2, la matrice de dotproduct et le gap
    simsw = add_col_row(simlsw)
    Hsw = np.zeros((len(seq1sw), len(seq2sw)))
    for i in range(1,len(seq1sw)-1):
        for j in range (1,len(seq2sw)-1):
            asso = Hsw[i-1, j-1] + (simsw[i,j] if seq1sw[i-1] == seq2sw[j-1] else - abs(simlsw[i,j]))
            dele = Hsw[i-1,j] - gap
            inse = Hsw[i, j-1] - gap
            Hsw[i,j] = max(asso,dele,inse, 0)
    return Hsw

def matrice_nw(seq1nw, seq2nw, simlnw, gapnw=0):
# Production la matrice de score avec la méthode de Needleman-Wunsch
# Input = sequences des protéine 1 et 2, la matrice de dotproduct et le gap
    simnw = add_col_row(simlnw)
    Hnw = np.zeros((len(seq1nw), len(seq2nw)))
    for i in range(1,len(seq1nw)-1):
        for j in range (1,len(seq2nw)-1):
            assonw = Hnw[i-1, j-1] + (simnw[i,j] if seq1nw[i-1] == seq2nw[j-1] else - abs(simnw[i,j]))
            delenw = Hnw[i-1,j] - gapnw
            insenw = Hnw[i, j-1] - gapnw
            Hnw[i,j] = max(assonw,delenw,insenw)
    return Hnw


def backtrak_sw(seq1sw_b, seq2sw_b, Msw_b):
# Réalisation du backtracking et de l'alignement avec la méthode de Smith-Waterman
# Input = sequences des protéine 1 et 2, matrice de score obtenue par la méthode de Smith-Waterman
    coord = np.max(Msw_b)
    c , d = np.where(Msw_b==coord)
    i = np.max(c)
    j = np.max(d)
    al1 = ''
    al2 = ''
    while Msw_b[i,j] > 0:
        if Msw_b[i-1,j-1] > Msw_b[i-1,j] and Msw_b[i-1,j-1] > Msw_b[i,j-1] or Msw_b[i-1,j-1] == Msw_b[i-1,j] or Msw_b[i-1,j-1] == Msw_b[i,j-1]:
            if i == 0 or j == 0:
                i = i-1
                j = j-1
            else :
                al1 = al1 + seq1sw_b[i]
                al2 = al2 + seq2sw_b[j]
                i = i-1
                j = j-1
        elif Msw_b[i-1,j] > Msw_b[i-1,j-1] and Msw_b[i-1,j] > Msw_b[i,j-1]:
            if i == 0 or j == 0:
                i = i-1
                j = j
            else:
                al1 = al1 + '-'
                al2 = al2 + seq2sw_b[j]
                i = i-1
                j = j
        elif Msw_b[i,j-1] > Msw_b[i-1,j-1] and Msw_b[i,j-1] > Msw_b[i-1,j]:
            if i == 0 or j == 0:
                i = i
                j = j-1
            else :
                al1 = al1 + seq1sw_b[i]
                al2 = al2 + '-'
                i = i
                j = j-1
    return al1, al2


def backtrak_nw(seq1sw_b, seq2sw_b, Msw_b):
# Réalisation du backtracking et de l'alignement avec la méthode de Needleman-Wunsch
# Input = sequences des protéine 1 et 2, matrice de score obtenue par la méthode de Needleman-Wunsch
    i = len(seq1sw_b)-1
    j = len(seq2sw_b)-1
    al1 = ''
    al2 = ''
    while i >=0 and j >=0:
        if Msw_b[i-1,j-1] > Msw_b[i-1,j] and Msw_b[i-1,j-1] > Msw_b[i,j-1] or Msw_b[i-1,j-1] == Msw_b[i-1,j] or Msw_b[i-1,j-1] == Msw_b[i,j-1]:
            if i == 0 or j == 0:
                i = i-1
                j = j-1
            else :
                al1 = al1 + seq1sw_b[i]
                al2 = al2 + seq2sw_b[j]
                i = i-1
                j = j-1
        elif Msw_b[i-1,j] > Msw_b[i-1,j-1] and Msw_b[i-1,j] > Msw_b[i,j-1]:
            if i == 0 or j == 0:
                i = i-1
                j = j
            else :
                al1 = al1 + '-'
                al2 = al2 + seq2sw_b[j]
                i = i-1
                j = j
        elif Msw_b[i,j-1] > Msw_b[i-1,j-1] and Msw_b[i,j-1] >= Msw_b[i-1,j]:
            if i == 0 or j == 0:
                i = i
                j = j-1
            else :
                al1 = al1 + seq1sw_b[i]
                al2 = al2 + '-'
                i = i
                j = j-1
    return al1, al2

def backtrak_sg(seq1sg_b, seq2sg_b, Msg_b):
# Réalisation du backtracking et de l'alignement Semi-global
# Input = sequences des protéine 1 et 2, matrice de score obtenue pour la méthode Semi-global
    coord = np.max(Msg_b[:,len(seq2sg_b)-1])
    c , d = np.where(Msg_b==coord)
    i = np.max(c)
    j = len(seq2sg_b)-1
    al1 = ''
    al2 = ''
    while i >= 0 or j >= 0 :
        if Msg_b[i-1,j-1] > Msg_b[i-1,j] and Msg_b[i-1,j-1] > Msg_b[i,j-1] or Msg_b[i-1,j-1] == Msg_b[i-1,j] or Msg_b[i-1,j-1] == Msg_b[i,j-1]:
            if i == 0 or j == 0:
                i = i-1
                j = j-1
            else :
                al1 = al1 + seq1sg_b[i]
                al2 = al2 + seq2sg_b[j]
                i = i-1
                j = j-1
        elif Msg_b[i-1,j] > Msg_b[i-1,j-1] and Msg_b[i-1,j] > Msg_b[i,j-1]:
            if i == 0 or j == 0:
                i = i-1
                j = j
            else:
                al1 = al1 + '-'
                al2 = al2 + seq2sg_b[j]
                i = i-1
                j = j
        elif Msg_b[i,j-1] > Msg_b[i-1,j-1] and Msg_b[i,j-1] > Msg_b[i-1,j]:
            if i == 0 or j == 0:
                i = i
                j = j-1
            else :
                al1 = al1 + seq1sg_b[i]
                al2 = al2 + '-'
                i = i
                j = j-1
    return al1, al2

def alignement (prot1_al, prot2_al, dotp_al):
# Lancement de l'alignement choisi par l'input
# Input = sequences des protéine 1 et 2, la matrice de dotproduct 
    if input_alignement == 'SW':
        result = matrice_sw(prot1_al, prot2_al, dotp_al)
        bt = backtrak_sw(prot1_al, prot2_al, result)
    elif input_alignement == 'NW':
        result2 = matrice_nw(prot1_al, prot2_al, dotp_al)       
        bt = backtrak_nw(prot1_al, prot2_al, result2)
    elif input_alignement == 'SG':
        result3 = matrice_nw(prot1_al, prot2_al, dotp_al)       
        bt = backtrak_sg(prot1_al, prot2_al, result3)
    return bt

def result_file(prot1_alt, prot2_alt, dotp_alt):
# Création du fichier de sortie avec les alignements, en inversant les alignement pour les méthodes SW et SG
# Input = sequences des protéine 1 et 2, la matrice de dotproduct 
    almt1 , almt2 = alignement (prot1_alt, prot2_alt, dotp_alt)
    almt1 = almt1[::-1]
    almt2 = almt2[::-1]
    new_file = open('alignement.txt', 'w')
    new_file.write(almt1)
    new_file.write("\n")
    new_file.write(almt2)
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

# Lecture des fichier d'embedding
emb1=np.loadtxt(input_emb1)
emb2=np.loadtxt(input_emb2)

# Lecture et mise en forme de séquences de protéines
with open(input_prot1) as f:
    pr1 = f.readlines()
    pro1 = pr1[1:]

with open(input_prot2) as f:
    pr2 = f.readlines()
    pro2 = pr2[1:]

# Transformation des séquences en des listes avec autant d'éléments qu'il y a d'AA dans la séquence
listAA1 = ''.join(pro1)
prot1 = list(listAA1)
listAA2 = ''.join(pro2)
prot2 = list(listAA2)

# Ajout d'un élement vide devant les séquences de protéines
prot1.insert(0,' ')
prot2.insert(0,' ')


# Calcul du dot product entre les 2 matrices d'embedding
# Retournement du 2ème array
emb2_t = emb2.T

# Réalisation du dot product
dotp = np.dot(emb1, emb2_t, out = None)

# Lancement de la fonction produisant l'alignement
#Input = sequences des protéine 1 et 2, la matrice de dotproduct 
result_file(prot1, prot2, dotp)
