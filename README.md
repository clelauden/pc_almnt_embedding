# pc_almnt_embedding

Ce programme a été conçu pour réaliser un alignement de deux protéines à partir de valeurs d'embedding en utilisant trois méthodes d'alignement différentes, Local (Smith-Waterman), Global (Needleman-Wunsch) et Semi-global.
Pour que le programme fonctionne les fichiers d'input doivent être dans le dossier principal pc_almnt_embedding
Les fichiers 6PF2K_1bif.t5emb, SKI_1shka.t5emb, 6PF2K_1BIF.fasta et SKI_1SHKA.fasta sont des exemples d'input pour faire fonctionner le programme. D'autres fichiers pouvant être utilisé en exemple sont fourni dans le dossier entree.
Des exemples d'alignement réalisé avec le programme sont également fourni : expl_al_6PF2K_SKI_NW.txt et expl_al_6PF2K_SKI_SW.txt
 ils ont été produit ont utilisant 6PF2K_1bif.t5emb, SKI_1shka.t5emb, 6PF2K_1BIF.fasta et SKI_1SHKA.fasta en input et respectivement les méthode Needleman-Wunsch et Smith-Waterman.

Suite de commande en bash a exécuter pour lancer le programme d'alignement

Télécharger le répertoire
	
	git clone https://github.com/clelauden/pc_almnt_embedding.git

Si conda n'est pas installé, installer conda 

	Install conda

Se placer dans le dossier

	cd pc_almnt_embedding

Créer l'environnement

	conda env create -f env_pc_al.yml

Activer l'environnement

	conda activate pc_env

Lancer le programme

	python pc_almnt.py

Entrer les fichier des embedding et des séquences de protéines (avec les extensions de fichiers), pour chaque protéine, il faut une embedding et un fichier fasta de la séquence 
Exemple : 

	6PF2K_1bif.t5emb
	SKI_1shka.t5emb
	6PF2K_1BIF.fasta
	SKI_1SHKA.fasta
	
Selectionner la méthode d'alignement

Pour le local (Smith-Waterman) : SW
Pour le global (Needleman-Wunsch) : NW
Pour le semi-global : SG

Le programme créer un fichier alignement.txt avec l'alignement dans le dossier pc_almnt_embedding




