# HMM-BreastCancer

Code source et ressources pour le projet de BIN 702 (Algorithmes pour la Bioinformatique), utilisant les modèles de Markov cachés (HMM) pour analyser les séquences d'expression génétique liées à la récidive du cancer du sein.

## Sujet
**Utilisation des modèles de Markov cachés pour l'analyse d'expression génétique dans le cancer**

L'objectif général de projet est d’implémenter et d’évaluer un algorithme basé sur les modèles de Markov cachés (HMM) dans le but d'identifier les séquences d’expression génétique associées à la récidive du cancer du sein, tout en optimisant la structure et l'entraînement des modèles pour améliorer la précision de prédiction de la progression tumorale. Les principales étapes incluent le nettoyage des données, la sélection des gènes, la construction des ensembles représentatifs, l'entraînement et l'évaluation des modèles HMM.

---

## Structure du dépôt

- **Scripts :**  
  Les scripts doivent être exécutés dans l'ordre suivant :  
  `clean_data.R` : Nettoyage des données d'expression génétique et préparation pour l'analyse.  
  `gene_set_selection.R` : Sélection des modules de gènes pertinents basés sur des critères statistiques et biologiques.  
  `gene_set_representative.R` : Détermination des représentants de module de gènes.
  `train_HMM.R` : Entraînement et évaluation des modèles de Markov cachés à partir des séquences génétiques nettoyées et sélectionnées.  

- **Dossier `data/` :**  
  Ce dossier contient les données nécessaires pour exécuter les scripts. Notez que, pour des raisons de confidentialité et de taille, ce dépôt GitHub ne contient pas l'ensemble complet des données. Veuillez télécharger les fichiers de données à partir du lien fourni ci-dessous.

---

## Données

Pour télécharger les données complètes nécessaires au projet, utilisez le lien suivant :  
[**Télécharger les données**](https://drive.google.com/drive/folders/1y0zpIOVX_JyP4txjIAnxw2NBD83u3PlL?usp=drive_link)  

Une fois téléchargées, placez le contenu (garder le nom des repertoires) dans le dossier `data/` à la racine du projet.

---

## Pré-requis

Chaque fichier R contient le code pour télécharger les packages nécéssaires ( princiapalement `BiocManager`,`HMM`, `hgu133a.db`)

