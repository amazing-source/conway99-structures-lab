# Journal du labo « structures »

## 25/09 nuit
- Ouverture du labo après la clôture de la transposition du contenu de C3 (c1_recherche/notes/113).
- Partage convenu avec C1 : C1 fait c1_experimental (graphe lointain d'une droite, réseaux lointains) ; ce labo fait
  D11 + D2 (lemmes de C1 sur de vrais graphes, balayage des tables).
- S1 : genus_scan.py (classe propre : d impair sans carré, premier à v). Pour 99, redonne exactement L22 (plafond exclu,
  α pair). Validé sur Petersen et T(7)…T(21) (réseaux exacts). Le seul écart, T(5), vient du harnais (complément mal
  apparié), pas de la théorie.
- En chemin, 12 ensembles « sans α admissible ». Diagnostic : tous tombent au niveau le plus simple, le discriminant
  mod 4 d'un réseau pair de déterminant impair → S5.
- S5 : disc4_scan.py → 36 exclusions pour v ≤ 1300, 185 pour v ≤ 5000. Contrôle contre la table de Brouwer téléchargée
  (brouwer_crosscheck.py) : 0 graphe existant exclu sur 1183 ; (289,54,1,12) déjà connu inexistant ; 35 « ? » exclus.
  Test décisif à v = 625 et 841 : NL_e connus tous à e pair. Rédigé dans NOTE_S5_DISCRIMINANT_MOD4.md.
