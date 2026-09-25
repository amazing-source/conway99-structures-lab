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
- S3 : s3_lemmas_real.py. L21, forme générale : S = A_G − rI + (1/c)(D + J/(c−3)), c = r + 1. Pour BvLS et srg9,
  S|_U = 0 et le rang local vaut rang(A_G − rI) ; 891/891 droites, contrôle négatif 100 %. L23 : le signe des blocs de
  rang plein égale le ε(V) du genre sur T(7), T(9), T(13), T(15). Le contrôle négatif (mauvais multiple de J) ne
  discrimine pas : à noter.
- S6 : c2_balance_check.py. Équilibre de C2 en forme générale, BBᵀ1_C = (k−14)1_C + (5/3)|C|1_N, vérifié sur BvLS.
  Relecture du lemme de la queue (stabilité de X_In(k), v et σv dans C_Out(k), conséquences) : tient. Envoyé à C2.
- S7 (aide à C2) : partage acté, C2 garde les mélanges presque ancrés, je prends le tout-tordu pur.
  k7_cycle_decomp.py : 121 215 décompositions étiquetées, 20 types.
  star_S_spectrum.py : W = intérieurs + étoile(k), 7680 étiquetages, spectre unique dans [−4, 3] ; les extrêmes sont les
  vecteurs triviaux θe_k + Ce_k. En tout-tordu |D_h(x)| ≤ 2 : budgets non saturés. Aucun levier local ; (7,7,7) est le
  cœur.
- S9 : 27.2 de C2 relu, il tient (équilibre croisé en 0 et en 3).
- S8 : zeros_budget.py. Deux contraintes ajoutées après que le contrôle Fano a d'abord donné 42, puis 70 :
  (a) φ_x(jk) = 0 si x est dans le triangle {i,j,k} ; (b) au plus un R_j-voisin par triangle {j,r,s} ∌ i (compte exact 1).
  Fano donne alors 84 = compte de C2. (3⁵,6) : 48 > 42, exclu ; (3⁴,4,5) : 32–38, survit.
  Le cœur est redéfini avec C2 comme « G_c connexe en tout bloc » ; le résidu tout-tordu compte 18 types.
- Consigne de l'utilisateur : « be extremely careful with errors … once you're done with your part report to your peers ».
  Pratique adoptée : contrôle positif (Fano = 84), solveur indépendant au niveau des étiquettes, recomptage
  indépendant de la classification, tests unitaires, et rien d'envoyé avant la fin des contrôles.
- S8 : preuve à la main de l'exclusion de (3⁵,6) (NOTE_S8_TOUT_TORDU.md), relue ligne à ligne. La borne vaut 48 et
  e(S₁,S₂) ≥ 2 vient d'un argument de degrés.
- Certificat d'Astra (mission géométrique) vérifié par code indépendant (NOTE_ASTRA_CERTIFICAT.md) : wᵀCw = −2 et
  qᵀAq − 3‖q‖² = 16. Ses 8 dimensions sont exactement le Π₋₄ de rang 8 du X12 de C1 pour la droite {x, 0, 1}.
