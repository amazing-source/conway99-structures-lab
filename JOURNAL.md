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
- S10 : doubles_budget3 (système local avec C8), calcul exhaustif. Les types à ≥ 3 triangles sont tous exclus, et une
  partie des types à 2 triangles. explain_blocks.py met au jour deux mécanismes : l'orbite morte (lemme S8.2, preuve à la
  main −4 = 0 ; une autre forme force φ ≥ 3) et le triangle forcé dans le 6-cycle (Fano, S8.3).
  Contrôle croisé par solveur indépendant au niveau des étiquettes : 670/670 clés d'accord.
  Résidu tout-tordu : 14 types (6 sans triangle, 5 à 1 triangle, 3 à 2 triangles, en partie).
- Nouvelle règle de l'utilisateur, transmise à C1 et C2, mémoire mise à jour : les mathématiques d'abord ; solveurs et
  calcul seulement en petits tests de 5 minutes au plus, s'ils sont primordiaux pour une vérification importante, et
  jamais comme appui.
- Familles sans chaîne à la main pour C2 : Fano, étoile (6,6), étoile B et étoile + transversale acceptées par C2.
  Étoile A envoyée (48 > 42). Triangle de triangles : format à convenir.
  Erreur de nom relevée par C2 : « quadrilatère » au lieu de « étoile + transversale ». Corrigée dans la note.
- Recalculs indépendants pour C1 : lemme du triangle (complément de Schur 3I − P_∂, cas tendu) et certificat 3 (Sturm,
  une racine positive ≈ 0,0543 ; admissibilité locale OK). Confirmés.
- Triangle de triangles (famille 4) à la main : 43 > 42, grâce au **lemme du complément** (CL : si k porte {k,c,d} avec
  c, d ∈ R_x, alors φ(ab) = 3 − N_ab sur la paire complémentaire). Il faut trois faits en plus : X1, X2, X3. C2 a tout
  relu et conclut « TIENT » ; CL devient un lemme général à son §26. Tous les types à au moins 3 triangles sont donc
  exclus à la main. Contrôle : verify_pasch_hand.py, avec un solveur propre en noms d'origine qui concorde sur 1340 clés.
- Corollaire S8.5 de CL : **deux triangles disjoints ⇒ impossible**, sans triangle propre, par un argument de parité au
  7ᵉ bloc. Il tue exactement les 5880 décompositions à deux triangles que le calcul exhaustif tuait d'emblée. Restent
  les 27 720 à deux triangles sécants, dont les bornes locales valent 26 à 34, et les types à 1 ou 0 triangle.
- S11 (25/09, soir) : **budget exact paire par paire** (PPE).
  - C'est 16.1 de C2 en tout-tordu. Chaque paire est doublée par exactement 2 orbites, a 12 uns et 6 zéros. Le total 42
    est donc exact, alors que mon calcul exhaustif n'utilisait que Σ ≤ 42.
  - Avec le **couplage** E[P][Q] = E[Q][P], il ferme les 8 classes à deux triangles sécants.
  - Nouveaux lemmes : le lemme T (R_x = triangle ∪ {r} ⇒ y_r = 1 − n_k − n_r, jamais de double) et le lemme A (deux
    triangles sécants ⇒ les paires de base se doublent mutuellement et sont saturées).
  - Preuves à la main, étapes locales vérifiées par solveur exact. Des tests PLNE/PL courts (moins de 5 min chacun,
    règle du 25/09) ont servi à localiser les arguments.
  - Limite : les types à 1 triangle restent faisables pour ces outils.
- Rapport d'Astra reçu de l'utilisateur (recherches indépendantes : 5 triangles disjoints reliés deux à deux par 2
  arêtes impossibles, contraintes mod 3, tests spectraux sur témoins à 51 sommets). Consigne de l'utilisateur : ne pas
  se réorienter dessus sauf si c'est utile. Rien n'a été envoyé à Astra.
- S13, PPE signé (accord de C2). μ = 2 par sommet (t,ε) donne des types aux doubles :
  - (L1) un y-double est de type r ;
  - (L2) les doubles incidents à un bloc ne sont jamais mélangés ; un double proche est toujours du type du bloc non
    partagé, ce qui vaut aussi dans les mélanges M1 ;
  - (P) une orbite de chaque type par paire.
  - Test court (règle du 25/09) sur les 3 échantillons à 1 triangle que le modèle non signé laissait passer : PLNE
    indécis en 45 s, et **relâchement linéaire faisable**. Aucun certificat linéaire : le signe au premier ordre ne
    suffit pas seul.
- S12, vérification demandée par C1 (notes/116). (P2) et sa diagonale, le critère lointain et le Fano tendu (table des
  g, κ, F5) sont refaits à la main et contrôlés exactement sur Paley 9 et BvLS : TIENT. Précision : la diagonale vaut
  16 − 2e(R) ; sur les vrais graphes, e(R) = 4 et elle vaut 8.
- S14 : (P2) de C1 transposé au cadre de C2.
  - La diagonale de (P2) redonne exactement 16.1/PPE. Les identités hors diagonale sont nouvelles : 105 couples de
    paires disjointes, qui donnent 16 − E[R][R′] par orbites, et 105 couples partageant un bloc, qui donnent
    12 − N/2.
  - La valeur val_Z(Q) vaut : φ/y locale si Q est disjointe de P_Z ; n_r si Q contient le bloc libre de Z ; 2[k = r]
    si Q contient son bloc d'appariement ; 0 si Q = P_Z.
  - Formes exactes contrôlées sur BvLS.
  - Test court : relâchement McCormick faisable sur les 3 échantillons à 1 triangle. La PLNE n'a pas pu tourner dans
    le budget, et je ne l'ai pas prolongée (règle du 25/09).
- C1 a corrigé notes/116 d'après S12 (commit b9fb15e dans c1_recherche) : la diagonale générale vaut 16 − 2e(R).
- Correction répercutée (C1, d02ac18). « Torsion maximale ⇒ e(R) = 0 » était faux : la torsion maximale ne porte que
  sur τ, et σ_c(v) = v̄ (coïncidences S de L18) peut créer une arête interne au bloc.
  - NOTE_S12 et le registre S12 sont corrigés. Au Fano tendu, e(R) = 0 tient par tension.
  - Rien ne change pour le tout-tordu (S8, S11, S14) : il n'utilise que n_s = 0, théorème de C2 pour les carrés
    « un bloc », avec e(R) = 2n_s.
