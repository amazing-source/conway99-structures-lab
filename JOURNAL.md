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
- Relecture pour C2 (sa note 30, a0900cc ; carnet §28, 0fed809) : **TIENT**.
  - Lemme 30.1 : une paire liée tordue en c donne f_a ≥ 8. Par l'équilibre en c (réunion fermée de 12 sommets, donc
    2 voisins par lointain), P n'a aucune racine et chaque autre Z en a une seule. A = {X₁, X₂, P} est indépendant,
    et un 6-cycle force P à 3 voisins.
  - Conséquence : H = H-β. Il n'y a aucun bloc α, par 28.5 (paire ancrée, dans les deux sens) ou par 30.1 (paire
    tordue).
- 26/09. L'utilisateur : « toi poursuis pour le moment ». Je comprends : pas de gros calcul de mon côté, je continue
  par les maths et des tests de 5 min au plus. C2 a obtenu dans sa propre session une recherche plafonnée à 3 h sur
  H-β.
- S15, modèle de configurations. Chaque orbite choisit (paire de voisins, solution locale exacte), ce qui fait 5 700
  à 6 100 colonnes par décomposition à 1 triangle. Toute valeur val_Z(Q) devient une constante de la colonne, donc
  (P2) devient linéaire, et le PL est bien plus serré que la linéarisation de McCormick.
  - Résultats (tests courts) :
    - témoin, triangle de triangles : infaisable ;
    - échantillon 0 (3,4,4,5,5) : faisable sans (P2), **PL infaisable avec (P2)** ;
    - échantillon 1 (3,4,4,4,6) : **PLNE infaisable avec (P2)** en 7 s ;
    - échantillon 2 : PL faisable avec (P2), PLNE indécise en 110 s. Les types en plus rendent le PL trop lent.
  - Certificat de Farkas pour l'échantillon 0 : 162 lignes (P2), 62 de couplage, 18 PPE, 48 de cycle, poids
    fractionnaires. La contradiction est globale, d'allure spectrale : pas de lemme court visible.
  - Lecture : le second ordre mord là où le premier ordre (PPE, couplage, signe) ne mordait pas. Une preuve à la main
    demandera un principe spectral par blocs (compressions de Rayleigh de C1 sur les 21 paires), pas une chasse aux
    doubles.
- S16 (tests courts). Coupes de sous-tours (tailles 3, 4, 5) ajoutées au modèle, suggestion de C2.
  - 15 classes non isomorphes à 1 triangle. Avec (P2) : 10 sont PL-infaisables. Les 5 plus symétriques sont
    PL-faisables mais **PLNE-infaisables** (30 à 45 s chacune).
  - **Toutes les classes à 1 triangle sont donc infaisables pour le solveur**, et faisables sans (P2).
  - Contrôle d'implémentation : sur une solution du premier ordre, la somme des 210 résidus (P2) vaut exactement 0
    (identité globale), et 178 identités sont violées.
  - Sans triangle (6 types, 12 600 colonnes) : PL avec (P2) indécis en 25 s.
  - Statut : indication forte, pas une preuve. Il faudra des certificats exacts (Farkas rationnels pour les 10
    classes PL, branchement avec certificats pour les 5), ce qui dépasse les tests de 5 min : **décision de
    l'utilisateur requise**.
  - C2 a transposé le modèle à H-β, avec la même conclusion côté (P2) : squelette circulant certifié infaisable, et
    552 classes en cours dans sa session, calcul plafonné autorisé.
- S17, calcul certifié (autorisé par l'utilisateur, plafond d'environ 1 h).
  - **Tout-tordu à 1 triangle : impossible.** 15 classes, 84 certificats exacts en entiers, 937 s, 258 Mo. Revérifié
    par `verify_core_certs.py` (code séparé, lignes reconstruites depuis leurs étiquettes, arbres complets). La
    couverture des 41 790 étiquetées est vérifiée.
  - **Sans triangle** : les 14 classes ont un PL faisable à la racine, même avec (P2) et les types. Le branchement
    est estimé à des heures ; arrêt dans le plafond, décision requise.
  - Erreur relevée dans mon message à l'utilisateur, corrigée aussitôt : 84 certificats et non 83 ; 5 classes à
    branchement, avec 74 certificats, et non 4 classes avec 69.
- Contre-vérification indépendante de H-β (C2), par échantillon, avec `verify_hb_independent.py` : code propre, écrit
  d'après la seule description de C2.
  - Couverture : 552 classes régénérées (7 680 squelettes), identiques aux enregistrements, multiplicités égales.
  - 8 certificats tirés au hasard : 8/8 valides.
  - Mon nombre de colonnes est exactement celui de C2 sur chaque enregistrement testé.
  - Les 552 complets coûteraient environ 4,5 h : non lancé, accord de l'utilisateur requis.
- Décision de l'utilisateur après S17 : **maths d'abord** pour les 14 classes sans triangle (tests de 5 min au plus).
- S18, vérification indépendante demandée par C1 : certificat X19 (`c1_experimental/data/x19_cert_s1_c1.json`),
  méthode par types avec (P2) transposée au cœur de C1. **Il TIENT** (`scripts/verify_c1_x19.py`, rationnels exacts).
  - Spéciaux symétriques.
  - Les types refaits par énumération des vrais graphes d'étiquettes coïncident avec ceux de C1 (84/84 paires,
    4 968 colonnes).
  - Les 735 seconds membres, recalculés, n'ont aucun écart.
  - La valeur de Farkas vaut 5,9916 > 0, égale en exact à celle de C1.
  - Portée : une configuration (τ, σ). Le « 19/31 » reste un échantillon.
- S19 (maths d'abord, selon la décision de l'utilisateur). Ingrédient nouveau : l'identité de C2 au niveau des orbites.
  - Énoncé : Σ_Z q_Z(U) q_Z(W) + (q_U(W) + q_W(U))/2 = 12[U = W] + 4 − Sh(U,W), pour les 903 couples d'orbites.
    C'est Q_out² + Q_out = 12I + 4J − Sh, re-dérivé via A² = 12I + 2J − A et les chemins de longueur 2 par les
    lointains. Contrôle de somme : 6552 = 6552.
  - (P2) par paires en est la compression, ce qui explique pourquoi (P2) ne suffisait pas.
  - Le modèle de configurations raffiné affecte chaque voisin de valeur 1 à l'une des deux orbites de sa paire (un
    double touche les deux, le carré va à l'orbite D(Z)). Cela fait 352 800 colonnes pour (7,7,7)×120, traitées par
    génération de colonnes avec énumération exacte des répartitions.
  - Test : **PL infaisable à la racine** pour (7,7,7)×120 (convergé en 11 tours, 270 s). Premier signe que le
    tout-tordu sans triangle peut tomber sans branchement.
- S19, suite :
  - Cadre vérifié contre la note 10 de C2 : σ échange les deux étiquettes de chaque bloc et z ≁ σz. Donc Sh = Σ b_i b_iᵀ
    et le terme diagonal est nul, comme codé.
  - Ajouts à `orbit_p2_colgen.py` : vérification exacte entière sur toutes les colonnes-orbites (lignes P2o doublées),
    arrêt dès que le certificat tient, sauvegarde avec étiquettes.
  - Nouveau vérificateur indépendant `verify_orbit_certs.py`.
  - zero0 (4,4,4,4,5) : certificat exact au tour 8 (phase 1 = 34,30, non convergé : inutile), en 252 s pour 785 Mo de
    RAM au plus. Le vérificateur indépendant retrouve exactement les mêmes bornes (écart 26,3 à l'échelle 10⁹).
  - Incident de format : les bornes P2o étaient stockées non doublées. Elles ont été détectées par le vérificateur
    (qui ne les utilise pas), puis corrigées dans le générateur et dans le fichier.
  - Messages : à C2 (identité d'orbites, réponses à ses 3 questions sur les mélanges, transmission du repère ajusté de
    C1) et à C1 (même objet que W14 ; raffinement sans symétrie par cellules antipodales).
- S20 : C1 a appliqué les cellules antipodales (S19) à son cadre sans symétrie. Résultats : 10 survivants sur 11
  infaisables, et le témoin de Fano tendu aussi.
  - Sa demande : vérifier le certificat du témoin. `scripts/verify_c1_cells.py` le fait avec du code propre, sans
    import de C1.
  - Résultat : spéciaux identiques et symétriques ; 4 spéciaux distincts par paire ; types = énumération propre (0
    écart, 0 doublon, 8 génériques chacun) ; découpages complets ; 29 652 lignes reconstruites depuis leur nom,
    identiques à A et b ; Farkas reconstruit = 75 381 256 (échelle 10⁶), égal à sa valeur.
  - VALIDE, rapporté à C1. Durée 90 s, en parallèle du lot d'orbites (calcul léger).
- Lot des 13 classes lancé (accord de l'utilisateur, « Lancer les 13 classes ») : `scripts/run_orbit_batch.sh`, une
  classe à la fois. zero13 est CERTIFIÉ au tour 5 (106 s, écart 21,4).
- S19, contrôle sur un vrai graphe (`scripts/orbit_identity_real.py`, `logs/orbit_identity_real.log`). Graphe : BvLS
  = srg(243,22,1,2), x = 0, σ(v) = −v, même cadre que C2 (un seul point fixe, les deux étiquettes de chaque bloc
  échangées).
  - Q² + Q = 20I + 4J − Sh est **exacte**, sur 110 orbites.
  - Q est symétrique, de diagonale nulle, avec exactement un 2 par ligne (partenaire D) et les autres entrées ≤ 1.
  - Spectre : 20, −2 (×10), 4 (×55), −5 (×44), comme prévu (x² + x = 20 hors de l'espace des blocs).
- 26/09 :
  - **Décision de l'utilisateur** : il transmet lui-même à Astra un résumé de nos résultats (C1, C2, structures).
    L'indépendance d'Astra prend fin par son choix. Je n'ai rien envoyé à Astra moi-même.
  - Partage des ressources avec C1, à la demande de l'utilisateur transmise par C1 (même machine, Core Ultra 5 125H,
    15,4 Go de RAM) :
    - 3,5 Go au plus chacun en simultané, tout sous run_capped ;
    - 2 fils au plus par tâche ;
    - un message avant tout lancement de plus de 1,5 Go.
  - Répartition : je vérifie au fil de l'eau les certificats « cellules » de C1 (format x19k). Mon lot d'orbites tourne
    sur 1 cœur, avec moins de 0,8 Go.
- **S19, fin du lot : TOUT-TORDU PUR ENTIÈREMENT EXCLU.**
  - Les 14 classes sans triangle sont certifiées à la racine, au tour 5 à 8, en 106 à 282 s par classe, avec moins de
    0,8 Go de RAM.
  - Écart low − high : de 0,36 (zero4) à 31,4 (zero11), en unités de y. L'arrêt anticipé s'arrête dès que le
    certificat tient, d'où certains petits écarts, tous exacts et positifs.
  - `verify_orbit_certs.py coverage` : 14/14 VALIDES, low et high identiques au générateur (vérifié classe par classe),
    couverture 30 360/30 360, ALL VERIFIED.
  - Bilan : ≥ 2 triangles à la main (S8, S11) ; 1 triangle, 84 certificats (S17, 41 790 décompositions) ; 0 triangle,
    14 certificats (S19, 30 360 décompositions).
  - Vérification des certificats cellules de C1 : lot 1 (s11_1, s21_0) VALIDE ; lot 2 en cours.
  - Deux incidents de lancement, corrigés, qui n'ont rien vérifié : `\$b` dans le script, puis la conversion de
    chemins MSYS (corrigée par MSYS_NO_PATHCONV=1).
- Certificats « cellules » de C1 : lot 2 (s11_2, s11_3, s21_1) VALIDE. Cela fait 5/5, et chaque valeur de Farkas
  reconstruite est égale à celle de C1 (la plus serrée : 5 286 630 / 10⁶). Moins de 0,1 Go et 140 à 180 s par
  certificat.
- Relecture du modèle paire × bloc de C1 (x19n) : lignes saines (labdeg, pairblk, pairpair, cons ; produits entre
  blocs omis). Le vérificateur `scripts/verify_c1_pairblock.py` est préparé, sur la convention de noms proposée à C1 ;
  il sera testé dès l'arrivée des fichiers.
- **Contre-vérification de C2** (son code, carnet §32) :
  - 14/14 certificats d'orbites VALIDES, low et high identiques ;
  - couverture 30 360 par réénumération propre et Σ 5040/|Aut| ;
  - identité redérivée.
  Le tout-tordu pur est exclu avec double vérification indépendante (mon vérificateur séparé + le code de C2).
- Certificats de C1, vérifiés par du code propre :
  - « cellules », lot 3 (s11_4 à s11_6, s21_2 à s21_6) : VALIDE, avec les valeurs de Farkas de C1 à l'identique. Cela
    fait **13/13** certificats cellules VALIDES (moins de 0,1 Go et 70 à 180 s chacun, `logs/verify_c1_cells_batch.log`) ;
  - **paire × bloc, seed21#1** (`scripts/verify_c1_pairblock.py`) : VALIDE. Les 78 492 lignes des 8 familles sont
    reconstruites depuis leur nom sans aucun écart ; les 124 477 colonnes s sont complètes (toutes les parties de
    taille t[R]) ; le Farkas reconstruit vaut 147 842 015 (échelle 10⁶), égal à celui de C1 ; 147 s, 177 Mo
    (`logs/vpb_seed21_1.log.out`).
- Obstruction de rang de C1 (cellcross au Fano tendu), **recalculée avec mon propre code**
  (`scripts/fano_cellcross_rank.py`, `logs/fano_cellcross_rank.log`), sans lire le script x25 de C1 :
  - données : droites de Fano lues sur τ (un plan de Fano, identique pour les deux étiquettes de chaque classe) ;
    spéciaux depuis τ et σ ;
  - comptes par bloc c_R(z) : résolus par les comptes de familles tendues 2/1/3 et les équations de fente 4. Solution
    unique, entière et ≥ 0 pour les 84 paires. Contrôles : Σ = 12, 0 arête interne, équilibre des blocs ;
  - systèmes cellcross découplés par bloc, 24 inconnues : rang(A) = 16 pour les 21 blocs, et rang[A|b] = 17
    exactement pour 02, 03, 05, 23, 24, 25, 34, 35, 36. **Identique à C1.**
  - Remarque : les 6 blocs qui contiennent la classe 1 sont tous cohérents.
- Autres vérifications : paire × bloc seed11#0 VALIDE (78 879 lignes, 0 écart, Farkas 168 240 155).
- Certificats « cellules » de C1, lot 4 (s11_7 à s11_9, s21_7 à s21_9) : VALIDE. **Total : 19/19 certificats cellules
  VALIDES**, soit les 19 configurations tuées au niveau cellules (graines 11 et 21). S'y ajoutent **2/2 certificats
  paire × bloc VALIDES** (seed21#1, seed11#0). Chaque valeur de Farkas reconstruite est égale à celle de C1.
- Noyau universel de C1, recalculé (paramétrisation des atterrissages permis de C1 : a ∈ C4, b ∈ {Φ(P,a), Φ(Q,a)}) :
  - les 168 vecteurs de comptes (21 paires de blocs × 8 atterrissages) sont bien définis (solution unique, entière,
    ≥ 0), et chaque vecteur du témoin est l'un d'eux ;
  - le noyau universel est de dimension 3 pour les 21 blocs ;
  - le membre de droite du témoin y est orthogonal partout.
  **Identique à C1.** L'obstruction vit donc dans les 2 dimensions propres au témoin (noyau à gauche de dimension 5).
- **Fonctionnelle Ψ_R de C1** (« piste à vérifier »), vérifiée avec mon propre code (section 6 du script) :
  - noyau effectif (atterrissages avec c_R > 0) de dimension 5 pour les 21 blocs ;
  - ⟨c(z;a,b), Ψ_R⟩ = 0 dès que c_R(z;a,b) > 0, sur les 21 blocs × tous les atterrissages ; ±6 sinon (336 + 336) ;
  - au témoin, ⟨rhs_R, Ψ_R⟩ ≠ 0 exactement sur les 9 blocs, avec les valeurs ±24 et −48. I_R est identique à la
    liste de C1, signes compris (même convention de signe) ;
  - forme close −6 Σ_{y∈R} ε(y)[θ(z_P(y)) + θ(z_Q(y))] égale à ⟨rhs_R, Ψ_R⟩ sur les 21 blocs ; tous les
    atterrissages du témoin sont permis.
  - Conséquence : **toute configuration de Fano tendu doit vérifier I_R = 0 pour les 21 blocs** (condition
    nécessaire exacte, sans PL). Le lemme « aucun sommet de Fano tendu » revient à l'incompatibilité de ces 21
    conditions avec le modèle.
- 26/09, consigne de l'utilisateur relayée par C1 : **arrêt des solveurs, retour aux mathématiques pures** (C1 et C2 ;
  j'applique aussi, rien ne tourne chez moi).
- C1 : le Ψ évitable par les retournements (480 réglages avec le σ du témoin satisfont les 21 conditions I_R = 0) ; la
  relaxation locale par bloc ne tue que les 9 blocs du témoin ; x22 indéterminé. La tâche « lemme Fano » proposée pour
  Astra est donc caduque (réponse négative), ce que j'ai signalé à l'utilisateur.
- S21 (maths pures) : preuve à la main de l'item 1 de C1 (NOTE_S21). Forme close des comptes génériques au Fano tendu
  (contrôle arithmétique 168/168), puis trois positions relatives.
- Idée « commutant » testée pour une tâche à Astra : dim{w : [A, L_w] = 0} = nombre de composantes de K sur Paley(9) (6)
  et BvLS (11). Rien de plus que les indicatrices de composantes ; piste close.
- Relecture du **lemme P-Ψ de C1** (condition θ(z_X(y)) + θ(z_Y(y)) = 0 par paire, au Fano tendu) : **TIENT**. Toutes
  les étapes sont refaites à la main avec le lemme 1 de NOTE_S21 (§4 ajouté).
  - σ-spécial dans le bloc {X,A}, cas b = W :
    somme/σ = −1 + 9 + (0 − 9 + 0) + (0·9 + 1·(−1) + 2·(−9) + 2·1) = −18.
  - σ-spécial dans le bloc {X,A}, cas b = π_Y(A) :
    somme/σ = −1 + 0 + (−1 − 9 + 0) + (9 + 0 − 18 + 2) = −18.
- C2 demande un lemme pour les carrés « dis ». Je lui ai transmis :
  - un dictionnaire par les spéciaux (dis ⇔ le partenaire n'est pas spécial ; aucune arête τ dans le carré) ;
  - une dichotomie exacte par paire, tirée de l'identité des sœurs : n_s = 0 ⇒ soit 2 doubles (1,1), soit 0 double
    et un seul (2,1).
  Je n'ai pas de borne sur n_dis.
- Question de C2 : pour un carré « dis », φ_u(R_D) vaut-il 2 ou 3 (profil (2,0) ou (2,1)) ? Son lemme 11.2 (= 2) n'a pas
  de preuve écrite.
  - Exemple réel : BvLS avec σ = −v n'a que des carrés « frère » (110/110). Aucun cas dis à observer
    (`scripts/dis_squares_real.py`, exploration de quelques secondes).
  - Les identités par paire (Σφ = 24, Σφ² = 32 − 2n_s) admettent les deux cas : (2,1) donne 9 + 4 + 19 = 32, (2,0)
    donne 4 + 4 + 8 + 16 = 32.
  - Dans le cas (2,1) : les 4 racines de R forment un stable ; voisins communs (w,w′) = {(k,γ), u} ; U est l'unique
    orbite qui voit W et W′ ; D(W′) ne voit pas W ; u a deux voisins supplémentaires v₁ ≠ v₂ hors de R.
  - Non tranché. J'attends de C2 les énoncés exacts de ses lemmes d'apex et d'équilibre (note 39).
- **Fano tendu, bilan (26/09)** : P-Ψ est un **filtre exact**, pas un lemme.
  - Sur 9 σ solutions de KBLDP, environ 40 % des droites sont vivantes. Estimation de C1 : environ 0,4⁷ ≈ 0,2 % des
    solutions de KBLDP passeraient les 7 droites, et elles sont très nombreuses. P-Ψ + KBLDP n'est donc
    probablement pas une obstruction complète.
  - Trancher demanderait un PLNE (KBLDP + retournements + P-Ψ), ce que l'utilisateur ne veut pas pour l'instant.
  - Vérifié à la main : le graphe τ des 12 paires d'une droite est 2-régulier, sans triangle (λ = 1), de longueurs de
    cycle multiples de 3, donc C6 + C6 ou C12.
  - Les deux sessions reviennent au cas général.
- **Nouvelle piste de C1, cas général** (dimension forcée), vérifiée à la main :
  - Sur K = ker ∂ (rang 70), B a pour valeurs propres 3 et −4, de multiplicités 40 et 30. En effet tr(B|_K) = 0, car
    BN = N(J − P − I) et ce bloc a pour spectre 12, 0 (×7) et −2 (×6).
  - K_cell (fonctions constantes sur les 42 cellules antipodales, dans ker ∂) est de rang 35 : 14 dans E₋₂ et 21
    antisymétriques.
  - D'où dim(K_cell ∩ V₃) ≥ 5 sur Q, et ≥ 7 modulo 7 (avec k ≤ 28 de C1). En tout sommet de toute srg(99,14,1,2),
    au moins 5 fonctions constantes sur les cellules sont vecteurs propres de B pour la valeur 3. C'est une
    contrainte bilinéaire (comptes × F), invisible aux PL par types.
  - Pour le test (b) au Fano tendu : M(δ) = [M_β | M_γ(δ)], de taille 84 × 35, doit vérifier rang ≤ 30.
- Test (b) de C1 au Fano tendu : **rang M_β = 8/14 sur Q et sur F₇** (`scripts/fano_mbeta_rank.py`). Il existe 6
  directions β ∈ E₋₂ constantes par paire, avec γ = 0, qui sont déjà vecteurs propres de B pour la valeur 3. La
  dimension forcée ≥ 5 n'impose donc rien aux δ au témoin : la piste est vide ici, sauf borne forcée > 6.
- Triades de C1 : les types (i), (iv), (m) et R1 sont vérifiés (multigraphe 2-régulier sur les 4 racines d'un
  bloc ; une cerise ne peut pas être doublée).
- Bijection λ : l'équilibre étiquette (w, Z) contient le compte T(w), sans le forcer à 1. La question reste ouverte.
- Fano tendu : T ≡ 1 est faux au témoin, dont le σ est solution de KBLDP. Sur 42 paires ordonnées, seules 14 ont le
  profil (1,1,1,1) (`scripts/fano_T_check.py`). La route du dérangement est close, et C1 met le Fano tendu en pause.
- Critère sans PL proposé à C2 pour les feuilles à bloc A (note 38 de C2).
  - Hypothèses : Σ = ∅, partenaires un bloc, a_Z ≤ 2.
  - La diagonale de Gram (Σa² = 24, Σa = 20) donne exactement 2 orbites hors A avec a_Z(x) = 2. Ce sont les
    partenaires D hors A (o_x) et les doubles hors A.
  - Avec la dichotomie (doubles totaux = 2 − n_s) : **#doubles dedans sur {A,x} = o_x − n_s**, et o_x ≥ n_s.
  - Extension à Σ ≠ ∅ : #dedans = o_x + D_xx − n_s, à reprendre dans le cadre de C2.
  - Information d'orbite, invisible au PL par paires. En attente du test de C2.
- Nouvelle porte de C1 (cas général, P = 0) : la géométrie interne I_x (140 droites internes) et la torsion aux
  voisins.
  - Relu :
    - t₁ + t₂ + t₃ + t₄ = 140 et t₂ + 2t₃ + 3t₄ = 84, car toutes les arêtes σ sont internes à P = 0 ;
    - t₃ = #{z : les deux σ-partenaires de z sont adjacents} − 3t₄, donc un degré de liberté de premier ordre.
  - Pont : une boucle τ^u_y équivaut à un prisme, qui équivaut à un triangle du graphe des carrés K d'Astra. Les
    conditions « sans boucle aux voisins » de C1 sont donc les instances locales de « K sans triangle » (P = 0), et
    Astra utilise P = 0 globalement.
  - Pas de parité locale gratuite en y fixé : l'union de deux couplages parfaits est faite de cycles pairs, et « sans
    boucle » exclut les digones.
- C2 : mon critère « dedans » est **vide** pour un bloc A complet. Tout est ancré, donc o_x = 0 ; par F1, les orbites de
  A ne sont jamais des doubles ; la diagonale PPE donne déjà 2 doubles hors A. C'est accepté.
- **C2 signale 14 survivants « faisables, A fixé »** au niveau paires :
  - M1 n = 6, jobs 0 et 2 (2 + 8) ;
  - mélange 2C₃ (4).
  Tous sont du cœur (G_c connexe partout, équilibre vide) ; fichier `c2_experimental/logs/lp_survivors.json`.
  - C'est le profil où le niveau orbite a mordu dans le tout-tordu.
  - Réponse à C2 : colonnes et lignes libres de `orbit_p2_colgen.py` propres au tout-tordu ; lignes P2o/C/PPE,
    pricing, exact_check et vérificateur génériques.
  - Interface donnée : (zi, fixed, v0, opts). Deux adaptations : voisin libre ambigu, et option (2,1) des partenaires
    dis. Je ne lance aucun PL.
- **Astra, nouveaux résultats** (Fano tendu, transmis par l'utilisateur) :
  - lemme « aucun revêtement régulier de degré 4 de K₇ n'a λ = 1 » : **relu, juste** (argument θ mod 2, puis tiroirs
    dans une classe d'ordre 2) ;
  - contrainte τ–σ « β_j = β_k = 1 ⟹ e_i = 0 » (contrôle d'Astra sur 3 072 cas) ;
  - profils d'atterrissage liés aux 24 couplages parfaits de Heawood (conséquence de (P2)).
  - Transmis à C1. Dossier n° 2 envoyé à l'utilisateur pour Astra : tâche 1 (revêtements non réguliers, S₄), tâche 2
    (la contrainte τ–σ tue-t-elle les 480 réglages ?), les questions de C1, et les données du témoin (τ, σ, les 480
    réglages).
- C1 :
  - identité automatique des branches ; e_i d'Astra = E(R_i) mod 2 (la sommation est vérifiée) ;
  - « où P-Ψ est automatique », **vérifié à la main** : m_B = 1 ⇒ π_B unique (6 placements), #MM = #CC + 6 − 2Σm,
    et dans le cas tout-CM τ est déterminé par σ.
  - Addendum envoyé à l'utilisateur pour Astra.
- **C2 : le survivant 0 (M1 n = 6, job 0) MEURT au niveau orbite**, avec un certificat exact de C2. Son propre code
  donne un écart de 77,09 sur 414 260 colonnes-orbites.
  - Contrôles de C2 : faisable sans les identités, et faisable avec les identités comprimées par paires. La mort vient
    donc de la partie orbite.
  - J'ai relu son expansion : (a) exclusions des voisins libres (F1, pivot j, orbites distinctes), (b) aucune paire
    lointaine ne contient D(Z), (c) Sh = 2 pour les sœurs. Tout est juste.
- **Lemme des différences de retournement de C1 : VÉRIFIÉ à la main.**
  - L'application de retour de T sur les 4 paires de AB est affine sur F₂², de partie linéaire M(d).
  - Poids 1 : M d'ordre 3 (x² + x + 1), un point fixe, donc un τ-triangle : inadmissible.
  - Poids 3 : M = I, donc C6+C6. Poids 0 ou 2 : transvection, donc C12.
  - Parties de Fano qui rencontrent chaque droite en 0, 2 ou 3 points : ∅, le complémentaire d'une droite, le
    complémentaire d'un point, ou le plan entier.
  - Addendum 2c envoyé à l'utilisateur pour Astra.
- **Rapport d'Astra sur le dossier n° 2** (transmis par l'utilisateur ; Astra n'avait pas reçu les addendums 2b et 2c) :
  - **Revêtements de degré 4 de K₇ avec λ = 1 : un seul**, non régulier, de monodromie S₄, avec exactement 8 prismes.
    L'unicité vient de deux énumérations exactes indépendantes d'Astra ; la construction et le comptage des prismes
    sont prouvés à la main, sur F₄.
  - Conséquence : la branche « factorisation de Heawood » du Fano tendu force P ≥ 24, donc elle est exclue à P = 0.
  - **Le σ témoin de C1 meurt par P-Ψ pour tout τ.** Certificat sur la paire {0.0, 4.0} : σ₅^{±}(0.0) ∈ classe 2 et
    σ₅^{±}(4.0) ∈ classe 3, sur la même droite (2,3,5). **Revérifié sur les données.**
  - La contrainte τ–σ d'Astra est redondante (conséquence de P-Ψ, de la réciprocité de σ et de la symétrie des
    comptes).
  - Réponses aux questions de C1 (β, U, d_i, r_i) transmises à C1. Astra retrouve indépendamment la typologie de C1
    par l'application de retour affine.
- **Dossier n° 3 pour Astra**, préparé avec C1 à la demande de l'utilisateur.
  - Partage : Astra prend le Fano tendu (tâche principale : le cas général ; tâche secondaire : une généralisation en
    contrainte). C1 prend le cas général (I_x à P = 0). « structures » vérifie et fait les ponts.
  - Contenu :
    - Q1 à Q4 de C1 (Q2, décisive : F5 + KBLDP + P-Ψ + admissibilité sur 7 droites a-t-il une solution ?) ;
    - le découpage des branches ;
    - deux outils de C1 **revérifiés à la main** : les coordonnées φ (C4(L) ≅ vecteurs pairs, F5 et P-Ψ par droite),
      et (A + I)g_R = 12ψ_R − 30ζ_R (deg_Ψ(c) = 15s(c)) ;
    - les 9 σ avec leur statut par droite (`scripts/fano_perline_status.py`) ;
    - les notes 120 à 122.
  - Statut par droite **recoupé** : pour le témoin, les comptes P-Ψ sont 0,0,2,1,0,0,1, comme Astra ; pour s41_5, 5
    droites P-Ψ-vivantes (mortes : 013 et 602), comme C1, et 4 avec l'admissibilité.
- **C2, niveau orbite**, avec son propre code et un second vérificateur indépendant de son côté :
  - 14/14 survivants du niveau paires tués (écarts de 73 à 97) ;
  - la feuille M1 n = 6 job 0 n° 17 tuée entière, A variable (écart 27,9, un PL) ;
  - l'orientation M1 n = 6 job 0 tuée entière, toutes structures D libres (écart 5,03, 2,2 M colonnes-orbites,
    862 s ; le niveau paires y avait mis 358 PL).
  - « Tout libre » : faisable, comme attendu. En cours : T fixé, orientation libre, sur l'hexagone.
  - Note : ce sont de nouveaux PL de C2, sur la consigne qu'il a reçue directement de l'utilisateur.
- Argument de C2 (« le voisin libre de pivot j est exclu par F1 + couplage ») : juste pour les solutions entières.
  En fractionnaire, le couplage se compense en moyenne (poids ½ sur des colonnes où D(F) = Z). La relaxation est donc
  plus faible mais valide. Règle explicite suggérée si un PL reste faisable.
- **Rapport d'Astra sur le dossier n° 3** (transmis à C1) :
  - **Q2 : oui.** Témoin explicite de rang 8 (σ_c(d_t) = ρ^t_c(d)_{1−t}, avec g = (0 1)(4 6), tous les retournements
    à 1). Il satisfait F5, P-Ψ, l'admissibilité (7 × C12) et toutes les identités de bloc (P2 prouvée par
    Q² + Q = 12I + 8J − 2DᵀD). Donc « P-Ψ + admissibilité » n'est **pas** un lemme de Fano complet.
  - Ce témoin, et toute sa famille (deux factorisations par parité), meurt par un certificat de bornes ordonnées au
    niveau paire × bloc : 6 > 5 en y = {0.0, 2.0}, R = {0,1}, pour tout τ admissible.
  - L'outil général des bornes ordonnées (pairblk + choix de g_y(B) sommets distincts) est **vérifié**.
  - **Théorème (vérifié à la main)** : un revêtement de K₇ avec λ = 1, μ ≤ 2 et sans prisme a un degré ≥ 5. Le degré 3
    est exclu parce que la multiplicité vaudrait 42/5. À P = 0, le graphe auxiliaire H du Fano tendu a des
    composantes de degrés 12, 5 + 7 ou 6 + 6.
  - Ouvert : un Fano tendu qui satisfait aussi toutes les bornes ordonnées.
- **Décision de l'utilisateur** (AskUserQuestion, 26/09) : arrêter l'ancien pipeline par paires de C2 et lancer sa
  campagne au niveau orbite sur M1 n = 5–6 et les mélanges (75 à 125 h CPU, sur la grosse machine de C2). Relayée à
  C2.
  - État de C2 (son message) : sont morts les ancrés, le tout-tordu, H / H-β, M1 n = 3 et 4, et l'orientation M1
    n = 6 job 0.
  - Restent en calcul M1 n = 5 et 6, et les 47 mélanges. Restent hors calcul les carrés disjoints et les 9 motifs
    Σ ≠ ∅ : c'est le vrai verrou.
- C1 (I_x) : holonomie (prismes sur T = points fixes de σ_T), H + 3P = 213 444, graphe signé des triangles
  A_Σ² = 13A_Σ + 68I (spectre 17 ×44 et −4 ×187, cohérent), ombres. Toutes ces identités sont automatiques : le
  plafond doit venir des génériques effectifs. C1 cherche des strates rigides forcées à P = 0.
- Leçon (mémoire) : pour savoir ce que fait une session, le lui demander ; les commits sont en retard sur le travail.
- **Correction des crédits (C1)** : ses identités sur I_x étaient déjà dans le corpus.
  - Graphe signé des triangles : notes/99 §4 et 109 §4. Ombres : notes/114 §1. Holonomie et prismes : F1, notes/99
    §4 bis. H + 3P = 213 444 : théorème de Reimbayev.
  - Seul point nouveau : à P = 0 et S = 0, H(x) = 12 936 en tout sommet. C'est automatique.
  - C1 a réécrit notes/123 avec les bons crédits.
- **(2,0) contre (2,1) pour un carré dis** (partage avec C2, ma part), premier passage à la main :
  - Rigidité aux blocs k et l : hormis w, σw et w′, seuls v₁ ∋ (k,−γ) et v₂ ∋ (l,δ) touchent k ou l. Tous les autres
    voisins racines de u évitent k et l.
  - Σ c_v(R) = 5 se ferme automatiquement dans les deux sous-cas : t ∈ D(W′), ou t, y₁, y₂ distincts.
  - Globalement : Σ doubles = 42 − Σn_s − 2N₂₁, cohérent avec la dichotomie.
  - Pas de contradiction. Proposé : classification locale de Q ∪ W′, à confier éventuellement à Astra.
- Campagne orbite de C2, premières réponses : beaucoup d'orientations restent faisables avec D entièrement libre
  (classe 2 job 6, papillon job 10, M1 n = 5 job 10, M1 n = 6 jobs 11 et 12). C2 compare deux découpages sur M1 n = 6
  job 12 : branchement sur le partenaire D de l'orbite la plus étalée, contre union par feuilles.
- Lignes proposées à C2 pour resserrer la racine :
  - équilibre étiquette au niveau orbite : Σ_{Z : c ∈ paire(Z)} q_Z(W) = 4 (210 lignes) ;
  - unicité du partenaire D en ligne : la masse des colonnes qui mettent 2 sur W vaut 1 (42 lignes) ;
  - la règle explicite « voisin libre en j de pivot ≠ j ».
  - C2 intègre aussi mes filtres de colonnes dis (rigidité aux blocs k et l).
- **Dossier n° 4 pour Astra**. L'utilisateur m'a laissé le choix (« comme tu veux ») : j'ai pris les **carrés disjoints**
  de C2, vrai verrou de la branche involution, hors de portée du calcul.
  - Tâche 1 : (2,0) contre (2,1), soit le lemme 11.2, soit une configuration locale cohérente.
  - Tâche 2 : un lemme structurel pour la famille dis (Σ = ∅).
  - Tâche 3, facultative : les 9 motifs Σ ≠ ∅.
  - Joints : les notes 00, 07, 19, 35 et 39 de C2, la dichotomie, et mon analyse de (2,1). C2 est prévenu.
- **Dossier n° 5 pour Astra** : la cible 7-adique de C1 sur le cas général (L22).
  - R_L, matrice 60 × 60 mod 7 de rang k − 22, nulle à k = 22.
  - Question : existe-t-il des ombres β_u, β_ū et un graphe lointain 8-régulier de rang ≤ 6 (ou 0) ? Un « non »
    uniforme donnerait k ≥ 30, contre L22.
  - Cas dégénéré à exclure à part : les 231 droites 7-dégénérées, forcé si k ≤ 20.
  - Joints : notes 99–102, 106, 107, 110, 114 et 123 de C1. À traiter après le dossier n° 4.
- Bilan I_x de C1 : la fermeture de Δ_x est essentielle (des géométries tordues au maximum existent hors srg, avec
  P(x) = 0) ; pas de certificat linéaire uniforme (relaxations faisables jusqu'à l'ordre 8) ; il faut du non convexe.
  C'est cohérent avec ma réponse à l'utilisateur sur la méthode.
- Estimation donnée à l'utilisateur : C2 environ 60 % (reste le verrou dis / Σ ≠ ∅) ; C1 moins de 10 % (une idée
  nouvelle manque probablement).
- **Rapport d'Astra sur le dossier n° 4 (carrés disjoints)**, transmis à C2 :
  - Outil : 𝒫 = 6I − 2Q − Sh + J avec **𝒫² = 14𝒫**, **revérifié à la main**. Il vient de Stage Q : Sh² = 10Sh + 8J,
    Q·Sh = 8J − 2Sh, QJ = 12J, ShJ = 24J. Donc 𝒫 ⪰ 0, et ses sous-matrices principales sont des filtres exacts.
  - Témoin local (2,1) sur 5 lignes d'orbites, vérifié par deux implémentations d'Astra : le lemme 11.2 ne peut pas se
    démontrer localement.
  - Filtres démontrés :
    - (2,1) ⇒ aucun apex dans T (noyau (5,5,2,4) de K, revérifié) ;
    - un apex intérieur force l'autre à l'extérieur (résidu de Schur 1, revérifié) ;
    - le partenaire d'un apex intérieur a 5 supports possibles.
  - Borne : **3h ≤ 2n_dis**.
  - Le lemme 11.2 reste ouvert, et la famille dis n'est pas fermée.
- Suivi du dossier n° 5 envoyé à l'utilisateur (il avait déjà transmis la v1). Contenu : la vérification de Schur par C1
  (6 tirages sur 6, rang Q_BB = 37 − d(L)) et la correction de l'heuristique (R_L = KΞKᵀ, Ξ symétrique 26 × 26, soit
  environ 7^(−210) ≈ 10^(−177) et non 7^(−2916) ; Ξ = 0 de codimension 351, arithmétique revérifiée).
- Projecteurs : P₋₄ = Π_K(3I − B)/7, avec Π_K = I − ∂ᵀ(∂∂ᵀ)⁺∂ qui commute avec B, est de degré 1. Ses mineurs principaux
  sont le lemme du triangle de C1 (notes/115), déjà exploité. Ma remarque « degré plus haut » était fausse avec Π_K
  explicite. Le §3 du suivi était imprécis ; signalé à l'utilisateur.
- Règle de l'utilisateur : une fois un dossier transmis à Astra, ne lui renvoyer que les nouveautés (fichier de suivi),
  jamais une version refaite du dossier.
- Suivi n° 2 du dossier n° 5 pour Astra : la version réelle du complément de Schur lointain de C1 (notes/124).
  - Énoncé : Φ_F = −7Π, avec Π le projecteur sur K₋₄(L), de dimension 7 + c.
  - Vérifié à la main : sur K₋₄(L), A_BF f_F = 0, d'où Φ_F f_F = −7 f_F ; les dimensions se referment à 60.
  - Conséquence : A_F = 3I + NᵀRN − 7Π, donc (3I + NᵀRN − A_F)/7 est un projecteur de rang 7 + c. Ses mineurs principaux
    donnent un filtre réel sur les lointains.
- C2 : **𝒫(Q̄) ⪰ 0 est impliqué par le PL au niveau orbite**, et je l'ai vérifié.
  - M₂ + Q̄ = R et M₂ − Q̄² = Σ Cov ⪰ 0, donc Q̄² + Q̄ ⪯ R.
  - Q̄ préserve span(Nb), et spec(Q̄) ⊂ [−4, 3] sur le complément.
  - L'outil d'Astra ne vaut donc que comme filtre exact sur des squelettes, et C2 garde les filtres dis.
  - C2 teste un modèle au niveau racine (signé, S² + S = 12I − CCᵀ) sur M1 n = 6 job 12, avec environ 11 M colonnes
    implicites.
- **Témoin local (2,1) d'Astra (dossier n° 4) : VÉRIFIÉ indépendamment** (`scripts/verify_astra_21_witness.py`, avec
  son propre code, sans solveur).
  - 10 racines centrales de degré 14, adjacences symétriques, λ/μ exacts entre elles et avec x et les 14 intérieurs.
  - 0 dépassement de budget sur les 4 851 paires.
  - Lignes d'orbites avec un seul 2 et une somme de 12, profil (2,1) : u ~ w, σw, w′ et u ≁ σw′.
  - Identité d'orbites exacte sur le bloc central 5 × 5.
  - Conséquence : le lemme 11.2 ne se démontre pas à ce niveau local.

- **Énoncé autonome du dossier n° 4 pour C1** (NOTE_S23), à la demande de C1, qui ne connaît pas le cadre de C2.
  - Contenu : définitions de C2 avec leurs preuves courtes (racines, orbites, Q, D, π_i, F1, identité d'orbites, budget
    de Gram, apex) ; lemme « A est dis et Q_UV = Q_WV = 0 » ; théorème 3h ≤ 2n_dis ; 𝒫² = 14𝒫 ; réduction à 5 supports.
  - **3h ≤ 2n_dis : refait à la main, TIENT.** J'ai simplifié le point « A est de type dis » : si i ∈ supp D(A), l'un
    de v, σv contient (i,α), donc vaut π_i(a) = u, contre Q_AU = 1.
  - **6.1 à 6.3 : refaits à la main, TIENNENT.** On utilise la forme close du complément de Schur sur (B,V) après
    K = 𝒫[U,W,A], avec det K = 16 et K⁻¹ = (1/16)[[21,19,16],[19,21,16],[16,16,16]]. Les tableaux sont recoupés en
    fractions exactes (`scripts/check_dossier4_compressions.py`, `logs/check_dossier4_compressions.log`) : 0 écart.
  - Restent exactement {p,q}, {p,r}, {q,r}, {j,r}, {l,r}, avec Q_BV = 0 forcé pour les deux derniers, comme chez Astra.
- **C1 relit NOTE_S23 en troisième regard : tout tient, sans écart.** Il a tout refait à la main et recalculé le §5 avec
  son propre script en fractions. Sa remarque : la clause « aucune orbite comptée n'est aussi apex choisi » s'obtient en
  appliquant « U n'est apex intérieur d'aucun carré » au carré (A, V) lui-même, quand il est compté. Cela s'applique
  tel quel.
- **Lemme (4,4,4) de C1 (G(L) à P = 0, S = 0) : TIENT.** Pour tout v, r_v est dans P_u ou dans P_ū, jamais un
  lointain. Refait à la main, avec les 8 configurations vérifiées une à une. Preuve en une ligne envoyée à C1 :
  - par (c), τ_ū = σM̄σ⁻¹, où M̄ est le couplage complémentaire de τ_u ;
  - r_v ∈ P_u ⟺ σ⁻¹(v) = τ_u(v), et r_v ∈ P_ū ⟺ σ⁻¹(v) = M̄(v) ;
  - σ⁻¹(v) est dans l'autre classe, qui vaut {τ_u(v), M̄(v)}.
  Les quatre r_v sont du même côté si et seulement si σ est une involution sur la partie.
- **Rapport d'Astra sur le dossier n° 5 (relayé par l'utilisateur) : entièrement correct** (NOTE_S24,
  `astra/RAPPORT_ASTRA_DOSSIER5.md`).
  - §1 : la boule tordue explicite est exclue à k = 22. Reproduit avec mon propre code
    (`scripts/verify_astra_d5_seven_shadows.py`) : τ = 1 928 719 172 295 ≡ 6, mêmes candidats, aucune complétion en
    27 626 nœuds. Le contrôle positif retrouve le témoin d'Astra.
  - §2 (R = 0 ⇒ SRG exacte), §3 (MR = 0 et rang ≤ 2 ⇒ SRG exacte) et §4 (un échange de ≤ 6 ombres force deux prismes) :
    refaits à la main.
  - Mon ajout : MR = 0 est automatique dans une vraie SRG. En effet Q² ≡ 0 mod 7, et la condition d'image de notes/102
    §8 s'applique. Donc k = 24 est aussi un problème de construction complet.
- **Modèle de colonnes au niveau racine de C2** (`root_colmodel.py`) : relecture en cours, à la demande de C2. Points
  prioritaires : l'équilibre 2/2 aux blocs lointains et le choix des racines m et π.
- **Relecture de `root_colmodel.py` de C2 (commit 8cb7fb0) : aucune erreur de modélisation.**
  - Équilibre 2/2 aux blocs lointains : TIENT. Les 4 racines proches, dont celles de D, sont bien retranchées avec
    leurs signes ; il y a exactement 8 lointaines, au plus une par orbite.
  - Choix de m et π : TIENT. Dans chaque orbite, le code prend la racine qui porte le bon intérieur signé, et
    (F_m, F_π) est ordonné. La règle typée « pas toutes deux dans {j,k} » est redondante avec l'équilibre. F_π peut
    être la sœur de F_m sur {j,c} avec c ≠ k : c'est légitime, car c'est un double (1,1).
  - Lignes symP, symM et I(U,W,±) : TIENNENT. C'est (B² + B) = 12I + 2J − CCᵀ lu dans les colonnes, avec un second
    membre doublé.
  - Tests exacts sur le job 0 de M1 n = 6 (le job 12 n'est pas dans le clone), sans PL, script
    `scripts/review_c2_root_colmodel_tests.py` :
    - T1 : 26 720 colonnes relues, 0 problème ;
    - T2 : far_sets identique à une énumération indépendante (44 clés) ;
    - T3 : entries(), 126/126 ;
    - T4 : scores() = Aᵀy, 96/96.
  - Plantage latent l. 74 : reshape(0, −1) sur un ensemble lointain vide.
  - Réserve : un EVEN = 1 infaisable ne prouverait pas un bug. En effet symP et symM restent séparées, et le 2/2 signé
    est plus fort que les sommes orbitales.
- **Suivi n° 3 du dossier n° 5 envoyé à l'utilisateur pour Astra** (`astra/SUIVI3_DOSSIER5.md`). C1 a choisi les tâches.
  - Contenu : le rapport d'Astra est vérifié, et MR = 0 est automatique, ce qui fait de k = 24 un problème de
    construction complet.
  - **Tâche principale (a)** : rendre uniforme l'obstruction de l'étoile à k = 22.
    - Motif exact de l'étoile de u (sans prisme par u) : 2 apex isolés, et un couplage parfait par a = 1 des 8 autres
      lointains, avec l'égalité Σ_{w≠v} a = 1 − (Bn_v)_u.
    - La somme en étoile est automatique mod 7 (C1, revérifié), donc l'obstruction est dans les valeurs par paire.
  - **Tâche secondaire (b)** : le cas où toutes les droites sont dégénérées. Fait de C1 : si G(L) a c composantes,
    alors rang_Q(A_B − 3I) = 38 − c et d(L) ≥ c. Contrôlé en calcul exact (`scripts/check_c1_disconnected_degenerate.py`)
    sur la boule d'Astra (d = c = 1) et sur des boules (4,4,4), (6,6), (4,8) (d = c).
  - (c), la strate (4,4,4), reste chez C1.
- **C2 a appliqué mes deux correctifs** : far_sets vide rend zeros((0,8)), et symP et symM sont fusionnées en mode EVEN.
  - « Pas la sœur » vise la sœur de Z. F_π sœur de F_m sur {j,c} est voulu : c'est un double proche (1,1).
  - **Résultat négatif de C2** : le système A-local signé reste faisable pour les 4 types de D_A (6-cycle, 4+2, 3+3,
    2+2+2). Il n'y a donc pas de lemme purement local « pas de bloc A », même avec les signes.
  - **Correction de ce que j'avais dit à l'utilisateur.** L'objectif de phase 1 du PL racine (job 12) est passé de 18,4
    (itération 29) à 8,4 (itération 33). La génération de colonnes n'a pas convergé, et un objectif positif en cours de
    route ne prouve rien. Seule une convergence avec un objectif > 0 et un certificat vérifié serait une infaisabilité.
- **notes/125 de C1 (Σ_L c₄(G(L)) = S + 3P et Σ_L c₃ = 2P ; lemme de Klein) : TIENT.**
  - Classement des 4-cycles refait à la main, ainsi que la bijection (prisme, barreau) ↔ cycle et Σ_L S(L) = 2S.
  - **Contrôle positif sur les vrais membres de la famille** (`famille/control_c1_note125.py`) :
    - H(2,3) : 12 = 2P et 18 = 3P ;
    - BvLS : 17 820 = 2·8 910 et 26 730 = 3·8 910, avec S = 0 calculé par la définition.
  - Réponse à la question de C1 : ses échantillons « S = 0 » n'imposaient que σ ≠ v̄. On y compte 1 à 13 événements
    σ = τ_u et 2 à 12 événements σ = τ_ū par configuration (`scripts/count_S_events_c1_samples.py`). Seul le témoin de
    Fano est à S = 0 au sens complet.
- **Programme « famille » (demande de l'utilisateur : étudier à fond les graphes à 9 et 243 sommets, chercher une règle
  qui échoue à 99, voire un nouveau membre).** Démarré dans `famille/`.
  - H(2,3) et BvLS construits (Cayley sur F₃ⁿ ; code de Golay vérifié : dimension 6, distance 5).
  - Invariants de C1 et C2 calculés sur les deux graphes (`famille/family_invariants.py`) :
    - l'holonomie de chaque droite est triviale ;
    - P = P_max = L(m − 1), soit 6 et 8 910 ;
    - 100 % de points de Hamming ;
    - rang₃(A − rI) = 4 et 67 ;
    - forme de Smith 3-adique de A − 4I (BvLS) : 1^67 3^34 9^10 0^132.
  - Les deux graphes sont à l'extrême « torsion nulle ». Or tous les lemmes de C1 (L8, L16…) excluent justement ce
    côté à 99. Les graphes connus ne renseignent donc pas sur le régime fortement tordu, où vivrait un Γ.
  - Observation arithmétique : l'écart des valeurs propres vaut t = √(8m − 7), et t divise 63. Les membres qui
    existent ont t ∈ {3, 9}, les trois cas ouverts ont t ∈ {7, 21, 63}, c'est-à-dire exactement 7 | t.
  - Recherche bibliographique en cours (agent) : unicité du srg(243,22,1,2) et méthode de preuve.
- **Rapport d'Astra « étoiles d'ombres » (réponse au suivi n° 3) : correct** (NOTE_S26,
  `astra/RAPPORT_ASTRA_DOSSIER5_SUITE.md`).
  - §1 : la famille cyclique est exclue à k = 22. Code indépendant (`scripts/star_shadows.py`,
    `scripts/verify_astra_d5_cyclic.py`) : les 6 triplets génériques à τ non carré sont retrouvés, avec les mêmes
    candidats. L'étoile complète échoue sur les 6.
  - §2 : à k = 24 avec τ carré, |J_bad| ≥ 8. Refait à la main, y compris l'indépendance de 5 colonnes et le cas de 7.
    Le lien « τ carré ⇒ R anisotrope » est à confirmer par C1 (signes de L23).
  - §3 : 66 composantes 4+4+4 admissibles, τ ∈ {8 100, 8 640} (carrés non nuls mod 7), d(L) = 3. C'est compatible
    avec le lemme de Klein de C1 : elles ont des événements S.
- **Calibration heuristique** (4,5 min, ce n'est pas une preuve ; `scripts/random_balls_star_calibration.py`,
  `logs/random_balls_star_calibration.log`).
  - Sur 92 boules aléatoires, 35 sont admissibles (K ≥ 0), génériques et à τ non carré. **L'étoile complète en A₀ les
    tue toutes les 35.** Avec la boule d'origine d'Astra et les 6 cycliques, cela fait 42 sur 42.
  - Contrôles du module (`logs/star_shadows_controls.log`) : le témoin partiel v = 1..7 est retrouvé, le négatif
    (1..6, 8) aussi.
  - Suggère une obstruction universelle à k = 22 ; tâche (a′) proposée à C1 pour Astra.
- **Revue de littérature (agent) : l'unicité de srg(243,22,1,2) est OUVERTE** (`famille/REVUE_LITTERATURE_9_243.md`,
  NOTE_S25 §7).
  - La « règle » vraie pour 9 et 243 et fausse à 99 est la propriété de grille (*) (P = P_max). Makhnev 1988 (Th. 2)
    prouve qu'aucun srg(99,14,1,2) ne la vérifie, et L16 de C1 en est une autre preuve.
  - La pièce manquante est la conjecture de Reimbayev (tout membre vérifie (*)). L'unicité de 243 bute sur le même lemme.
- **Suivi n° 4 du dossier n° 5 envoyé à l'utilisateur pour Astra** (`astra/SUIVI4_DOSSIER5.md`). Contenu :
  - rapport vérifié ;
  - lien de signe démontré par C1 : χ(disc R_L) = (−1)^{k/2} χ(τ) ;
  - calibration 42/42 ;
  - tâche (a′) : trouver pourquoi l'étoile échoue toujours, avec les pistes α (noyau des 2 apex) et β (potentiel) de
    C1 ;
  - tâche (b′) : l'étoile avec j mauvais membres à k = 24.
- C2 demande l'énumération du transfert de coin (parité le long des traînées) : à faire.
- **Initiative « lemme unique par 243 » : réorientée.**
  - Le seul forçage connu de la propriété de grille serait arithmétique. À 243, c'est 3-adique avec 3 | v ; à 99, c'est
    7-adique avec 7 ∤ v. Les situations diffèrent (Q² ≢ 0 mod 3 à 243, alors que Q² ≡ 0 mod 7 à 99). Le local seul est
    exclu (faux mondes), et l'ordre ≤ 8 est faisable (notes/91). Donc aucun transfert direct.
  - Travail recentré sur le plafond 7-adique : l'obstruction de l'étoile à k = 22.
- **Ablation de l'étoile complète à k = 22** (`scripts/star_ablation_one.py`, `scripts/star_ablation_ball.py`,
  `logs/star_ablation_astra.log`, `logs/star_ablation_balls.log`).
  - Sur la boule d'Astra, l'étoile reste morte sans λ/μ contre la boule, sans les zéros de K et sans les
    recouvrements. Sans g(n,n) = 6, la recherche est trop grosse (> 110 s).
  - **Avec les seules contraintes « de Green »**, l'étoile est morte sur 4 boules sur 4 (Astra, cyclique (2,4,1),
    2 aléatoires), en environ 90 s chacune. Ces contraintes sont : forme des ombres, g(n,n) ≡ 6, g ∈ {2,3} par paire,
    et motif de l'étoile 4K₂ + 2K₁ avec capacités égales.
  - **L'obstruction est donc algébrique** : la forme de Green H⁺ mod 7 de la boule ne représente pas le Gram
    3I + 3J − (4K₂ + 2K₁) par des ombres de forme admise.
