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

