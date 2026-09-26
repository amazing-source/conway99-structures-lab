# Wagner : fermeture extérieure et calcul distant

26 septembre 2026. Recherche conduite seul. Les calculs de cette reprise sont exécutés sur SageMathCell ; aucun processus Python de calcul n'est lancé sur le PC de l'utilisateur. Les fichiers sont des sorties de cette conversation.

**Résultat nouveau, vérifié indépendamment : à P=0, tout support Wagner de norme carrée 8 utilise les sept blocs.** Les 25 classes à six blocs sont exclues par des certificats entiers. Les 120 classes à sept blocs restent ouvertes. Ce résultat ne clôt pas la branche involution.

## Cadre et portée

On suppose un srg(99,14,1,2), une involution fixant x, et P=0. On étudie seulement un vecteur b de L_T de norme carrée 8, dans le cas Wagner restant après l'exclusion du cube. Les réductions antérieures sont dans ../attaque_conjectures/RAPPORT.md.

Il y a huit représentants positifs p_i, leurs images q_i, et B=J-I-A(W_8), avec indices modulo 8. Ainsi B_ij=1 exactement lorsque i-j appartient à {2,3,5,6}. Le support induit a pour adjacency [[0,B],[B,0]]. Les p_i sont un stable, de même que les q_i.

Les paires d'indices dont les lignes de B ont un seul voisin commun sont les arêtes de H=C4+C4 :

```
(0,2), (0,6), (1,3), (1,7), (2,4), (3,5), (4,6), (5,7).
```

## Lemme des six ou sept blocs — démontré

Pour un bloc c, r_c(i) vaut +1, -1 ou 0 selon que p_i contient c+, c- ou aucun intérieur de ce bloc. Cb=0 donne autant de +1 que de -1. Soit h_c leur nombre commun.

Deux indices positifs dans cette ligne partagent déjà c+. Comme ils ont au plus deux voisins communs dans le graphe complet, leurs lignes de B doivent avoir au plus un voisin commun. Ils forment donc une clique de H. H est sans triangle : h_c<=2.

Si h_c=2, les deux indices positifs et les deux négatifs sont des arêtes disjointes E,F de H. Après une symétrie de Wagner, E={0,2}.

- F={4,6} donne un prisme induit : les triangles {c+,p_0,q_6} et {c-,q_2,p_4} ont exactement les trois barreaux c+c-, p_0q_2, q_6p_4. Il est exclu par P=0.
- Parmi les quatre arêtes de la composante impaire de H, F={3,5} ou {5,7} donne deux voisins communs sur une arête intérieur-racine, contre lambda=1.
- Restent F={1,3} et F={1,7}.

Ainsi chaque bloc à h_c=2 porte quatre indices consécutifs du cycle de Wagner, avec signes alternés, à un retournement de bloc près.

Deux tels blocs ne peuvent partager deux indices : sur leur intersection, les signes alternés ont un rapport constant. Deux intérieurs appropriés des blocs auraient alors x et au moins deux racines choisies comme voisins communs, contre mu=2.

Deux intervalles cycliques de quatre indices ont une intersection de taille au plus 1 seulement si leurs débuts sont à distance cyclique au moins 3. Trois débuts imposeraient trois écarts d'au moins 3 sur un cycle de longueur 8. Il y a donc au plus deux blocs à h_c=2.

Chaque p_i utilise deux blocs : somme h_c=8. Si a est le nombre de blocs actifs et l le nombre de blocs à h_c=2, alors a+l=8. Avec a<=7 et l<=2, on obtient exactement :

| Blocs actifs | Blocs à quatre représentants | Blocs à deux représentants | Blocs inutilisés |
|---|---:|---:|---:|
| 7 | 1 | 6 | 0 |
| 6 | 2 | 4 | 1 |

Ce lemme ferme uniformément les profils à cinq blocs ou moins. La fermeture des six blocs, donnée ci-dessous, laisse seulement sept blocs.

## Classification finie des noyaux

Un noyau comprend x, les quatorze intérieurs et les seize racines choisies : 31 sommets. On vérifie seulement les bornes locales lambda<=1, mu<=2 et l'absence de prisme induit.

Le calcul exact trouve 28 lignes de bloc à deux représentants et 8 lignes à quatre représentants. Il y a 1 968 ensembles de lignes admissibles, modulo permutation et retournement des blocs, avec les indices de Wagner fixés : 240 à six blocs et 1 728 à sept blocs. Une réduction supplémentaire par translations et réflexions de W_8 donne 145 représentants : 25 à six blocs, 120 à sept blocs. Ce sont des noyaux, pas des SRG.

## Comptes extérieurs : correction des unités

Pour un sommet extérieur z, t_z est son nombre de voisins positifs, égal à son nombre de voisins négatifs. Les 83 sommets extérieurs ont les nombres (n_2,n_1,n_0)=(8,64,11). En orbites, cela donne quatre doubles, trente-deux simples, cinq zéros et x.

Une fois les intérieurs retirés, les 34 orbites de racines extérieures ont les profils :

| Blocs actifs | t=2 | t=1 | t=0 |
|---|---:|---:|---:|
| 7 | 3 | 26 | 5 |
| 6 | 2 | 28 | 4 |

Le partenaire D(U_i) est une orbite extérieure dont les deux sommets voient à la fois p_i et q_i. Un simple voisin commun à deux orbites différentes n'est pas leur partenaire D. L'argument D du rapport Deep Research ne réfute donc pas à lui seul le noyau partiel.

## Premier noyau : certificat exact de non-prolongeabilité

Lignes du noyau, avec un septième bloc nul :

```
 1  0  0  0  0  0 -1  0
 1 -1  1 -1  0  0  0  0
 0  1  0  0  0 -1  0  0
 0  0  1  0 -1  0  0  0
 0  0  0  1  0  0  0 -1
 0  0  0  0  1 -1  1 -1
```

La relaxation décrit seulement les voisins des 34 orbites extérieures dans le noyau. Elle impose les degrés des racines choisies et les comptes exacts de voisins communs des paires du noyau. Les arêtes entre racines extérieures restent inconnues.

Sur ses 949 colonnes admissibles, un certificat entier à 37 termes vérifie A^T y>=0 et b^T y=-6. Une combinaison d'équations égale donc -6, tandis que chaque colonne extérieure a une contribution non négative. Ce noyau précis ne peut être prolongé en SRG.

Une seconde construction, par graphes à 33 sommets représentés en masques de bits, retrouve exactement les 949 colonnes. Leurs contributions au certificat sont 0 pour 860 colonnes et 2 pour 89 colonnes. La vérification entière est donc indépendante du solveur et de la première construction des colonnes.

## Fermeture des 25 classes à six blocs — certificat fini vérifié indépendamment

La campagne donne pour chacune des 25 classes un vecteur entier y tel que A^T y>=0 et b^T y<0. Ici A est la matrice des contributions des voisinages extérieurs possibles et b celle des besoins exacts du noyau. Toute complétion fournirait des choix non négatifs x avec Ax=b. On aurait alors b^T y=x^T A^T y>=0, contradiction. Le solveur sert à trouver y ; sa réponse numérique ne sert pas de preuve.

L'exhaustivité des voisinages extérieurs utilise l'équation propre hors support : un représentant z voit autant de p_i que de q_i. Ce nombre commun vaut 0, 1 ou 2. Trois p_i voisins de z imposeraient trois paires dont les lignes de B ont au plus un voisin commun, donc un triangle dans H=C4+C4, impossible. Pour chaque orbite extérieure, on énumère tous les deux ensembles de voisins de ces tailles, avec les deux intérieurs de son étiquette. Les bornes lambda/mu déjà présentes dans le graphe à 33 sommets éliminent uniquement des profils impossibles. Les arêtes entre sommets extérieurs ne sont pas prescrites.

La seconde implémentation reconstruit séparément la couverture des noyaux à six blocs. Elle fixe un premier grand intervalle à {0,1,2,3}. À réflexion près, le deuxième commence en 3 ou en 4. Elle remplit ensuite les déficits d'incidence des huit représentants par quatre petits blocs, en appariant les occurrences restantes. Elle vérifie les bornes lambda/mu et l'absence de prisme, puis réduit par les translations/réflexions. Elle retrouve exactement 25 classes. L'ensemble de ces classes est égal à l'ensemble des 25 modèles certifiés, sans doublon.

Pour chacune, cette même implémentation régénère tous les voisinages extérieurs admissibles, vérifie le second membre du certificat et contrôle chaque coefficient en arithmétique entière. **Les 25 certificats passent.** Il ne s'agit donc pas de l'extrapolation d'un premier noyau.

Exemple particulièrement court : pour le modèle [0,4,16,21,32,35], un certificat à trois termes combine les besoins de voisins communs des paires (9,18), (18,20), (18,30), avec coefficients -1,-1,+1. Leurs besoins sont chacun 1, donc le second membre vaut -1. Sur les 993 colonnes possibles, la contribution est 0 pour 971 colonnes et 1 pour 22 colonnes. Les sommets sont numérotés x=0, c+=2c+1, c-=2c+2, p_i=15+i, q_i=23+i ; les identifiants des lignes de bloc sont définis dans les programmes.

**Conclusion.** Le lemme humain laisse six ou sept blocs. Les certificats excluent tous les six blocs. Tout support Wagner considéré doit donc utiliser les sept blocs : un bloc porte quatre représentants, les six autres en portent deux.

## Données et reproduction

- [Certificats des 25 classes](certificats_six_blocs.json), avec leurs trois lots d'origine et [manifest](MANIFEST.json).
- [Programme de découverte distant](campagne_distante.sage) : classification et recherche de certificats ; ses deux modèles par défaut ne constituent pas une campagne complète.
- [Vérificateur indépendant sans solveur](verifier_famille.py) : couverture des six blocs, génération indépendante des colonnes et vérification entière de tous les certificats.
- [Sorties des vérifications](sorties_verification.txt).

L'empreinte SHA-256 du tableau consolidé, après tri récursif des clés des objets et sérialisation JSON compacte, est :

```
b91eb6f9e16267042cfe18fb93a826b6201a204e41ee120df73c5d38fe20ed41
```

Cette empreinte a été calculée sur les données vérifiées à distance et retrouvée sur les fichiers sauvegardés. Le vérificateur a été exécuté sur SageMathCell, pas sur le PC de l'utilisateur.

## Limites

Même une exclusion complète de Wagner à norme 8 ne tuerait pas g^-=12 : il resterait les normes 9 et 10 que permet la borne de Blichfeldt. La branche involution et la cible « involution impose un prisme » restent ouvertes.

Les contrôles positifs ont été exécutés exactement à distance : pour m=2, A²+A=2I+2J, S²+S=2I et rang_3(S-I)=2 ; pour m=11, A²+A=20I+2J, S²+S=20I et rang_3(S-4I)=35. Le nouveau lemme emploie spécifiquement sept blocs et B=J-I-A(W_8), donc il ne s'applique pas à ces graphes de valeurs propres différentes. Ces contrôles vérifient le cadre général ; ils ne remplacent pas la couverture et les certificats propres à Wagner.
