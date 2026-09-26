# Attaque des conjectures prioritaires : supports courts, ombres et raccords

26 septembre 2026. Mathématiques et petits calculs exacts, sans solveur externe.

## Résultat obtenu

La conjecture « les vecteurs courts forcent un cube interdit » est remplacée par deux résultats démontrés :

1. Un vecteur non nul de (L_T\) de norme carrée au plus 10 a des coefficients non nuls égaux à ±1. Après choix des représentants, ses sommets positifs forment un stable. Son support complet est biparti et 4-régulier.
2. À norme carrée 8, il existe exactement deux formes nécessaires : le double biparti du complément du cube et celui du complément du graphe de Wagner. **La première force un prisme dans le cadre à sept blocs.** Elle est donc exclue à (P=0\).

La seconde forme possède un noyau partiel explicitement étiqueté sur sept blocs, sans prisme, passant les bornes locales de λ/μ. Sa complétion au graphe entier reste inconnue. Les normes 9 et 10 restent également à traiter.

Les attaques des deux autres conjectures ont précisé leurs limites : les projections binaires seules ne donnent pas l'obstruction d'ombres espérée ; les images de deux cartes ont des intersections forcées qu'une hypothèse de transversalité générique ignore.

**Aucune involution n'est exclue par ce rapport.** Les valeurs encore ouvertes de (g^-) restent 12, 14, 16, 18 et 20. La nouvelle exclusion porte sur une famille de supports courts, sous (P=0\).

## 1. Un support court est forcément biparti

### Hypothèses et conventions

On suppose un vrai srg(99,14,1,2) avec l'involution du dossier 6. Soit (b\in L_T\), non nul : ses sept coordonnées de bloc sont nulles, et ses coordonnées de racine satisfont (Tb=-4b\), (Cb=0\).

Le minimum 8 a déjà été démontré dans `../dossier6_involution/PREUVE_RESEAU.md`. Supposons ici (\|b\|^2\le10\).

Si une coordonnée vaut (a\), avec (|a|\ge2\), l'équation de sa ligne impose

\[
4|a|\le\sum_{j\ne i}|b_j|\le\sum_{j\ne i}b_j^2\le10-a^2.
\]

C'est impossible dès (|a|=2\), et a fortiori au-delà. Toutes les coordonnées non nulles valent donc ±1. Leur nombre (m\) est 8, 9 ou 10.

Choisissons dans chaque orbite le représentant (p_i\) de coefficient +1, et posons (q_i=\sigma p_i\). Soient

\[
A_{ij}=[p_i\sim p_j],\qquad B_{ij}=[p_i\sim q_j].
\]

Les matrices sont symétriques ; leurs diagonales sont nulles. Le support complet a l'adjacence

\[
\begin{pmatrix}A&B\\B&A\end{pmatrix},
\qquad B\mathbf1=A\mathbf1+4\mathbf1.
\]

Notons (e\) le nombre d'arêtes de (A\), et (d_i\) ses degrés. Le vecteur entier du graphe complet correspondant vaut +1 sur les (p_i\), −1 sur les (q_i\), et 0 ailleurs. C'est un vecteur propre de valeur −4.

### Deux comptes extérieurs exacts

Pour un sommet extérieur (z\), posons

\[
t_z=|N(z)\cap\{p_i\}|=|N(z)\cap\{q_i\}|.
\]

L'égalité vient de l'équation propre en (z\). Le sommet fixe (x\) a (t_x=0\), puisque les (p_i,q_i\) sont des racines. L'involution groupe tous les autres sommets extérieurs par deux et donne (t_{\sigma z}=t_z\).

Le nombre d'incidences extérieures du côté positif est

\[
s:=\sum_z t_z=10m-4e.
\]

En comptant les voisins communs extérieurs de toutes les paires positives, λ/μ donnent

\[
\begin{aligned}
D:=\sum_z\binom{t_z}{2}
&=m(m-1)-e-\sum_i\binom{d_i}{2}-\sum_i\binom{d_i+4}{2}\\
&=m(m-7)-7e-\sum_i d_i^2.
\end{aligned}
\]

Ces deux identités utilisent la fermeture du graphe à 99 sommets, pas seulement le noyau induit.

### Première conséquence : (e\) est pair

Dans (D\), chaque paire extérieure d'involution apporte deux fois la même quantité, et le sommet fixe apporte 0. Donc (D\) est pair.

Or (m(m-7)\) est pair, et (\sum d_i^2\equiv\sum d_i=2e\equiv0\pmod2\). La formule de (D\) impose

\[
\boxed{e\equiv0\pmod2.}
\]

### Seconde conséquence : (e=0\) aux normes 8, 9 et 10

On a (\sum d_i^2\ge\sum d_i=2e\), d'où

\[
D\le m(m-7)-9e.
\]

Pour (m=8\), tout (e\ge1\) donnerait (D<0\). Ainsi (e=0\).

Pour (m=9,10\), si (e>0\), la parité impose (e\ge2\). Pour tout entier (t\ge0\),

\[
t-\binom t2\le[,t>0,].
\]

Le nombre de sommets extérieurs actifs est donc au moins

\[
s-D=m(17-m)+3e+\sum d_i^2\ge m(17-m)+5e.
\]

Mais il y a au plus (98-2m\) sommets extérieurs actifs, puisque (x\) est inactif. À (m=9,e\ge2\), le minorant vaut au moins 82, contre 80 disponibles. À (m=10,e\ge2\), il vaut au moins 80, contre 78 disponibles. Contradiction.

**Conclusion démontrée :** (A=0\) et (B\mathbf1=4\mathbf1\) aux trois petites normes. Il reste à classifier des matrices d'adjacence symétriques 4-régulières, plutôt que des supports arbitraires.

## 2. Norme 8 : classification et prisme forcé

### Les deux formes

À (m=8\), deux lignes distinctes de (B\) ont au plus deux 1 communs, puisque les (p_i\) sont non adjacents. Posons (C=J-I-B\). C'est l'adjacence d'un graphe cubique à huit sommets, et

\[
B^2=I+2C+C^2.
\]

Une arête de (C\) appartenant à un triangle donnerait une entrée hors diagonale de (B^2\) au moins 3. Donc (C\) est sans triangle.

Un graphe cubique sans triangle à huit sommets a un quadrangle : sinon, les deux premiers niveaux autour d'un sommet contiendraient (1+3+6=10\) sommets distincts. Chaque sommet du quadrangle a un voisin extérieur.

Les quatre voisins extérieurs sont distincts. En effet, un sommet extérieur ne peut toucher plus de deux sommets du quadrangle, et deux doivent être opposés. Les répartitions non bijectives des quatre attaches sont (2,1,1,0) ou (2,2,0,0). Dans le premier cas, le sommet sans attache est adjacent aux trois autres extérieurs ; les deux extérieurs à une attache doivent aussi être adjacents entre eux, créant un triangle. Dans le second cas, les deux extérieurs sans attache sont adjacents et ont un autre voisin commun, créant encore un triangle.

Le graphe consiste donc en deux quadrangles et un couplage entre eux. Les 24 permutations du couplage donnent deux classes : le cube (Q_3\), ou le graphe de Wagner (W_8\), cycle de longueur 8 avec ses quatre cordes opposées.

Une seconde implémentation énumère directement **les 19 355 graphes cubiques étiquetés** d'ordre 8. Parmi eux, les 3 360 sans triangle se répartissent exactement en 840 cubes et 2 520 Wagner. Cela contrôle la classification sans reprendre la réduction aux deux quadrangles.

### La forme cubique force un prisme

Dans le cas (C=Q_3\), les seules paires de lignes distinctes de (B\) ayant moins de deux 1 communs sont les quatre paires antipodales du cube ; elles ont zéro commun.

Fixons un bloc (c\). Soient (P_c\) les indices des (p_i\) contenant ((c,+)\), et (M_c\) ceux contenant ((c,-)\). La relation (Cb=0\) donne (|P_c|=|M_c|\).

Deux indices de (P_c\) ont déjà l'intérieur ((c,+)\) comme voisin commun. Leur nombre de communs parmi les (q_i\) doit donc être au plus 1. Ils sont nécessairement antipodaux dans (Q_3\). Ainsi (|P_c|\le2\), et de même pour (M_c\).

Les huit racines sélectionnées portent seize incidences aux blocs. Si chaque bloc avait (|P_c|\le1\), les sept blocs en porteraient au plus quatorze. Il existe donc un bloc tel que

\[
P_c=\{i,i'\},\qquad M_c=\{j,j'\},
\]

avec deux paires antipodales distinctes. Entre ces deux paires, (B\) est un couplage. Choisissons les noms pour avoir (B_{ij}=B_{i'j'}=1\).

On obtient les deux triangles

\[
\{(c,+),p_i,q_j\},\qquad
\{(c,-),q_{i'},p_{j'}\}.
\]

Leurs trois barreaux sont

\[
(c,+)(c,-),\qquad p_iq_{i'},\qquad q_jp_{j'}.
\]

Les deux dernières arêtes existent parce que (B\) relie chaque paire antipodale. Les autres arêtes croisées sont absentes : elles relient deux positifs, deux négatifs, ou une racine à un intérieur qu'elle ne contient pas. Les six sommets induisent donc exactement un prisme triangulaire.

**Résultat démontré et contrôlé par six certificats explicites : à (P=0\), un vecteur de norme carrée 8 dans (L_T\) ne peut avoir la forme cubique.**

## 3. Wagner : un témoin partiel qui empêche une fausse conclusion

Pour la seconde forme, indexons (W_8\) par 0,…,7, avec les arêtes (i(i+1)\) modulo 8 et (i(i+4)\). Posons (B=J-I-A(W_8)\).

Voici huit représentants de racine, en classes et signes ; le bloc 6 n'est pas utilisé par ce support :

| indice | racine (p_i\) |
|---:|---|
| 0 | (0,+)(1,+) |
| 1 | (0,−)(2,+) |
| 2 | (0,+)(3,+) |
| 3 | (0,−)(4,+) |
| 4 | (2,−)(5,+) |
| 5 | (1,−)(3,−) |
| 6 | (1,+)(4,−) |
| 7 | (1,−)(5,−) |

Inclure (x\), les quatorze intérieurs, ces huit (p_i\) et leurs huit images (q_i\). Mettre les arêtes usuelles du moulin, les incidences des racines aux intérieurs, et (p_i\sim q_j\) exactement si (B_{ij}=1\). Il n'y a pas d'arête entre positifs ni entre négatifs.

Ce noyau induit à 31 sommets a :

- une involution fixant (x\) ;
- un vecteur propre de valeur −4, égal à ±1 sur les seize racines et nul ailleurs ;
- 19 triangles et aucun prisme ;
- au plus un commun sur chaque arête et au plus deux sur chaque non-arête.

Ces propriétés sont reproduites par un programme indépendant utilisant des ensembles de voisins, sans importer le générateur NumPy.

**Ce témoin n'est pas un srg(99,14,1,2).** Beaucoup de demandes λ/μ doivent encore être remplies par les 68 racines restantes. Il interdit seulement de déclarer Wagner impossible sur ce noyau, ou d'annoncer que tout support court serait déjà cubique.

### Ce qui reste pour exclure (g^-=12\)

La borne de Blichfeldt sur (L_T\) à (g^-=12\) impose un vecteur de norme carrée au plus 10. À (P=0\), les seuls cas encore à exclure sont désormais :

1. norme 8, forme Wagner ;
2. norme 9, double biparti d'un graphe 4-régulier à neuf sommets ;
3. norme 10, double biparti d'un graphe 4-régulier à dix sommets.

Dans les deux derniers cas, chaque paire de lignes de (B\) a au plus deux 1 communs. L'étiquetage, les communs extérieurs et les apex restent nécessaires. Ce rapport ne les a pas classifiés.

## 4. Ombres : la conjecture 33 demande réellement le Green cubique

Pour le Gram complet attendu de dix ombres, posons

\[
K=3I+3J-A(4K_2+2K_1)\pmod7.
\]

Il est inversible : son déterminant vaut 2 modulo 7. Si (d_i=n_i-n_i'\) sont les quatre différences de membres appariés, alors

\[
g(d_i,d_j)=\delta_{ij}\pmod7.
\]

Chacun des seize choix de signes dans (\sum_i\varepsilon_i d_i\) a donc norme 4 modulo 7.

J'ai testé l'espoir que **la seule forme binaire à deux points par couche** forcerait néanmoins un choix isotrope. Ce mécanisme est insuffisant : on peut construire dix telles ombres et une forme symétrique (G\), de rang 35, telle que

\[
G\mathbf1=0,\qquad M^TGM=K.
\]

Les dix ombres sont indépendantes, ainsi que leurs projections sur (1^\perp\). On complète ces projections en une base et on prescrit la forme (K\) sur leur span, puis une forme non dégénérée sur son complément. Le programme fournit la matrice complète 36×36 et vérifie les seize normes, toutes égales à 4.

**Portée exacte :** cela écarte une preuve utilisant seulement les projections binaires, le rang et les sommes nulles. Cette (G\) n'est pas annoncée comme le pseudo-inverse d'un laplacien cubique. L'obstruction uniforme de l'étoile complète reste ouverte ; elle doit exploiter cette condition supplémentaire, ou d'autres contraintes véritables sur les ombres.

## 5. Deux cartes : une borne que la transversalité doit respecter

Cette attaque vise la piste 25. Utilisons la forme globale β de C1, sur (G\) de dimension (14+k\). Écrivons

\[
W_x=\langle g_v:v\in B_1(x)\rangle,\qquad G_x=W_x^{\perp_\beta}.
\]

Pour deux centres adjacents (x,a\), leurs boules ont en commun la droite ({x,a,b}\). Le Gram de ses trois générateurs est nul sur la diagonale et vaut 4 hors diagonale ; son déterminant est (2\cdot4^3\equiv2\pmod7\), donc il est non dégénéré de rang 3.

Ce span est contenu dans (W_x\cap W_a\), et

\[
G_x+G_a\subseteq(W_x\cap W_a)^{\perp_\beta},
\qquad
\boxed{\dim(G_x+G_a)\le k+11\le39.}
\]

Pour deux centres non adjacents, leurs deux communs sont non adjacents : s'ils étaient adjacents, leur arête aurait deux communs. Le Gram des deux générateurs vaut (
\begin{pmatrix}0&3\\3&0\end{pmatrix}
\), de rang 2. On obtient

\[
\boxed{\dim(G_x+G_a)\le k+12\le40.}
\]

Même les trois cartes centrées sur une même droite restent orthogonales au span de cette droite, donc leur somme a dimension au plus 39.

**Conséquence :** deux cartes ne peuvent être traitées comme deux sous-espaces génériquement transverses engendrant tout le radical global à (k=28\). Une attaque par raccord doit exploiter des restrictions sur les arêtes effectives, ou au moins trois centres ne portant pas la même droite. Ces bornes ne réfutent pas la possibilité qu'un raccord forcé viole précisément 39 ou 40 ; aucun tel raccord n'est démontré ici.

## 6. Vérifications et fichiers

| fichier | contenu |
|---|---|
| `support_huit.py` | matrices nécessaires, six prismes cubiques, recherche d'un noyau Wagner |
| `support_huit_resultats.json` | matrices, étiquettes et certificats exacts |
| `verifier_independant.py` | énumération indépendante des cubiques et contrôle par ensembles de voisins |
| `verification_independante.json` | 19 355 cubiques ; 3 360 sans triangle ; noyau Wagner vérifié |
| `contre_test_ombres.py` | construction exacte de la forme abstraite sur les ombres |
| `contre_test_ombres.json` | dix supports, matrice 36×36, seize normes |

Les contrôles existants (m=2\) et (m=11\) ont été relancés avec `../dossier6_involution/verifier_dossier6.py` : identités SRG, involution, identités impaires et vecteurs universels restent validés. Ils ne sont pas exclus. Le lemme de support court ci-dessus utilise spécifiquement le degré 14, la valeur propre −4 et les 99 sommets ; il n'est pas appliqué à leurs paramètres différents.

## 7. Priorité suivante

Le prochain sous-problème concret est de **compléter ou tuer le support Wagner avec les communs extérieurs et le couplage D**. Le noyau publié donne un point de départ exact ; il faut examiner toutes les étiquettes admissibles avant toute exclusion uniforme de cette forme. En parallèle, les normes 9 et 10 sont maintenant réduites à des supports positifs stables.

La conjecture de l'étoile complète reste prometteuse, mais le raccourci proposé par une somme signée de différences ne suffit pas. Les déterminants 49 et (7^4\) ne sont pas exclus : aucun nouveau minimum dual n'est prétendu.
