# Dossier 6 — Involution : minimum, genre et réduction à un réseau de déterminant 49

26 septembre 2026. Dossier initial et complément de C1 traités ensemble. Recherche manuelle ; aucun SAT, CP, PLNE ni solveur externe. Aucune consultation des dépôts ou documents préexistants du projet.

**Résultat principal démontré et relu indépendamment :**

\[
\boxed{\min\Lambda_{-4}^{-}\ge8,\qquad g^-\ge12.}
\]

Les cas **g⁻ = 8 et g⁻ = 10 sont exclus** par une preuve combinatoire, un calcul exact du déterminant du réseau racines, puis une borne de Blichfeldt dont la démonstration est rappelée ci-dessous. Les certificats numériques sont des inégalités entre entiers, vérifiées indépendamment.

**La cible principale n'est pas atteinte : g⁻ = 18 et g⁻ = 20 restent ouverts.** De plus, les résultats fournis ne forcent pas l'existence d'une droite-bloc générique. Sans cette hypothèse, exclure seulement 18 et 20 ne suffirait pas à exclure l'involution.

Deux autres résultats démontrés :

- Tout réseau entier contenant le repère issu de Γ a un minimum au moins égal à 2 ; plus précisément, \(\min N^*\ge8/7\).
- Tous les rangs restants se ramènent à une **contenance** dans un réseau entier de rang 22, de déterminant 49 et de discriminant 7-adique anisotrope de dimension 2. La condition de génération intégrale n'est conservée que lorsque g⁻ = 20.

Ces résultats n'utilisent pas P = 0 : ils valent pour toute involution du cadre du dossier.

## 1. Notations et audit des prémisses

Pour éviter les collisions de notation, on écrit

\[
S=\begin{pmatrix}-I_7&C\\ C^\top&T\end{pmatrix},\qquad
CC^\top=12I_7,\quad CT=0,
\]

où C est l'incidence signée des 42 racines aux 7 blocs, et T le bloc racines × racines. On réserve N au **réseau** de Gram \(G=3I-S\), et on pose \(B_7=S-3I\pmod7\).

Le Gram G est semi-défini positif sur 49 générateurs, avec radical de rang 27. Le réseau N est son quotient par le radical entier ; il est défini positif de rang 22. Ainsi « réseau de Gram 3I−S » ne signifie pas une matrice inversible de taille 49.

Les assertions du dossier suivantes tiennent, avec les précisions indiquées :

| Assertion | Audit |
|---|---|
| \(\det\Lambda_3^-=\det\Lambda_{-4}^-=7^{g^-}\) | Tient : recollement des deux réseaux primitifs dans le réseau standard unimodulaire. |
| g⁻ pair et g⁻ ≠ 22 par l'oddité | Tient. Le nombre M dans \(497-4M\) compte les arêtes négatives **entre racines**, pas toutes les entrées négatives de S. |
| Vecteurs universels z_c et w_c | Tient ; les signes d'incidence doivent être gardés si les représentants sont fixés globalement. |
| \(\det N=7^{22-g^-}\) | Tient : \(N\simeq\sqrt7(\Lambda_{-4}^-)^*\). |
| Les 42 vecteurs racines engendrent N intégralement | Tient, mais les seules relations \(3u_c+\sum C_{co}u_o=0\) ne suffisent pas ; voir ci-dessous. |
| Exclure 18 et 20 suffit à exclure l'involution | Conditionnel : cela suffit si au moins un Q_c est inversible. Cette existence n'est pas démontrée dans les données. |

**Raccord pour la génération intégrale.** Les relations donnent d'abord que l'indice du sous-réseau engendré par les racines divise une puissance de 3. Modulo 3, S est semi-simple, car \(S^2+S=0\). La projection de \(\ker S\) sur les sept coordonnées blocs est surjective : son annihilateur serait un vecteur \((a,0)\in\operatorname{im}S\), donc \(S(a,0)=-(a,0)\), ce qui impose \(C^\top a=0\). Or Cᵀ est injective modulo 3 : les deux colonnes de chaque support donnent \(a_i+a_j=a_i-a_j=0\). La saturation et le projecteur 3-entier relèvent cette surjectivité au réseau propre entier. Avec les relations dont les coordonnées blocs sont \(3e_c\), l'indice est donc 1.

La formule de genre donne bien

\[
\chi(\operatorname{disc}\beta^-)=(-1)^{g^-/2}.
\]

Ainsi **β⁻ est scindée**. Cela ne contredit pas le discriminant non scindé du réseau N décrit au §4 : il s'agit de deux formes différentes.

## 2. Lemme uniforme de minimum : aucun vecteur propre entier de norme carrée ≤ 7

**Énoncé démontré.** Si \(0\ne f\in\mathbf Z^{49}\) et \(Sf=-4f\), alors \(f^\top f\ge8\).

Supposons \(f^\top f\le7\). Une coordonnée de valeur absolue 2 laisserait aux autres coordonnées une somme des valeurs absolues au plus égale à 3. L'équation de cette ligne demande pourtant une contribution hors diagonale de valeur absolue 8 sur une racine, ou 6 sur un bloc. C'est impossible. Toutes les coordonnées non nulles valent donc ±1 ; notons m leur nombre, avec m ≤ 7.

En changeant le représentant de chaque orbite sélectionnée, rendons ces m coordonnées égales à +1. Appelons p_i les représentants et q_i = σp_i leurs images. Ces deux ensembles de sommets sont disjoints.

Soit A₊ l'adjacence induite sur les p_i, et soit

\[
D_{ij}=[p_i\sim q_j].
\]

La matrice D est symétrique, à entrées 0 ou 1. Sa diagonale vaut 1 exactement pour les orbites qui sont des blocs. L'équation propre devient

\[
D\mathbf1=A_+\mathbf1+4\mathbf1.
\]

Chaque degré d_j de D est donc au moins 4. Deux p_i distincts ont au plus deux voisins communs parmi les q_j, par λ = 1 et μ = 2. En comptant ces voisins communs par colonnes,

\[
6m\le\sum_j\binom{d_j}{2}\le2\binom m2=m(m-1).
\]

Il faut m ≥ 7, donc m = 7. Toutes les inégalités sont des égalités : D est 4-régulière, A₊ = 0, et deux lignes distinctes ont exactement deux 1 communs. Par conséquent

\[
D^2=2I+2J.
\]

La valeur propre sur les constantes est 4. Les six autres valeurs propres sont ±√2. L'intégralité du polynôme caractéristique impose trois occurrences de chaque signe ; donc tr D = 4. Exactement **quatre** des sept orbites sélectionnées seraient des blocs.

Écrivons maintenant f = (a,b), selon blocs et racines. Les lignes blocs donnent

\[
Cb=-3a.
\]

Si h blocs sont sélectionnés, ils demandent au moins 3h incidences de racines sélectionnées. Mais les 7−h racines sélectionnées n'en portent que 2(7−h). Donc

\[
3h\le2(7-h),\qquad h\le2.
\]

Cela contredit h = 4. ∎

**Portée précise.** Cette preuve emploie les adjacences réelles entre représentants et images, avec la borne sur leurs voisins communs. Elle n'est pas une conséquence de la seule identité signée \(S^2+S=12I\). Elle ne nécessite ni la classification des carrés D ni l'absence de prismes.

## 3. Exclusion complète des rangs g⁻ = 8 et g⁻ = 10

### 3.1 Première application : le réseau complet de rang 22

Le réseau \(\Lambda=\Lambda_{-4}^-\) est de rang 22 et de déterminant \(7^{g^-}\). La borne de Blichfeldt donne

\[
\min\Lambda\le\frac2\pi(12!)^{1/11}7^{g^-/22}.
\]

Pour g⁻ = 8, le membre de droite vaut environ 7,948955, strictement moins que le minimum 8 du §2.

**Certificat exact.** En utilisant seulement \(\pi>157/50\), l'inégalité stricte requise est impliquée par

\[
12!\,7^4\,25^{11}<314^{11}.
\]

Les deux entiers sont

```text
gauche = 2742011169433593750000000000
droite = 2925675324340604967987390464
ecart  =  183664154907011217987390464
```

Ils ont été recalculés indépendamment en Python et en JavaScript BigInt. Avec la parité et le plafond déjà démontrés dans le dossier :

\[
\boxed{g^-\in\{10,12,14,16,18,20\}.}
\]

**Démonstration de la borne utilisée.** Normalisons un réseau de dimension n pour que sa distance minimale soit √2. Posons \(F(x)=\max(0,1-\|x\|^2)\). Pour les t points du réseau à distance < 1 d'un même x, notons y_i leurs différences à x. L'identité

\[
t\sum_i\|y_i\|^2-\Big\|\sum_i y_i\Big\|^2
=\sum_{i<j}\|y_i-y_j\|^2\ge t(t-1)
\]

entraîne \(\sum_i F(y_i)\le1\). L'intégration sur un domaine fondamental donne

\[
\operatorname{covol}\Lambda\ge\int F
=\frac{\pi^{n/2}}{\Gamma(n/2+2)}.
\]

En rétablissant le minimum et en prenant n = 22, on obtient exactement la borne ci-dessus. Aucune classification de réseaux n'intervient.

Sur le réseau complet de rang 22, la même borne vaut environ 9,487 à g⁻ = 10 et ne suffit pas. Le passage au sous-réseau de rang 15 permet de renforcer le résultat.

### 3.2 Déterminant exact du sous-réseau racines

Posons \(L_T=\{r\in\mathbf Z^{42}:Cr=0,\ Tr=-4r\}\), de rang 15. Il s'identifie aux vecteurs de \(\Lambda=\Lambda_{-4}^-(S)\) dont les sept coordonnées blocs sont nulles. Son minimum est donc lui aussi au moins 8.

**Lemme démontré et relu indépendamment :**

\[
\boxed{\det L_T=2^{12}7^{g^--7}.}
\]

Voici le raccord entier, qui est essentiel. Soient W le réseau engendré par les sept \(w_c=(-4e_c,C^\top e_c)\), et P l'image de Λ par projection sur les coordonnées blocs. On a

\[
\operatorname{Gram}W=28I_7,\qquad
4\mathbf Z^7\subseteq P,\qquad
\Lambda/(W\perp L_T)\simeq P/4\mathbf Z^7.
\]

Il reste à calculer P, et non simplement son rang.

Modulo 2, le projecteur \((3I-S)/7\) est entier et Λ se réduit sur \(\ker S\), puisque \(S^2=S\). L'annihilateur de la projection de \(\ker S\) sur les blocs est constitué des a tels que \((a,0)\in\operatorname{im}S\). L'équation \(S(a,0)=(a,0)\) équivaut à \(C^\top a=0\). Or, modulo 2, les colonnes de C sont les \(e_i+e_j\), donc ce noyau est la droite des constantes. Ainsi **P modulo 2 est l'hyperplan de somme paire**.

Le groupe \(\mathbf Z^7/P\) est donc cyclique d'ordre 2 ou 4 : il a un seul générateur modulo 2 et son exposant divise 4. Mais

\[
\frac12\sum_cw_c\in\Lambda
\]

est entier, car chaque colonne de C a somme paire ; sa projection est \(-2\mathbf1\). Si le quotient avait ordre 4, P serait le noyau d'un caractère vers Z/4 dont tous les coefficients sont impairs, puisque sa réduction modulo 2 est la somme des coordonnées. Ce caractère prendrait la valeur \(-2\cdot7\equiv2\) sur \(-2\mathbf1\), contradiction.

Donc \([\mathbf Z^7:P]=2\) et \([P:4\mathbf Z^7]=2^{13}\). L'identité des déterminants donne

\[
7^{g^-}=\frac{28^7\det L_T}{(2^{13})^2},
\]

d'où la formule annoncée. ∎

### 3.3 Deuxième application de Blichfeldt : g⁻ = 10 est impossible

À g⁻ = 10, le réseau L_T a rang 15 et déterminant

\[
2^{12}7^3=1\,404\,928.
\]

La même preuve de Blichfeldt donne

\[
\min L_T\le\frac2\pi\Gamma(19/2)^{2/15}
\big(2^{12}7^3\big)^{1/15}
\simeq7{,}7733005734<8,
\]

ce qui contredit le §2. On a utilisé \(\Gamma(19/2)=17!!\sqrt\pi/2^9\).

**Certificat entièrement entier.** Avec \(17!!=34\,459\,425\) et \(\pi>157/50\), il suffit de vérifier

\[
(17!!)^2\,7^3\,50^{14}<2^{36}157^{14}.
\]

```text
gauche = 248593765970705795288085937500000000000000
droite = 379904972191464137167101541730810607960064
ecart  = 131311206220758341879015604230810607960064
```

Ce certificat a été recalculé indépendamment ; la seconde voie utilise l'expression équivalente en 18! et 9!. La même borne est croissante en g⁻ et exclut aussi g⁻ = 8.

**Conclusion finale :**

\[
\boxed{g^-\in\{12,14,16,18,20\}.}
\]

À g⁻ = 12, la borne en rang 15 vaut environ 10,076 et reste compatible avec le minimum 8. Aucune exclusion de ce rang, de 18 ou de 20 n'en résulte.

## 4. Réduction uniforme au déterminant 49, avec une condition de minimum sur le dual

Soient u_a les 49 générateurs du réseau N. Puisque \(G^2=7G\), ils forment un repère serré :

\[
\sum_a(v,u_a)^2=7\|v\|^2,
\qquad S((v,u_a))_a=-4((v,u_a))_a.
\]

Pour \(v\in N^*\), tous ces produits sont entiers. Le §2 donne donc

\[
\boxed{\min N^*\ge8/7.}
\]

Tout sur-réseau **entier** L contenant N est contenu dans N* et a des normes entières. Il vérifie donc \(\min L\ge2\).

Passons à la réduction arithmétique. Sur Z₇, écrivons

\[
\Lambda_{-4}^-=U\perp7V,
\qquad \operatorname{rang}U=h=22-g^-.
\]

Alors

\[
N_7\simeq7U^*\perp V^*.
\]

Le groupe discriminant de N est \(\mathbf F_7^h\), avec une forme de caractère

\[
\chi(\det U)=(-1)^{g^-/2}.
\]

Le caractère d'une forme scindée en dimension h serait

\[
(-1)^{h/2}=(-1)^{11-g^-/2}=-(-1)^{g^-/2}.
\]

La forme discriminante est donc **non scindée**, d'indice de Witt \((h-2)/2\). Le relèvement d'un sous-espace totalement isotrope maximal fournit un sur-réseau entier N_max tel que

\[
\boxed{\det N_{\max}=49,\qquad
[N_{\max}:N]=7^{(20-g^-)/2}.}
\]

Son discriminant est anisotrope de dimension 2. En outre

\[
\min N_{\max}\ge2,
\qquad \min N_{\max}^*\ge8/7,
\]

car \(N_{\max}^*\subseteq N^*\).

**Ce que cette réduction permet réellement.** Pour exclure toutes les involutions par la voie des réseaux, il suffirait d'exclure la **présence du repère**, avec son incidence et sa réalisation combinatoire, dans tous ces réseaux de déterminant 49. Il ne suffirait pas d'exclure les réseaux de déterminant 49 **engendrés** par le repère : lorsque g⁻ < 20, les générateurs engendrent N, de strict indice dans N_max.

En particulier, le cas g⁻ = 18 est ramené à un sous-réseau d'indice 7 d'un tel N_max.

## 5. Un candidat concret de déterminant 49 éliminé par la nouvelle borne

Ce contrôle montre à la fois la limite du genre seul et l'utilité du minimum du dual.

Considérons le réseau entier unimodulaire impair

\[
D_{12}^+=D_{12}\cup\big((1/2)^{12}+D_{12}\big),
\qquad D_{12}=\{z\in\mathbf Z^{12}:\sum z_i\text{ pair}\}.
\]

Le vecteur

\[
a=(3/2,3/2,1/2,\ldots,1/2)
\]

est primitif et de norme carrée 7. Son orthogonal entier K a rang 11, déterminant 7 et minimum 2. Ainsi

\[
L_0=K\perp K
\]

a rang 22, déterminant 49, minimum 2 et discriminant anisotrope \(\mathbf F_7^2\). K possède 594 vecteurs de norme 3 : le comptage exact des signes donne 45 + 504 + 45. L₀ en possède donc 1 188. Cela ne construit pas le repère du dossier.

Mais avec \(b=e_1+e_2\), la projection

\[
b-\frac37a\in K^*
\]

a pour norme

\[
2-\frac97=\frac57<\frac87.
\]

**Conclusion démontrée : aucun réseau L₀ de cette construction ne peut contenir le repère provenant d'une involution de Γ.** Le genre, le déterminant, le minimum 2 et l'abondance de vecteurs de norme 3 ne le rejetaient pas ; la borne du dual le rejette.

Cela exclut ce candidat explicite, pas tous les réseaux de déterminant 49.

## 6. Vecteurs de norme 2 : deux profils seulement

Si v est de norme carrée 2 dans un sur-réseau entier contenant N, les produits \(f_a=(v,u_a)\) sont entiers et satisfont \(f^\top f=14\).

Les équations propres et λ/μ excluent les profils comportant deux ou trois coordonnées de valeur absolue 2. Elles excluent aussi qu'une unique coordonnée double soit une racine. Il ne reste que :

- quatorze produits de valeur absolue 1 ;
- un produit de valeur absolue 2 sur un **bloc** et dix produits de valeur absolue 1.

Dans le second cas, après changement de signe, \((v,u_c)=2\), et

\[
u_c=v+v',\quad \|v\|^2=\|v'\|^2=2,\quad(v,v')=0.
\]

Donc un des sept vecteurs de norme 4 se décompose en deux racines orthogonales. Les preuves détaillées, relues indépendamment, sont dans [PREUVE_RESEAU.md](PREUVE_RESEAU.md). Aucun des deux profils restants n'est exclu.

## 7. Complément de C1 : ce qui renforce réellement la cible

### 7.1 Généricité et rangs encore ouverts

La formule transmise

\[
\operatorname{rang}_7 B_{7,I_cI_c}=6+\operatorname{rang}_7 Q_c
\]

implique qu'un seul Q_c inversible force g⁻ ≥ 18. Avec le genre, il reste alors effectivement **18 ou 20**. Je n'ai pas réexécuté les 3 000 boules de C1 : ce contrôle reste attribué à C1.

Si les deux rangs sont exclus, la conséquence inconditionnelle serait que **les sept Q_c doivent être singuliers**, et non encore qu'il n'existe aucune involution. Le cas entièrement dégénéré doit être fermé séparément.

Le couplage avec la partie paire donne également :

- k = 28 force \((g^-,g^+)=(20,22)\), sans hypothèse de généricité ;
- exclure g⁻ = 20 forcerait k ≤ 26 dans la branche involution ;
- exclure aussi g⁻ = 18 forcerait k ≤ 24, mais laisserait les rangs inférieurs de la partie impaire dans le cas dégénéré.

### 7.2 Résidus : nilpotence vérifiée algébriquement

Toutes les égalités suivantes sont **modulo 7**. Écrivons

\[
B_7=\begin{pmatrix}P&X\\X^\top&Y\end{pmatrix},
\quad B_7^2=0,\quad X=PC.
\]

Le symbole C de ce paragraphe désigne une solution de X = PC, distincte de l'incidence du §1. Posons

\[
R=Y-C^\top PC.
\]

Le premier bloc de \(B_7^2=0\) donne \(P^2+PCC^\top P=0\). Le deuxième donne

\[
XR=-(P^2+XX^\top)C=0.
\]

Puis, en écrivant \(Y=R+C^\top PC\), les termes croisés s'annulent et le dernier bloc donne

\[
0=Y^2+X^\top X
=R^2+C^\top(PCC^\top P+P^2)C=R^2.
\]

La congruence de Schur fournit aussi \(\operatorname{rang}B_7=\operatorname{rang}P+\operatorname{rang}R\). Ainsi les identités \(XR=0\), \(R^2=0\) et, pour Q_c inversible, \(\operatorname{rang}R=g^--18\), sont compatibles et démontrées sous la condition d'image énoncée.

Le signe affine encore ce résidu : à g⁻ = 20, **det Q_c carré donne un résidu de rang 2 anisotrope** ; det Q_c non carré donne un résidu scindé. Dans le premier cas, une entrée diagonale nulle de R force toute sa ligne à être nulle. C'est la même propriété de rang 2 qui était utile au dossier 5 ; elle ne constitue pas, seule, une nouvelle contradiction.

### 7.3 Précision indispensable sur le critère « 7-isotrope »

Il faut distinguer le réseau propre complet sur 49 coordonnées et le sous-réseau

\[
L_T=\{r\in\mathbf Z^{42}:Cr=0,\ Tr=-4r\}.
\]

L'annulation de la somme des **seules coordonnées racines** pour tout vecteur propre complet est fausse : le vecteur universel

\[
w_c=(-4e_c,C^\top e_c)
\]

a une somme de carrés sur les racines égale à \(12\equiv5\pmod7\). Sa norme complète est bien \(16+12=28\equiv0\).

Pour L_T, la formulation correcte est : « toutes les normes sont nulles modulo 7 ». Par polarisation, cela équivaut à un **Gram entièrement nul modulo 7**, donc de rang 0.

Sur Z₇, on a la décomposition orthogonale

\[
\Lambda_{-4}^-(S)=\langle w_1,\ldots,w_7\rangle\perp L_T,
\qquad \operatorname{Gram}(w_c)=28I_7.
\]

En effet, si f = (a,b), on retranche la combinaison des w_c qui annule les coordonnées blocs ; 4 est une unité dans Z₇. Il s'ensuit

\[
\boxed{\operatorname{rang}_7\operatorname{Gram}(L_T)=22-g^-.}
\]

Ce rang vaut donc **2 à g⁻ = 20**, et **4 à g⁻ = 18**. Prouver le rang ≤ 1 suffirait à contredire le genre ; prouver l'annulation de toutes les normes serait plus fort. Mais le dossier ne donne pas encore de mécanisme qui impose l'une de ces deux propriétés. Le critère est une reformulation de l'obstruction recherchée.

Enfin, \(g^-=\operatorname{rang}_7(T-3I)\) tient : la décomposition \(\mathbf F_7^{42}=\operatorname{im}C^\top\perp\ker C\), avec \(CC^\top=12I\), permet de comparer directement les deux rangs.

## 8. Contrôles positifs et fichiers reproductibles

Le programme [verifier_dossier6.py](verifier_dossier6.py) reconstruit S depuis les deux graphes certifiés produits auparavant dans cette conversation. Il vérifie l'automorphisme, l'identité SRG, l'identité impaire, les vecteurs universels, les rangs modulaires et une diagonalisation par congruence sur Z/3⁶Z.

| Contrôle | Taille de S | Valeurs propres | Rang₃(S−rI) | Jordan du Gram N en 3 | det N |
|---|---:|---|---:|---|---|
| H(2,3), m = 2 | 4 | 1, −2 | 2 | \(1^2\) | 1 |
| BvLS, m = 11 | 121 | 4, −5 | 35 | \(1^{35}3^{10}9^{10}\) | \(3^{30}\) |

Les déterminants des réseaux propres négatifs sont respectivement \(3^2\) et \(3^{80}\). À m = 11, l'écart propre vaut 9 : on ne remplace donc pas le recollement par un groupe élémentaire de dimension égale au seul rang modulo 3.

**Aucun contrôle positif n'est exclu.** En particulier, les vecteurs racines du contrôle m = 2 ont norme 1 dans son réseau N. La borne 8/7 du §4 ne lui est pas appliquée : elle utilise les constantes −4 et 7 propres à 99.

Fichiers :

- [controles_exacts.json](controles_exacts.json) : sortie complète des contrôles et certificat entier ;
- [PREUVE_RESEAU.md](PREUVE_RESEAU.md) : preuves détaillées de minimum et de profils ;
- [AUDIT_GENRE.md](AUDIT_GENRE.md) : relecture indépendante du genre et des nouveaux arguments ;
- [VERIFICATION_INDEPENDANTE_G10.md](VERIFICATION_INDEPENDANTE_G10.md) : deuxième preuve du déterminant de L_T et de l'exclusion de g⁻ = 10 ;
- [ETAT_REPRISE.md](ETAT_REPRISE.md) : état exact pour une reprise manuelle.

## 9. Ce qui reste à démontrer

Le dossier produit deux exclusions uniformes supplémentaires, **g⁻ = 8 et g⁻ = 10**, et des restrictions qui éliminent aussi un candidat de réseau de déterminant 49. Il ne tue pas g⁻ = 18 ou 20 et ne prouve pas la présence d'un prisme.

La réduction la plus forte obtenue ici est désormais : un vrai Γ avec involution impose un réseau entier N_max de rang 22 et déterminant 49, de discriminant anisotrope, avec minimum du dual ≥ 8/7, contenant le repère et sa réalisation combinatoire. À g⁻ = 18, ce repère engendre un sous-réseau d'indice 7. Une classification fondée seulement sur le déterminant et sur les générateurs de norme 3 perdrait cette distinction.

Pour une attaque par les étoiles, la cible de C1 reste pertinente lorsque Q_c est inversible : résidu nul à g⁻ = 18, ou rang 2 à g⁻ = 20. Le cas où les sept Q_c sont singuliers reste distinct et ouvert.

**FIN DU RAPPORT — aucune exclusion complète de l'involution annoncée.**
