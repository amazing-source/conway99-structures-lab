# Rapport d’Astra, suite du dossier n° 5 (réponse au suivi n° 3), relayé par l’utilisateur le 26/09/2026

Texte d’Astra, copié tel quel. Relecture de « structures » : JOURNAL et NOTE_S26.

Astra — Étoiles d’ombres : exclusion d’une famille cyclique et contrainte à k=24

1. Une famille entière exclue à k=22

Statut : réduction démontrée à la main ; deux cas résiduels exclus par des énumérations exactes indépendantes.

Famille étudiée

Dans chacune des trois parties A,B,C, les douze sommets sont
(i,t)\in\mathbb Z_6\times\{0,1\}.
Leur indice local est i+6t.

Les couplages internes sont
(i,0)\sim(i+\delta,1),
\qquad(\delta_A,\delta_B,\delta_C)=(0,r,s),
avec r,s\in\{1,\ldots,5\}, r\ne s.

Les arêtes entre parties sont
A_{i,t}\sim B_{i,t},\qquad
A_{i,t}\sim C_{i,t},\qquad
B_{i,t}\sim C_{i+q,t},
où q\in\{1,\ldots,5\}.

Cela représente 100 triplets de paramètres étiquetés. Ces boules sont sans triangle, et leurs sous-graphes sur deux parties n’ont aucun cycle de longueur quatre. Elles passent donc ces premières conditions nécessaires à P=0.

Cette forme cyclique n’est pas forcée par P=0. Le résultat porte sur cette famille précise.

Réduction par le discriminant

Notons H=3I-B. Comme les six racines de l’unité appartiennent à \mathbb F_7, la décomposition de Fourier donne
\[
\tau(G)\equiv
2\det H(-1)\det H(2)^2\det H(3)^2\pmod7.
\]

Le facteur 2 vient des valeurs propres non nulles du bloc trivial :
\[
\frac{2\cdot3^2\cdot5^2}{36}\equiv2\pmod7.
\]

Posons
T_q=
\begin{pmatrix}
0&1&1\\
1&0&(-1)^q\\
1&(-1)^q&0
\end{pmatrix},
\qquad
D=\operatorname{diag}(1,(-1)^r,(-1)^s).
Alors
\det H(-1)=
\det(3I-T_q-D)\det(3I-T_q+D).

Les deux déterminants d’ordre trois donnent : Parités de r,s

q pair

q impair

0,0

0

6

0,1, 1,0, 1,1

1

2

Les valeurs du tableau sont celles de \det H(-1), modulo 7.

À k=22, une droite générique doit avoir \tau(G) non carré par L23. Le tableau force donc
\{r,s\}=\{2,4\},\qquad q\in\{1,3,5\}.

Les changements de coordonnées ramènent ces six paramètres à deux représentants :
(r,s,q)=(2,4,1),\qquad(2,4,3).

Leurs nombres d’arbres couvrants sont :
\[
\tau_1=3\,096\,828\,960\,300\equiv5\pmod7,
\]
\[
\tau_3=2\,905\,653\,667\,500\equiv3\pmod7.
\]

Ils sont donc bien génériques et échappent à L23 seul.

Exclusion par les ombres d’une étoile

On fixe u=A_0. Pour chaque non-partenaire v, on énumère les ombres n_v dont la projection sur A est \{0,v\}.

Chaque candidat vérifie :

* deux non-partenaires dans chaque partie ;
* Bn\le 2\mathbf1-n ;
* aucune paire de son support interdite par MM^\top ;
* \(n^\top H^+n=6\pmod7\).

Pour deux candidats, on impose
g(n,n')\in\{2,3\},\qquad
a(n,n')=3-g(n,n'),
n^\top n'\le2-a(n,n'),
et, dans l’étoile sélectionnée,
\deg_a(n_v)\le1-(Bn_v)_u.

Ce dernier test exprime la capacité du couplage dans N(u).

Sept groupes suffisent déjà à obtenir une contradiction :

Représentant

Indices v

Candidats par groupe

q=1

1,2,3,4,5,7,11

381,336,354,336,381,362,362

q=3

1,2,3,4,5,7,9

375,348,405,348,375,392,351

Cas

Recherche par masques, ordre adaptatif

Recherche indépendante, ordre fixe

q=1

39 574 nœuds, aucune solution

50 646 nœuds, aucune solution

q=3

41 826 nœuds, aucune solution

52 872 nœuds, aucune solution

Les recherches sont exhaustives, sans plafond atteint. J’ai également refait séparément les deux étoiles complètes : mêmes candidats, aucune solution, respectivement 29 925 et 30 225 nœuds.

L’égalité
\sum_v n_v=Ke_u
n’a pas été nécessaire : une relaxation de l’étoile complète est déjà impossible.

En revanche, toutes les sous-étoiles de six groupes parmi les dix passent cette relaxation, pour chacun des deux représentants. L’obstruction vient donc réellement du raccord entre plusieurs ombres.

Conclusion. Aucun membre générique de cette famille cyclique ne peut apparaître à k=22. Les membres dégénérés et les autres familles de boules restent ouverts.

2. À k=24, au moins huit ombres doivent sortir de la norme 6

Statut : démontré ; relecture indépendante effectuée.

Supposons la droite générique, k=24, et \tau(G) carré modulo 7. Posons
Z=M^\top H^+M,\qquad
R=A_F-3I-3J+Z.

Alors R a rang 2, et sa partie non dégénérée est anisotrope sur \mathbb F_7. On peut écrire
R=uu^\top+vv^\top,
avec u,v indépendants.

Puisque -1 n’est pas carré modulo 7,
R_{yy}=0\iff u_y=v_y=0.
Or
R_{yy}=Z_{yy}-6.

Ainsi :

Toute ombre de norme 6 a une ligne de R entièrement nulle. Toutes ses adjacences lointaines sont donc déjà forcées par les ombres, comme à k=22.

Notons
J_{\mathrm{bad}}=\{y:Z_{yy}\ne6\}.
L’identité MR=0 donne Mu=Mv=0. Les deux vecteurs étant portés par J_{\mathrm{bad}},
\dim\ker M_{J_{\mathrm{bad}}}\ge2.

Lemme d’indépendance des petites familles d’ombres

Cinq colonnes quelconques de M sont indépendantes sur \mathbb F_7.

En effet, leur projection sur une partie est l’incidence non orientée de cinq arêtes distinctes au plus. Une dépendance à ce niveau est portée par un C_4, avec coefficients alternés +1,-1,+1,-1.

Si elle était une dépendance de M, on aurait un échange
n_1+n_2=n_3+n_4.
Le passage du modulo 7 aux entiers est valide : chaque entrée de la différence appartient à [-2,2].

Dans chaque partie, une ombre positive rencontre chacune des deux ombres négatives en un sommet. Cela donne trois voisins communs, contre \mu=2.

Pourquoi sept mauvaises ombres ne suffisent pas

Six colonnes ont donc une nullité au plus 1.

Avec sept mauvaises ombres, leur matrice aurait rang au plus 5. Or l’incidence non orientée de sept arêtes simples a rang au moins 5 ; si son rang vaut 5, son graphe contient un C_4.

Les rangs de la projection et de la matrice complète seraient donc égaux à 5, et leurs noyaux seraient égaux. La relation alternée de ce C_4 remonterait dans M, donnant l’échange interdit précédent.

Par conséquent,
\boxed{|J_{\mathrm{bad}}|\ge8.}

Portée. C’est une contrainte sur les ombres, valable à k=24 lorsque \tau(G) est carré. Elle ne donne pas encore de contradiction. Lorsque \tau(G) est non carré, le résidu de rang deux est scindé : une diagonale nulle ne force plus une ligne nulle, et cette preuve ne s’applique pas.

3. Droites dégénérées : classification des composantes de taille 4+4+4

Statut : classification démontrée ; déterminants vérifiés exactement.

Considérons une composante de G(L) qui rencontre chacune des trois parties en quatre sommets.

Après normalisation, deux bijections entre parties sont l’identité et la troisième est \pi\in S_4. Posons
a=(01)(23),\qquad b=(02)(13),\qquad c=(03)(12).

L’absence de triangle et de prisme passant par la droite laisse exactement trois types, à renommage des sommets et permutation des parties près :

Type

Couplages internes

Images de \pi, dans l’ordre 0,1,2,3

\tau

I

a,b,b

(1,3,0,2)

8 100

II

a,b,b

(2,3,1,0)

8 640

III

a,b,c

(1,0,3,2)

8 640

La classification utilise l’action S_4\to S_3 sur les trois couplages. Une énumération indépendante des 3^3\times9=243 données normalisées retrouve les 66 données admissibles et ces trois types.

Comme
\[
8100\equiv1,\qquad8640\equiv2\pmod7,
\]
chaque composante contribue exactement une dimension au noyau de son laplacien modulo 7.

Donc, pour une boule composée de trois telles composantes,
\boxed{d(L)=3,\qquad k\ge20.}

Cela distingue nettement deux mécanismes : ici, la dégénérescence vient uniquement de la déconnexion, sans torsion supplémentaire du groupe critique en 7.

Aucun patron déconnecté n’est exclu par ce résultat.

4. Contrôles et questions restantes

Les calculs ont été exécutés en mémoire, sans solveur externe ni écriture de fichier. Les déterminants ont été recalculés en entiers par Bareiss.

Sur les contrôles positifs m=2 et m=11, j’ai revérifié l’identité SRG et les hypothèses combinatoires du lemme d’indépendance des ombres : projections bijectives sur les paires non partenaires et intersections d’ombres au plus 2. Aucun conflit. Les conclusions spécifiques à k=22 ou 24 pour les paramètres 99 ne leur ont pas été appliquées.

Restent indécis :

* l’exclusion de k=22 pour toutes les boules génériques ;
* une obstruction à k=24 couvrant les deux signes du discriminant ;
* l’impossibilité que les 231 droites soient toutes dégénérées.

La réduction cyclique est une exclusion de famille ; elle ne constitue pas encore le plafond général recherché.
