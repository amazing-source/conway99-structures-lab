# Rapport d’Astra (26/09) : réponse au suivi n° 4 du dossier 5, et audit du critère S5

Texte relayé par l’utilisateur, copié tel quel (sorties jointes comprises).

# Retransmission intégrale — Astra, 26 septembre 2026

# MESSAGE 1 — Suivi n° 4 du dossier n° 5

## 1. Les 18 conditions des deux apex ne suffisent pas

**VÉRIFIÉ EXACTEMENT, par des implémentations indépendantes.** Il existe dix ombres qui satisfont toutes les conditions proposées sur les deux apex, les dix normes et la somme entière de l'étoile. Elles ne satisfont pas le Gram complet.

On utilise la boule cyclique du précédent rapport, de paramètres
\[
(r,s,q)=(2,4,1).
\]
Pour rendre le témoin autonome : chaque partie \(A,B,C\) est indexée par \((i,t)\in\mathbb Z_6\times\{0,1\}\), d'indice local \(i+6t\). Les couplages internes sont \((i,0)\sim(i+\delta,1)\), avec \(\delta_A=0,\delta_B=2,\delta_C=4\). Les arêtes croisées sont
\[
A_{i,t}\sim B_{i,t},\qquad A_{i,t}\sim C_{i,t},\qquad B_{i,t}\sim C_{i+1,t}.
\]
Les indices globaux sont \(A=0,\ldots,11\), \(B=12,\ldots,23\), \(C=24,\ldots,35\). On fixe \(u=A_0\).

Notations : \(B\) désigne ici la matrice d'adjacence de la boule, \(H=3I-B\),
\[
H^+=(H+J)^{-1}-J\pmod7,\qquad
K=12I-B+2J-B^2-D.
\]
Dans cette formule, \(J=J_{36}\), \(D=\operatorname{diag}(J_{12},J_{12},J_{12})\), et \(g(n,n')=n^\top H^+n'\).

| Groupe \(v\), ombre sur \(A=\{0,v\}\) | Support complet |
|---:|---|
| 1 | 0, 1, 21, 23, 30, 33 |
| 2 | 0, 2, 13, 19, 32, 33 |
| 3 | 0, 3, 13, 16, 29, 31 |
| 4 | 0, 4, 12, 23, 27, 35 |
| 5 | 0, 5, 15, 18, 27, 28 |
| 7 | 0, 7, 15, 22, 25, 26 |
| 8 | 0, 8, 14, 19, 31, 34 |
| 9 | 0, 9, 14, 17, 32, 35 |
| 10 | 0, 10, 20, 22, 28, 29 |
| 11 | 0, 11, 16, 21, 24, 26 |

Les deux apex sont les groupes **4** et **11** : ils contiennent respectivement \(u_B=12\) et \(u_C=24\).

Le vérificateur reconstruit les matrices et contrôle :

- deux non-partenaires dans chaque partie ;
- \(Bn\le2\mathbf1-n\) pour chaque ombre ;
- \(g(n,n)=6\) pour les dix ombres ;
- produit \(g=3\) de chacun des apex avec chacune des neuf autres ombres ;
- somme exacte \(\sum n=Ke_u\), en entiers ;
- intersections de deux ombres au plus égales à deux ;
- borne entière \(M_\star M_\star^\top\le K\), entrée par entrée.

Les capacités \(1-(Bn)_u\), dans l'ordre du tableau, sont
\[
(1,1,1,0,1,1,1,1,1,0).
\]

Le premier programme a trouvé le témoin après 525 nœuds. Une reconstruction indépendante en Python entier, puis le vérificateur joint, ont confirmé ses propriétés.

**Pourquoi cela ne donne pas une étoile complète.** Par exemple,
\[
g(n_1,n_2)=4,
\]
alors qu'il faudrait \(2\) ou \(3\). Les relations entre les huit ombres non apex échouent. Le témoin ne fournit ni un graphe local complet, ni une matrice globale \(M\), ni un SRG.

**Conséquence précise.** Une contradiction uniforme ne peut pas être tirée des seules 18 conditions proposées, même renforcées des dix normes, de la somme entière et des bornes entières ci-dessus. Elle doit utiliser d'autres relations, notamment celles entre les huit ombres non apex.

Programme : [verifier_temoin_apex.py](verifier_temoin_apex.py). Données et Gram calculé : [temoin_apex_verifie.json](temoin_apex_verifie.json).

## 2. Sommes et formulation par potentiels

**DÉMONTRÉ, mais automatique.** La formule du brief est correcte :
\[
\sum_{v\ne w}g(n_v,n_w)=5+(Bn_w)_u\pmod7.
\]
Avec le terme diagonal, elle devient \(4+(Bn_w)_u\).

En particulier, sous la somme exacte de l'étoile, la condition entre les deux apex est déjà impliquée par leurs normes et leurs seize produits avec les huit autres ombres. Pour un apex \(a\) et l'autre apex \(b\),
\[
5=g(a,Ke_u)=6+8\cdot3+g(a,b),
\]
d'où \(g(a,b)=3\).

La formulation
\[
x_w=H^+n_w,\qquad Hx_w=n_w+\mathbf1
\]
est également correcte. Elle ne rend pas le problème strictement local autour d'une arête : \(H^+\) dépend de toute la boule et les six positions du support de \(n_w\) interviennent.

Le témoin du §1 passe ces formulations. Elles ne fournissent donc pas, seules, l'invariant uniforme recherché.

## 3. À \(k=24\), préciser le relâchement des mauvaises ombres

**DÉMONTRÉ ET RELU INDÉPENDAMMENT.** Dans la branche \(\tau(G)\) carré, écrivons
\[
R=uu^\top+vv^\top,\qquad Mu=Mv=0,
\]
et \(t_y=(u_y,v_y)\), \(h_y=R_{yy}=g(n_y,n_y)-6\). L'anisotropie donne
\[
h_y=0\iff t_y=0.
\]

Une ombre bonne a donc **toute sa ligne de \(R\) nulle**. Pour chaque paire distincte ayant au moins une ombre bonne, on doit conserver
\[
g(n_y,n_z)=3-A_{yz}\in\{2,3\}.
\]
Seules les paires de deux mauvaises ombres peuvent échapper à cette règle.

Pour chaque étoile \(S(a)\), l'équation \(M(u,v)=0\) impose
\[
\sum_{y\in S(a)\cap J_{\rm bad}}t_y=0.
\]
Ainsi le nombre \(j_a\) de mauvaises ombres dans une étoile vérifie
\[
\boxed{j_a=0\ \text{ou}\ j_a\ge2.}
\]
S'il vaut exactement deux, pour les ombres \(y,z\) concernées :
\[
t_z=-t_y,\quad h_z=h_y\ne0,\quad R_{yz}=-h_y,
\]
donc
\[
\boxed{g(n_y,n_z)=3-h_y-A_{yz}.}
\]

Enfin,
\[
\sum_{a\in N(L)}j_a=6|J_{\rm bad}|.
\]
Des minorants \(\ell_a\) étoile par étoile donnent donc
\[
|J_{\rm bad}|\ge
\left\lceil\frac{\sum_a\ell_a}{6}\right\rceil.
\]
La borne \(6j_{\min}\) exige bien un minorant valable sur **les 36 étoiles**.

**INDÉCIS.** Aucun \(j_{\min}>0\) uniforme n'est établi ici ; aucune recherche exhaustive du minimum avec mauvaises ombres n'est annoncée. Ces contraintes précisent le test à effectuer, sans remplacer le calcul demandé. La borne globale précédente \(|J_{\rm bad}|\ge8\) reste valable.

## 4. Portée et contrôles

Le retour de C1/« structures » sur les composantes à quatre sommets par partie est intégré : elles peuvent subsister à \(P=0,S>0\), mais sont exclues à \(P=S=0\). La classification précédente ne doit donc pas être présentée comme une branche ouverte dans ce dernier sous-cas. Cette exclusion est un résultat transmis par C1, non une nouvelle preuve de ce rapport.

Les contrôles positifs \(m=2,11\) valident les propriétés universelles des ombres utilisées dans le rapport précédent. Le présent témoin est contrôlé directement et n'est pas un critère d'exclusion de ces graphes. Les conditions spécifiques à \(99\), \(k=22\) ou \(24\), ne leur sont pas appliquées.

**L'exclusion uniforme de \(k=22\) reste ouverte.** Les 42 boules éliminées, rapportées par C1/« structures », constituent une calibration et non une preuve universelle.

## Annexe de retransmission : Gram complet et valeurs interdites

Ordre commun des lignes et colonnes : `1, 2, 3, 4, 5, 7, 8, 9, 10, 11`. Toutes les valeurs sont modulo 7.

```text
6 4 5 3 2 2 4 6 4 3
4 6 2 3 3 0 3 1 0 3
5 2 6 3 1 6 5 4 4 3
3 3 3 6 3 3 3 3 3 3
2 3 1 3 6 3 6 6 6 3
2 0 6 3 3 6 3 6 0 3
4 3 5 3 6 3 6 2 4 3
6 1 4 3 6 6 2 6 2 3
4 0 4 3 6 0 4 2 6 3
3 3 3 3 3 3 3 3 3 6
```

Les valeurs hors diagonale qui ne sont ni 2 ni 3 sont, sous la forme `(v,w,g)` :

```text
(1, 2, 4)
(1, 3, 5)
(1, 8, 4)
(1, 9, 6)
(1, 10, 4)
(2, 7, 0)
(2, 9, 1)
(2, 10, 0)
(3, 5, 1)
(3, 7, 6)
(3, 8, 5)
(3, 9, 4)
(3, 10, 4)
(5, 8, 6)
(5, 9, 6)
(5, 10, 6)
(7, 9, 6)
(7, 10, 0)
(8, 10, 4)
```

Ces 19 valeurs interdisent la complétion du Gram de ce témoin. Les lignes des deux apex, groupes 4 et 11, passent pourtant toutes leurs conditions.

### Retour reçu lors de la retransmission

C1 rapporte que quatre boules échouent encore au test complet après suppression des bornes lambda/mu contre la boule, des zéros de K et des bornes de recouvrement. Ce retour est compatible avec le témoin ci-dessus : les seules conditions des apex passent, tandis que le couplage complet des huit autres ombres échoue. Le problème résiduel est un problème algébrique sur des ombres binaires et un motif de couplage, pas une obstruction de rang du Gram. Cette observation sur quatre boules ne constitue pas encore une preuve uniforme ; je ne la présente pas comme une nouvelle vérification exécutée ici.

**FIN DU MESSAGE 1.**

---

# MESSAGE 2 — Audit adverse du critère « table de Brouwer »

## 1. Verdict sur la preuve

**DÉMONTRÉ : le critère proposé tient.** Aucun pas faux ni contre-exemple n'a été trouvé. La répartition des facteurs premiers de l'étape 4 demande toutefois d'écrire l'argument ci-dessous.

| Étape | Verdict |
|---|---|
| 1. Réseau propre pair | **Tient.** Pour \(\theta\ne k\), tout \(x\in\Lambda_\theta\) a somme nulle ; \(x^\top x\equiv\sum x_i=0\pmod2\). |
| 2. Projecteur 2-entier et unimodularité | **Tient, y compris pour \(v\) pair.** L'hypothèse assure l'intégralité 2-adique de \(E_\theta\). La primitivité donne \(E_\theta\mathbb Z^v=\Lambda_\theta^\#\), donc l'image locale est autoduale. |
| 3. Déterminant modulo 4 | **Tient.** Ajouter explicitement que \(m_\theta\) est pair, ce qui rend \(m_\theta/2\) bien défini. Chaque plan pair unimodulaire a déterminant \(-1\bmod4\). |
| 4. Facteurs \(q\mid v,\ q\nmid d\) | **Tient**, avec la justification de l'attribution donnée ci-dessous. |

Pour l'égalité \(E_\theta\mathbb Z^v=\Lambda_\theta^\#\), une forme linéaire entière sur le sous-réseau primitif \(\Lambda_\theta\) s'étend à \(\mathbb Z^v\). Elle est représentée par un vecteur entier ; sa projection représente la forme sur \(\Lambda_\theta\). Cela prouve la surjectivité, et pas seulement une inclusion.

### Justification du « tout ou rien »

Soient \(q\mid v,\ q\nmid d\), \(t=v_q(v)\), et
\[
W=\mathbb Z_q^v\cap\mathbf1^\perp.
\]
Les projecteurs sur \(W\) sont \(q\)-entiers, donc
\[
W=\Lambda_r\otimes\mathbb Z_q\ \perp\ \Lambda_s\otimes\mathbb Z_q.
\]
Son groupe discriminant est cyclique d'ordre \(q^t\). Une somme directe de deux groupes \(q\)-primaires non triviaux n'étant pas cyclique, toute cette partie appartient à un seul réseau propre.

Sur ce groupe discriminant, \(A\) agit par \(k\), en utilisant les représentants \(e_i-\mathbf1/v\). Sur la partie issue de \(\Lambda_\theta\), il agit par \(\theta\). Il faut donc
\[
k-\theta\equiv0\pmod{q^t}.
\]
Or
\[
(k-r)(k-s)=\mu v,\qquad (k-r)-(k-s)=-d.
\]
Exactement un facteur est divisible par \(q\), et sa valuation est au moins \(t\). On obtient exactement l'attribution annoncée :
\[
q^t\mid\det\Lambda_\theta
\iff v_q(k-\theta')<t.
\]

Le corollaire \(NL_e(n)\) tient également : \(n\equiv1\bmod4\), \(m_{e-n}=e(n+1)\equiv2\bmod4\), tandis que tous les facteurs premiers possibles du déterminant sont \(1\bmod4\). Contradiction avec le signe imposé.

## 2. Graphes réels et réseaux saturés

**VÉRIFIÉ EXACTEMENT.** Les graphes ont été construits directement, leurs paramètres vérifiés, puis les noyaux entiers saturés calculés par transformations unimodulaires. Les Gram et leurs déterminants ont été calculés par Bareiss.

| Graphe | \(\theta\) | Rang | Déterminant | Hypothèse 2-adique |
|---|---:|---:|---:|---|
| Petersen | \(1\) | 5 | \(162=2\cdot3^4\) | Fausse |
| Petersen | \(-2\) | 4 | \(405=5\cdot3^4\) | Vraie |
| \(T(7)\) | \(3\) | 6 | \(109375=7\cdot5^6\) | Vraie |
| \(T(7)\) | \(-2\) | 14 | \(46875=3\cdot5^6\) | Vraie |

Toutes les conclusions applicables passent. Petersen vérifie notamment que l'hypothèse lorsque \(v\) est pair ne doit pas être supprimée.

**Seconde vérification, analytique.** Pour \(T(n)\), \(n\) impair, l'incidence \(C\) sommet–arête de \(K_n\) donne \(A=C^\top C-2I\) et \(CC^\top=(n-2)I+J\). Les réseaux saturés ont les déterminants
\[
\det\Lambda_{n-4}=n(n-2)^{n-1},
\qquad
\det\Lambda_{-2}=\frac{n-1}{2}(n-2)^{n-1}.
\]
Le premier vient de \(C^\top A_{n-1}\), le second de \(\ker_{\mathbb Z}C\), sachant que l'image entière de \(C\) a indice deux. Le complément de \(T(5)\) fournit Petersen. Ces formules reproduisent indépendamment les quatre déterminants.

Programme : [verifier_reseaux.py](verifier_reseaux.py). Bases saturées, Gram et résultats : [reseaux_resultat.json](reseaux_resultat.json).

## 3. Balayage indépendant et comparaison de la liste

**VÉRIFIÉ EXACTEMENT PAR DEUX IMPLÉMENTATIONS.**

Les 26 pages de la [table de Brouwer, jusqu'à 1 300 sommets](https://aeb.win.tue.nl/graphs/srg/srgtab.html) ont été consultées intégralement, avec récupération des fins de page initialement tronquées. Les statuts ont été lus dans la colonne d'existence : « + », « ! » et nombres suivis de « ! ». Les couleurs HTML n'ont pas été utilisées.

Résultats, en comptant une seule orientation par complément lorsque \(k\le(v-1)/2\) :

| Contrôle | Résultat |
|---|---:|
| Lignes extraites | 4 538 |
| Lignes marquées existantes | 1 183 |
| Lignes existantes exclues | **0** |
| Paramètres « ? » exclus, à complément près | **35** |
| Éléments de la liste fournie non retrouvés | **0** |
| Paramètres « ? » exclus supplémentaires | **0** |

Le critère retrouve aussi **deux** ensembles déjà marqués inexistants :
\[
(154,51,8,21),\qquad(289,54,1,12).
\]
Le premier n'était pas cité dans le brief comme contrôle ; ce n'est pas une exclusion supplémentaire parmi les « ? ».

Une génération indépendante à partir de
\[
k(k-\lambda-1)=(v-k-1)\mu
\]
et de l'intégralité des valeurs propres et des multiplicités retrouve les mêmes 37 exclusions canoniques. Il s'agit d'une génération de paramètres **arithmétiquement admissibles**, non d'une décision de leur réalisabilité.

Un second calcul, en JavaScript et utilisant explicitement les valuations, retrouve les mêmes 35 paramètres « ? » et aucune ligne existante exclue.

- [audit_parametres.py](audit_parametres.py)
- [audit_parametres_independant.js](audit_parametres_independant.js)
- [Données extraites de la table](brouwer_rows_2026-09-26.json)
- [Résultats, liste complète et différences](audit_parametres_resultats.json)

## 4. Antériorité : ne pas annoncer plus que ce qui est établi

**Lemme de réseau classique.** Jacques Martinet, *On parity classes*, notes de Lausanne, mai 2003, corollaire 6.4, donne la congruence de signature d'un réseau pair à déterminant impair. Dans le cas positif, c'est exactement
\[
(-1)^{m/2}\det\Lambda\equiv1\pmod4.
\]
Cette partie n'est pas nouvelle. [Texte primaire, p. 14](https://jamartin.perso.math.cnrs.fr/Othertexts/paritylat.pdf#page=14).

**Cas \((289,54,1,12)\) déjà exclu.** Bondarenko–Radchenko, *J. Combin. Theory Ser. B* 103 (2013), 521–531, DOI 10.1016/j.jctb.2013.05.005, théorème 1, excluent explicitement ces paramètres dans leur classification d'une famille à \(\lambda=1\). Leur preuve utilise des symétries et une structure sur \(\mathbb F_3\). C'est un contrôle, pas une nouvelle exclusion. [Article primaire](https://arxiv.org/pdf/1201.0383).

**Antériorités proches, de portée différente.** Le théorème 3 de Koyama, *On Applications of the Bruck-Ryser-Chowla Theorem* (1991), fournit une condition pour les SRG d'ordre impair ; avec valeurs propres entières, l'équation considérée admet déjà une solution triviale parce que son coefficient est \((r-s)^2\). [Article primaire](https://teapot.lib.ocha.ac.jp/record/34740/files/KJ00004470758.pdf).
Bagchi emploie les discriminants et invariants p-adiques des espaces propres SRG, mais son théorème 4.2 vise des designs quasi symétriques ayant le SRG comme graphe de blocs, ce qui est une hypothèse supplémentaire. [Article primaire](https://arxiv.org/pdf/2102.05576).

**INDÉCIS : nouveauté du critère paramétrique exact et du corollaire \(NL_e(n)\).** Ils n'ont pas été retrouvés dans cette recherche bibliographique bornée. L'absence de référence retrouvée ne constitue pas une preuve de nouveauté. Le statut « ? » dans la table ne garantit pas non plus l'absence d'une publication ailleurs.

## 5. Reproduction

Les scripts n'utilisent aucun SAT, CP, PLNE ou solveur externe. Le vérificateur des ombres utilise NumPy avec arithmétique entière modulo 7 ; les vérifications de réseaux et de paramètres utilisent la bibliothèque standard Python.

Depuis ce dossier :

~~~text
python -B verifier_temoin_apex.py
python -B verifier_reseaux.py
python -B audit_parametres.py
node audit_parametres_independant.js
~~~

Le second audit des paramètres a également été exécuté directement en JavaScript dans l'environnement de travail. Tous les résultats d'impossibilité ici sont limités à leur portée annoncée : aucune résolution de srg(99,14,1,2) n'est revendiquée.

## Annexe de retransmission : liste des 35 paramètres retrouvés

```text
(169, 42, 5, 12)
(169, 70, 27, 30)
(209, 52, 15, 12)
(235, 52, 9, 12)
(289, 90, 23, 30)
(289, 126, 53, 56)
(375, 154, 53, 70)
(375, 170, 85, 70)
(391, 140, 39, 56)
(430, 39, 8, 3)
(441, 190, 89, 76)
(494, 85, 12, 15)
(507, 44, 1, 4)
(616, 205, 90, 57)
(625, 130, 15, 30)
(625, 182, 45, 56)
(625, 234, 83, 90)
(625, 286, 129, 132)
(630, 68, 1, 8)
(690, 318, 121, 168)
(721, 220, 69, 66)
(750, 210, 55, 60)
(750, 214, 63, 60)
(780, 369, 158, 189)
(841, 150, 11, 30)
(841, 210, 41, 56)
(841, 270, 79, 90)
(841, 330, 125, 132)
(841, 390, 179, 182)
(871, 348, 137, 140)
(875, 184, 33, 40)
(931, 180, 29, 36)
(1026, 451, 212, 187)
(1148, 481, 210, 195)
(1248, 203, 22, 35)
```

La comparaison à la liste du brief est exacte : aucun manquant, aucun ajout. Les compléments sont également exclus ; ils ne sont pas répétés dans cette liste.

**FIN DU MESSAGE 2.**

**FIN DE LA TRANSMISSION INTÉGRALE.**
SORTIES INTEGRALES ENREGISTREES — aucune nouvelle recherche lancée.

=== temoin_apex_verifie.json ===
{
  "supports": [
    [
      0,
      1,
      21,
      23,
      30,
      33
    ],
    [
      0,
      2,
      13,
      19,
      32,
      33
    ],
    [
      0,
      3,
      13,
      16,
      29,
      31
    ],
    [
      0,
      4,
      12,
      23,
      27,
      35
    ],
    [
      0,
      5,
      15,
      18,
      27,
      28
    ],
    [
      0,
      7,
      15,
      22,
      25,
      26
    ],
    [
      0,
      8,
      14,
      19,
      31,
      34
    ],
    [
      0,
      9,
      14,
      17,
      32,
      35
    ],
    [
      0,
      10,
      20,
      22,
      28,
      29
    ],
    [
      0,
      11,
      16,
      21,
      24,
      26
    ]
  ],
  "groups": [
    1,
    2,
    3,
    4,
    5,
    7,
    8,
    9,
    10,
    11
  ],
  "apex_groups": [
    4,
    11
  ],
  "caps": [
    1,
    1,
    1,
    0,
    1,
    1,
    1,
    1,
    1,
    0
  ],
  "green_gram": [
    [
      6,
      4,
      5,
      3,
      2,
      2,
      4,
      6,
      4,
      3
    ],
    [
      4,
      6,
      2,
      3,
      3,
      0,
      3,
      1,
      0,
      3
    ],
    [
      5,
      2,
      6,
      3,
      1,
      6,
      5,
      4,
      4,
      3
    ],
    [
      3,
      3,
      3,
      6,
      3,
      3,
      3,
      3,
      3,
      3
    ],
    [
      2,
      3,
      1,
      3,
      6,
      3,
      6,
      6,
      6,
      3
    ],
    [
      2,
      0,
      6,
      3,
      3,
      6,
      3,
      6,
      0,
      3
    ],
    [
      4,
      3,
      5,
      3,
      6,
      3,
      6,
      2,
      4,
      3
    ],
    [
      6,
      1,
      4,
      3,
      6,
      6,
      2,
      6,
      2,
      3
    ],
    [
      4,
      0,
      4,
      3,
      6,
      0,
      4,
      2,
      6,
      3
    ],
    [
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      6
    ]
  ],
  "all_norms_6": true,
  "exact_shadow_sum": true,
  "gram_integer_bound": true,
  "overlap_at_most_2": true,
  "structural_matching_groups": [
    [
      1,
      3
    ],
    [
      2,
      5
    ],
    [
      7,
      8
    ],
    [
      9,
      10
    ]
  ],
  "failed_pair": {
    "groups": [
      1,
      2
    ],
    "g": 4
  },
  "status": "counterexample_to_apex_relaxation_only"
}

=== reseaux_resultat.json ===
{
  "status": "all exact assertions passed",
  "method": "unimodular integer kernel, Bareiss, independent incidence formula",
  "graphs": [
    {
      "graph": "Petersen",
      "vertices": [
        [
          0,
          1
        ],
        [
          0,
          2
        ],
        [
          0,
          3
        ],
        [
          0,
          4
        ],
        [
          1,
          2
        ],
        [
          1,
          3
        ],
        [
          1,
          4
        ],
        [
          2,
          3
        ],
        [
          2,
          4
        ],
        [
          3,
          4
        ]
      ],
      "srg_parameters": [
        10,
        3,
        0,
        1
      ],
      "checks": [
        {
          "eigenvalue": 1,
          "lattice_rank": 5,
          "equation_rank": 5,
          "unimodular_transform_determinant": 1,
          "determinant": 162,
          "determinant_factorization": {
            "2": 1,
            "3": 4
          },
          "analytic_determinant": 162,
          "mod4_hypothesis": false,
          "signed_determinant_mod4_if_applicable": null,
          "primes_dividing_v_not_d": [
            {
              "q": 2,
              "v_q_v": 1,
              "assigned": true
            },
            {
              "q": 5,
              "v_q_v": 1,
              "assigned": false
            }
          ],
          "basis_rows": [
            [
              0,
              0,
              -1,
              0,
              1
            ],
            [
              -1,
              2,
              1,
              0,
              -1
            ],
            [
              1,
              -1,
              0,
              0,
              0
            ],
            [
              0,
              -1,
              0,
              0,
              0
            ],
            [
              1,
              -2,
              0,
              0,
              1
            ],
            [
              -1,
              1,
              1,
              1,
              0
            ],
            [
              0,
              1,
              0,
              -1,
              -2
            ],
            [
              0,
              0,
              -1,
              -1,
              -1
            ],
            [
              0,
              0,
              0,
              1,
              1
            ],
            [
              0,
              0,
              0,
              0,
              1
            ]
          ],
          "gram": [
            [
              4,
              -6,
              -2,
              -1,
              2
            ],
            [
              -6,
              12,
              3,
              0,
              -6
            ],
            [
              -2,
              3,
              4,
              2,
              -1
            ],
            [
              -1,
              0,
              2,
              4,
              4
            ],
            [
              2,
              -6,
              -1,
              4,
              10
            ]
          ]
        },
        {
          "eigenvalue": -2,
          "lattice_rank": 4,
          "equation_rank": 6,
          "unimodular_transform_determinant": 1,
          "determinant": 405,
          "determinant_factorization": {
            "3": 4,
            "5": 1
          },
          "analytic_determinant": 405,
          "mod4_hypothesis": true,
          "signed_determinant_mod4_if_applicable": 1,
          "primes_dividing_v_not_d": [
            {
              "q": 2,
              "v_q_v": 1,
              "assigned": false
            },
            {
              "q": 5,
              "v_q_v": 1,
              "assigned": true
            }
          ],
          "basis_rows": [
            [
              0,
              -1,
              0,
              -2
            ],
            [
              -1,
              0,
              1,
              0
            ],
            [
              -1,
              0,
              0,
              0
            ],
            [
              -1,
              -2,
              2,
              -1
            ],
            [
              1,
              1,
              -1,
              0
            ],
            [
              1,
              1,
              -2,
              0
            ],
            [
              1,
              -1,
              0,
              -1
            ],
            [
              0,
              2,
              -1,
              2
            ],
            [
              0,
              0,
              1,
              1
            ],
            [
              0,
              0,
              0,
              1
            ]
          ],
          "gram": [
            [
              6,
              3,
              -6,
              0
            ],
            [
              3,
              12,
              -9,
              9
            ],
            [
              -6,
              -9,
              12,
              -3
            ],
            [
              0,
              9,
              -3,
              12
            ]
          ]
        }
      ]
    },
    {
      "graph": "T(7)",
      "vertices": [
        [
          0,
          1
        ],
        [
          0,
          2
        ],
        [
          0,
          3
        ],
        [
          0,
          4
        ],
        [
          0,
          5
        ],
        [
          0,
          6
        ],
        [
          1,
          2
        ],
        [
          1,
          3
        ],
        [
          1,
          4
        ],
        [
          1,
          5
        ],
        [
          1,
          6
        ],
        [
          2,
          3
        ],
        [
          2,
          4
        ],
        [
          2,
          5
        ],
        [
          2,
          6
        ],
        [
          3,
          4
        ],
        [
          3,
          5
        ],
        [
          3,
          6
        ],
        [
          4,
          5
        ],
        [
          4,
          6
        ],
        [
          5,
          6
        ]
      ],
      "srg_parameters": [
        21,
        10,
        5,
        4
      ],
      "checks": [
        {
          "eigenvalue": 3,
          "lattice_rank": 6,
          "equation_rank": 15,
          "unimodular_transform_determinant": 1,
          "determinant": 109375,
          "determinant_factorization": {
            "5": 6,
            "7": 1
          },
          "analytic_determinant": 109375,
          "mod4_hypothesis": true,
          "signed_determinant_mod4_if_applicable": 1,
          "primes_dividing_v_not_d": [
            {
              "q": 3,
              "v_q_v": 1,
              "assigned": false
            },
            {
              "q": 7,
              "v_q_v": 1,
              "assigned": true
            }
          ],
          "basis_rows": [
            [
              0,
              -1,
              -5,
              2,
              -2,
              1
            ],
            [
              -1,
              -13,
              -69,
              26,
              -26,
              13
            ],
            [
              -1,
              -14,
              -76,
              29,
              -30,
              13
            ],
            [
              -1,
              -14,
              -75,
              30,
              -29,
              14
            ],
            [
              -1,
              -14,
              -75,
              30,
              -30,
              15
            ],
            [
              -1,
              -14,
              -75,
              28,
              -28,
              14
            ],
            [
              1,
              14,
              76,
              -30,
              30,
              -14
            ],
            [
              1,
              13,
              69,
              -27,
              26,
              -14
            ],
            [
              1,
              13,
              70,
              -26,
              27,
              -13
            ],
            [
              1,
              13,
              70,
              -26,
              26,
              -12
            ],
            [
              1,
              13,
              70,
              -28,
              28,
              -13
            ],
            [
              0,
              1,
              5,
              -3,
              2,
              -2
            ],
            [
              0,
              1,
              6,
              -2,
              3,
              -1
            ],
            [
              0,
              1,
              6,
              -2,
              2,
              0
            ],
            [
              0,
              1,
              6,
              -4,
              4,
              -1
            ],
            [
              0,
              0,
              -1,
              1,
              -1,
              -1
            ],
            [
              0,
              0,
              -1,
              1,
              -2,
              0
            ],
            [
              0,
              0,
              -1,
              -1,
              0,
              -1
            ],
            [
              0,
              0,
              0,
              2,
              -1,
              1
            ],
            [
              0,
              0,
              0,
              0,
              1,
              0
            ],
            [
              0,
              0,
              0,
              0,
              0,
              1
            ]
          ],
          "gram": [
            [
              10,
              135,
              725,
              -280,
              280,
              -135
            ],
            [
              135,
              1830,
              9830,
              -3800,
              3800,
              -1830
            ],
            [
              725,
              9830,
              52810,
              -20415,
              20420,
              -9825
            ],
            [
              -280,
              -3800,
              -20415,
              7910,
              -7905,
              3805
            ],
            [
              280,
              3800,
              20420,
              -7905,
              7910,
              -3800
            ],
            [
              -135,
              -1830,
              -9825,
              3805,
              -3800,
              1840
            ]
          ]
        },
        {
          "eigenvalue": -2,
          "lattice_rank": 14,
          "equation_rank": 7,
          "unimodular_transform_determinant": 1,
          "determinant": 46875,
          "determinant_factorization": {
            "3": 1,
            "5": 6
          },
          "analytic_determinant": 46875,
          "mod4_hypothesis": true,
          "signed_determinant_mod4_if_applicable": 1,
          "primes_dividing_v_not_d": [
            {
              "q": 3,
              "v_q_v": 1,
              "assigned": true
            },
            {
              "q": 7,
              "v_q_v": 1,
              "assigned": false
            }
          ],
          "basis_rows": [
            [
              0,
              0,
              0,
              0,
              1,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0
            ],
            [
              1,
              0,
              0,
              0,
              1,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0
            ],
            [
              -1,
              1,
              2,
              0,
              1,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0
            ],
            [
              0,
              -1,
              -1,
              0,
              -1,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0
            ],
            [
              0,
              0,
              -1,
              1,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0
            ],
            [
              0,
              0,
              0,
              -1,
              -2,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0
            ],
            [
              -1,
              0,
              0,
              0,
              -2,
              0,
              0,
              0,
              1,
              0,
              0,
              0,
              0,
              0
            ],
            [
              1,
              -1,
              -2,
              0,
              -2,
              1,
              0,
              0,
              1,
              0,
              0,
              0,
              0,
              0
            ],
            [
              0,
              1,
              1,
              0,
              1,
              -1,
              1,
              0,
              0,
              0,
              0,
              0,
              0,
              0
            ],
            [
              0,
              0,
              1,
              -1,
              0,
              0,
              -1,
              1,
              0,
              0,
              0,
              0,
              0,
              0
            ],
            [
              0,
              0,
              0,
              1,
              2,
              0,
              0,
              -1,
              -2,
              0,
              0,
              0,
              0,
              0
            ],
            [
              0,
              0,
              0,
              0,
              1,
              -1,
              0,
              0,
              -2,
              0,
              0,
              1,
              0,
              0
            ],
            [
              0,
              0,
              0,
              0,
              0,
              1,
              -1,
              0,
              -1,
              1,
              0,
              1,
              0,
              0
            ],
            [
              0,
              0,
              0,
              0,
              0,
              0,
              1,
              -1,
              0,
              -1,
              1,
              0,
              0,
              0
            ],
            [
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              1,
              2,
              0,
              -1,
              -2,
              0,
              0
            ],
            [
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              1,
              -1,
              0,
              -2,
              0,
              1
            ],
            [
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              1,
              -1,
              -1,
              1,
              1
            ],
            [
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              1,
              2,
              -1,
              -2
            ],
            [
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              1,
              -1,
              -2
            ],
            [
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              1,
              1
            ],
            [
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              1
            ]
          ],
          "gram": [
            [
              4,
              -2,
              -4,
              0,
              0,
              1,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0
            ],
            [
              -2,
              4,
              6,
              0,
              5,
              -2,
              1,
              0,
              -1,
              0,
              0,
              0,
              0,
              0
            ],
            [
              -4,
              6,
              12,
              -2,
              8,
              -3,
              0,
              1,
              -2,
              0,
              0,
              0,
              0,
              0
            ],
            [
              0,
              0,
              -2,
              4,
              4,
              0,
              1,
              -2,
              -2,
              0,
              0,
              0,
              0,
              0
            ],
            [
              0,
              5,
              8,
              4,
              22,
              -4,
              1,
              -2,
              -10,
              0,
              0,
              1,
              0,
              0
            ],
            [
              1,
              -2,
              -3,
              0,
              -4,
              4,
              -2,
              0,
              2,
              1,
              0,
              0,
              0,
              0
            ],
            [
              0,
              1,
              0,
              1,
              1,
              -2,
              4,
              -2,
              1,
              -2,
              1,
              -1,
              0,
              0
            ],
            [
              0,
              0,
              1,
              -2,
              -2,
              0,
              -2,
              4,
              4,
              1,
              -2,
              -2,
              0,
              0
            ],
            [
              0,
              -1,
              -2,
              -2,
              -10,
              2,
              1,
              4,
              16,
              -2,
              -2,
              -9,
              0,
              1
            ],
            [
              0,
              0,
              0,
              0,
              0,
              1,
              -2,
              1,
              -2,
              4,
              -2,
              2,
              1,
              0
            ],
            [
              0,
              0,
              0,
              0,
              0,
              0,
              1,
              -2,
              -2,
              -2,
              4,
              5,
              -2,
              -3
            ],
            [
              0,
              0,
              0,
              0,
              1,
              0,
              -1,
              -2,
              -9,
              2,
              5,
              16,
              -4,
              -9
            ],
            [
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              1,
              -2,
              -4,
              4,
              6
            ],
            [
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              0,
              1,
              0,
              -3,
              -9,
              6,
              12
            ]
          ]
        }
      ]
    }
  ]
}


=== audit_parametres_resultats.json ===
{
  "rows": 4538,
  "existing_rows": 1183,
  "status_counts": {
    "!": 190,
    "2!": 2,
    "-": 466,
    "15!": 1,
    "10!": 2,
    "4!": 2,
    "41!": 1,
    "3854!": 2,
    "180!": 2,
    "32548!": 2,
    "+": 973,
    "28!": 2,
    "78!": 2,
    "167!": 2,
    "?": 2889
  },
  "excluded_existing": [],
  "unknown_canonical": [
    {
      "parameters": [
        169,
        42,
        5,
        12
      ],
      "status": "?",
      "fail": [
        {
          "theta": 3,
          "m": 126,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -10,
          "m": 42,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 3
    },
    {
      "parameters": [
        169,
        70,
        27,
        30
      ],
      "status": "?",
      "fail": [
        {
          "theta": 5,
          "m": 98,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -8,
          "m": 70,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 3
    },
    {
      "parameters": [
        209,
        52,
        15,
        12
      ],
      "status": "?",
      "fail": [
        {
          "theta": 8,
          "m": 76,
          "P": 11,
          "reason": "determinant_mod4"
        },
        {
          "theta": -5,
          "m": 132,
          "P": 19,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 4
    },
    {
      "parameters": [
        235,
        52,
        9,
        12
      ],
      "status": "?",
      "fail": [
        {
          "theta": 5,
          "m": 140,
          "P": 47,
          "reason": "determinant_mod4"
        },
        {
          "theta": -8,
          "m": 94,
          "P": 5,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 4
    },
    {
      "parameters": [
        289,
        90,
        23,
        30
      ],
      "status": "?",
      "fail": [
        {
          "theta": 5,
          "m": 198,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -12,
          "m": 90,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 5
    },
    {
      "parameters": [
        289,
        126,
        53,
        56
      ],
      "status": "?",
      "fail": [
        {
          "theta": 7,
          "m": 162,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -10,
          "m": 126,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 5
    },
    {
      "parameters": [
        375,
        154,
        53,
        70
      ],
      "status": "?",
      "fail": [
        {
          "theta": 4,
          "m": 308,
          "P": 3,
          "reason": "determinant_mod4"
        },
        {
          "theta": -21,
          "m": 66,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 7
    },
    {
      "parameters": [
        375,
        170,
        85,
        70
      ],
      "status": "?",
      "fail": [
        {
          "theta": 20,
          "m": 68,
          "P": 3,
          "reason": "determinant_mod4"
        },
        {
          "theta": -5,
          "m": 306,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 7
    },
    {
      "parameters": [
        391,
        140,
        39,
        56
      ],
      "status": "?",
      "fail": [
        {
          "theta": 4,
          "m": 322,
          "P": 17,
          "reason": "determinant_mod4"
        },
        {
          "theta": -21,
          "m": 68,
          "P": 23,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 7
    },
    {
      "parameters": [
        430,
        39,
        8,
        3
      ],
      "status": "?",
      "fail": [
        {
          "theta": -4,
          "m": 300,
          "P": 43,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 8
    },
    {
      "parameters": [
        441,
        190,
        89,
        76
      ],
      "status": "?",
      "fail": [
        {
          "theta": 19,
          "m": 98,
          "P": 9,
          "reason": "determinant_mod4"
        },
        {
          "theta": -6,
          "m": 342,
          "P": 49,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 8
    },
    {
      "parameters": [
        494,
        85,
        12,
        15
      ],
      "status": "?",
      "fail": [
        {
          "theta": -10,
          "m": 208,
          "P": 19,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 9
    },
    {
      "parameters": [
        507,
        44,
        1,
        4
      ],
      "status": "?",
      "fail": [
        {
          "theta": 5,
          "m": 308,
          "P": 3,
          "reason": "determinant_mod4"
        },
        {
          "theta": -8,
          "m": 198,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 10
    },
    {
      "parameters": [
        616,
        205,
        90,
        57
      ],
      "status": "?",
      "fail": [
        {
          "theta": -4,
          "m": 560,
          "P": 11,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 12
    },
    {
      "parameters": [
        625,
        130,
        15,
        30
      ],
      "status": "?",
      "fail": [
        {
          "theta": 5,
          "m": 494,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -20,
          "m": 130,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 12
    },
    {
      "parameters": [
        625,
        182,
        45,
        56
      ],
      "status": "?",
      "fail": [
        {
          "theta": 7,
          "m": 442,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -18,
          "m": 182,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 12
    },
    {
      "parameters": [
        625,
        234,
        83,
        90
      ],
      "status": "?",
      "fail": [
        {
          "theta": 9,
          "m": 390,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -16,
          "m": 234,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 12
    },
    {
      "parameters": [
        625,
        286,
        129,
        132
      ],
      "status": "?",
      "fail": [
        {
          "theta": 11,
          "m": 338,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -14,
          "m": 286,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 12
    },
    {
      "parameters": [
        630,
        68,
        1,
        8
      ],
      "status": "?",
      "fail": [
        {
          "theta": 5,
          "m": 440,
          "P": 63,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 12
    },
    {
      "parameters": [
        690,
        318,
        121,
        168
      ],
      "status": "?",
      "fail": [
        {
          "theta": 3,
          "m": 644,
          "P": 15,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 13
    },
    {
      "parameters": [
        721,
        220,
        69,
        66
      ],
      "status": "?",
      "fail": [
        {
          "theta": 14,
          "m": 308,
          "P": 103,
          "reason": "determinant_mod4"
        },
        {
          "theta": -11,
          "m": 412,
          "P": 7,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 14
    },
    {
      "parameters": [
        750,
        210,
        55,
        60
      ],
      "status": "?",
      "fail": [
        {
          "theta": -15,
          "m": 308,
          "P": 3,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 14
    },
    {
      "parameters": [
        750,
        214,
        63,
        60
      ],
      "status": "?",
      "fail": [
        {
          "theta": -11,
          "m": 428,
          "P": 3,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 14
    },
    {
      "parameters": [
        780,
        369,
        158,
        189
      ],
      "status": "?",
      "fail": [
        {
          "theta": -36,
          "m": 104,
          "P": 15,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 15
    },
    {
      "parameters": [
        841,
        150,
        11,
        30
      ],
      "status": "?",
      "fail": [
        {
          "theta": 5,
          "m": 690,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -24,
          "m": 150,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 16
    },
    {
      "parameters": [
        841,
        210,
        41,
        56
      ],
      "status": "?",
      "fail": [
        {
          "theta": 7,
          "m": 630,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -22,
          "m": 210,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 16
    },
    {
      "parameters": [
        841,
        270,
        79,
        90
      ],
      "status": "?",
      "fail": [
        {
          "theta": 9,
          "m": 570,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -20,
          "m": 270,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 16
    },
    {
      "parameters": [
        841,
        330,
        125,
        132
      ],
      "status": "?",
      "fail": [
        {
          "theta": 11,
          "m": 510,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -18,
          "m": 330,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 16
    },
    {
      "parameters": [
        841,
        390,
        179,
        182
      ],
      "status": "?",
      "fail": [
        {
          "theta": 13,
          "m": 450,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -16,
          "m": 390,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 16
    },
    {
      "parameters": [
        871,
        348,
        137,
        140
      ],
      "status": "?",
      "fail": [
        {
          "theta": 13,
          "m": 468,
          "P": 67,
          "reason": "determinant_mod4"
        },
        {
          "theta": -16,
          "m": 402,
          "P": 13,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 17
    },
    {
      "parameters": [
        875,
        184,
        33,
        40
      ],
      "status": "?",
      "fail": [
        {
          "theta": 9,
          "m": 552,
          "P": 7,
          "reason": "determinant_mod4"
        },
        {
          "theta": -16,
          "m": 322,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 17
    },
    {
      "parameters": [
        931,
        180,
        29,
        36
      ],
      "status": "?",
      "fail": [
        {
          "theta": 9,
          "m": 588,
          "P": 19,
          "reason": "determinant_mod4"
        },
        {
          "theta": -16,
          "m": 342,
          "P": 49,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 18
    },
    {
      "parameters": [
        1026,
        451,
        212,
        187
      ],
      "status": "?",
      "fail": [
        {
          "theta": -8,
          "m": 836,
          "P": 27,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 20
    },
    {
      "parameters": [
        1148,
        481,
        210,
        195
      ],
      "status": "?",
      "fail": [
        {
          "theta": 26,
          "m": 328,
          "P": 7,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 22
    },
    {
      "parameters": [
        1248,
        203,
        22,
        35
      ],
      "status": "?",
      "fail": [
        {
          "theta": 8,
          "m": 896,
          "P": 39,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 24
    }
  ],
  "missing_claimed": [],
  "extra_unknown": [],
  "excluded_known_nonexistent": [
    {
      "parameters": [
        154,
        51,
        8,
        21
      ],
      "status": "-",
      "fail": [
        {
          "theta": 2,
          "m": 132,
          "P": 7,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 3
    },
    {
      "parameters": [
        289,
        54,
        1,
        12
      ],
      "status": "-",
      "fail": [
        {
          "theta": 3,
          "m": 234,
          "P": 1,
          "reason": "determinant_mod4"
        },
        {
          "theta": -14,
          "m": 54,
          "P": 1,
          "reason": "determinant_mod4"
        }
      ],
      "source_page": 5
    }
  ],
  "arithmetic_parameter_exclusions": 37,
  "table_exclusions_absent_from_independent_enumeration": [],
  "arithmetic_exclusions_not_in_table": []
}

=== FIN DES SORTIES ===
