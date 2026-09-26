# Rapport d’Astra, dossier n° 5 (relayé par l’utilisateur le 26/09/2026)

Texte d’Astra, copié tel quel (formules LaTeX d’origine). Relecture de « structures » : voir JOURNAL et NOTE_S24.

> J’ai une exclusion nouvelle à \(k=22\), pour une boule explicitement tordue qui échappe à L23. Elle est vérifiée par deux méthodes indépendantes. Le cas uniforme \(k\le28\) reste ouvert.

Dossier 5 — Une obstruction par sept ombres et deux critères exacts
1. Résultat concret : une boule exclue à \(k=22\)
Statut : vérifié exactement, indépendamment.
On considère les 36 sommets \(A_i,B_i,C_i\), avec \(i\in\mathbb Z/12\mathbb Z\). Les arêtes entre parties sont
\[
A_iB_i,\qquad A_iC_i,\qquad B_iC_{i+1}.
\]Les trois couplages internes sont donnés par les involutions suivantes :
a = [9,4,7,8,1,6,5,2,3,0,11,10]
b = [7,2,1,9,11,8,10,0,5,3,6,4]
c = [7,3,10,1,11,9,8,0,6,5,2,4]
Ce graphe \(G\) vérifie :
- il est cubique et sans triangle ;
- chacune des trois réunions de deux parties induit un \(C_{24}\) ;
- il ne contient donc aucun des carrés donnant un prisme de barreau dans \(L\) ;
- \(\operatorname{rang}_7(3I-B)=35\) ;
- son nombre d’arbres couvrants vérifie\[
  \boxed{\tau(G)\equiv6\pmod7.}
  \]
Ainsi, L23 n’exclut pas \(k=22\) pour cette boule, puisque \(6\) est un non-carré modulo \(7\).
Les sept ombres qui donnent la contradiction
Fixons \(u=A_0\). Les bijections d’ombre imposent notamment sept sommets lointains \(y_v\) dont les voisins dans la partie \(A\) sont
\[
\{A_0,A_v\},\qquad v\in\{1,2,3,4,5,6,8\}.
\]Écrivons \(n_v\in\{0,1\}^{36}\) pour leur ombre complète. Posons
\[
H=3I-B,\quad
g(n,n')=n^\top H^+n',\quad
K=12I-B+2J-B^2-D.
\]Pour chaque \(v\), on énumère les ombres satisfaisant :
1. deux points non partenaires dans chacune des trois parties ;
2. \(Bn\le2\mathbf1-n\), par les bornes \(\lambda,\mu\) ;
3. aucune paire du support de \(n\) sur une entrée nulle de \(K=MM^\top\) ;
4. \(g(n,n)=6\), nécessaire à \(R=0\).
Les nombres de candidats sont :
\[
\begin{array}{c|rrrrrrr}
v&1&2&3&4&5&6&8\\ \hline
\#&350&356&338&356&307&294&345
\end{array}
\]Pour deux ombres retenues, \(R=0\) impose
\[
g(n,n')\in\{2,3\},\qquad
a(n,n')=3-g(n,n')\in\{0,1\}.
\]Il faut aussi
\[
n^\top n'\le2-a(n,n').
\]Le raccord décisif utilise l’étoile de \(u\). Tous les \(y_v\) sont voisins de \(u\). Pour chacun d’eux,
\[
\boxed{\sum_{w\ne v}a(n_v,n_w)\le1-(Bn_v)_u.}
\]En effet, l’arête \(uy_v\) a exactement un voisin commun ; \((Bn_v)_u\) de ces voisins sont déjà dans \(G\).
Aucun choix simultané d’une ombre dans chacun des sept groupes ne satisfait ces conditions.
Vérification complète	Nœuds explorés
Masques de bits, ordre dynamique	27 626
Tables explicites, ordre fixe	41 014


J’ai également reconstruit séparément le graphe, le Green, les candidats et la seconde recherche : mêmes candidats, même réfutation.
Contrôles et portée
Les groupes \(v=1,\ldots,7\), avec \(7\) à la place de \(8\), admettent un choix compatible. Voici un témoin, avec les indices \(A_i=i\), \(B_i=12+i\), \(C_i=24+i\) :
1 : 0 1 15 19 31 34
2 : 0 2 17 18 32 35
3 : 0 3 15 18 26 29
4 : 0 4 18 20 27 33
5 : 0 5 15 20 24 35
6 : 0 6 20 21 25 26
7 : 0 7 16 22 25 33
Ce contrôle positif a été revérifié indépendamment. Il concerne ces contraintes partielles, pas une matrice d’ombres complète.
Conclusion exacte : toute SRG contenant cette boule a \(k\ne22\), donc \(k\ge24\). Toutes les bijections d’ombre admissibles de cette boule sont couvertes par la réfutation.
Cela ne traite pas encore \(k=24,26,28\), ni les autres boules.
2. À \(k=22\), le problème d’ombres est déjà un problème de construction complet
Statut : démontré, relu indépendamment.
Les hypothèses suivantes sont indispensables :
- \(B\) est la boule cubique prescrite, avec \(d(L)=1\) ;
- \(M\) est binaire, avec deux voisins par partie et dix par ligne ;
- l’identité entière\[
  MM^\top=12I-B+2J-B^2-D
  \]est satisfaite ;
- \(A_F\) est binaire, symétrique, sans boucle et 8-régulière.
Alors
\[
\boxed{R=0\pmod7\quad\Longrightarrow\quad
\text{le graphe reconstruit est une SRG exacte en entiers}.}
\]Preuve
Posons \(Z=M^\top H^+M\). Modulo \(7\),
\[
MM^\top=-H^2-D+2J.
\]Les colonnes de \(H^+M\) ont somme nulle dans chacune des trois parties. On obtient
\[
MZ=-HM,\qquad Z^2=J-M^\top M,\qquad ZJ=0.
\]Si \(R=0\), alors \(A_F=3I+3J-Z\), d’où
\[
(B+I)M+MA_F\equiv2J,
\]\[
A_F^2+A_F\equiv12I+2J-M^\top M.
\]Considérons les résidus entiers
\[
F=(B+I)M+MA_F-2J,
\]\[
E=A_F^2+A_F-12I-2J+M^\top M.
\]Chaque entrée de \(F\), et chaque entrée hors diagonale de \(E\), est un multiple de \(7\) supérieur ou égal à \(-2\). Elle est donc non négative. La diagonale de \(E\) est nulle.
Or leurs sommes de lignes valent
\[
F\mathbf1=(40+80-120)\mathbf1=0,
\]\[
E\mathbf1=(64+8-12-120+60)\mathbf1=0.
\]Ainsi \(F=E=0\). Les autres blocs de l’identité SRG sont déjà imposés par la boule. ∎
Conséquence pratique : un témoin complet \(R=0\), avec ces hypothèses, construirait réellement le graphe de Conway. Il ne serait pas seulement un survivant d’une relaxation.
Les seules conditions individuelles sur \(g\) ne suffisent pas : le Gram exact et la régularité exacte doivent rester dans le modèle.
3. Le même phénomène vaut au rang résiduel \(2\)
Statut : démontré, vérifié indépendamment.
Sous les mêmes hypothèses,
\[
\boxed{MR=0,\quad \operatorname{rang}_7R\le2
\quad\Longrightarrow\quad\text{SRG exacte}.}
\]En effet, \(MR=0\) donne
\[
E\equiv R^2.
\]La diagonale de \(E\) étant nulle,
\[
\operatorname{diag}(R^2)=0.
\]Par ailleurs, le moment déjà connu donne \(\operatorname{tr}R=0\).
Écrivons \(R=UCU^\top\), avec \(U\) de plein rang \(r\le2\), et posons \(T=CU^\top U\). Alors
\[
\operatorname{tr}T=\operatorname{tr}(T^2)=0.
\]Cayley–Hamilton en dimension au plus deux donne \(T^2=0\), puis
\[
\operatorname{rang}(R^2)\le1.
\]Une matrice symétrique non nulle de rang un ne peut avoir toute sa diagonale nulle. Donc \(R^2=0\), et le relèvement du §2 s’applique.
Cette extension s’arrête ici pour l’argument algébrique. Au rang quatre, le contre-exemple
\[
R=
\begin{pmatrix}
0&1&4&2\\
1&0&2&4\\
4&2&0&1\\
2&4&1&0
\end{pmatrix}
\oplus uu^\top,\qquad u=(1,2,4)^\top
\]vérifie, sur \(\mathbf F_7\),
\[
\operatorname{rang}R=4,\quad
R\mathbf1=0,\quad
\operatorname{tr}R=0,\quad
\operatorname{diag}(R^2)=0,
\]mais \(R^2\ne0\).
Ce contre-exemple est algébrique ; il ne prétend pas provenir d’ombres admissibles.
4. Un petit échange d’ombres force deux prismes à \(k=22\)
Statut : démontré, relu indépendamment.
Supposons \(k=22\) et une droite générique. Si deux ensembles disjoints de lointains \(S,T\) vérifient
\[
M\mathbf1_S=M\mathbf1_T,\qquad
1\le |S|=|T|\le6,
\]alors nécessairement :
\[
\boxed{|S|=|T|=6,\quad
S,T\text{ induisent deux prismes, sans arête entre eux}.}
\]En particulier, à \(P=0\), aucun tel échange de six ombres au plus n’est permis.
Preuve
Posons \(z=\mathbf1_S-\mathbf1_T\), prolongé par zéro au graphe entier. Puisque \(Mz=0\) et \(R=0\),
\[
(A-3I)z\equiv0\pmod7.
\]Comme chaque signe apparaît au plus six fois, hors du support \(Az=0\) exactement. Sur le support, les seules valeurs possibles sont \(3z\) et \(-4z\).
La décomposition spectrale de \(z\) sépare donc son support en deux vecteurs à coefficients \(0,\pm1\), de valeurs propres \(3\) et \(-4\).
Un tel vecteur non nul de valeur propre \(-4\), avec \(q\) points de chaque signe, exige au moins quatre voisins de signe opposé par point. En comptant les voisins communs entre deux points de même signe,
\[
q\binom42\le2\binom q2,
\]donc \(q\ge7\). Cette composante est impossible ici.
Ainsi \(Az=3z\). Les graphes induits sur \(S\) et \(T\) ont degré minimal au moins trois. Sur au plus six sommets, les bornes \(\lambda\le1,\mu\le2\) forcent exactement le prisme triangulaire à six sommets. Son degré interne est trois, donc aucune arête ne relie \(S\) à \(T\). ∎
Limite : je n’ai pas démontré que les bijections d’ombre doivent contenir un tel échange.
5. Ce qui reste ouvert
- L’exclusion uniforme des complétions de rang résiduel \(4\) ou \(6\).
- L’exclusion uniforme de \(k=22\), au-delà de la boule explicite ci-dessus.
- L’existence forcée d’une droite générique.
- Le cas où les 231 droites sont toutes 7-dégénérées.
Les contrôles positifs \(m=2\) et \(m=11\) conservent leurs identités SRG exactes. Aucun argument probabiliste n’intervient dans les résultats ci-dessus.
Annexe — Vérificateur de l’exclusion par sept ombres
Calcul entier, avec NumPy ; aucune écriture de fichier ni aucun solveur. Le résultat attendu est None, avec 41014 nœuds.
import itertools as it
import numpy as np

fs = [
    [9,4,7,8,1,6,5,2,3,0,11,10],
    [7,2,1,9,11,8,10,0,5,3,6,4],
    [7,3,10,1,11,9,8,0,6,5,2,4],
]
B = np.zeros((36,36), dtype=np.int64)
for t,f in enumerate(fs):
    for i,j in enumerate(f):
        B[12*t+i,12*t+j] = 1
for i in range(12):
    for v,w in [(i,12+i), (i,24+i),
                (12+i,24+(i+1)%12)]:
        B[v,w] = B[w,v] = 1

def inverse7(A):
    n = len(A)
    X = np.column_stack([A % 7, np.eye(n, dtype=np.int64)])
    for j in range(n):
        i = j + int(np.flatnonzero(X[j:,j])[0])
        X[[i,j]] = X[[j,i]]
        X[j] = X[j] * pow(int(X[j,j]), -1, 7) % 7
        for i in range(n):
            if i != j:
                X[i] = (X[i] - X[i,j]*X[j]) % 7
    return X[:,n:]

I = np.eye(36, dtype=np.int64)
J = np.ones((36,36), dtype=np.int64)
H = 3*I-B
Green = (inverse7(H+J)-J) % 7
assert np.array_equal(H @ Green % 7, (I-J) % 7)

K = 12*I-B+2*J-B@B
for t in range(3):
    K[12*t:12*t+12,12*t:12*t+12] -= 1

pairs = []
for t,f in enumerate(fs):
    pairs.append([
        (12*t+i,12*t+j)
        for i in range(12) for j in range(i+1,12)
        if f[i] != j
    ])

groups = [1,2,3,4,5,6,8]
choices = {}
for v in groups:
    pool = []
    for eb in pairs[1]:
        for ec in pairs[2]:
            support = (0,v)+eb+ec
            n = np.zeros(36, dtype=np.int64)
            n[list(support)] = 1
            if np.any(B@n > 2-n):
                continue
            if any(K[i,j] == 0
                   for i,j in it.combinations(support,2)):
                continue
            if int(n@Green@n) % 7 == 6:
                pool.append(n)
    choices[v] = np.array(pool)

assert [len(choices[v]) for v in groups] == [
    350,356,338,356,307,294,345
]

adj, compatible = {}, {}
capacity = {v: 1-ns@B[0] for v,ns in choices.items()}
for v in groups:
    for w in groups:
        if v == w:
            continue
        g = choices[v]@Green@choices[w].T % 7
        overlap = choices[v]@choices[w].T
        adj[v,w] = (g == 2)
        compatible[v,w] = (
            ((g == 2) | (g == 3))
            & (overlap <= 2-adj[v,w])
        )

nodes = 0
def search(selected, domains):
    global nodes
    nodes += 1
    if not domains:
        return selected
    v = min(domains)
    for x in domains[v]:
        chosen = selected+[(v,int(x))]

        def degree(j,y):
            return sum(bool(adj[j,k][y,z])
                       for k,z in chosen if j != k)

        if any(degree(j,y) > int(capacity[j][y])
               for j,y in chosen):
            continue

        next_domains = {}
        possible = True
        for w,ys in domains.items():
            if w == v:
                continue
            keep = np.ones(len(ys), dtype=bool)
            deg = np.zeros(len(ys), dtype=np.int64)
            for j,y in chosen:
                keep &= compatible[j,w][y,ys]
                deg += adj[j,w][y,ys]
                if degree(j,y) == int(capacity[j][y]):
                    keep &= ~adj[j,w][y,ys]
            keep &= deg <= capacity[w][ys]
            remaining = ys[keep]
            if not len(remaining):
                possible = False
                break
            next_domains[w] = remaining

        if possible:
            result = search(chosen, next_domains)
            if result is not None:
                return result
    return None

answer = search([], {
    v: np.arange(len(choices[v])) for v in groups
})
print(answer, nodes)
assert answer is None
assert nodes == 41014
