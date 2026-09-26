# Revue de littérature (26/09) : srg(v,2m,1,2) / PQ(2,m−1,2), les graphes à 9 et 243 sommets, et ce qu'ils disent de 99

Faite par un agent de recherche pour « structures », en lecture seule et depuis les sources primaires sauf mention
« non vérifié ». Le dépôt GitHub conway-99-research n'a pas été utilisé. Des copies locales (scans et PDF arXiv) sont
dans le bloc-notes de la session (`scratchpad/litrev/`) ; elles ne sont pas versées dans ce dépôt.

## Q1. Unicité de srg(243,22,1,2) : OUVERTE
- **Table de Brouwer** (aeb.win.tue.nl/graphs/srg/srgtab201-250.html) : la ligne 243 porte « + » (exemple connu) et non
  « ! » (unique).
- **Brouwer–Van Maldeghem, *Strongly Regular Graphs* (CUP 2022).** Le §10.55 dit seulement « unique rank 3 strongly
  regular graph with parameters (243,22,1,2) ». On y lit aussi : « The only known examples have v = 9 or v = 243 ».
- **Soe Soe Zaw** (arXiv:1907.02800, 2020) : « its uniqueness … remain undecided ».
- **Unicités partielles prouvées :**
  - parmi les graphes de rang 3 (via la classification des groupes de rang 3) ;
  - parmi les graphes de Cayley abéliens (Arasu–Jungnickel–Ma–Pott, JCTA 67 (1994) 116–125) ;
  - parmi les graphes de translation (Calderbank 1982 et Tzanakis–Wolfskill, BvM Th. 7.1.1) : seul le dual du Golay
    ternaire [11,5]₃, de poids 6 et 9, convient.

## Q2. Caractérisations
- **Makhnev, « О сильно регулярных графах с λ=1 », Mat. Zametki 44(5) (1988) 667–672.** Scan lu en entier.
  - La condition (*) dit : deux triangles reliés par au moins deux arêtes le sont par exactement trois. Pour μ = 2,
    cela équivaut à « toute paire de droites sécantes est dans une grille 3×3 », à n₃ = 0 chez Reimbayev, et à P = P_max.
  - **Th. 2 : aucun srg(99,14,1,2) ne vérifie (*).** La preuve utilise m − 6 = 1, donc elle est propre à m = 7 : le
    sous-graphe Λ₀ serait un srg(33,12,1,6), infaisable (f = 180/7).
- **Vérification de l'agent (non publiée).** Dans BvLS, chaque classe de parallélisme (81 droites) est le
  srg(81,20,1,6) de Brouwer–Haemers. Pour m quelconque, l'analogue srg((2m²+1)/3, 2m−2, 1, 6) n'a des multiplicités
  entières que pour m ∈ {2, 11}.
- **Aucun théorème publié du type « PQ(2,t,2) + grilles ⟹ quotient de H(t+1,3) ».** Pour le diamètre ≥ 3, il existe
  des résultats voisins :
  - Egawa 1981 : caractérisation de H(n,q) ;
  - van Dam–Koolen–Tanaka, Th. 12.1 : c_i = i pour i ≤ 3 et a₂ = 2a₁ ⟹ partition complètement régulière de
    H(k/(a₁+1), a₁+2).
- **Bondarenko–Radchenko** (JCTB 2013) : ne couvre ni 99 ni 243.

## Q3. La famille
- **6273 et 494019** : rien au-delà de la faisabilité et de l'exclusion des graphes de Cayley abéliens.
- **Automorphismes de 99 :**
  - Wilbrink 1984 : pas d'automorphisme d'ordre 11, donc pas de sommet-transitivité ;
  - Makhnev–Minakova 2004 : restrictions sur les points fixes ;
  - Behbahani–Lam 2011 ;
  - Crnković–Maksimović 2020 ;
  - Cesarz–Woldar 2025 (arXiv:2308.02978, sans ordinateur) ;
  - Ishida 2026 (arXiv:2606.29183).
  - Combinaison faite par l'agent : **Aut ∈ {1, Z₂, Z₃}**.
- **Groupe critique** (Ducey et al., arXiv:1910.07686) : il est entièrement déterminé par les paramètres.

## Q4. p-rangs (Brouwer–van Eijl, J. Algebraic Combin. 1 (1992) 329–346)
- BvLS : rk₃(A − sI) = 67. « structures » retrouve 67 (`family_invariants.py`).
- Grille 3×3 : rk₃ = 4.
- Rien sur 99. Seul p = 7 y est intéressant, avec rk₇(A − 3I) ≤ 45.

## Q5. Littérature récente sur 99
- **Reimbayev** (arXiv:2409.10620, 2508.03377) : nombre d'hexagones ≥ … avec égalité si et seulement si (*). **Il
  conjecture l'égalité**, ce qui, avec Makhnev, entraînerait l'inexistence de 99.
- **Keramatipour** (arXiv:2604.23037, SAT) : redémontre Makhnev et conjecture qu'il n'y a aucune grille 3×3 à 99.
- **Lou–Murin** (PRIMES 2014) : pas de « grille moins une arête » induite.

## Pourquoi l'unicité de 243 est ouverte (mécanisme)
- La seule méthode éprouvée est celle de Brouwer pour le quasi-hexagone à 729 points : des grilles (Shult–Yanushka),
  puis des cubes H(3,3), puis une structure de F₃-espace, puis le code de Golay.
- Pour un srg de diamètre 2, **rien ne force les grilles**. Une preuve de l'unicité de BvLS devrait établir (*), puis la
  condition de cube, puis revêtir par H(11,3).
- **La première étape, (*), est exactement ce que Makhnev montre impossible à 99.**

## Lecture de « structures »
- La règle que cherche l'utilisateur (vraie pour 9 et 243, impossible à 99) existe depuis 1988 : c'est (*), et
  Makhnev a prouvé qu'elle échoue à 99.
- Le lemme L16 de C1 (Krein, P ≤ 1151 < P_max) en est une autre preuve.
- **La pièce manquante est la conjecture de Reimbayev** (« tout membre de la famille vérifie (*) »). Elle
  entraînerait à la fois l'inexistence de 99 et le premier pas de l'unicité de 243.
- Inversement, **un srg(243,22,1,2) sans (*) serait un nouveau graphe**, le premier membre « tordu » de la famille.
- L'unicité de 243 et l'existence de 99 butent donc sur le même lemme.
