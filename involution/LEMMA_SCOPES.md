# Portée des lemmes des modèles locaux (structures, 27/09, pour le barreau 1 de l'échelle)

Chaque ligne : ce que le lemme utilise, s'il reste valide en présence de carrés **dis** (n_dis > 0), et en présence de paires
**frères** (Σ ≠ ∅). « Oui » = démontré sans l'hypothèse, relu. Référence du code : `sat_local2.py` (encodage 2) et
`sat_local.py` (encodage 1).

| lemme (flag) | contenu | utilise | dis présents | Σ ≠ ∅ | relu par |
|---|---|---|---|---|---|
| système local (i)–(iii) | Φ bijection sur les paires non-μ ; K₂ = 2J − I − P_π − P_μ − P_{τπ} ; K₃ = 2J − P_π − P_τ − P_{τμ} − P_{μτ} | voisins communs avec sommets de correction (λ = 1, μ = 2), en tout bloc | **oui** | **oui** | structures + C1 (dérivations à la main indépendantes) ; BvLS dégénéré |
| L1, S1, S2 (C1, ponctuels) | aucune racine ne voit (c+,a) et (c−,a) ; S1/S2 sur τ | λ = 1, μ = 2 | **oui** | **oui** | C1, structures |
| `--cuts` | bloc ancré NON frère p, Sp = {p±, μp±} : a₄ = 0, a₃ ≤ 1 − [digone] | (i)–(iii) + couverture ; les blocs frères sont exclus (couverture 0) | **oui** | **oui** | C1 (relu), structures |
| `--digon` | digone ancré non frère {p,q} : k_Z = 2 − \|blocs(z) ∩ {p,q}\| | forme égalité du précédent | **oui** | **oui** | C1, structures |
| `--path` | bloc chemin ancré non frère : exactement un excès | forme égalité | **oui** | **oui** | C1, structures |
| `--twist` | bloc tordu, S₃ = {p*, p°, q*} : identité de membre droit 3 | (i)–(iii) + couverture | **oui** | **oui** | C1 (2812 blocs tordus), structures |
| `--ident[=s]` | identité générale : a₃ + 3a₄ + Σ_{Q∋p}(k−1)⁺ + #{Q∌p : k = 0} = 20 − 10\|S\| + Λ(S) + C_p(S) | (i)–(iii) + couverture | **oui** | **oui** | C1 (balayage des 3308 classes), structures (BvLS : 1200/1200) |
| **SEC** (`--sec`) | si v ∈ Φ(z) ∩ Φ(σz), alors z rencontre bloc(v) | la σ-paire est le partenaire D de l'orbite étoile de v (Lemme 1) ; son carré ne doit pas être dis | **NON** : faux en c dès qu'un carré dis contient une orbite de l'étoile de c | **oui** : les carrés frères ne sont pas dis | structures, C1, C2 |
| L1∘SEC (C1) | aucune sécante ne contient π(centre) | SEC + L1 | **non** (comme SEC) | **oui** | C1, structures |
| `--fullcyc` | bloc plein : équilibre (v) sur les D-cycles de longueur 3–4 | lemme d'équilibre de C2 (note 23 §8) | **oui** (C2) | **oui** (C2) | C2 |
| `--balance` | pour toute orbite S de ⟨π, μ, τ⟩ sur les étiquettes et tout z lointain : \|z ∩ S\| + \|Φ(z) ∩ S\| + \|Φ(σz) ∩ S\| = \|S\|/2 (généralise --digon, --fullcyc et la rigidité des blocs frères) | (i)–(iii) : Σ_z n = 30\|S\|, Σ_z n² = 15\|S\|², donc variance nulle (C2, note 23 §8, redérivé par structures) | **oui** | **oui** | C2, structures |
| deux blocs : couplage, CROSS, `--full` | toute adjacence décrite deux fois ; λ/μ exacts | définitions + λ/μ | **oui** | **oui** | C1 (couplage, CROSS) ; BvLS exact |
| deux blocs : réduction à 60 instances | c de type minimal, d = bloc 0 ; t_c ≤ t_d | T = K₇ (toutes les paires tordues) | — | — | structures ; C1 (parité des pivots) |
| lex-leader (C1) | brisure de symétrie sous Stab_B₆(μ, τ) | correcte pour UNSAT, hors déclaration (C2 exporte --nolex) | — | — | C1, C2 |

**Réponse à la question de C2.** SEC survit à Σ ≠ ∅ tant que n_dis = 0. En un bloc c, SEC ne demande que ceci : le carré D
de chaque orbite étoile NON pivot de c n'est pas dis. Les carrés frères ne concernent que des orbites pivots (le partenaire
est dans l'étoile), ils ne gênent donc pas. D'où :
- **n_dis = 0 (Σ quelconque) : SEC est valide en tout bloc** ; il sert à un énoncé intermédiaire « n_dis = 0 ».
- Avec des carrés dis : SEC reste valide en c si aucun carré dis n'a d'orbite passant par c. C'est une restriction locale,
  inutilisable sans connaître les carrés dis.

**Contrôle des membres droits codés en dur** (`lemma_scope_check.py`) : sur les 156 classes à ≥ 6 pivots (all_ge6.json de C2),
chaque Sp ancré non frère a un slack égal à 1 − [digone], et chaque S₃ tordu un slack égal à 3, selon l'identité générale.
**0 écart.** --cuts, --digon, --path et --twist sont donc sûrs sur les 156, frères compris.

**Blocs frères.** Pour un bloc frère p (μp+ = p−, pivots), S = {p+, p−} a un slack 0 : c'est une rigidité. Toute σ-paire
qui ne passe pas par p contient une étiquette de p, et aucune σ-paire qui passe par p n'en contient. `--ident=1` l'inclut
automatiquement : les 132 blocs frères des 77 classes non pertinentes l'ont, 0 manquant.
