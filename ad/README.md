# ad/ : niveau « paires » pour les squelettes de C2 (involution, Σ = ∅)

Inconnues φ_W(P) : le nombre de voisins d'un représentant de l'orbite W parmi les 4 racines de la paire P. Tout ce qui est
imposé est démontré :
- F1 ;
- zéros de pivot : si W pivote en t, φ(t·tiers) = 2 et 0 ailleurs en t ;
- plafonds libres : φ_W(tm) ≤ nombre d'orbites de tm qui ne pivotent pas en t ;
- équilibre (note 23 §8) sur chaque classe D fermée de G_k et sur son complément ;
- partenaire dis : φ ∈ {2,3} ;
- n_s ne peut valoir 1 que si les deux orbites de la paire sont dis (39.3) ;
- symétrie e(P,P′) ;
- dichotomie, identique au budget μ 16.1 : #{Z hors partenaires D : φ_Z(P) = 2} + 2·#{φ_Z(P) = 3} = 2 − n_s(P).

Fichiers :
- `pair_level.py` : modèle AD, structures énumérées ;
- `pair_level_relaxed.py` : modèle AD, avec les structures libres (bloc A, degré 4, graphe dis) en variables de partenaire ;
- `general_model.py` : squelette complet quelconque (ancrées, tordues, mixtes, dd), qui lit le format de C2 (`from_c2_key`) ;
- `ad_enum_dump.py`, `ad_enum_dump2.py` : copie de `ad_enum_all.py` de C2, avec export des classes ;
- `run_all.py`, `run_relaxed.py` : balayages AD ;
- `single_orbit_proof.py`, `proofs_ad2.py` : contradictions à une racine, lisibles à la main (AD2) ;
- `test_mix.py` : feuilles de mélanges de C2 (A énuméré).

Résultats : voir `NOTE_S28_AD_NIVEAU_PAIRES.md` à la racine du labo.
