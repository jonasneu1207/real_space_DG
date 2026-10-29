# Separabler 2D-QTBM-Solver fuer die RTD

Eigenstaendiger MATLAB-Einstieg: `main2D_rtd_qtbm`. Der Transport,
die Kontaktbesetzung, Energieintegration, Dichteprojektion und der
Schroedinger-Poisson-Zyklus sind implementiert. Bestehende DG-Dateien
werden nicht veraendert oder als Transportsolver aufgerufen.

**Validierungsgrenze:** Fester-Potential-Transport und die Kopplung sind
getestet. Eine durchgehend konvergierte selbstkonsistente IV-Kurve bis
0.6 V ist noch nicht bestaetigt. Im Resonanzbereich stagniert die aeussere
Iteration teilweise; die betroffenen IV-Werte werden als NaN markiert.
Konkrete Ergebnisse und verwendete Toleranzen stehen in [VALIDATION.md](VALIDATION.md).

## Start

Nur diesen Ordner zum MATLAB-Pfad hinzufuegen; benoetigte gemeinsame
Ordner werden innerhalb der Funktionen voruebergehend eingebunden.

```matlab
addpath('RTD/QTBM_2D')
eq = main2D_rtd_qtbm(struct('bias',0)); % selbstkonsistentes Gleichgewicht

fixed = main2D_rtd_qtbm(struct( ...
    'simulation','fixed-potential','bias',0:0.02:0.6));

iv = main2D_rtd_qtbm(struct( ...
    'simulation','self-consistent','bias',0:0.02:0.6));

plot(iv.sweep.bias, iv.sweep.currentDensity)
xlabel('Spannung / V'); ylabel('Stromdichte / (A/m^2)')
```

`simulation='prepare'` liefert weiterhin Geometrie, Kontaktmoden,
Startpotential und Poisson-System ohne Transportrechnung.
Es gibt keine automatischen Dateiausgaben, Plots, Workspace-Loeschung
oder Aenderung der DG-Konfiguration.

## Modell und Vergleich mit DG

- `H = Hx + Ey + Ez`: uniforme FD-Diskretisierung in X, konstante effektive
  Massen, Y-homogenes Material und Potential.
- Periodische Y-Moden mit Null-Blochphase: `ky = 2*pi*n/W`,
  `n=-modeOrder:modeOrder`. Die vorhandene Funktion
  `get_PeriodicContactModes_full2D` liefert Energien und Normierung.
  Kein kontinuierliches ky-Integral; kz wird kontinuierlich integriert.
- Kontaktneutralitaet aus derselben mx/mz-Modensummenformel wie DG.
  Fermi-Dirac-Verteilung, Spin-Faktor 2, zusaetzlich `mat.deg_factor`.
  Bei Dirichlet werden die Fermipegel an den verschobenen Bandkanten neu
  bestimmt. Bei Neumann werden sie einmal kalibriert und danach als
  muSource=mu0, muDrain=mu0-Bias festgehalten.
- Exakte retardierte Kontakt-Selbstenergien des diskreten Hx.
  Injektion, Zustandsdichte und Bondstrom verwenden dieselbe Dispersion.
  Es gibt keinen CAP und keine kuenstliche imaginaere Energieverbreiterung.
- `buildDevice2D_rtd`, `eval_Poisson2D_rtd` und
  `solve_Poisson2D_rtd` werden wiederverwendet: Dirichlet 0/Bias oder Neumann in X,
  periodisch in Y. Potentialenergie `Vc-phi` in eV.
- Die nodale QTBM-Dichte wird stueckweise linear rekonstruiert und exakt
  ueber die nodalen Poisson-Kontrollvolumina gemittelt. Halbe X-Randzellen
  und der doppelte periodische Y-Endpunkt werden korrekt beruecksichtigt.
- Der aeussere Zyklus benutzt die gemeinsame approximative Newton-
  Poisson-Aktualisierung, aber eigene Iteration, Mischung und Historie.
  Anderson-Mischung nutzt mehrere vorherige Updates. Schrittbegrenzung
  und Ruecknahme verschlechternder Probeschritte sichern Extrapolationen ab.
  Optional aktualisiert Broyden die inverse Ableitung des physikalischen
  Residuals. Sein Startwert stammt aus demselben Poisson-Operator.
  Er prueft ungemischte Potentialkorrektur, Dichteaenderung und das
  physikalische Poisson-Residuum mit der aktuellen Quantendichte.

Die Randbehandlung des DG-Transportoperators und die offenen QTBM-Kontakte
sind unterschiedliche Diskretisierungen. Gleiche Geometrie, Besetzung und
Poisson-Behandlung erzwingen daher noch keine identischen IV-Kurven.
Die kontinuierliche Kontaktneutralitaet und die diskrete X-Dispersion
unterscheiden sich um einen mit dx verschwindenden Diskretisierungsfehler.

Die Option `poissonXBoundary='neumann'` aktiviert feldfreie X-Raender bei
weiterhin periodischem Y. Newton ist dabei erforderlich; die globale
Ladungsbilanz wird zusaetzlich geprueft. Konfiguration und physikalische
Details fuer beide Solver: [Poisson-Randbedingungen](../POISSON_BOUNDARIES.md).

## Numerische Kontrollen

Adaptives Gauss-Legendre-Paar 8/16, mit Intervallteilung an den Kontakt-
Bandkanten und sin^2-Transformation gegen die Schwellensingularitaet.
`energyStepEV` begrenzt die anfaengliche Intervallbreite; bei zu grosser
geschaetzter Abweichung werden Intervalle halbiert. Bei ueberschrittenem
Budget oder fehlender Konvergenz wird kein Transportresultat akzeptiert.
Die kz-Integration wird gegen doppelte Quadraturordnung geprueft.

Die Fehlerabschaetzung kann extrem schmale, von beiden Regeln verfehlte
Resonanzen nicht ausschliessen. Fuer quantitative IV-Vergleiche
`energyStepEV` halbieren, Toleranz verschaerfen sowie X-Aufloesung und
Modenzahl getrennt auf Konvergenz pruefen. Die Abbruchfehler der Energie-
quadratur sagen nichts ueber die X-Diskretisierung aus.

Bondstroeme an allen X-Verbindungen werden mit dem Landauer-Strom
verglichen. Dichte, Strom und Potential des zurueckgegebenen Endzustands
stammen stets aus derselben Transportauswertung. Nicht konvergierte
Poisson-Punkte stoppen die Fortsetzung; ihr IV-Eintrag ist NaN, ihr
letzter Zustand bleibt unter `points` zur Diagnose erhalten.

### Echte gebundene Zustaende

Zustaende unterhalb beider Kontaktkontinua werden mit den evaneszenten
Kontakt-Selbstenergien bestimmt. Ihre Normierung enthaelt die ausserhalb
des Rechengebiets liegenden Auslaeufer. Im Gleichgewicht gilt die gemeinsame
Fermi-Besetzung; die diskreten Zustaende werden der Streuzustandsdichte
hinzugefuegt und tragen keinen Strom.

Bei unterschiedlichen Fermipegeln legt das kohärente stationaere
Streuproblem allein ihre Besetzung nicht fest. Der Standard
`boundStateOccupation='error'` meldet diesen Fall. Fuer gezielte
Sensitivitaetsrechnungen koennen explizit `'source'`, `'drain'` oder
`'average'` gesetzt werden. Diese Werte sind physikalische Annahmen,
keine numerischen Konvergenzhilfen; die verwendete Annahme und Energien
stehen unter `transport.boundStates`.

## Optionen und Ergebnisse

| Option | Standard | Bedeutung |
|---|---:|---|
| simulation | self-consistent | prepare / fixed-potential / self-consistent |
| bias | 0:0.02:0.6 | Spannungsfolge in V, Reihenfolge bleibt erhalten |
| initialPhi | [] | optionales Startpotential [Nx,Ny] in V fuer den ersten selbstkonsistenten Punkt |
| Ly / dy | 3 / 0.25 | nm, Breite und Material-/Poisson-Y-Raster |
| modeOrder | 2 | periodische Moden -2 bis +2 |
| energyStepEV | 0.001 | maximale anfaengliche Energieintervallbreite |
| energyRelativeTolerance | 2e-4 | integrierte Dichte-/Strom-Fehlerabschaetzung |
| energyTailKT | 24 | oberes Energieende oberhalb mu/Bandkanten |
| kzPoints | 64 | Gauss-Punkte des positiven kz-Halbraums |
| maxEnergyRefinements | 12 | maximale adaptive Halbierungsstufen |
| maxEnergyEvaluations | 500000 | Energiebudget pro Transportauswertung |
| currentTolerance | 1e-6 | relativer Bond-/Landauer-Stromvergleich |
| poissonMethod | newton-raphson | alternativ direct |
| poissonXBoundary | dirichlet | alternativ neumann; Y bleibt periodisch |
| chargeNeutralityTolerance | 1e-3 | globale relative Ladungsabweichung, nur bei Neumann als Abbruchkriterium |
| couplingMethod | anderson | alternativ broyden, aeusserer nichtlinearer Loeser |
| maxIterations / mixing | 150 / 0.3 | aeussere Iteration |
| andersonDepth / maxPotentialStepV | 5 / 0.03 | Mischhistorie und Schrittbegrenzung; Tiefe 0 schaltet Anderson aus |
| potentialTolerance | 5e-3 | ungemischte Potentialkorrektur in V |
| densityTolerance | 5e-3 | relative Dichteaenderung |
| poissonResidualTolerance | 1 | physikalisches Poisson-Residuum |
| boundStateOccupation | error | Besetzungsannahme nur bei echten gebundenen Zustaenden unter Bias |
| verbose | true | Iterationsausgabe |

`sweep.currentDensity`: A/m^2. `sweep.current` und `transport.j`: A/m,
ueber eine Y-Periode integriert (wie die DG-Ausgabe). Positive Stroeme
entsprechen positivem Elektronenfluss von Source nach Drain, wie im
bestehenden DG-Code. Es wird keine z-Dicke angenommen.

Pro Bias: `points{k}.phi` [V], `potentialEV` [eV],
`transport.n` [m^-3, Punktwerte], `nPoisson` [m^-3, Volumenmittel],
`transport.jx` [A/m^2, X-Verbindungen], `transport.spectrum`
(Energie/Transmission), `transport.reservoir` (Fermipegel),
`transport.diagnostics`, `poissonInfo`, `history`, `converged`.

Fuer den direkten Vergleich an einem gespeicherten DG-Potential:

```matlab
opt = initParams_qtbm2D(struct('modeOrder',2,'verbose',false));
contacts = get_contacts_qtbm2D(mat, Vc_minus_phi, opt);
transport = solve_transport_qtbm2D(mat, Vc_minus_phi, contacts, opt);
```

Dabei muss `mat` dieselben Y-homogenen Felder und konstanten Massen wie
das verwendete RTD-Device enthalten. Allgemeines V(X,Y), Modenmischung,
variable Massen, andere Blochphasen und inelastische Streuung sind nicht
Teil dieses separablen Vergleichsfalls.

## Funktionsaufteilung und Tests

- Einstieg: main2D_rtd_qtbm, run_RTD_qtbm2D, initParams_qtbm2D.
- Transport: build_hamiltonian_qtbm2D, get_lead_self_energy_qtbm2D,
  solve_scattering_qtbm2D, solve_transport_qtbm2D.
- Besetzung: get_contacts_qtbm2D, get_occupation_qtbm2D,
  contact_supply_qtbm2D, get_bound_states_qtbm2D.
- Quadratur/Observablen: gauss_rule_qtbm2D, adaptive_energy_qtbm2D,
  integrate_observables_qtbm2D.
- Kopplung: project_density_qtbm2D, prepare_poisson_qtbm2D,
  poisson_step_qtbm2D, solve_self_consistent_qtbm2D.
- Hilfen: bias_ramp_qtbm2D, collect_sweep_qtbm2D.

```matlab
results = runtests('RTD/QTBM_2D/tests/qtbm2DTest.m');
assertSuccess(results)
```

Die Tests vergleichen den freien Kanal mit T=1 und der analytischen
diskreten Zustandsdichte, pruefen die Kontaktneutralitaet, Energie-
verfeinerung, Eichinvarianz, einen analytisch loesbaren gebundenen Zustand,
Ladungserhaltung, Poisson-Randbedingungen und RTD-Selbstkonsistenz.
