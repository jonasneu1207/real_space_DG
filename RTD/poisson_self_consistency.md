# RTD: Poisson-Kopplung mit periodischem Y

Stand: 17.10.2026.

## Aufruf

```matlab
options = struct;
options.simulation = 'self-consistent';
options.bias = [0, 0.01, 0.02];
options.selfConsistent = struct( ...
    'poissonMethod','newton-raphson', ... % alternativ: 'direct'
    'maxIterations',100, ...
    'mixing',0.1, ...
    'potentialTolerance',1e-5, ...       % Volt, UNGEMISCHTER Schritt
    'densityTolerance',1e-4, ...         % relative Dichteaenderung
    'poissonResidualTolerance',1e-5);    % tatsaechliches Poisson-Residuum
options.dg.full2D_solver = 'gmres';
options.dg.full2D_solverTol = 1e-8;
options.dg.full2D_gpu = false;           % auf GPU-Rechner optional true
[results,mat] = run_RTD_full2D(options);
```

Das Beispiel verwendet die bestehenden Device- und Gitterdefaults; es ist
keine Zusicherung, dass diese bei jedem Bias bereits konvergiert/aufgeloest
sind. Bei RTDs sind Gitter-, Toleranz- und Bias-Schritt-Studien notwendig.
Die bisherigen Warmstarts innerhalb der Schleife und zwischen erfolgreich
konvergierten Biaspunkten bleiben aktiv. Ohne `simulation` bleibt Flatband
der Default. Nur im selbstkonsistenten Runner ist die Standard-Transport-
toleranz jetzt `1e-8` statt `1e-3`; explizite `options.dg`-Werte haben Vorrang.

Direkt mit einem vorbereiteten Materialobjekt:

```matlab
mat.dg.params.full2D_solver = 'gmres';
sc = struct('poissonMethod','direct','mixing',0.1,'maxIterations',100);
results = solve2D_rtd_self_consistent('DG',mat,sc);
```

Wenn `poissonMethod` nicht gesetzt ist, wird `mat.poisson.opt.solve` verwendet,
falls vorhanden, sonst `direct`. Die alte Schreibweise `newton-rhapson` ist
weiterhin akzeptiert. Ein expliziter `poissonStep`-Callback ersetzt den
eingebauten Schritt; dessen Physik liegt dann in Verantwortung des Aufrufers.

## Physikalische Konvention und Review der vorhandenen Funktionen

Die gesuchte Variable ist das elektrostatische Potential `phi` in Volt.
Fuer vollstaendig ionisierte Dotierung und ausschliesslich Elektronen gilt

\[
 \nabla\!\cdot(\epsilon_0\epsilon_r\nabla\phi)
   = q\,(n-N_D+N_A), \qquad E_c[\mathrm{eV}]=E_{c,0}[\mathrm{eV}]-\phi[\mathrm V].
\]

Dies entspricht der ueblichen Konvention `-div(epsilon grad phi)=rho_charge`;
positive Donatorladung erzeugt bei geerdeten Kontakten positives Potential,
Elektronenueberschuss negatives Potential. Siehe [QTCAD: Poisson und
Bandkanten](https://docs.nanoacademic.com/qtcad/theory_spin_fem/poisson/).

| Bestandteil | Befund / Umsetzung |
| --- | --- |
| Duennbesetzter Fuenfpunktoperator | Beibehalten; Spaltenreihenfolge `ix+(iy-1)*Nx` |
| Harmonische Permittivitaet an Flaechen | Beibehalten, fuer stetigen normalen Verschiebungsfluss ohne Flaechenladung |
| Laengen in nm / SI-Operator | Beibehalten: Faktoren `epsilon0*1e18/dx^2` bzw. `dy^2` |
| Bisheriger Quellterm | Vorzeichen korrigiert passend zu `Vc-phi`; vorher entsprach er nicht diesem elektrostatischen Potential |
| Dotierung | `Nd/Na` sind m^-3, `N_d/N_a` cm^-3. Nur Legacy-Felder werden mit 1e6 umgerechnet; doppelt vorhandene Felder muessen uebereinstimmen |
| Neumann-Endpunkte | Halbe Knoten-Kontrollvolumina beruecksichtigt (Faktor 2 zum vorhandenen Nachbarn) |
| Newton-Jacobian | Lokale Boltzmann-Ableitung beibehalten, Konvention korrigiert und als Naeherung kenntlich gemacht |
| Funktionsnamen | Stimmen mit den Dateinamen ueberein; die entsprechenden manuellen Korrekturen im MATLAB-Ordner bleiben erhalten |

Es werden keine Loecher, unvollstaendige Ionisation, Austausch-Korrelation
oder Ladungen an Materialgrenzen hinzugefuegt. Permittivitaetsspruenge liegen
im nodalen Stencil effektiv zwischen den benachbarten Materialpunkten;
die genaue Grenzflaechenlage muss mit der Gitterverfeinerung untersucht werden.

## Randbedingungen und Ladungsuebertragung

An den vollen X-Kontaktkanten gilt `phi(source)=0`, `phi(drain)=bias`.
Das setzt die bisherige RTD-Referenz mit gleichen Kontaktmaterialien voraus;
ungleiche Kontakt-Arbeitsfunktionen sind nicht automatisch enthalten.

Der Builder enthaelt bereits beide Endpunkte `y=0` und `y=W`. Fuer periodisches
Y werden daher nur `Ny-1` unabhaengige Knoten verwendet. Die erste Zeile ist
mit den Nachbarn `iy=2` und `iy=Ny-1` gekoppelt; die letzte Spalte wird durch
`phi(:,end)-phi(:,1)=0` gebunden. Intern wird diese Duplikatgleichung eliminiert.
Weder Breite noch Ladung werden doppelt gezaehlt. Materialdaten und Randwerte
muessen an den identifizierten Endpunkten uebereinstimmen.

`eval_Poisson2D_rtd` folgt standardmaessig der periodischen DG-Y-Option.
Standalone laesst sich `mat.poisson.yBoundary='periodic'` oder `'neumann'`
setzen; der gekoppelte Runner lehnt widerspruechliche Y-Periodizitaeten ab.
Die gepolsterten alten Knotentabellen bleiben erhalten. Neu ist Knotentyp 3
fuer die doppelte letzte Y-Spalte. Der Loeser versteht gepolsterte und
ungepolsterte Tabellen.

Die Darstellung fuer Plots `transport.n` wird **nicht** als Poisson-Ladung
interpoliert. Stattdessen integriert `project_DG_density_RTD` die rohen
DG-Polynome `transport.nDG` elementweise in beiden Richtungen ueber die
Poisson-Kontrollvolumina. Doppelte DG-Knoten bleiben einseitige Elementwerte.
Am periodischen Rand werden die beiden halben Kontrollvolumina vereinigt.
Die resultierende Volumendichte erhaelt das Integral des diskreten DG-Felds
(bis auf Rundung). Die Dotierung wird als Material-Kontrollvolumenwert benutzt.

`results.poisson.chargeTransfer` enthaelt Integral und Erhaltungsfehler,
`results.poisson.electronDensity` die tatsaechlich verwendete Poisson-Dichte.
Diese ist eine Kontrollvolumen-Dichte, keine neue Punktabtastung.
Negative Kontrollvolumen-Dichten werden **nicht abgeschnitten**, sondern mit
`RTD:NegativeElectronDensity` abgewiesen. Sie koennen einen unzureichend
aufgeloesten bzw. nicht positivitaetserhaltenden Transportlauf anzeigen.
Poisson benoetigt m^-3; die Kz-Integration des Kontaktreservoirs muss aktiv sein.

## Was direct und newton-raphson hier bedeuten

`direct` loest die lineare Poisson-Gleichung mit festgehaltener aktueller
Quantentransportdichte. Danach mischt die aeussere Schleife das Potential.

`newton-raphson` fuehrt pro Transportiteration **einen approximativen
Newton-Schritt** aus. Fuer

\[
 F(\phi)=M\phi+q(N_D-N_A-n),\quad
 J_{\mathrm{approx}}=M-\operatorname{diag}(q^2 n/(k_B T))
\]

wird `J*deltaPhi=-F` mit homogenen Korrekturen an den Dirichlet-Kontakten
geloest. Das setzt lokal `dn/dphi=q*n/(kB*T)` als Stabilisierung an. Die
tatsaechliche Dichte wird anschliessend wieder durch DG berechnet, nicht
durch eine Boltzmann-Dichte ersetzt. Insbesondere in einer entarteten RTD
nahe einer Resonanz ist dies **kein exakter Newton-Raphson-Schritt fuer das
vollstaendige gekoppelte Quantensystem**. Solche approximativen
Quantentransport/Poisson-Iterationen sind von einem vollstaendigen Newton-
Verfahren zu unterscheiden; siehe [Pourfath: self-consistent iteration
scheme](https://iue.tuwien.ac.at/phd/pourfath/node89.html).

Beide Varianten haben denselben angestrebten diskreten Fixpunkt. Dass sie
zu demselben Zustand konvergieren, ist bei mehreren moeglichen RTD-Zweigen
nicht garantiert. Eine Modifikation des physikalischen Ladungsmodells findet
durch den approximativen Jacobian allein nicht statt.

Die eingebaute Kopplung verwendet ausschliesslich `selfConsistent.mixing`;
`mat.poisson.opt.alpha` wird dort auf 1 gesetzt, um doppelte Daempfung zu
vermeiden. Beim Standalone-Aufruf bleibt `opt.alpha` verfuegbar.

## Konvergenz, Diagnostik und verbleibende Grenzen

Aeussere Konvergenz erfordert gleichzeitig:

1. konvergierten Transport mit echtem GMRES-Residuum,
2. kleinen ungemischten Potentialschritt,
3. kleine Dichteaenderung,
4. kleines Poisson-Residuum **am aktuellen Potential mit der aktuellen
   DG-Dichte**, nicht nur das Residuum der linearen Newton-Korrektur.

Das relative Poisson-Residuum ist die Maximumsnorm des freien Residualvektors,
geteilt durch das Maximum der Normen von Elektronenladungsdichte,
Netto-Dotierungsladung und `M*phi` (kleiner Floor `q`). Es wird vor und nach
dem Poisson-Schritt ausgegeben. `physicalRelresAfter` haelt die Eingabedichte
fest und ist daher nach einer Newton-Korrektur noch kein neues gekoppeltes
Residuum. `info.converged` des Einzelschritts bedeutet nur, dass dessen
lineares System geloest wurde. Massgeblich fuer den gesamten Lauf ist
`results.selfConsistency.converged`.

Historie: `results.selfConsistency.history{biasIndex}` mit
`poissonPhysicalRelres`, `potentialDefect`, `densityChange` und den
GMRES-Warmstartdaten. Bei Iterationslimit bleibt das letzte zusammengehoerige
Paar aus Potential und Transport erhalten; der naechste Bias wird nicht
begonnen. Nach einem fehlgeschlagenen Transport werden keine veralteten
Poisson-Diagnosen als aktuell ausgegeben.

Poisson und die Ladungsuebertragung laufen auf der CPU. Der Transport kann
weiter auf der GPU laufen; seine Observablen liegen ohnehin gesammelt auf
der CPU. Es wurde kein GPU-Poisson-Loeser implementiert.

Validierung und offene Punkte: analytische Ladungs-/Bias-Loesungen,
periodische Fourier-Loesung mit Gitterkonvergenz, Fluss am Dielektrikumssprung,
konservative DG-Projektion, Newton-Vorzeichen und Kopplungstests. Der kleine
gekoppelte CPU-Test verwendet eine schwach dotierte RTD mit niedrigen
Barrieren, nicht die Produktionsparameter. Grobe Laeufe mit den originalen
Barrieren erzeugten bereits vor der Poisson-Kopplung negative Dichten.
Die volle RTD benoetigt daher weiter Aufloesungs-, Positivitaets-, Strom-
erhaltungs- und Bias-Konvergenztests; GPU-Ausfuehrung bleibt separat zu testen.

### Messergebnisse der kleinen CPU-Validierung

`rtdPoissonTest`: 24 Tests bestanden (23 gemeinsam, der ergaenzte Runner-
Abbruchtest separat). Der schwach dotierte Low-Barrier-Test
konvergierte mit beiden Methoden nach 21 aeusseren Iterationen:

| Verfahren | Relatives gekoppeltes Poisson-Residuum |
| --- | --- |
| direct | 8.728e-7 |
| newton-raphson | 9.341e-7 |

Die maximale Potentialdifferenz zwischen beiden Ergebnissen betrug
6.667e-12 V. Dies ist eine Funktions-/Kopplungsvalidierung, kein Nachweis
der Konvergenz oder physikalischen Genauigkeit der Produktions-RTD.

Zusaetzlich: 155 bestehende CPU-Regressionspruefungen bestanden, keine
Fehlschlaege; 6 GPU-Tests mangels unterstuetzter GPU uebersprungen.
Geprueft wurden Warmstarts, Schleifenabbruch/Biasfortsetzung, periodische und
Hard-Wall-Raender, variable Masse, Kontaktmoden und der Vorkonditionierer.
Der zusaetzliche Runner-Test prueft, dass sehr kleine Potentialmischung und
lockere Schritt-/Dichtetoleranzen ein zu grosses echtes Poisson-Residuum
nicht als Konvergenz durchlassen. MATLAB Code Analyzer: keine Befunde in
den neun hinzugefuegten/geaenderten MATLAB-Dateien dieses Poisson-Schritts.
