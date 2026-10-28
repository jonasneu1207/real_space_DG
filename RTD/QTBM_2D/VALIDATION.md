# QTBM-Validierung, 28.10.2026

Device: `n-RTD-constM-comp`, dx=0.25 nm, W=3 nm, dy=0.25 nm,
periodische Moden n=-2:2, T=300 K. Unveraenderte gemeinsame
Device-/Poisson-Funktionen und unveraenderter DG-Code.

## Physikalische und numerische Tests

Abschliessender Lauf: **14 bestanden, 0 fehlgeschlagen, 0 unvollstaendig**.
MATLAB Code Analyzer: **21 Dateien, 0 Hinweise**. Der separate Main-Einstieg
wurde zusaetzlich im Prepare-Modus ausgefuehrt. Die gruenen Tests bestaetigen
die unten genannten Eigenschaften, keine vollstaendige IV-Konvergenz.

`tests/qtbm2DTest.m` prueft:

- freien Kanal mit T=1 und analytischer diskreter lokaler Zustandsdichte;
- retardierte und evaneszente Kontakte sowie beide Bandkanten;
- neutrale Modenbesetzung, homogene Dichte und Gleichgewichtsstrom;
- Bond-/Landauer-Strom, Dichtepositivitaet und Energieverfeinerung;
- Eichinvarianz bei gemeinsamer Energieverschiebung;
- analytische Energie eines lokal gebundenen Zustands und dessen Besetzung;
- konservative Uebergabe an Poisson und periodische/Dirichlet-Randwerte;
- Zurueckweisung nichtseparabler Potentiale, Modenenergien bei W=10 nm;
- Energiebudget und kontrollierten Abbruch eines nicht konvergierten Sweeps;
- selbstkonsistentes RTD-Gleichgewicht mit Anderson und Broyden.

Die kontinuumsbasierte Kontaktneutralitaet ergibt mu-Ec = 0.00898286 eV.
Die diskrete X-Zustandsdichte weicht im freien Kanal um weniger als 0.1 %
von der Kontinuumsdichte ab; dieser Gitterfehler verschwindet erst mit dx.
Bei festem Potential und 0.3 V stimmen Bond- und Landauer-Strom bis auf
Rundungsfehler ueberein. Halbierte Energieintervalle, vierfach strengere
Quadraturtoleranz und doppelte kz-Quadraturordnung werden im Test verglichen.

## Selbstkonsistenz: bestaetigte Punkte und verbleibende Grenze

Mit Potentialtoleranz 1e-5 V, relativer Dichteaenderung 5e-4 und
Poisson-Residualtoleranz 1e-3 wurden folgende Punkte erreicht. Dies sind
einzelne Testpunkte bzw. Fortsetzungen, kein vollstaendiger 20-mV-Sweep.

| Bias / V | j / (A/m^2), gerundet |
|---:|---:|
| 0 | 0 |
| 0.10 | 1.3927e9 |
| 0.20 | 4.8525e9 |
| 0.22 | 5.8712e9 |
| 0.24 | 6.9716e9 |
| 0.26 | 8.1191e9 |

Bei 0.28 V stagnierten strengere Testlaeufe trotz Schrittverkleinerung und
alternativer Broyden-Kopplung bei etwa 0.6 % relativem Poisson-Residuum
und etwa 1.3e-4 V ungemischter Potentialkorrektur. Diese Werte erfuellen
die strengen Standardkriterien nicht; sie sind keine akzeptierten IV-Punkte.
Eine zehnfach feinere Energieaufteilung mit Quadraturtoleranz 1e-5 und
16 Verfeinerungsstufen veraenderte am zuletzt untersuchten Zustand die
Dichte um 5.8e-8 relativ und den Strom um 7e-12 relativ. Das ist eine
Kontrolle an diesem Zustand, kein allgemeiner Nachweis fuer alle Resonanzen.

Ein separater Vergleichslauf mit den bestehenden DG-Abbruchkriterien
(Potential 5e-4 V, Dichteaenderung 5e-4, Poisson-Residuum 8e-2) erreichte
aus einem bereits berechneten Startpotential:

| Bias / V | j / (A/m^2), gerundet | Status |
|---:|---:|---|
| 0.28 | 9.2547e9 | konvergiert gemaess DG-Kriterien |
| 0.30 | 1.0231e10 | konvergiert gemaess DG-Kriterien |
| 0.32 | NaN | nach 80 Iterationen nicht konvergiert |

Diese lockereren Vergleichskriterien wurden explizit fuer den Test gesetzt.
Die QTBM-Standards bleiben strenger. Ab 0.32 V wurde der Sweep nicht
fortgesetzt; eine Peaklage oder vollstaendige selbstkonsistente Referenzkurve
kann aus diesen Ergebnissen noch nicht verlaesslich angegeben werden.
Die weitergehende Untersuchung der nichtlinearen Konvergenz im Resonanzbereich
ist offen; die ausgefuehrten Tests identifizieren noch keine eindeutige Ursache.

Beispiel zum bewussten Abgleich der Abbruchbedingungen mit DG:

```matlab
result = main2D_rtd_qtbm(struct('bias',0:0.02:0.6, ...
    'potentialTolerance',5e-4,'densityTolerance',5e-4, ...
    'poissonResidualTolerance',8e-2));
% Immer result.sweep.converged pruefen. Keine automatische Fortsetzung
% ueber einen nicht konvergierten Punkt hinweg.
```

Der Vergleich an einem identischen festgehaltenen DG-Potential ist bereits
moeglich und trennt Transportunterschiede von dieser aeusseren Konvergenzfrage.
