# Note prof. Bellan 29/03/2021 - Esp. 5 e 6
(Valerio)

### Restrizioni dei fit caratteristica in uscita: ammesse
A basse correnti e basse tensioni il modello del diodo non è completamente corretto (resistenza di giunzione + andamento esponenziale). Si prendono i parametri della curva interpolata con restrizione e si calcola il chi^2 incrementale aggiungendo i punti uno alla volta. <p>

### Fit delle caratteristiche in ingresso
Le caratteristiche in ingresso si possono fittare con l'equazione di Schockley. Non è proprio il caso di plottare 18 serie, ne bastano 3 <p>

### Medie delle Va
Se nel fit di Early le Va che mediamo non sono compatibili tra loro si possono comunque mediare, ma gli errori devono essere la semidifferenza tra valore medio e valore più esterno.
-> Fare dei test di compatibilità tra i valori. `ANOVA` va benissimo.

Il `betaF` non va fittato.

Le Vce si possono considerare tutte compatibili perché non siamo in regione di saturazione e quindi la caratteristica varia poco.

Nei MultiGraphs usare dei markers tutti diversi, linee e colori uguali.

