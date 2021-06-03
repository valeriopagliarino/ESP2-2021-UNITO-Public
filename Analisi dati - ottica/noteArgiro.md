# Note per le relazioni e le analisi dati di Ottica:

- L'analisi dati dovrà essere fatta "online" cioè contemporanea alla presa dati e caricata insieme ai dati sperimentali, quindi preparata in anticipo rispetto all'esperienza di laboratorio
- Riportare nell'abstract anche il risultato della misura
- I test statistici si possono evitare quando la compatibilità è entro 1 o 2 sigma
- Il limite per la consegna dei dati sperimentali e delle analisi preliminari è estero da "24 ore" a "qualche giorno" da quello che ha detto Argirò
- La consegna va sempre fatta attraverso la cartella condivisa GDrive
- A questo punto forse conviene non utilizzare più Excel, ma conviene scrivere direttamente in lab dentro a files CSV in modo da poter eseguire un'analisi in ROOT già in laboratorio.
- A differenza del modulo di elettrotecnica, in modo più simile ad ESP I, ci saranno diversi esperimenti di misure ripetute. Utilizzare `TH1D` per gli istogrammi e test di Student per il confronto tra le medie. Attenzione alle ipotesi di "uguale varianza" che devono essere soddisfatte.

