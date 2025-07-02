# Domande e Risposte su Intelligenza Artificiale

## 1. Differenza tra Machine Learning e Deep Learning

* **Machine Learning (ML)**: Algoritmi che apprendono da dati per fare previsioni o decisioni.
* **Deep Learning (DL)**: Sottoinsieme di ML che usa reti neurali con molti strati (deep) per modelli complessi.

## 2. Storia dell’Intelligenza Artificiale

* **Origini**: Anni ‘50, primi tentativi di creare macchine intelligenti.
* **Inverni dell’IA**: Periodi di scarso interesse (anni ’70 e ’90).
* **Rinascita moderna**: Dal 2011 grazie a DL, big data e GPU.

## 3. Paradigma ML vs Programmazione Tradizionale

* **Tradizionale**: Regole scritte a mano.
* **ML**: Il modello impara le regole dai dati.

## 4. Preparazione dei Dati

* **Training Set**: Per addestrare il modello.
* **Validation Set**: Per scegliere i parametri migliori.
* **Test Set**: Per valutare la performance finale.

## 5. Task del Machine Learning

* **Classificazione**: Assegnare etichette (es: spam/non-spam).
* **Regressione**: Prevedere valori numerici (es: prezzo di una casa).
* **Clustering**: Raggruppare dati simili senza etichette (es: segmentazione clienti).

## 6. Supervisionati vs Non Supervisionati

* **Supervisionati**: Usano dati etichettati.
* **Non supervisionati**: Scoprono strutture nei dati non etichettati.

## 7. Neurone Artificiale vs Biologico

* **Neurone biologico**: Riceve segnali, li elabora e li trasmette.
* **Neurone artificiale**: Simula il biologico con input, pesi, funzione di attivazione e output.

## 8. Funzioni di Attivazione

* **RELU**, **Sigmoid**, **Tanh**: Aggiungono non-linearità al modello.

## 9. Percettrone a soglia e limiti

* **Percettrone**: Modello semplice per classificazione lineare.
* **Limite**: Non gestisce problemi non lineari.

## 10. Struttura di una Rete Neurale

* **Input Layer**: Riceve i dati.
* **Hidden Layer**: Elabora i dati (possono essere molti).
* **Output Layer**: Produce il risultato finale.

## 11. Tipi di Reti

* **Feedforward**: Flusso in avanti senza cicli.
* **Ricorrenti (RNN)**: Con memoria, utili per sequenze.

## 12. MultiLayer Perceptron (MLP)

* **MLP**: Rete neurale con più layer nascosti. Apprende rappresentazioni complesse.

## 13. Forward e Backward Propagation

* **Forward Propagation**: Calcola l'output della rete.
* **Backward Propagation**: Aggiusta i pesi usando il gradiente dell’errore.

## 14. Iperparametri di una Rete

* Esempi: learning rate, numero layer, numero neuroni, epoche, batch size.

## 15. Funzione Costo (Loss Function)

* **Regressione**: Errore quadratico medio (MSE).
* **Classificazione**: Cross-entropy (log loss).

## 16. Reti Neurali Convoluzionali (CNN)

* **Parte convoluzionale**: Estrae feature (es: contorni).
* **Pooling**: Riduce la dimensionalità.
* **Parte fully-connected**: Classifica le feature estratte.

## 17. Formula Backpropagation (esempio semplice)

* In una rete MLP con 1 input, 2 nodi nascosti, 1 output:

  $$
  w := w - \eta \cdot \frac{\partial C}{\partial w}
  $$

  dove $\eta$ è il learning rate.

## 18. Ottimizzazione (Gradient Descent)

* **Batch GD**: Usa tutti i dati per aggiornare i pesi (preciso ma lento).
* **SGD**: Un dato alla volta (più veloce ma instabile).
* **Mini-batch GD**: Compromesso tra batch e SGD (il più usato).

## 19. Convergenza del Gradient Descent

* Serve passo piccolo e funzione convessa per garantire convergenza al minimo globale.

## 20. Non Convessità della Funzione Costo

* Le funzioni di costo reali non sono convexe → rischio di minimi locali.

## 21. Learning Rate

* **Troppo grande**: Diverge.
* **Troppo piccolo**: Converge lentamente.

## 22. Momentum

* Aggiunge “inerzia” all’ottimizzazione per superare ostacoli locali.
* Formula:

  $$
  v = \beta v - \eta \nabla C,\quad w = w + v
  $$

## 23. Learning Rate Scheduling

* **Step decay**: Riduce a intervalli regolari.
* **Decadimento esponenziale**: Calo continuo.
* **Time decay**: Calo in base al tempo/epoche.

## 24. Ottimizzatori Adattivi

* **Adagrad**: Riduce il passo per parametri frequenti.
* **RMSProp**: Corregge Adagrad per evitare blocchi.
* **Adadelta**: Variante senza learning rate fisso.
* **Adam**: Combina Momentum e RMSProp → molto usato.
* **Formula Adam** (semplificata):

  $$
  m_t = \beta_1 m_{t-1} + (1-\beta_1)\nabla C  
  $$

  $$
  v_t = \beta_2 v_{t-1} + (1-\beta_2)(\nabla C)^2  
  $$

  $$
  w = w - \eta \cdot \frac{m_t}{\sqrt{v_t} + \epsilon}
  $$
