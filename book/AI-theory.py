#########
# Descrivi in cosa consiste la fase di forward-propagation e back-propagation nell'ambito del training di una rete neurale,
# soffermandoti sull'importanza dell'algoritmo di backpropagation per il calcolo delle derivate parziali della funzione 
# costo rispetto ai pesi di tutti i layer
#---
# Forward-propagation
# È la fase in cui i dati di input passano attraverso la rete neurale, strato per strato, per calcolare l'output finale. 
# Ogni neurone riceve un input, lo moltiplica per i pesi, aggiunge un bias e applica una funzione di attivazione (es. ReLU, sigmoid).

# Backward-propagation
# È la fase in cui si correggono i pesi della rete usando l’errore tra output previsto y^ e reale y
# Si calcolano le derivate della funzione costo rispetto ai pesi usando la regola della catena (derivata composta), 
# propagando l’errore dall’output verso l’input.


#########
# Learning rate adattivo per ogni peso (durante il processo di ottimizzazione) : Adagrad, RMSProp, Adadelta, Adam.
#---
# Durante l’addestramento, il learning rate determina quanto vengono modificati i pesi ad ogni passo. 
# I metodi adattivi regolano automaticamente il learning rate per ogni peso, rendendo l’ottimizzazione più stabile e veloce.

# Adagrad -> Adatta il learning rate a ogni parametro in base alla somma dei quadrati dei gradienti passati.
# vantaggi: utile per dati sparsi, svantaggi: il learning rate può diventare troppo piccolo

# RMSProp -> Modifica Adagrad mantenendo una media mobile dei gradienti.
# vantaggi: risolve il problema del learning rate che si annulla.

# Adadelta -> Variante di RMSProp che non richiede di fissare a priori il learning rate.
# vantaggi: più robusto non ce il bisogno di scegliere iperparametri

# Adam -> 
# vantaggi: adattivo per ogni peso è ottimo pre reti profonde.


#########
# Spiega in modo dettagliato come il learning rate influenza la convergenza di una rete neurale durante il training. 
# Quali sono le conseguenze di un learning rate troppo alto o troppo basso nel training di una rete neurale. 
# Descrivi le principali strategie di  aggiornamentodel  learning rate durante il training di una rete neurale.
#---
# Il learning rate (tasso di apprendimento) è un parametro che controlla quanto i pesi vengono aggiornati a ogni passo della discesa del gradiente.

# Troppo alto -> il modello oscilla intorno al minimo o diverge, quindi non converge, oppure lo fa in modo instabile
# Troppo basso -> il modello converge molto lentamente, potrebbe fermarsi in un minimo locale o impiegare molto tempo a trovare la soluzione


#########
# Descrivi in dettaglio l'algoritmo di discesa del gradiente con momento. Quali sono le motivazioni che hanno portato alla sua 
# introduzione rispetto alla discesa del gradiente standard? Fornisci la formula matematica dell'aggiornamento dei pesi in 
# questo algoritmo e spiega il ruolo del termine di momento
#---
# problema del gradiente con momento: in presenza di valli strette nella superficie della funzione costo, l'algoritmo oscilla (zig-zag) e converge lentamente 
# Motivazioni del momento, è che sommando una parte della direzione del gradiente alla direzione attule, ci si ritrova per 
# accelerare lungo le direzioni corrette e rendere più "morbidi" le oscillazioni

# aggiornamento pesi: W(t + 1) = w(t) + v(t) -> i termini dentro le parentesi sono gli step, non dei parametri, quindi W(t) è il corrente e W(t + 1) il successivo
# il ruolo del termine del momento, controlla quanta parte della direzione precendete viene mantenuta


#########
# Quale è il ruolo del learning rate nella formula di aggiornamento dei pesi mediante gradient descent. Aggiornamento del 
# learning rate programmato (learning rate scheduling) : step decay, decadimento esponenziale, decadimento dipendente dal tempo.
#---
# Il learning rate è un parametro fondamentale nel metodo di discesa del gradiente. 
# Indica quanto grandi sono i passi con cui i pesi della rete vengono aggiornati in ogni iterazione.

# Troppo grande, si rischia di saltare il minimo della funzione di costo e di non far convergere il modello.
# Troppo piccolo, la convergenza sarà molto lenta e il training richiederà troppo tempo.

# Per migliorare la convergenza, il learning rate può essere modificato durante il training secondo alcune strategie
# step decay -> si riduce il learning rate a intervalli fissi
# decadimento esponenziale ->  il learning rate diminuisce gradualmente e in modo esponenziale con il numero di iterazioni
# decadimento dipentende dal tempo -> il learning rate si aggiorna diminuendo in funzione del tempo o delle iterazioni


#########
# Ricavare la formula di aggiornamento dei pesi mediante algoritmo di backpropagation nel caso di una rete MLP formata da un nodo di input, 2 
# layer nascosti ciascuno dei quali costituito da un solo neurone ed un nodo di output.
#---
# Abbiamo 1 nodo di input x, 2 layer nascosti e 1 nodo di output y
# input -> layer 1 -> layer 2 -> output

# Forward pass -> uscita = f(peso * input + bias) [f funzione tipo signoide o ReLU]
# Backward pass -> parte dall'output e torna indietro corrente i pesi, per farlo usa W(nuovo) = W(vecchio) - n * gradiente
# n in questo caso è la velocità di apprendimento

#########
# Descrivi e confronta le seguenti varianti del metodo di discesa del gradiente utilizzate per ottimizzare la loss function durante 
# l'addestramento di una rete neurale per un task di regressione:
#    Discesa del Gradiente (Batch Gradient Descent)
#    Discesa del Gradiente Stocastico (Stochastic Gradient Descent - SGD)
#    Discesa del Gradiente a Mini-Batch (Mini-Batch Gradient Descent)
#---
# discesa del gradiente, calcola il gradiente della funzione costo usando l'intero sateset a ogni iterazione, aggioran i pesi ogni volta per epoca
# vantaggi: discesa stabile e converge in modo regolare, svantaggi: lenta se il dataset è molto grande 

# discesa del gradiente stoccastico, aggiorna i pesi dopo ogni singolo esempio, è molto rumorosa, ogni passo può andare in direzioni diverse
# vantaggi: molto veloce nei dateset grandi, può uscire più facilemte dai monimi locali, svantaggi: convergenza meno stabile, potrebbe impiegare più iterazioni per raggiungere il minimo

# discesa del gradiente e mini-batch, divide il dataset in mini-batch, ognuno viene usato per aggiornare i pesi
# vantaggi: più stabile del precendete, strutta bene l'hardware, è il più usato, svantaggio: va scelta bene la dimensione del mini-batch, può influenzare la qualità dell'addestramento
