This repository is an (adapted) copy of the official web site of the open source book entitled *"Teoria degli errori e fondamenti di statistica (introduzione alla fisica sperimentale)"*, written by Prof. Maurizio Loreti.
The official web site is hosted at the URL http://wwwcdf.pd.infn.it/labo/INDEX.html but it is unclear for how long it is going to be available.
With the idea of preserving the work of Prof. Loreti, I saved a copy of the repository in a more durable location.
Here I also report the landing page of the website, appropriately adapted for this new location.

--_Cristiano Fontana_

# Teoria degli errori e fondamenti di statistica
_Maurizio Loreti_

## Contenuto della directory

Questa directory contiene, in forma direttamente stampabile (Portable Document Format), il libro dal titolo **"Teoria degli errori e fondamenti di statistica (introduzione alla fisica sperimentale)"**, scritto negli anni che vanno dal 1987 al 2005 da me (Maurizio Loreti).
La subdirectory [source](./source) ospita inoltre i file sorgente a partire dai quali sono stati generati sia il libro che le tabelle e le illustrazioni in esso contenute.

Tutte queste tabelle ed alcune delle illustrazioni sono state, a loro volta, prodotte da programmi in linguaggio **C** anch'essi presenti sotto questa URL e che si basano sulla libreria **GSL** (acronimo di _GNU Scientific Library_); quest'ultima è una collezione di funzioni scientifiche rilasciate dalla [Free Software Foundation](https://www.fsf.org/)).
Per maggiori informazioni sulla GSL, fate riferimento a [questa URL](http://www.gnu.org/software/gsl/).

La subdirectory [builds](./builds) contiene l'intero libro in un unico file in formato PDF, in due versioni: una versione standard col testo in nero ed una versione con testo ed illustrazioni in marrone scuro (invece che in nero), sempre su fondo bianco.
A parer mio la versione in marrone è più riposante se si vuole leggere il libro sul display di un computer.

## Sorgente

Come anticipato, la subdirectory [source](./source) contiene il sorgente LaTeX usato per comporre il libro e tutti i file accessori.
In questa subdirectory, in particolare, è presente un [Makefile](./source/Makefile) che è servito (tra le altre cose) a generare il PDF.
Le versioni dei singoli file sorgente sono elencate nel file [VERSION](./source/VERSION).

## Stato

Il testo è stato _congelato_ alla data del **1 marzo 2005**; da questo momento in avanti non intendo ampliare il testo, ma solo correggere gli errori che mi verranno segnalati e compiere modifiche "minori": ad esempio è stata aggiunta, in data 13 aprile 2005, una appendice contenente il testo inglese e la traduzione in italiano della licenza GNU GPL sotto la quale (vedi in proposito il prossimo paragrafo) il libro viene distribuito.

## Licenza

Una versione ridotta del libro che si trova in questa directory, e sempre intitolato **"Teoria degli errori e fondamenti di statistica (introduzione alla fisica sperimentale)**, è stata pubblicata dalla casa editrice Decibel-Zanichelli nel gennaio 1998 (ISBN 88-08-09785-4).
Riporto qui di seguito parte di alcuni messaggi di posta elettronica scambiati tra me e Giorgio Villella, titolare della Decibel; il secondo chiede esplicitamente informazioni sulla disponibilità di copie per gli studenti dell'Anno Accademico 2005-2006.

> Date: Tue, 1 Mar 2005 10:56:42 +0100
> From: Giorgio Villella <villella AT tin DOT it>
> To: maurizio DOT loreti AT pd DOT infn DOT it
> Subject: Re: 9785 loreti teoria degli errori
>
> ... ... ...
> Il 27/01/2005 erano presenti nel magazzino centrale della
> Zanichelli 91 copie del suo libro e altre 10 erano nei
> magazzini periferici.  Viste le vendite annuali, da queste
> cifre si può dedurre che molto probabilmente il libro andrà
> esaurito verso la fine di questo anno.
> ... ... ...
> quando il volume sarà esaurito non verrà più ristampato e
> verrà tolto dal catalogo. Gliene darò immediatamente
> comunicazione ufficiale e il contratto tra lei e me sarà
> subito sciolto e lei sarà libero di ristamparlo come vuole e
> con chi vuole o di metterlo in rete, senza aspettare un anno.
> ... ... ...
> Cordiali saluti,
> Giorgio Villella

--

>  Date: Thu, 23 Jun 2005 07:28:43 +0200 (CEST)
>  From: Maurizio Loreti <loreti AT pd DOT infn DOT it>
>  To: Giorgio Villella <villella AT tin DOT it>
>  Subject: Re: 9785 loreti teoria degli errori
>
>  ... ... ...
>  Il dipartimento mi chiede di definire con esattezza le
>  caratteristiche del mio corso di laboratorio, che inizierà nel
>  secondo trimestre del prossimo anno accademico
>  (orientativamente il 9 gennaio 2006); ed io mi trovo in
>  conseguenza nella necessità di dover scegliere tra l'indicare
>  nel bollettino "libro di testo pubblicato dalla
>  Decibel-Zanichelli" oppure "libro di testo messo
>  gratuitamente a disposizione degli studenti".
>
>  Saprebbe, dopo aver controllato le giacenze attuali,
>  consigliarmi in proposito?
>  ... ... ...
> Maurizio Loreti

--

>  Date: Mon, 27 Jun 2005 14:56:59 +0200
>  From: Giorgio Villella <villella AT tin DOT it>
>  To: maurizio DOT loreti AT pd DOT infn DOT it
>  Subject: Re: 9785 loreti teoria degli errori
> 
>  Un accurato conteggio ha stabilito che ci sono in tutto
>  attualmente 140 copie del libro.  Se ne vendono circa 135
>  copie all'anno.  Direi che almeno fino al 30 giugno 2006 resta
>  in catalogo.
> 
>  Appena saranno rimaste poche unità sarà mia cura avvertirla.
> 
>  Giorgio Villella

--

>  Date: Mon, 23 Jan 2006 16:24:13 +0100
>  From: Villella <villella DOT giorgio AT alice DOT it>
>  To: Maurizio Loreti <maurizio DOT loreti AT pd DOT infn DOT it>
>  Subject: Re: 9785 loreti teoria degli errori
>
>  Le comunico ufficialmente che il suo libro è stato eliminato
>  dal catalogo della Zanichelli e che pertanto può distribuire
>  agli studenti gratuitamente la sua versione o anche farne una
>  nuova edizione con un altro editore.
>
>  Cordiali saluti, Giorgio Villella

Dal momento che le copie in possesso della Decibel-Zanichelli sono state esaurite, è volontà precisa dell'autore che il libro che si trova in questa directory diventi **liberamente ridistribuibile** sotto i termini della [GNU General Public License](http://www.gnu.org/licenses/licenses.html) o **GPL**, così come pubblicata dalla [Free Software Foundation](http://www.gnu.org/fsf/fsf.html) (**FSF**).

In breve, la licenza GPL permette che il libro possa essere stampato, modificato e ridistribuito da chiunque ed in qualunque forma, però con alcune restrizioni: ad esempio, la distribuzione deve comprendere i sorgenti (in modo da non impedire a nessun altro di compiere a sua volta modifiche) o, quanto meno, indicare dove essi possano essere reperiti; e (se l'opera è stata modificata) questo fatto deve essere esplicitamente menzionato ed un elenco delle modifiche deve essere incluso nella distribuzione.
In particolare la licenza GPL implica che **il nome dell'autore, Maurizio Loreti, debba essere menzionato esplicitamente**.

Una versione testo della GPL può essere trovata nel file [LICENSE](./LICENSE).

Ultime modifiche:
* `03-Mar-2005`: prima stesura
* `10-Mar-2005`: aggiunta una figura nel cap. 12; cambiato il `Makefile`.
* `13-Apr-2005`: aggiunto il testo della licenza GPL (nell'Appendice F), ed una _Prefazione alla sesta edizione_
* `05-Sep-2005`: aggiunta nel Makefile la path degli include file e delle librerie della GSL; e, nel preambolo, sottolineata la licenza di distribuzione.
* `20-Jan-2006`: corretto un errore nel Makefile (le versioni 2-up non erano corrette) ed aggiunta la versione "marrone".
* `27-Jan-2006`: aggiunti Source.tar.gz e Source.tar.bz2.
* `01-Feb-2006`: aggiunto un contatore web.
* `20-Feb-2006`: errata pag. 102, terza riga (1.7 diventa 2.35; grazie a Enrico Mazzaglia).
* `05-May-2006`: aggiunto, in seconda di copertina, un pointer alla URL.
* `16-May-2006`: corretto un refuso (pag. vii).
* `06-Nov-2006`: aggiunta una precisazione a pagina 194; aggiunto Makefile.html.
* `01-Jan-2008`: sono in pensione!  Arrivederci; un grazie a tutti i miei colleghi, ed un monte di auguri a tutti quelli che sono stati miei studenti o che hanno utilizzato (o che utilizzeranno in futuro) questo libro.
