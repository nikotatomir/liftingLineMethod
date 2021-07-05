program je napravljen sa 3 klase

Prva klasa se zove "LLTgeometryClass.py" i racuna geomterijske parametre potrebne za proracun (npr. evaluacijske tocke, koordinate podkovicastih vrtloga, itd.)

Druga klasa se zove "LLTsolutionClass.py" gdje se nalazi sam proracun. Tu se racuna cirkulacija podkovicastih vrtloga, te panelne sile i momente, aerodinamicki koeficijenti
i distribucija sila duz krila. Ova klasa je vezana za "LLTgeometryClass", sto znaci da koristi atribute izracunate u "LLTgeometryClass"

Treca klasa se zove "LLTpostProcessingClass.py". Ona daje rezultate proracuna u obliku textfile-a i slika. Ona je vezana za "LLTgeometryClass" i "LLTsolutionClass".

Za pokrenuti kod, koristi se "LLT_run.py"

U sustini, za proracun jedino se trebaju zadati sljedeci parametri u fajlu "LLT_run.py":

PRIMJER:

b= 10.
wakeLenFactor = 40. 
Lam = 0. 
beta = 10. 
cRoot = 1. 
cTip = 1. 
Nspan = int(200) 
rho = 1.
freestreamVel = 1. 
alphaRange = np.linspace(5,5,1) 
cRef = 1. 
HSbreak = 0.25 
typeSpacing = "uniform" # moze biti "uniform" ili "cosine"
typeEvalPt = "center" # moze biti "center" ili "glauert". No ako je typeSpacing "uniform" onda typeEvalPt MORA biti "center"

Oni se nalaze odmah na pocetku koda, i objasnjenje parametra je dano. Kod je moze vrtiti za jedan napadni kut ili za raspon kuteva. 
Na primjer, ako se stavi da je alphaRange = np.linspace(5,10,1), kod ce napraviti proracun za napadni kut od 5 stupnjeva.
Ako se stavi da je alphaRange = np.linspace(0,10,11), kod ce napraviti proracun za sljedeci raspon kuteva [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

Za kraj, ako se proracun izvrti za raspon kuteva, u folderu gdje se nalaze python kodovi ce se pojaviti 2 fajla
1. textfile sa integralnim aerodinamickim velicinama u zavisnosti od kuta
2. graf sa integralnim aerodinamickim velicinama u zavisnosti od kuta

ako se proracun izvrti za jedan kut, u folderu gdje se nalaze python kodovi ce se pojaviti 3 fajla
1. textfile sa aerodinamickim velicinama duz krila + integralne velicine za propisani napadni kut 
2. graf sa aerodinamickim velicinama duz krila za propisani napadni kut
3. graf sa geometrijom krila i podkovicastih vrtloga za propisani napadni kut
