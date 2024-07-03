# SOR iteracija za razpršene matrike

V tej nalogi sem implementiral reševanje linearnega sistema enačb z uporabo SOR iteracije za razpršene matrike. Na začetku sem definiral tip RazprsenaMatrika, ki diagonalno dominantno razpršeno matriko zaradi prostorskih zahtev hrani v dveh matrikah. Za ta tip sem nato definiral metode, ki jih bom uporabljal v funkciji za izračun SOR iteracije, in sicer metode za indeksiranje (getindex, setindex, firstindex, lastindex) in množenje matrike z desne z vektorjem. Nato sem definiral funkcijo, ki rešuje sistem Ax=b s SOR iteracijo, kjer je A diagonalno dominantna razpršena matrika. Ta funkcija vrne rešitev sistema in število iteracij, ki so bile potrebne za dosego rešitve.

Na koncu sem še to metodo uporabil za vožitev grafa v ravnino. Definiral sem funkcijo, ki za dani graf fiksira vsako drugo vozlišče in nato uporabi funkcijo za SOR iteracijo, ki konstruira vožitev grafa v ravnino, tako da ostala vozlišča zavzamejo ravnovesno lego med fiksiranimi vozlišči.

Za naš primer smo vzeli preprost graf, sestavljen iz 5 vozlišč, kjer so povezave v grafu naslednje: (1, 2), (2, 3), (3, 4), (4, 1), (1, 3), (1, 5), (4, 5). Za ta graf smo preizkusili nekaj vrednosti za parameter omega in pogledali, koliko iteracij je potrebnih, da konvergira SOR. Na sliki lahko vidimo odvisnost hitrosti konvergence od izbire parametra omega:
