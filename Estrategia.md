Estratègia actual

- Actualitzar informació de ronda
- Buscar el millor oponent
   - Generar solucions aleatòries
   - Simular contra un dummy
   - Comparar els resultats per quedar-nos amb la millor
   - Generar noves solucions mitjançant un procés evolutiu
   - (repetir els dos passos anteriors fins passat un 20% del temps)
- Buscar la millor solució
   - Agafar l'anterior millor solució i corre-la una posició
   - Generar solucions aleatòries
   - Simular contra un dummy i contra l'oponent trobat
   - Comparar els resultats per quedar-nos amb la millor (maxmin)
   - Generar noves solucions mitjançant un procés evolutiu
   - (repetir els dos passos anteriors fins passat el temps)
- Printar el primer moviment de la millor solució

-----

Estratègia objectiu

- Actualitzar informació de ronda
- Generar múltiples oponents aleatòris
   - Generar solucions aleatòries
   - Simular contra un dummy
   - Comparar els resultats per quedar-nos amb la millor
   - (repetir per tantes solucions finals com volguem generar)
- Buscar les millors solució per cada oponent que comparteixin el primer moviment
   - Generar solucions aleatòries
       - Cada solució té un primer moviment comú i un seguit de moviments específics per a cada oponent
   - Simular contra un dummy i contra tots els oponents
   - Comparar els resultats per quedar-nos amb la millor
       - Si el primer moviment és igual, quedar-nos per cada oponent amb els moviments que donen millor resultat
       - Si el primer moviment és diferent, quedar-nos amb la millor sencera (maxmin)
   - Generar noves solucions mitjançant un procés evolutiu
   - (repetir els dos passos anteriors fins passat el temps)
- Printar el primer moviment de la millor solució
