"# vitber---prosjekt-3" 

Må gjøres

Oppgave 1

  a - Her burde vi også implementere kode for å plotte partikkelens bane slik at vi kan gjenbruke koden senere. Også viktig at funksjonene tar inn t-verdier til tross for at det ikke er nødvendig i Vw-ligningen for å sikre gjenbruk som en mulighet. Her er ikke analyse en veldig stor nødvendighet siden vi ikke har mye å sammenligne med.
  EQUATION 2
  EULER'S METHOD 
	1. Plot (log-log) final global error for each h (10+)
	2. Find max h with error less than 10 metersz

  b - Viktig her også å følge samme retningslinjer som 1a. Her burde vi sette et større fokus på analyse enn tidligere, det er nok en stor del av oppgaven. Den deloppgaven i oppgave 1 som er mest relevant for løsning av oppgaver 2 og 3.
  EQUATION 2
  TRAPEZOID METHOD
	1. Plot (log-log) final global error for each h (10+)
	2. Find max h with error less than 10 metersz
  
  c - Stor oppgave siden vi også må implementere den analytiske løsningen så vi har noe å sjekke opp mot. Dette er i mye større grad en "self-contained" oppgave siden eq1 har lite å gjøre med oppgaver 2 og 3. Analyse blir derfor viktig! Burde understreke at plot skiller seg fra 1a både i diff.ligning (eq2 vs eq1) og løser(Euler vs ETM).
  EQUATION 1
  TRAPEZOID METHOD
	1. Calculate trajectory for 10 timesteps from hundreds to thousands
	2. Plot final global error for each h
	3. Plot trajectory for adequate h
	4. Compare to 1a. Why much smaller step needed here?
  
  d
  VARIABLE TIMESTEP
	1. Apply to scenario from 1c 
	2. Plot trajectory
	3. Plot development of timestep over time
	4. Discuss the result (for every task really)
  
Oppgave 2

  a
	1. X0 given, transport 10 days
	2. Plot trajectory in xy coordinates using different start dates (1. -10. feb)
	3. Compare with others, should be exactly same
	
  b
	1. Plot trajectory on MAP
	2. Plot trajectory for 3 additional positions
  
Oppgave 3

  a - OBS! Her må alt være i numpy arrays og vi kjører vektorisering for å få en effektiv kode
	1. Randomly place Np particles i intiial square (vil ha Np = 100000 ved endelig kjøring)
	2. Plot position of all particles at t = 0, 2, 4, 6, 8, 10
  
  b
	1. Define a grid with 800 meters x 800 meters cells with same coordinates as map.
	2. Calculate concentration by counting the number of particles in each cell
	3. Plot the concentration grid at t = 0, 2, 4, 6, 8, 10. 
	4. Her angir vi fargestyrke proporsjonalt med concentration?
