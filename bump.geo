Mesh.MshFileVersion = 2.2;

// Definizione lunghezze caratteristiche

L  = 1;
h  = 0.3;
lc = 0.005;   // Target mesh size

// Definizione punti: entità elementare di gmsh, univocamente definito tramite 3 
//					  coordinate (X,Y,Z) ed lc

Point(1) = {0,0,0,lc};
Point(2) = {L,0,0,lc};
Point(3) = {L,h,0,lc};
Point(4) = {0,h,0,lc};
Point(5) = {L/2,h/10,0,lc};
Point(6) = {L/2-L/8,0,0,lc};
Point(7) = {L/2+L/8,0,0,lc};

// Definizione linee: si ottengono dall'unione dei punti
// Ex: Line(1) = {2,3} -> linea 1 che va dal punto 2 al punto 3

Line(1)   = {2,3};
Line(2)   = {3,4};
Line(3)   = {4,1};
Spline(4) = {1,6,5,7,2};   // Indice (4) perché le spline condividono la stessa 
                           // numerazione delle linee

// Definizione superficie: si crea un loop dalle curve che poi può essere utilizzato 
// per definire la superficie da meshare
// Le curve devono essere coordinate: linee chiuse con verso orario o antiorario
// Nel caso di dominio "bucato" si dovranno definire più line loop

Line Loop(1)     = {1,2,3,4}; 
Plane Surface(1) = {1};  

// Costruzione mesh strutturata

// Comando Transfinite: suddivide la curva considerata in N elementi (sovrascrive la dimensione lc degli elementi)
// Progression: crea elementi la cui dimensione aumenta progressivamente secondo
// il rateo indicato a partire dall'inizio della curva
//		Ex: Progression 1.1 -> l'elemento successivo sarà il 10% più grande 
//          del precedente
// Bump: simile alla funzione Progression, però la crescita parte dal centro
//       della linea (?)

Transfinite Line(1)  = 120; \Using Progression 1;
Transfinite Line(-3) = 120; \Using Progression 1;
Transfinite Line (2) = 400; \Using Bump 1;
Transfinite Line (4) = 400; \Using Bump 1;

// Transfinite Surface: rende la griglia strutturata tramite un algoritmo di 
//                      connessione dei nodi 

Transfinite Surface(1) = {1,2,3,4};

// Recombine Surface: rende la mesh da triangolare a quadrata

//Recombine Surface(1);

// Condizioni al contorno
// (99): parete solida
// (2): ingresso
// (1): uscita

Physical Surface(1) = {1};

Physical Line(99) = {4,2};  // wall
Physical Line(2)  = {3};    // inlet
Physical Line(3)  = {1};    // outlet





