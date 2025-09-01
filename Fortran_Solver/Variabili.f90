Module Variabili
use strutture
Implicit none
save

		Integer:: k,kinf,kout,kfinal ! variabili intere per gestire i passi temporali e la frequenza di output
		real(4) :: gam,ga,gb,gc,gd,ge,gf,gg,gh,gi,gj !rapporto dei valori specifici e variabili correlate
		real(4)::time,dt!tempo e passo temporale
		real(4)::CFL !Numero di CFL
		real(4)::periodo !passo lungo y in caso di problema periodico
		character(100)::mesh_file ! Nome del file mesh in formato Gmsh 2
		real(4)::ttotin,ptotin,alpha,pexit,machin ! Condizioni al contorno
		integer::nnodi,ninterf,nele,nele_interni,nele_bordi,nentity,ndim,neqs ! numero delle varie entità geometriche, dimensioni ed equazioni
		real(4),dimension(4)::norm2_residuals,norminf_residuals
		
		
		type(tipo_nodo),dimension(:), allocatable  :: nodo  ! array dei nodi 
		type(tipo_interfaccia),dimension(:), allocatable  :: interf !array delle interfacce
		type(tipo_elemento),dimension(:), allocatable  :: ele_all !array che contiene tutti gli elementi (sia quelli 2D che quelli 1D)
		type(tipo_elemento_solido),dimension(:), allocatable  :: ele ! array che contiene gli elmenti 2D
		type(tipo_bordo),dimension(:), allocatable  :: ele_bordi ! array che contiene gli elementi 1D usati sul bordo
		
		type(tipo_entity),dimension(6)  :: entita !array con le varie entità fisiche

		integer,dimension(:),allocatable::nodi_vis,nodi_agg_vis !array indici necessari per output in formato Tecplot
		integer::nnodi_vis,index_perio1

End module





