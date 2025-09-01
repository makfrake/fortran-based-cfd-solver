module strutture
implicit none


type tipo_elemento !tipo di elemento generico (può essere sia 2D che 1D)
    integer:: nnodi,tipo_ele
    integer,dimension(:),pointer::nodi
    integer::entity
end type tipo_elemento


type tipo_nodo
   real(4),dimension(2):: x   !coordiante del nodo
   integer:: nnghbrs,neles 
   integer,dimension(:),allocatable :: nghbr    !list of neighbors
   integer,dimension(:),allocatable :: ele      !list of elements
end type tipo_nodo
! ****************************************************


type tipo_interfaccia
    integer:: e1, e2  !elementi associati(e1 è quello di riferimento per la normale uscente)
    integer::edge_e1,edge_e2 !Lato corrispondente degli elementi separati dall''interfaccia
    integer::nodo1,nodo2 !nodi agli estremi dell'interfaccia
    real(4),dimension(2) ::normal,x0 !normale e centro dell'interfaccia
    real(4)::length !lunghezza dell'interfaccia
       
    integer::entity !0 all''interno, >0 sui bordi
    real(4),dimension(4)::f !flussi che attraversano l'interfaccia
    integer::int_perio

end type tipo_interfaccia


type tipo_elemento_solido
    integer:: nnodi,tipo_ele,nlati !numero nodi, tipo elemento e numero lati
    integer,dimension(:),allocatable::nodi    ! lista dei nodi dell'elemento
    integer,dimension(:,:),allocatable::edge !lista dei lati dell'elemento
    integer:: nnghbrs  !numero dei vicini dell'elemento
    integer,   dimension(:), allocatable :: nghbr    !lista dei vicini dell'elemento
    integer,dimension(:),allocatable::interfaces ! lista delle interfacce dell'elemento
    
    real(4)::area   !area dell'elemento
    real(4),dimension(2)::x0  !coordinate baricentro dell'elemento
    real(4),dimension(4)::ucons !vettore delle grandezze conservative: rho*e,rho,rho*u,rho*v
    real(4)::u,v,a,P,T,S !grandezze primitive
    real(4),dimension(4)::d_dt !vettore delle derivate temporali delle grandezze conservative


end type tipo_elemento_solido
! ****************************************************
type tipo_bordo
    integer:: nnodi,tipo_ele
    integer,dimension(:),pointer::nodi
    integer::entity
end type tipo_bordo
! ****************************************************
type tipo_entity
    integer:: indx,tipo_ele
    character(20)::name
    integer:: nmembers
    integer,   dimension(:), allocatable :: members
    
end type tipo_entity

! ****************************************************





end module

