subroutine processa_elementi
use variabili
implicit none
integer::i,j,j2,j3,j4,vtx1,vtx2,vtx3,vtx4,temp,iele,edge_ele,numnodi,numnodik,nodo1,nodo2,inte,contatore,nfs,nisw
real(4)::deltax,deltay,lato,vetx,vety,omp_get_wtime,xmax,ymax,xmin,ymin,L,nx,ny
real(4),dimension(2)::csi,xnodo
real(4),dimension(:,:),allocatable::temp_vet,temp_vet2
logical::file_exists
write(*,*)' '
write(*,*)'GRID PRE-PROCESSING...'


! ********************************************************
! Ordina lati elemento
write(*,*)'Re-order the edges of the element...'

do i=1,nele_interni
    if(ele(i)%tipo_ele.eq.2)then
        allocate(ele(i)%edge(3,2))
        ele(i)%edge(1,1)=ele(i)%nodi(1)
        ele(i)%edge(1,2)=ele(i)%nodi(2)

        ele(i)%edge(2,1)=ele(i)%nodi(2)
        ele(i)%edge(2,2)=ele(i)%nodi(3)

        ele(i)%edge(3,1)=ele(i)%nodi(3)
        ele(i)%edge(3,2)=ele(i)%nodi(1)
    elseif(ele(i)%tipo_ele.eq.3)then
        allocate(ele(i)%edge(4,2))
        ele(i)%edge(1,1)=ele(i)%nodi(1)
        ele(i)%edge(1,2)=ele(i)%nodi(2)

        ele(i)%edge(2,1)=ele(i)%nodi(2)
        ele(i)%edge(2,2)=ele(i)%nodi(3)

        ele(i)%edge(3,1)=ele(i)%nodi(3)
        ele(i)%edge(3,2)=ele(i)%nodi(4)

        ele(i)%edge(4,1)=ele(i)%nodi(4)
        ele(i)%edge(4,2)=ele(i)%nodi(1)
    end if

end do


! ********************************************************
! Calcola baricentro elemento i-esimo
write(*,*)'Compute element''s center of mass ...'

do i=1,nele_interni

	call compute_center_ele(i)
	

    
end do


! Calcola area elemento i-esimo
write(*,*)'Compute element''s center of mass ...'
do i=1,nele_interni

	call compute_area_ele(i)
	

    
end do


! ********************************************************
!Scopri quanti elementi condividono il nodo i
write(*,*)'Find how many elements share node i...'


    do i=1,nnodi
       nodo(i)%neles=0

       do j=1,nele_interni
           do k=1,ele(j)%nnodi
           if(ele(j)%nodi(k).eq.i) nodo(i)%neles=nodo(i)%neles+1
           end do
       end do
      allocate(nodo(i)%ele(nodo(i)%neles))
    end do


    ! ********************************************************
    !Scopri quali sono gli elementi che condividono il nodo i
    write(*,*)'Find which elements share node i...'

    do i=1,nnodi
       nodo(i)%neles=0

       do j=1,nele_interni
           do k=1,ele(j)%nnodi
           if(ele(j)%nodi(k).eq.i) then
                nodo(i)%neles=nodo(i)%neles+1
                nodo(i)%ele(nodo(i)%neles)=j
           end if
           end do
       end do

    end do

    ! ********************************************************

    




! ********************************************************
!Alloca il vettore per i vicini dell'elemento i-esimo
write(*,*)'Allocate the vector for the neighbours of element i...'
do i=1,nele_interni
    if(ele(i)%tipo_ele.eq.2) then
                                ele(i)%nnghbrs=3
    elseif(ele(i)%tipo_ele.eq.3) then
                                ele(i)%nnghbrs=4
    end if


    ele(i)%nlati=ele(i)%nnghbrs


    allocate(ele(i)%nghbr(ele(i)%nnghbrs))

    ele(i)%nghbr(:)=0   !Tutti i vicini inizializzati a ZERO: resteranno zero le interfacce!
end do
! ********************************************************

!Cerca vicini dell'elemento i-esimo
write(*,*)'Find neighbours of element i...'

!            3
!        4-------3                   3
!        |       |                  / \
!     4  |       |  2             3/   \2
!        |       |                /     \
!        1-------2               1-------2
!            1                       1



    do i=1,nele_interni

    numnodi=min(ele(i)%nnghbrs,ele(i)%nnodi)


        do j=1,numnodi
        if(j.lt.numnodi)  j2=j+1
        if(j.eq.numnodi)  j2=1

        ! j e j2 identificano i nodi all'estremo di un lato
        !ora bisogna cercare quali sono gli altri elementi che possiedono sia j che j2

        do k=1,nele_interni
        if(k.ne.i) then
        numnodik=min(ele(k)%nnghbrs,ele(k)%nnodi)

                    do j3=1,numnodik
                        if(ele(k)%nodi(j3).eq.ele(i)%nodi(j)) then !Ho trovato j
                        do j4=1,numnodik !Cerco j2
                            if(ele(k)%nodi(j4).eq.ele(i)%nodi(j2)) then !Ho trovato j2

                                if(j.eq.1) then
                                    ele(i)%nghbr(1)=k
                                elseif(j.eq.2) then
                                    ele(i)%nghbr(2)=k
                                elseif(j.eq.3) then
                                    ele(i)%nghbr(3)=k
                                elseif(j.eq.4) then
                                    ele(i)%nghbr(4)=k
                                end if


                            end if
                        end do
                        end if
                    end do
        end if


        end do

        end do


    end do
   


 

! **************************************************************
! Scopri quante interfacce ci sono (comprese quelle sui bordi)
write(*,*)'Find the total number of interfaces (boundaries included)...'

ninterf = 0
do i = 1, nele_interni


 if (ele(i)%nnodi.eq.3) then

    if ( ele(i)%nghbr(3) > i  .or. ele(i)%nghbr(3)==0 ) then
     ninterf = ninterf + 1
    endif

    if ( ele(i)%nghbr(1) > i .or. ele(i)%nghbr(1)==0 ) then
     ninterf = ninterf + 1
    endif

    if ( ele(i)%nghbr(2) > i .or. ele(i)%nghbr(2)==0 ) then
     ninterf = ninterf + 1
    endif

elseif (ele(i)%nnodi.eq.4) then

    if ( ele(i)%nghbr(3) > i .or. ele(i)%nghbr(3) ==0 ) then
     ninterf = ninterf + 1
    endif

    if ( ele(i)%nghbr(4) > i .or. ele(i)%nghbr(4) ==0 ) then
     ninterf = ninterf + 1
    endif

    if ( ele(i)%nghbr(1) > i .or. ele(i)%nghbr(1) ==0 ) then
     ninterf = ninterf + 1
    endif

    if ( ele(i)%nghbr(2) > i .or. ele(i)%nghbr(2) ==0 ) then
     ninterf = ninterf + 1
    endif

   endif

  end do

  write(*,*)'NUMBER OF ELEMENTS          = ',nele
  write(*,*)'NUMBER OF INTERNAL ELEMENTS = ',nele_interni
  write(*,*)'NUMBER OF BOUNDARY ELEMENTS = ',nele_bordi
  write(*,*)'NUMBER OF INTERFACES        = ',ninterf
  write(*,*)'NUMBER OF NODES             = ',nnodi


! **********************************************************

  allocate(interf(ninterf))
  interf(:)%entity=0


! **********************************************************
! Scopri elementi e nodi sull'interfaccia i-esima, e1 ed e2 provvisori
write(*,*)'Find nodes and elements related to interface i, e1 and e2 temp...'

  interf(:)%e1=0
  interf(:)%e2=0


ninterf = 0
do i = 1, nele_interni

   vtx1 = ele(i)%nodi(1)
   vtx2 = ele(i)%nodi(2)
   vtx3 = ele(i)%nodi(3)


 if (ele(i)%nnodi.eq.3) then

    if ( ele(i)%nghbr(3) > i  .or. ele(i)%nghbr(3)==0 ) then
     ninterf = ninterf + 1
     interf(ninterf)%e1=i
     interf(ninterf)%e2=ele(i)%nghbr(3)
     interf(ninterf)%nodo1=vtx1
     interf(ninterf)%nodo2=vtx3
    endif

    if ( ele(i)%nghbr(1) > i .or. ele(i)%nghbr(1)==0 ) then
     ninterf = ninterf + 1
     interf(ninterf)%e1=i
     interf(ninterf)%e2=ele(i)%nghbr(1)
     interf(ninterf)%nodo1=vtx1
     interf(ninterf)%nodo2=vtx2
    endif

    if ( ele(i)%nghbr(2) > i .or. ele(i)%nghbr(2)==0 ) then
     ninterf = ninterf + 1
     interf(ninterf)%e1=i
     interf(ninterf)%e2=ele(i)%nghbr(2)
     interf(ninterf)%nodo1=vtx2
     interf(ninterf)%nodo2=vtx3
    endif

   elseif (ele(i)%nnodi.eq.4) then
    vtx4 = ele(i)%nodi(4)

    if ( ele(i)%nghbr(3) > i .or. ele(i)%nghbr(3) ==0 ) then
     ninterf = ninterf + 1
     interf(ninterf)%e1=i
     interf(ninterf)%e2=ele(i)%nghbr(3)
     interf(ninterf)%nodo1=vtx3
     interf(ninterf)%nodo2=vtx4
    endif

    if ( ele(i)%nghbr(4) > i .or. ele(i)%nghbr(4) ==0 ) then
     ninterf = ninterf + 1
     interf(ninterf)%e1=i
     interf(ninterf)%e2=ele(i)%nghbr(4)
     interf(ninterf)%nodo1=vtx1
     interf(ninterf)%nodo2=vtx4
    endif

    if ( ele(i)%nghbr(1) > i .or. ele(i)%nghbr(1) ==0 ) then
     ninterf = ninterf + 1
     interf(ninterf)%e1=i
     interf(ninterf)%e2=ele(i)%nghbr(1)
     interf(ninterf)%nodo1=vtx1
     interf(ninterf)%nodo2=vtx2
    endif

    if ( ele(i)%nghbr(2) > i .or. ele(i)%nghbr(2) ==0 ) then
     ninterf = ninterf + 1
     interf(ninterf)%e1=i
     interf(ninterf)%e2=ele(i)%nghbr(2)
     interf(ninterf)%nodo1=vtx2
     interf(ninterf)%nodo2=vtx3
    endif

   endif


  end do


! **********************************************************



! Scopri corrispondenza interfaccia lato elemento provvisori
write(*,*)'Find corrispondence between interface and edge of the elements (temp)'

do i=1,ninterf

    !Elemento e1
    if(interf(i)%e1.ne.0)then
    do j=1,ele(interf(i)%e1)%nnghbrs

        if (((ele(interf(i)%e1)%edge(j,1).eq.interf(i)%nodo1).and.((ele(interf(i)%e1)%edge(j,2).eq.interf(i)%nodo2))).or.&
        &( (ele(interf(i)%e1)%edge(j,1).eq.interf(i)%nodo2).and.((ele(interf(i)%e1)%edge(j,2).eq.interf(i)%nodo1)))) then
        interf(i)%edge_e1=j
        end if

    end do
    end if

    !Elemento e2
    if(interf(i)%e2.ne.0)then
    do j=1,ele(interf(i)%e2)%nnghbrs

        if (((ele(interf(i)%e2)%edge(j,1).eq.interf(i)%nodo1).and.((ele(interf(i)%e2)%edge(j,2).eq.interf(i)%nodo2))).or.&
        &( (ele(interf(i)%e2)%edge(j,1).eq.interf(i)%nodo2).and.((ele(interf(i)%e2)%edge(j,2).eq.interf(i)%nodo1)))) then
        interf(i)%edge_e2=j
        end if

    end do
    end if

end do





! ********************************************************
!  Calcola normali punti di quadratura interfaccia i-esima e correggi e1 ed e2

do i=1,ninterf


    interf(i)%x0=0.5*(nodo(interf(i)%nodo1)%x+nodo(interf(i)%nodo2)%x)
    


   


        call compute_normal_inte(i)
 



        if(interf(i)%e1.ne.0.and.interf(i)%e2.ne.0) then
            vetx=ele(interf(i)%e1)%x0(1)-interf(i)%x0(1)
            vety=ele(interf(i)%e1)%x0(2)-interf(i)%x0(2)

            if((vetx*interf(i)%normal(1)+vety*interf(i)%normal(2)).gt.0) then
                temp=interf(i)%e1
                interf(i)%e1=interf(i)%e2
                interf(i)%e2=temp
                temp=interf(i)%edge_e1
                interf(i)%edge_e1=interf(i)%edge_e2
                interf(i)%edge_e2=temp
            end if

        elseif(interf(i)%e1.ne.0) then
            vetx=ele(interf(i)%e1)%x0(1)-interf(i)%x0(1)
            vety=ele(interf(i)%e1)%x0(2)-interf(i)%x0(2)

            if((vetx*interf(i)%normal(1)+vety*interf(i)%normal(2)).gt.0) then
                interf(i)%normal=-interf(i)%normal
                
            end if

        elseif(interf(i)%e1.eq.0) then
            vetx=ele(interf(i)%e2)%x0(1)-interf(i)%x0(1)
            vety=ele(interf(i)%e2)%x0(2)-interf(i)%x0(2)

            if((vetx*interf(i)%normal(1)+vety*interf(i)%normal(2)).lt.0) then
                interf(i)%normal=-interf(i)%normal
                interf(i)%e1=interf(i)%e2
                interf(i)%edge_e1=interf(i)%edge_e2
            end if
        end if





end do




! **********************************************************

! Correggi corrispondenza interfaccia lato elemento

do i=1,ninterf

!Elemento e1
if(interf(i)%e1.ne.0)then
do j=1,ele(interf(i)%e1)%nnghbrs

    if (((ele(interf(i)%e1)%edge(j,1).eq.interf(i)%nodo1).and.((ele(interf(i)%e1)%edge(j,2).eq.interf(i)%nodo2))).or.&
    &( (ele(interf(i)%e1)%edge(j,1).eq.interf(i)%nodo2).and.((ele(interf(i)%e1)%edge(j,2).eq.interf(i)%nodo1)))) then
    interf(i)%edge_e1=j
    end if

end do
end if

!Elemento e2
if(interf(i)%e2.ne.0)then
do j=1,ele(interf(i)%e2)%nnghbrs

    if (((ele(interf(i)%e2)%edge(j,1).eq.interf(i)%nodo1).and.((ele(interf(i)%e2)%edge(j,2).eq.interf(i)%nodo2))).or.&
    &( (ele(interf(i)%e2)%edge(j,1).eq.interf(i)%nodo2).and.((ele(interf(i)%e2)%edge(j,2).eq.interf(i)%nodo1)))) then
    interf(i)%edge_e2=j
    end if

end do
end if

end do


! **********************************************************

! Scopri se l'interfaccia i-esima appartiene alle BC

do i=1,ninterf

 do j=1,nele_bordi

    if((interf(i)%nodo1.eq.ele_bordi(j)%nodi(1).and.interf(i)%nodo2.eq.ele_bordi(j)%nodi(2)).or.(interf(i)%nodo1.eq.ele_bordi(j)%nodi(2).and.interf(i)%nodo2.eq.ele_bordi(j)%nodi(1))) then
        interf(i)%entity=ele_bordi(j)%entity
        
    end if



 end do


end do





! ***************************************************************
! Associa BC alle interfacce sui bordi


do j=1,nentity

    if(entita(j)%indx.eq.99) then
    !if(entita(j)%name.eq.("isw")) then
        nisw=0
        do i=1,ninterf
            if(interf(i)%entity.eq.entita(j)%indx)then
                nisw=nisw+1
            end if
        end do
    entita(j)%nmembers=nisw

    allocate(entita(j)%members(entita(j)%nmembers))
    nisw=0
    do i=1,ninterf
        if(interf(i)%entity.eq.entita(j)%indx)then
           nisw=nisw+1
           entita(j)%members(nisw)=i
        end if
    end do


    end if

! ****************************************************



!    if(entita(j)%name.eq.("asw")) then
!        nisw=0
!        do i=1,ninterf
!            if(interf(i)%entity.eq.entita(j)%indx)then
!                nisw=nisw+1
!            end if
!        end do
!    entita(j)%nmembers=nisw

!    allocate(entita(j)%members(entita(j)%nmembers))
!    nisw=0
!    do i=1,ninterf
!        if(interf(i)%entity.eq.entita(j)%indx)then
!           nisw=nisw+1
!           entita(j)%members(nisw)=i
!        end if
!    end do


!    end if

!! ****************************************************



!    if(entita(j)%name.eq.("aasw")) then
!        nisw=0
!        do i=1,ninterf
!            if(interf(i)%entity.eq.entita(j)%indx)then
!                nisw=nisw+1
!            end if
!        end do
!    entita(j)%nmembers=nisw

!    allocate(entita(j)%members(entita(j)%nmembers))
!    nisw=0
!    do i=1,ninterf
!        if(interf(i)%entity.eq.entita(j)%indx)then
!           nisw=nisw+1
!           entita(j)%members(nisw)=i
!        end if
!    end do


!    end if

!! ***********************************************

!   if(entita(j)%name.eq.("fs")) then
!        nfs=0
!        do i=1,ninterf
!            if(interf(i)%entity.eq.entita(j)%indx)then
!                nfs=nfs+1
!            end if
!        end do
!    entita(j)%nmembers=nfs

!    allocate(entita(j)%members(entita(j)%nmembers))
!    nfs=0
!    do i=1,ninterf
!        if(interf(i)%entity.eq.entita(j)%indx)then
!           nfs=nfs+1
!           entita(j)%members(nfs)=i
!        end if
!    end do


!    end if

! ***********************************************

   !if(entita(j)%name.eq.("insub")) then
   if(entita(j)%indx.eq.2) then
        nfs=0
        do i=1,ninterf
            if(interf(i)%entity.eq.entita(j)%indx)then
                nfs=nfs+1
            end if
        end do
    entita(j)%nmembers=nfs

    allocate(entita(j)%members(entita(j)%nmembers))
    nfs=0
    do i=1,ninterf
        if(interf(i)%entity.eq.entita(j)%indx)then
           nfs=nfs+1
           entita(j)%members(nfs)=i
        end if
    end do


    end if

! ***********************************************

   if(entita(j)%indx.eq.40) then
   !if(entita(j)%name.eq.("perio1")) then
        nfs=0
        do i=1,ninterf
            if(interf(i)%entity.eq.entita(j)%indx)then
                nfs=nfs+1
            end if
        end do
    entita(j)%nmembers=nfs

    allocate(entita(j)%members(entita(j)%nmembers))
    nfs=0
    do i=1,ninterf
        if(interf(i)%entity.eq.entita(j)%indx)then
           nfs=nfs+1
           entita(j)%members(nfs)=i
        end if
    end do


    end if

! ***********************************************

   if(entita(j)%indx.eq.50) then
   !if(entita(j)%name.eq.("perio2")) then
        nfs=0
        do i=1,ninterf
            if(interf(i)%entity.eq.entita(j)%indx)then
                nfs=nfs+1
            end if
        end do
    entita(j)%nmembers=nfs

    allocate(entita(j)%members(entita(j)%nmembers))
    nfs=0
    do i=1,ninterf
        if(interf(i)%entity.eq.entita(j)%indx)then
           nfs=nfs+1
           entita(j)%members(nfs)=i
        end if
    end do


    end if

! ***********************************************
 
 

! ***********************************************

   if(entita(j)%indx.eq.3) then
   !if(entita(j)%name.eq.("outsub")) then
        nfs=0
        do i=1,ninterf
            if(interf(i)%entity.eq.entita(j)%indx)then
                nfs=nfs+1
            end if
        end do
    entita(j)%nmembers=nfs

    allocate(entita(j)%members(entita(j)%nmembers))
    nfs=0
    do i=1,ninterf
        if(interf(i)%entity.eq.entita(j)%indx)then
           nfs=nfs+1
           entita(j)%members(nfs)=i
        end if
    end do


    end if


if(entita(j)%name.ne.'inside')write(*,*)entita(j)%name,' = ',entita(j)%nmembers,' members'
end do


write(*,*)'Find ele(i)%interfaces(:)...'




    do i=1,nele_interni
        allocate(ele(i)%interfaces(ele(i)%nnghbrs))
        do j=1,ele(i)%nnghbrs
            nodo1=ele(i)%edge(j,1)
            nodo2=ele(i)%edge(j,2)
            do j2=1,ninterf
                if((nodo1.eq.interf(j2)%nodo1.and.nodo2.eq.interf(j2)%nodo2).or.&
                &(nodo1.eq.interf(j2)%nodo2.and.nodo2.eq.interf(j2)%nodo1))then
                    ele(i)%interfaces(j)=j2
                    if(interf(j2)%e2.eq.(i)) ele(i)%interfaces(j)=-j2
                    exit
                end if
            end do
        end do

    end do

    



write(*,*)'Find ele(i)%interfaces(:)...  OK'


do i=1,ninterf
	call compute_length_inte(i)
end do


write(*,*)'nentity =',nentity

    periodo=-1.d0
    do i=1,nentity
    write(*,*)'entita(i)%indx = ',entita(i)%indx
        if(entita(i)%indx.eq.40) then
        !if(entita(i)%name.eq."perio1") then
            index_perio1=entita(i)%indx
            call prepare_periodicity
            exit
        end if
    end do


write(*,*)'Grid pre-processing........COMPLETED'
write(*,*)' '


call prepar_tecplot

end subroutine







subroutine prepare_periodicity
use Variabili
implicit none
integer::i,j,perio1,perio2,intsopra,intsotto,l,m,mm,elesotto_int,elesotto,ent_outsub,inte
real(8)::ymax,ymin

ymin=huge(1.d0)
ymax=-huge(1.d0)


do i=1,nentity
    !if(entita(i)%name.eq.'perio1') perio1=i
    !if(entita(i)%name.eq.'perio2') perio2=i
    if(entita(i)%indx.eq.40) perio1=i
    if(entita(i)%indx.eq.50) perio2=i    
end do


do i=1,entita(perio1)%nmembers
    intsopra=entita(perio1)%members(i)

 

    do j=1,entita(perio2)%nmembers

        intsotto=entita(perio2)%members(j)


        if((abs(nodo(interf(intsopra)%nodo1)%x(1)-nodo(interf(intsotto)%nodo1)%x(1)).lt.1.d-8.and.abs(nodo(interf(intsopra)%nodo2)%x(1)-nodo(interf(intsotto)%nodo2)%x(1)).lt.1.d-8).or.&
        &(abs(nodo(interf(intsopra)%nodo1)%x(1)-nodo(interf(intsotto)%nodo2)%x(1)).lt.1.d-8.and.abs(nodo(interf(intsopra)%nodo2)%x(1)-nodo(interf(intsotto)%nodo1)%x(1)).lt.1.d-8))then


            do l=1,ele(interf(intsopra)%e1)%nnghbrs
                if(ele(interf(intsopra)%e1)%nghbr(l).eq.0) then
                    ele(interf(intsopra)%e1)%nghbr(l)=interf(intsotto)%e1
                    exit
                end if
            end do

            do l=1,ele(interf(intsotto)%e1)%nnghbrs
                if(ele(interf(intsotto)%e1)%nghbr(l).eq.0) then
                    ele(interf(intsotto)%e1)%nghbr(l)=interf(intsopra)%e1
                    exit
                end if
            end do

            interf(intsopra)%e2=-interf(intsotto)%e1
            interf(intsopra)%edge_e2=interf(intsotto)%edge_e1


            if(abs((interf(intsopra)%x0(1)-interf(intsotto)%x0(1))).lt.1.d-8)then
                interf(intsopra)%int_perio=intsotto
            else
                interf(intsopra)%int_perio=-intsotto
            end if


            !Vai a mettere intsopra (negativa) nell'interfaccia di sotto dell'elemento di sotto

                do mm=1,ele(interf(intsotto)%e1)%nnghbrs
                    
                    if((ele(interf(intsotto)%e1)%interfaces(mm)).eq.intsotto) then
                        elesotto_int=mm
                        exit
                    end if
                end do


 
            ele(interf(intsotto)%e1)%interfaces(elesotto_int)=-intsopra








        end if


    end do

end do

!pause


! ***********************************************************
!Calcola periodo della mesh


periodo=-1.

do i=1,nentity
    !if(entita(i)%name.eq.'outlet') ent_outsub=i
    if(entita(i)%indx.eq.3) ent_outsub=i
end do


do i=1,entita(ent_outsub)%nmembers
    inte=entita(ent_outsub)%members(i)

    if(nodo(interf(inte)%nodo1)%x(2).gt.ymax) ymax=nodo(interf(inte)%nodo1)%x(2)
    if(nodo(interf(inte)%nodo2)%x(2).gt.ymax) ymax=nodo(interf(inte)%nodo2)%x(2)

    if(nodo(interf(inte)%nodo1)%x(2).lt.ymin) ymin=nodo(interf(inte)%nodo1)%x(2)
    if(nodo(interf(inte)%nodo2)%x(2).lt.ymin) ymin=nodo(interf(inte)%nodo2)%x(2)
end do

periodo=ymax-ymin

write(*,*)'Computed period = ',periodo

 

end subroutine

