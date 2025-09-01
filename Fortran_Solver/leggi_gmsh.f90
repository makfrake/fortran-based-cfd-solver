subroutine leggi_gmsh
use strutture
use variabili
implicit none

integer::i,j,ntags,jj,tipo_ele,temp
integer,dimension(20)::tempvet
integer,dimension(11)::leggiint
character::car
character(300)::stringa
logical::trovato
real(4)::machinf2,temp_real

 

100 format(A)



write(*,*)'Opening mesh file: ',mesh_file



open(unit=1,file=mesh_file)

read(1,*)car
read(1,*)car
read(1,*)car
read(1,*)car

read(1,*)nnodi

allocate(nodo(nnodi))

do i=1,nnodi
    read(1,*)j,nodo(i)%x(1),nodo(i)%x(2),temp_real
end do

read(1,*)car
write(*,*)'Reading nodes...               OK'

! *********************************************************
! Lettura elementi
read(1,*)car
read(1,*)nele

nentity=0

nele_bordi=0
nele_interni=0
allocate(ele_all(nele))

do i=1,nele

    read(1,100)stringa

    read(stringa,*)j,tipo_ele,ntags

    if(tipo_ele.eq.1) then        !SEGMENTO
        nele_bordi=nele_bordi+1
        ele_all(i)%nnodi=2
        allocate(ele_all(i)%nodi(ele_all(i)%nnodi))    
    elseif(tipo_ele.eq.2) then     !TRIANGOLO LINEARE
        nele_interni=nele_interni+1
        ele_all(i)%nnodi=3
        allocate(ele_all(i)%nodi(ele_all(i)%nnodi))
    elseif(tipo_ele.eq.3) then     !QUADRILATERO LINEARE
        nele_interni=nele_interni+1
        ele_all(i)%nnodi=4
        allocate(ele_all(i)%nodi(ele_all(i)%nnodi))
    else
        write(*,*)'Element',i,'UNKNOWN!!!'
        stop
    end if


   read(stringa,*)j,ele_all(i)%tipo_ele,ntags,ele_all(i)%entity,leggiint(2:ntags),ele_all(i)%nodi(1:ele_all(i)%nnodi)
    


end do


! Finora ho messo tutti gli elementi in ele_al: ora divido gli interni dai bordi


allocate(ele(nele_interni))
allocate(ele_bordi(nele_bordi))


nele_interni=0
nele_bordi=0

do i=1,nele

    if(ele_all(i)%tipo_ele.eq.1) then
        nele_bordi=nele_bordi+1
        ele_bordi(nele_bordi)%nnodi=ele_all(i)%nnodi
        ele_bordi(nele_bordi)%entity=ele_all(i)%entity
        trovato=.false.
        do j=1,nentity
			if(entita(j)%indx.eq.ele_all(i)%entity)then
				trovato=.true.
				exit
			end if
        end do
        
        if(trovato.eqv..false.)then
			nentity=nentity+1
			entita(nentity)%indx=ele_all(i)%entity
		end if
        
        
        allocate(ele_bordi(nele_bordi)%nodi(ele_bordi(nele_bordi)%nnodi))
        ele_bordi(nele_bordi)%nodi=ele_all(i)%nodi
        ele_bordi(nele_bordi)%tipo_ele=ele_all(i)%tipo_ele
    elseif(ele_all(i)%tipo_ele.eq.2.or.ele_all(i)%tipo_ele.eq.3) then
        nele_interni=nele_interni+1
        ele(nele_interni)%nnodi=ele_all(i)%nnodi
        allocate(ele(nele_interni)%nodi(ele(nele_interni)%nnodi))
        ele(nele_interni)%nodi=ele_all(i)%nodi
        ele(nele_interni)%tipo_ele=ele_all(i)%tipo_ele
        if(ele_all(i)%tipo_ele.eq.2)then
                                ele(nele_interni)%nnghbrs=3
        elseif(ele_all(i)%tipo_ele.eq.3)then
                                ele(nele_interni)%nnghbrs=4
        end if
    end if

end do


read(1,*)car
write(*,*)'Reading elements ...           OK'

   




   
end subroutine


