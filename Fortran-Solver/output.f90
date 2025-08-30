subroutine write_tecplot_file(k_file,traslay)
 use variabili
 implicit none
 character(80):: datafile_tec,kcar
 integer :: i,os,j,k_file,iposj2,ndofsele,iele,jj
 real(4)::rho_node,rhou_node,rhov_node,e_node,u_node,v_node,p_node,M_node,S_node,muSA_node,wall_dist_node,phii,&
 phiorton,mua_node,u6_node,u7_node,u8_node,somma,sensore_residui,vort_thick_node,traslay
 real(4),dimension(2)::csi,xy,xyloc
logical::trovato




 	if(k_file.lt.10) then
	write(kcar,'(i1)') k_file
	else if(k_file.lt.100) then
	write(kcar,'(i2)') k_file
	else if(k_file.lt.1000) then
	write(kcar,'(i3)') k_file
	else if(k_file.lt.10000) then
	write(kcar,'(i4)') k_file
	else if(k_file.lt.100000) then
	write(kcar,'(i5)') k_file
	else if(k_file.lt.1000000) then
	write(kcar,'(i6)') k_file
	else if(k_file.lt.10000000) then
	write(kcar,'(i7)') k_file
	else if(k_file.lt.100000000) then
	write(kcar,'(i8)') k_file
	else if(k_file.lt.1000000000) then
	write(kcar,'(i9)') k_file
	end if


      iposj2=index(kcar,' ')

	datafile_tec=kcar(1:iposj2-1)//'.plt'




!--------------------------------------------------------------------------------
 open(unit=1, file=datafile_tec)

 write(1,*) 'title = "grid"'


 write(1,'(a80)') 'variables = "x","y","rho","P","u","v","M","S"'

 write(1,*) 'zone n=',nnodi,' e =', nele_interni,' et=quadrilateral, f=fepoint'

!--------------------------------------------------------------------------------
! Nodal quantities: x, y, rho, u, v, p, Mach number



   do i = 1, nnodi

    rho_node=0.
    rhou_node=0.
    rhov_node=0.
    e_node=0.
    muSA_node=0.
    wall_dist_node=0.
    mua_node=0.
    vort_thick_node=0.
    u6_node=0.
    u7_node=0.
    u8_node=0.
    somma=0.



    trovato=.FALSE.
    do j=1,nnodi_vis
    if(i.eq.nodi_vis(j)) then
        trovato=.TRUE.
        exit
    end if
    end do

    if(trovato)then
        !do j=1,1
        do j=1,nodo(i)%neles
            if((nodo(i)%ele(j).ne.0))then
            
  

                
                
                rho_node=rho_node+ele(nodo(i)%ele(j))%ucons(2)
                rhou_node=rhou_node+ele(nodo(i)%ele(j))%ucons(3)
                rhov_node=rhov_node+ele(nodo(i)%ele(j))%ucons(4)
                e_node=e_node+ele(nodo(i)%ele(j))%ucons(1)
                
                
                somma=somma+1.d0
            

            end if
        end do

            rho_node=rho_node/somma
            rhou_node=rhou_node/somma
            rhov_node=rhov_node/somma
            e_node=e_node/somma


            u_node=rhou_node/rho_node
            v_node=rhov_node/rho_node
            P_node=(gam-1.)*(e_node-0.5*rho_node*((rhou_node/rho_node)**2+(rhov_node/rho_node)**2))
            M_node=sqrt(u_node**2+v_node**2)/sqrt(gam*P_node/rho_node)
            S_node=gam*log(P_node/rho_node)-(gam-1.)*log(p_node)

    else

            rho_node=rho_node+ele(iele)%ucons(2)
            rhou_node=rhou_node+ele(iele)%ucons(3)
            rhov_node=rhov_node+ele(iele)%ucons(4)
            e_node=e_node+ele(iele)%ucons(1)
            
            u_node=rhou_node/rho_node
            v_node=rhov_node/rho_node
            P_node=(gam-1.)*(e_node-0.5*rho_node*((rhou_node/rho_node)**2+(rhov_node/rho_node)**2))
            M_node=sqrt(u_node**2+v_node**2)/sqrt(gam*P_node/rho_node)
            S_node=gam*log(P_node/rho_node)-(gam-1.)*log(p_node)


    end if



            write(1,*) nodo(i)%x(1), nodo(i)%x(2)+traslay,rho_node,P_node,u_node,v_node,M_node,S_node







   end do

!--------------------------------------------------------------------------------
! Both quad and tria elements in quad format:

 do i = 1, nele_interni
   

      if (ele(i)%nnodi.eq.3) then

       write(1,*) ele(i)%nodi(1), ele(i)%nodi(2), ele(i)%nodi(3), ele(i)%nodi(3)

      elseif (ele(i)%nnodi.eq.4) then

       write(1,*) ele(i)%nodi(1), ele(i)%nodi(2), ele(i)%nodi(3), ele(i)%nodi(4)

      else

       !Impossible
       write(*,*) " Error in ele%vtx data... Stop..: i,nele_interni,ele(i)%nvtx=",i,nele_interni,ele(i)%nnodi
       stop

      endif
  
 end do

!--------------------------------------------------------------------------------
 close(1)
 end subroutine write_tecplot_file
 
 
 
 
 
subroutine prepar_tecplot
use Variabili
implicit none
integer::i,j
logical::trovato
integer,dimension(:),allocatable::temp

nnodi_vis=0

if(allocated(nodi_vis)) then
    deallocate(nodi_vis)
    deallocate(nodi_agg_vis)
end if

allocate(nodi_vis(nnodi))
allocate(nodi_agg_vis(5*nele_interni))

do i=1,ninterf

  

        if(nnodi_vis.gt.0)then

            TROVATO=.FALSE.
            do j=1,nnodi_vis
                if(interf(i)%nodo1.eq.nodi_vis(j)) TROVATO=.TRUE.
            end do

            if(TROVATO.eqv..FALSE.) then
                nnodi_vis=nnodi_vis+1
                nodi_vis(nnodi_vis)=interf(i)%nodo1
            end if

            TROVATO=.FALSE.
            do j=1,nnodi_vis
                if(interf(i)%nodo2.eq.nodi_vis(j)) TROVATO=.TRUE.
            end do

            if(TROVATO.eqv..FALSE.) then
                nnodi_vis=nnodi_vis+1
                nodi_vis(nnodi_vis)=interf(i)%nodo2
            end if

        else
            nnodi_vis=nnodi_vis+1
            nodi_vis(nnodi_vis)=interf(i)%nodo1
            nnodi_vis=nnodi_vis+1
            nodi_vis(nnodi_vis)=interf(i)%nodo2
        end if


    

end do

call bubble(nodi_vis,nnodi_vis)

allocate(temp(nnodi_vis))

temp=nodi_vis(1:nnodi_vis)

do i=1,nnodi_vis
    nodi_vis(i)=temp(nnodi_vis-i+1)
end do

deallocate(temp)

write(*,*)'nnodi_vis = ',nnodi_vis

end subroutine











 

















SUBROUTINE Bubble (X,N)
! X e' l'array da ordinare in verso decrescente, di dimensione N
IMPLICIT NONE
INTEGER :: N, J, I, JMAX, TEMP
INTEGER :: X(N)
JMAX=N-1
spike1: DO I=1,N-1
spike2: DO J=1,JMAX
IF(X(J).GT.X(J+1)) GO TO 100
TEMP=X(J)
X(J)=X(J+1)
X(J+1)=TEMP
100 END DO spike2
JMAX=JMAX-1
END DO spike1
RETURN
END subroutine


