subroutine compute_fluxes
use variabili
implicit none
integer::i,e1,e2,j

do i=1,nele_interni
	ele(i)%d_dt(:)=0.
end do

do i=1,ninterf

	e1=interf(i)%e1
	e2=abs(interf(i)%e2)

	if(e1*e2.ne.0) then

!       call compute_internal_flux(i)
        call compute_internal_flux_Roe(i)
	else

		call compute_boundary_flux(i)

	end if

	do j=1,4
		ele(e1)%d_dt(j)=ele(e1)%d_dt(j)-interf(i)%F(j)*interf(i)%length/ele(e1)%area
	end do


	if(e2.gt.0)then
		do j=1,4
			ele(e2)%d_dt(j)=ele(e2)%d_dt(j)+interf(i)%F(j)*interf(i)%length/ele(e2)%area
		end do
	end if

end do



end subroutine

subroutine compute_internal_flux(i)
use variabili
implicit none
integer::i
real(4),dimension(4)::uconsl,uconsr,fluxes
real(4)::utildeA,utildeB,vtildeA,vtildeB,pa,pb,Ta,Tb,ua,ub,va,vb,lambda_max



! In questa subroutine si calcolano i flussi all'interfaccia usando il solutore di Lax-Friedrichs locale o di Roe

uconsl=ele(interf(i)%e1)%ucons(:)
uconsr=ele(abs(interf(i)%e2))%ucons(:)

!calcolo grandezze primitive

ua= uconsl(3)/uconsl(2)
ub= uconsr(3)/uconsr(2)

va= uconsl(4)/uconsl(2)
vb= uconsr(4)/uconsr(2)

pa= (gam-1.)*(uconsl(1)-0.5*uconsl(2)*(ua**2+va**2))
pb= (gam-1.)*(uconsr(1)-0.5*uconsr(2)*(ub**2+vb**2))

Ta= pa/uconsl(2)
Tb= pb/uconsr(2)

! proietto velocità in direzione normale e tangenziale
utildeA= interf(i)%normal(1)*ua + interf(i)%normal(2)*va
utildeB= interf(i)%normal(1)*ub + interf(i)%normal(2)*vb

vtildeA= -interf(i)%normal(2)*ua + interf(i)%normal(1)*va
vtildeB= -interf(i)%normal(2)*ub + interf(i)%normal(1)*vb


! calcolo massima velocità di propagazione
lambda_max= max(abs(utildeA)+sqrt(gam*Ta),abs(utildeB)+sqrt(gam*Tb))


! Flussi nel sistema ruotato con metodo di Lax-Friedrichs locale
fluxes(1)= 0.5 * (utildeA*(pa+uconsl(1)) + utildeB*(pb+uconsr(1))) - 0.5 * lambda_max * (uconsr(1)-uconsl(1))
fluxes(2)= 0.5 * (utildeA*uconsl(2)+utildeB*uconsr(2)) - 0.5 * lambda_max * (uconsr(2)-uconsl(2))
fluxes(3)= 0.5 * (pa+uconsl(2)*utildeA**2 + pb+uconsr(2)*utildeB**2) - 0.5 * lambda_max * (uconsr(2)*utildeB-uconsl(2)*utildeA)
fluxes(4)= 0.5 * (uconsl(2)*utildeA*vtildeA + uconsr(2)*utildeB*vtildeB) - 0.5 * lambda_max * (uconsr(2)*vtildeB-uconsl(2)*vtildeA)



!trasformazione dei flussi dal sistema locale a quello globale
interf(i)%F(1)=fluxes(1)
interf(i)%F(2)=fluxes(2)
interf(i)%F(3)= interf(i)%normal(1)*fluxes(3) - interf(i)%normal(2)*fluxes(4)
interf(i)%F(4)= interf(i)%normal(2)*fluxes(3) + interf(i)%normal(1)*fluxes(4)





end subroutine


subroutine compute_internal_flux_Roe(i)
use variabili
implicit none
integer::i
real,dimension(4)::ul,ur,fluxes



ul(1)=ele(interf(i)%e1)%ucons(2)
ul(2)=ele(interf(i)%e1)%ucons(3)
ul(3)=ele(interf(i)%e1)%ucons(4)
ul(4)=ele(interf(i)%e1)%ucons(1)

ur(1)=ele(abs(interf(i)%e2))%ucons(2)
ur(2)=ele(abs(interf(i)%e2))%ucons(3)
ur(3)=ele(abs(interf(i)%e2))%ucons(4)
ur(4)=ele(abs(interf(i)%e2))%ucons(1)


call Roe(uL, uR, interf(i)%normal(1),interf(i)%normal(2),fluxes)

interf(i)%F(1)=fluxes(4)
interf(i)%F(2)=fluxes(1)
interf(i)%F(3)=fluxes(2)
interf(i)%F(4)=fluxes(3)


end subroutine



subroutine Roe(uL, uR, nx, ny,fluxes)
 real :: uL(4), uR(4) !  Input: conservative variables rho*[1, u, v, E]
 real :: nx, ny       !  Input: face normal vector, [nx, ny] (Left-to-Right)
 real :: fluxes(4)       ! Output: Roe flux function (upwind)
!Local constants
 real :: gamma                          ! Ratio of specific heat.
 real :: zero, fifth, half, one, two    ! Numbers
!Local variables
 real :: tx, ty       ! Tangent vector (perpendicular to the face normal)
 real :: vxL, vxR, vyL, vyR             ! Velocity components.
 real :: rhoL, rhoR, pL, pR             ! Primitive variables.
 real :: vnL, vnR, vtL, vtR             ! Normal and tangent velocities
 real :: aL, aR, HL, HR                 ! Speeds of sound.
 real :: RT,rho,vx,vy,H,a,vn, vt        ! Roe-averages
 real :: drho,dvx,dvy,dvn,dvt,dp,dV(4)  ! Wave strenghs
 real :: ws(4),dws(4), Rv(4,4)          ! Wave speeds and right-eigevectors
 real :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
 integer :: i, j

!Constants.
     gamma = 1.4
      zero = 0.0
     fifth = 0.2
      half = 0.5
       one = 1.0
       two = 2.0

!Tangent vector (Do you like it? Actually, Roe flux can be implemented
! without any tangent vector. See "I do like CFD, VOL.1" for details.)
  tx = -ny
  ty = nx

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
     vxL = uL(2)/uL(1)
     vyL = uL(3)/uL(1)
     vnL = vxL*nx+vyL*ny
     vtL = vxL*tx+vyL*ty
      pL = (gamma-one)*( uL(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(4) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
     vxR = uR(2)/uR(1)
     vyR = uR(3)/uR(1)
     vnR = vxR*nx+vyR*ny
     vtR = vxR*tx+vyR*ty
      pR = (gamma-one)*( uR(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(4) + pR ) / rhoR

!First compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL
    vx = (vxL+RT*vxR)/(one+RT)
    vy = (vyL+RT*vyR)/(one+RT)
     H = ( HL+RT* HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*(vx*vx+vy*vy)) )
    vn = vx*nx+vy*ny
    vt = vx*tx+vy*ty

!Wave Strengths
   drho = rhoR - rhoL
     dp =   pR - pL
    dvn =  vnR - vnL
    dvt =  vtR - vtL

  dV(1) = (dp - rho*a*dvn )/(two*a*a)
  dV(2) = rho*dvt/a
  dV(3) =  drho - dp/(a*a)
  dV(4) = (dp + rho*a*dvn )/(two*a*a)

!Wave Speed
  ws(1) = abs(vn-a)
  ws(2) = abs(vn)
  ws(3) = abs(vn)
  ws(4) = abs(vn+a)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
! only for the nonlinear fields.
  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(4) = fifth
   if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )

!Right Eigenvectors
  Rv(1,1) = one
  Rv(2,1) = vx - a*nx
  Rv(3,1) = vy - a*ny
  Rv(4,1) =  H - vn*a

  Rv(1,2) = zero
  Rv(2,2) = a*tx
  Rv(3,2) = a*ty
  Rv(4,2) = vt*a

  Rv(1,3) = one
  Rv(2,3) = vx
  Rv(3,3) = vy
  Rv(4,3) = half*(vx*vx+vy*vy)

  Rv(1,4) = one
  Rv(2,4) = vx + a*nx
  Rv(3,4) = vy + a*ny
  Rv(4,4) =  H + vn*a

!Dissipation Term
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
   end do
  end do

!Compute the flux.
  fL(1) = rhoL*vnL
  fL(2) = rhoL*vnL * vxL + pL*nx
  fL(3) = rhoL*vnL * vyL + pL*ny
  fL(4) = rhoL*vnL *  HL

  fR(1) = rhoR*vnR
  fR(2) = rhoR*vnR * vxR + pR*nx
  fR(3) = rhoR*vnR * vyR + pR*ny
  fR(4) = rhoR*vnR *  HR

  fluxes = half * (fL + fR - diss)

 end subroutine
