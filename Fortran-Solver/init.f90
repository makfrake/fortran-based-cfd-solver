subroutine init
use variabili
implicit none
integer::i
real(4)::rho,u,v,p,t,alphaloc

	gam=1.4 !rapporto dei calori specifici

    ! variabili derivate dal rapporto dei calori specifici
    GA=gam/(gam-1.)
    GB=1./(gam-1.)
    GC=(gam+1.)/(gam-1.)
    GD=(gam-1.)/2.
    GE=(gam+1.)/2.
    GF=sqrt(gam)
    GG=2./(gam-1.)
    GH=(gam+1.)/(2.*gam)
    GI=(gam-1.)/(2.*gam)
    GJ=(gam-1.)/gam

    do i = 1,nele_interni

        ele(i)%P = ptotin / (1+0.5*(gam-1.)*machin**2)**(gam/(gam-1.))
        ele(i)%T = ttotin / (1+0.5*(gam-1.)*machin**2)
        ele(i)%a = sqrt(gam*ele(i)%T)

        if (ele(i)%x0(1).lt.0) then
            alphaloc = alpha
        elseif (ele(i)%x0(1).lt.1) then
            alphaloc = alpha + (-60.*((4.*atan(1.0))/(180.))-alpha)/1.0*(ele(i)%x0(1))
        else
                alphaloc = -60.*((4.*atan(1.0))/(180.))
        end if

        ele(i)%u = machin*ele(i)%a*cos(alphaloc)
        ele(i)%v = machin*ele(i)%a*sin(alphaloc)
        ele(i)%s = gam*log(ele(i)%T)-(gam-1.)*log(ele(i)%P)

        ! Vettore grandezze conservative
        ele(i)%ucons(2) = ele(i)%P/ele(i)%T
        ele(i)%ucons(1) = ele(i)%P/(gam-1.) + 0.5*ele(i)%ucons(2) * (ele(i)%u**2 + ele(i)%v**2)
        ele(i)%ucons(3) = ele(i)%ucons(2)*ele(i)%u
        ele(i)%ucons(4) = ele(i)%ucons(2)*ele(i)%v

    end do

end subroutine
