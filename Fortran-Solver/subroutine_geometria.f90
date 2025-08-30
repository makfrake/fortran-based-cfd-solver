subroutine compute_center_ele(i)
use variabili
implicit none
integer::i,j

 ! in questa subroutine si calcola il baricentro dell'elemento facendo la media aritmetica delle coordinate dei noti dell'ele

    ele(i)%x0=0

    do j=1,ele(i)%nnodi
        ele(i)%x0 = ele(i)%x0 + nodo(ele(i)%nodi(j))%x
    end do

    ele(i)%x0 = ele(i)%x0 / ele(i)%nnodi

end subroutine



subroutine compute_area_ele(i)
use variabili
implicit none
integer::i,j
real::a,b,c,d,f,At1,At2,p1,p2
real,dimension(3) :: x,y

!Calcolare l'area dell'elemento:
! * se è un triangolo con lati di lunghezza a,b,c e semiperimetro p1=0.5*(a+b+c) usare la formula di Erone At=sqrt(p1*(p1-a)*(p1-b)*(p1-c))
! * se è un quadrilatero calcolare la sua area come somma dei due triangoli che lo compongono, usando la formula di Erone su ciascuno

if (ele(i)%nnodi.eq.3) then
    a = sqrt( (nodo(ele(i)%nodi(2))%x(1)-nodo(ele(i)%nodi(1))%x(1))**2 + (nodo(ele(i)%nodi(2))%x(2)-nodo(ele(i)%nodi(1))%x(2))**2 )
    b = sqrt( (nodo(ele(i)%nodi(3))%x(1)-nodo(ele(i)%nodi(2))%x(1))**2 + (nodo(ele(i)%nodi(3))%x(2)-nodo(ele(i)%nodi(2))%x(2))**2 )
    c = sqrt( (nodo(ele(i)%nodi(3))%x(1)-nodo(ele(i)%nodi(1))%x(1))**2 + (nodo(ele(i)%nodi(3))%x(2)-nodo(ele(i)%nodi(1))%x(2))**2 )

    p1 = (a+b+c) * 0.5

    At1=sqrt(p1*(p1-a)*(p1-b)*(p1-c))

    ele(i)%area = At1

end if

if (ele(i)%nnodi.eq.4) then
    a = sqrt( (nodo(ele(i)%nodi(2))%x(1)-nodo(ele(i)%nodi(1))%x(1))**2 + (nodo(ele(i)%nodi(2))%x(2)-nodo(ele(i)%nodi(1))%x(2))**2 )
    b = sqrt( (nodo(ele(i)%nodi(3))%x(1)-nodo(ele(i)%nodi(2))%x(1))**2 + (nodo(ele(i)%nodi(3))%x(2)-nodo(ele(i)%nodi(2))%x(2))**2 )
    c = sqrt( (nodo(ele(i)%nodi(3))%x(1)-nodo(ele(i)%nodi(4))%x(1))**2 + (nodo(ele(i)%nodi(3))%x(2)-nodo(ele(i)%nodi(4))%x(2))**2 )
    d = sqrt( (nodo(ele(i)%nodi(4))%x(1)-nodo(ele(i)%nodi(1))%x(1))**2 + (nodo(ele(i)%nodi(4))%x(2)-nodo(ele(i)%nodi(1))%x(2))**2 )
    f = sqrt( (nodo(ele(i)%nodi(3))%x(1)-nodo(ele(i)%nodi(1))%x(1))**2 + (nodo(ele(i)%nodi(3))%x(2)-nodo(ele(i)%nodi(1))%x(2))**2 )


    p1 = (a+b+f) * 0.5
    p2 = (c+d+f) * 0.5

    At1=sqrt(p1*(p1-a)*(p1-b)*(p1-f))
    At2=sqrt(p2*(p2-c)*(p2-d)*(p2-f))

    ele(i)%area = At1 + At2

end if




!    do j=1,ele(i)%nnodi
!        x(j) = nodo(ele(i)%nodi(j))%x(1)
!       y(j) = nodo(ele(i)%nodi(j))%x(2)
!    end do
!
!    if(ele(i)%nnodi.le.3) then
!            a = sqrt((x(1)-x(2))**2+(y(1)-y(2))**2)
!           b = sqrt((x(2)-x(3))**2+(y(2)-y(3))**2)
!            c = sqrt((x(1)-x(3))**2+(y(1)-y(3))**2)
!            p1 =  (a+b+c) * 0.5
!            p2 =  0
!            At1 =  sqrt(p1*(p1-a)*(p1-b)*(p1-c))
!            ele(i)%area = At1
!    else
!            a = sqrt((x(1)-x(2))**2+(y(1)-y(2))**2)
!            b = sqrt((x(2)-x(3))**2+(y(2)-y(3))**2)
!            c = sqrt((x(4)-x(3))**2+(y(4)-y(3))**2)
!            d = sqrt((x(4)-x(1))**2+(y(4)-y(1))**2)
!            f = sqrt((x(1)-x(3))**2+(y(1)-y(3))**2)
!            p1 =  (a+b+f) * 0.5
!            p2 =  (c+d+f) * 0.5
!            At1 =  sqrt(p1*(p-a)*(p1-b)*(p1-f))
!            At2 = sqrt(p2*(p2-c)*(p2-d)*(p2-f))
!            AA = At1 + At2
!            ele(i)%area = AA
!    end if


end subroutine

subroutine compute_length_inte(i)
use variabili
implicit none
integer::i

! Calcolare la lunghezzas dell'interfaccia usando Pitagora sulle coordinate dei nodi 1 e 2 dell'interfaccia

interf(i)%length = sqrt((nodo(interf(i)%nodo2)%x(1)-nodo(interf(i)%nodo1)%x(1))**2 + (nodo(interf(i)%nodo2)%x(2)-nodo(interf(i)%nodo1)%x(2))**2)


end subroutine


subroutine compute_normal_inte(i)
use variabili
implicit none
integer::i
real::tx,ty

!in questa subroutine si calcola la normale all'interfaccia i. calcolare il versore tangente (tx,ty) e poi calcolare la normale ad esso perpendicolare nx=-ty; ny=tx

tx = (nodo(interf(i)%nodo2)%x(1)-nodo(interf(i)%nodo1)%x(1))
ty = (nodo(interf(i)%nodo2)%x(2)-nodo(interf(i)%nodo1)%x(2))

interf(i)%normal(1) = - ty/sqrt(tx**2+ty**2)
interf(i)%normal(2) = tx/sqrt(tx**2+ty**2)



end subroutine

