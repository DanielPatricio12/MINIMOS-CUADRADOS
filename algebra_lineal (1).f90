module algebra_lineal
implicit none

contains

    subroutine metodo_gauss(AB, X, n)
        implicit none
        
        real(10):: AB(n,n+1), X(n), AUX(n+1), m
        integer:: n, i, k, pos(1)
           
        do k=1, n-1
            do i=k+1, n
                do while (AB(k, k)==0) 
                    pos= maxloc(abs(AB(k+1:n,k)))
                if (maxval(AB(:, k))==0) then
                    write(*,*) ' Sistema incompatible'
                    stop
                    else
                    AUX= AB(k, :)
                    AB(k, :) = AB(k+pos(1), :)
                    AB(k+pos(1), :)= AUX
                    end if 
                end do
                m=AB(i, k)/AB(k, k)
                AB(i, :) = AB(i, :) - m*AB(k, :)
            end do
        end do

        call sustitucion_abajo_arriba (AB, X, n)

    end subroutine

    subroutine metodo_LU(A, B, X, n)
        implicit none
        real(10):: A(n, n), B(n), X(n), U(n, n), L(n, n), P(n, n), Y(n), BA(n, n+1), AB(n, n+1), BAUX(n)
        integer:: n

        call factorizacion_LU (A, n, L, U, P)
        BAUX= matmul(P, B)
        BA(:, 1) = BAUX 
        BA(:, 2:n+1) = L
        call sustitucion_arriba_abajo(BA, Y, n)
        AB(:, 1:n)= U
        AB(1:n, n+1)= Y
        call sustitucion_abajo_arriba(AB, X, n)


    end subroutine

    subroutine factorizacion_LU (A, n, L, U, P)
        real(10):: A(n, n), U(n, n), L(n, n), P(n, n), AUX(n), n_aux, m
        integer:: n, i, k, pos(1)

        U=A
        L=0.d0
        P = 0.d0
        do i =1,n
            L(i,i) = 1.d0
            P(i,i) = 1.d0
        enddo

        do k=1, n
            do while (U(k, k)==0)
                pos= maxloc(abs(U(k+1:n,k)))
                if (maxval(A(:, k))==0) then
                    write(*,*) ' Sistema incompatible'
                    stop
                else
                    AUX= U(k, :)
                    U(k, :) = U(k+ pos(1), :)
                    U(k+ pos(1), :)= AUX

                    AUX= P(k, :)
                    P(k, :) = P(k+ pos(1), :)
                    P(k+ pos(1), :)= AUX

                    do i=1, k-1
                        n_aux= L(k, i)
                        L(k, i)= L(k+pos(1),i)
                        L(k+pos(1),i) = n_aux
                    end do
                    
                end if
            end do

            do i = k+1, n
                m = U(i,k)/U(k,k) 
                U(i,:) = U(i,:) - m*U(k,:)
                L(i,k) = m
            end do
        end do

    end subroutine

    subroutine sustitucion_abajo_arriba (AB, X, n)
        implicit none
        real(10):: AB(n, n+1), X(n), m
        integer:: k, n, i

         do k= n, 1, -1
            m = 0.d0
            do i=k+1, n
                m = m + AB(k, i)*X(i)
            end do
            X(k)=(AB(k, n+1)-m)/AB(k,k)
        end do
        
    end subroutine
    
    subroutine sustitucion_arriba_abajo (BA, X, n)
        implicit none
        real(10):: BA(n, n+1), X(n), m
        integer:: k, n, i

        do k= 1, n
            m = 0.d0
            do i=2, k
                m = m + BA(k, i)*x(i-1)
            end do
            X(k)=(BA(k, 1)-m)/BA(k,k+1)
        end do

    end subroutine


    subroutine minimos_cuadrados(a,m,b,n, grado)
        implicit none
        real(10)::a(n, 2), m(grado, grado), b(grado)
        integer::n, i, j, grado
        do i=1, grado
            do, j=1, grado
                m(i,j)=a(i,1)**(i+j-2)
                b(i)=a(i,1)**(i-1)*a(i,2)
            end do
        end do
    end subroutine
end module