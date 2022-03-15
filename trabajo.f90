program practica_cuadrados
    use algebra_lineal
    implicit none
    !Crear variables (nr= n real y ni= n entero)
    real(10), allocatable,dimension(:, :) :: A, m, mb
    real(10), allocatable,dimension(:) :: b, x
    integer::i, ni, grado
    real:: nr(1)

    !Abrimos archivo
    open(10, file='datos.dat')
    
    !Dar valores
    print*, 'Escribe el grado del polinomio'
    read*, grado
    read(10, *) nr(1)
    ni=int(nr(1))
    allocate(a(ni,2), m(grado, grado), mb(grado, grado+1), b(grado), x(grado))

    do i=1, ni
        read(10,*) A(i,:)
    end do

    call minimos_cuadrados(a,m,b,ni,grado)
    
    print*, '-----------m---------'
    do i=1, grado
        write(*, fmt='(*(e25.18, 1x))') m(i, :)
    end do

    print*, '-----------b---------'
    do i=1, grado
        write(*, fmt='(*(e25.18, 1x))') b(i)
    end do
   
    mb(:, 1:grado)=m
    mb(:, grado+1)=b
   
    print*, '-----------mb---------'
    do i=1, grado
        write(*, fmt='(*(e25.18, 1x))') mb(i, :)
    end do
    
    call metodo_gauss (mb, x, grado)

    print*, '-----------x---------'
    do i=1, grado
        write(*, fmt='(2(e25.18, 1x))') x(i)
    end do

    write(10, *)'TITLE= "TESTPLOT"'
    write(10, *)'VARIABLES="x", "y"'
    


end program practica_cuadrados