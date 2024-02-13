module jacobiDiag 

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         December, 2020         +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Module of the Jacobi method
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    use types
    implicit none

    contains 

    subroutine jacobi (d, nAtom, vec)

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         December, 2020         +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Jacobi method to diaginalize a symmetric square matrix
      !(Method not optimized for big matrixes)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                                                                                                                                  
        integer(ik), intent(in)                               :: nAtom
        real(rk), allocatable, dimension(:,:), intent(inout)  :: d, vec

        integer(ik)                                           :: i, j, iter, converg
        integer(ik)                                           :: p, q
        real(rk)                                              :: cont, threshold 
        real(rk)                                              :: bigger, rot, pi
        real(rk), allocatable, dimension(:,:)                 :: s

        threshold = 1E-10

        pi = 2.0_rk * asin(1.0_rk) !define pi

        allocate (s(nAtom,nAtom))

        !Initialize matrix that will store the eigenvectors
        vec = 0.0_rk
        do i = 1, nAtom
            vec(i,i) = 1.0_rk
        enddo

        !Loop to reach convergence
        do iter = 1, 9999999 !while loop not used to avoid the possible excesive cyclation

            !Initialize S
            s = 0.0_rk
            do i = 1, nAtom
                s(i,i) = 1.0_rk
            enddo

            !Largest off-diagonal element (magnitude) of d
            bigger = 0.0_rk
            do i = 1, nAtom-1
                do j = i+1, nAtom
                    if ( abs(d(i,j)) > abs(bigger) ) then
                        p = i
                        q = j
                        bigger = d(i,j)
                    endif
                enddo
            enddo

            !Angle of rotation
            if (d(p,p) == d(q,q)) then
                rot = pi/4.0_rk
            else
                rot = 2.0_rk * d(p,q)
                rot = rot / (d(q,q)-d(p,p))
                rot = atan(rot)
                rot = rot / 2.0_rk
            endif

            !Current rotation matrix
            s(p,p) = cos(rot)
            s(q,q) = cos(rot)
            s(p,q) = sin(rot)
            s(q,p) = -sin(rot)

            !(S^(-1))DS 
            d = matmul(transpose(s),matmul(d,s))

            !Multiplication of all S matrixes
            vec = matmul(vec,s)

            !Check convergence 
            !|d_ii| >> sum_(i /= j) |d_ij| 
            converg = 0_ik !counter for convenrgence
            do i=1, nAtom
                cont=0.0_rk
                do j = 1, nAtom
                    if (i .ne. j) then
                        cont = cont + abs(d(i,j))  
                    endif
                enddo
                if ( cont/(d(i,i)) .lt. threshold ) then
                    converg = converg + 1_ik
                endif
            enddo

            if (converg .ge. nAtom) exit

        enddo

        deallocate(s)

    end subroutine jacobi                    

end module jacobiDiag
