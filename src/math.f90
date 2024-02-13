module math

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         December, 2020         +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Contains some mathematical subroutines
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                                                                                                                                  
    use types
    implicit none

    contains

    subroutine lowToHigh (vec1, n)

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         December, 2022         +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Sort the values of the vector vec1 from lower to higher 
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                                                                                                                                  
        integer(ik), intent(in)                               :: n
        real(rk), allocatable, dimension(:), intent(inout)    :: vec1

        integer(ik)                                           :: i,j, cont
        real(rk)                                              :: bigger, vmax, vmin
        real(rk), allocatable, dimension(:)                   :: vec2


        allocate(vec2(n))

        !Find mavimum value of vec1 to latter add it to the lowes 
        vmax = maxval(vec1(:))

        do i = 1, n

            !Find lower value in vec1
            j=minloc(vec1(:),n)

            !Copy lower value of vec1 into vec2
            vec2(i) = vec1(j)

            !sum value to avoid counting it again
            vec1(j) = vmax + vmax

        enddo

        !Copy the sorted vector into the first one to return it
        vec1= vec2

        !Clean memory
        deallocate(vec2)

    end subroutine lowToHigh                    

    subroutine normVec(vec, n)

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         December, 2020         +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Normalizes the columns of the matrix vec of size (n,n), treating each of them as a vector
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        real(rk), allocatable, dimension(:,:), intent(inout)        :: vec
        integer(ik), intent(in)                                     :: n

        integer(ik)                                                 :: i, j, k
        real(rk)                                                    :: integ

        do k = 1, n

            !Sum over the coefficients squared
            integ = dot_product(vec(:,k),vec(:,k))

            if (integ /= 1.0_rk) then
                !divide each coefficient by the normalization constant
                vec(:,k) = vec(:,k)/sqrt(integ)
            end if

        enddo

    end subroutine 

end module math
