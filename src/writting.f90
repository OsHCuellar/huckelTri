module writting 

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         January, 2022          +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Module that contains subroutines that writes outputs
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    use types
    implicit none

    contains

    subroutine wrtHamMat(nAtom,mat)                                           
                                                                     
      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         January, 2022          +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Writes the Hamiltonian matrix 
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        integer(ik), intent(in)                             :: nAtom
        real(rk), allocatable, dimension(:,:), intent(in)   :: mat

        integer(ik)                             :: i, j                     
       
        !File where is going to be written the Hamiltonian matrix  
        open(unit=200, file='output/HamMat.dat', action='write')

        !write the Hamiltonian matrix
        write(200,*) '#********** H A M I L T O N I A N    M A T R I X **********'          
        do i = 1, nAtom                                                                                          
            write(200,*) (mat(i,j), j = 1, nAtom)                                                    
        enddo                                                                                                              

        close(200)
    
    end subroutine wrtHamMat                      

    subroutine enLvl(ener, nAtom)

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         January, 2022          +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Writes the sorted energy levels 
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        integer(ik), intent(in)                               :: nAtom
        real(rk), allocatable, dimension(:), intent(in)       :: ener

        integer(ik)                               :: i

        !File where are going to be written the energy levels
        open(unit=200, file='output/EnergyLvl.dat', action='write')

        !write the energy levels
        write(200,*) '#  Level         Energy '
        do i = 1, nAtom
            write(200,*)i, ener(i)
        enddo                            

        close(200)

    endsubroutine enLvl

    subroutine eigenvectors(nAtom,mat,evec)

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         January, 2022          +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Writes the Coefficients for each energy in different files 
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        integer(ik), intent(in)                               :: nAtom
        real(rk), allocatable, dimension(:,:), intent(in)     ::  mat, evec

        integer(ik)                               :: i,j,k
        character(len=50)                         :: fileName

        !For each energy energy level (as much of atoms because each atoms contributes with one orbital)

        !To write the name of the output file in the valid format (becasue of the number of characters)
        if (nAtom .lt. 10) then 

            !___________________________________________________________________________________________
            do i = 1, nAtom

                !Output file name
                write(fileName,'(A32,I1,A4)')'output/Eigenvectors/Eigenvectors',i,'.dat'

                open(unit=200, file=fileName, action='write')

                    !write eigenvector (coefficients) for the corresponding energy level
                    write(200,*)'#E = ', mat(i,i)
                    do k = 1, nAtom
                        write(200,*) evec(k,i)
                    enddo
                    write(200,*)''

                close(200)

            enddo
            !___________________________________________________________________________________________

        else if (nAtom .lt. 100) then

            !___________________________________________________________________________________________
            do i = 1, 9

                !Output file name
                write(fileName,'(A32,I1,A4)')'output/Eigenvectors/Eigenvectors',i,'.dat'

                open(unit=200, file=fileName, action='write')

                    !write eigenvector (coefficients) for the corresponding energy level
                    write(200,*)'#E = ', mat(i,i)
                    do k = 1, nAtom
                        write(200,*) evec(k,i)
                    enddo
                    write(200,*)''

                close(200)

            enddo

            !-------------------------------------------------------------------------------------------

            do i = 10, nAtom

                !Output file name
                write(fileName,'(A32,I2,A4)')'output/Eigenvectors/Eigenvectors',i,'.dat'

                open(unit=200, file=fileName, action='write')

                    !write eigenvector (coefficients) for the corresponding energy level
                    write(200,*)'#E = ', mat(i,i)
                    do k = 1, nAtom
                        write(200,*) evec(k,i)
                    enddo
                    write(200,*)''

                close(200)

            enddo
            !___________________________________________________________________________________________

        else if (nAtom .lt. 1000) then

            !___________________________________________________________________________________________
            do i = 1, 9

                !Output file name
                write(fileName,'(A32,I1,A4)')'output/Eigenvectors/Eigenvectors',i,'.dat'

                open(unit=200, file=fileName, action='write')

                    !write eigenvector (coefficients) for the corresponding energy level
                    write(200,*)'#E = ', mat(i,i)
                    do k = 1, nAtom
                        write(200,*) evec(k,i)
                    enddo
                    write(200,*)''

                close(200)

            enddo

            !-------------------------------------------------------------------------------------------

            do i = 10, 99

                !Output file name
                write(fileName,'(A32,I2,A4)')'output/Eigenvectors/Eigenvectors',i,'.dat'

                open(unit=200, file=fileName, action='write')

                    !write eigenvector (coefficients) for the corresponding energy level
                    write(200,*)'#E = ', mat(i,i)
                    do k = 1, nAtom
                        write(200,*) evec(k,i)
                    enddo
                    write(200,*)''

                close(200)

            enddo

            !-------------------------------------------------------------------------------------------
            do i = 100, nAtom

                !Output file name
                write(fileName,'(A32,I3,A4)')'output/Eigenvectors/Eigenvectors',i,'.dat'

                open(unit=200, file=fileName, action='write')

                    !write eigenvector (coefficients) for the corresponding energy level
                    write(200,*)'#E = ', mat(i,i)
                    do k = 1, nAtom
                        write(200,*) evec(k,i)
                    enddo
                    write(200,*)''

                close(200)

            enddo

            !___________________________________________________________________________________________

        else

            write(*,*) 'W A R N I N G !!!'
            write(*,*) 'The eigenvectors were not written'
            write(*,*) 'If there are more than 999 atoms and want its eigenvalues, you should modify the writting module'  

        endif

    endsubroutine eigenvectors

end module writting
