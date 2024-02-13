program huckel

       ! ++++++++++++++++++++++++++++++++++                                         
       ! +    Osvaldo Hernandez Cuellar   +                                         
       ! +         February, 2022         +                                         
       ! ++++++++++++++++++++++++++++++++++                                         
                                                                                    
       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<                                  
       ! H U C K E L     P R O G R A M
       !
       !Using the Huckel method, creates the Hamiltonian matrix for the structure defined in the input files,
       !diagonalize the matrix, and return the eigenvalues and eigenvectors (energies and coefficients)
       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                                  

    use types
    use start
    use math
    use jacobiDiag
    use dsyevDiag
    use writting

    integer(ik)                               :: i
    integer(ik)                               :: nAtom, nBond, nBeta, nAlpha
    integer(ik)                               :: diagSub
    integer(ik), allocatable, dimension(:,:)  :: conBond
    integer(ik), allocatable, dimension(:)    :: posPrevNew
    real(rk), allocatable, dimension(:,:)     :: mat, evec,u
    real(rk), allocatable, dimension(:)       :: conBeta, conAlpha
    real(rk), allocatable, dimension(:)       :: alphas, betas
    real(rk), allocatable, dimension(:)       :: ener
    logical                                   :: debug


    !READ THE INPUT
    call input(nAtom,nBond,nAlpha, nBeta,conBond,conAlpha,conBeta, debug, diagSub)

    !MAKE THE HAMILTONIAN MATRIX
    call matrixMake(mat,nAtom, conAlpha, nBond, conBond, conBeta)

    !WRITE THE HAMILTONIAN MATRIX IF IT IS DEBUGGING
    if(debug) call wrtHamMat(nAtom,mat)

    !ALLOCATE THE MATRIX THATH WILL CONTAIN THE EIGENVECTORS
    allocate(evec(nAtom,nAtom))

    !CHOOSE THE METHOD TO DIAGONALIZE THE HAMILTONIAN MATRIX 
    !(Use subdsyevv preferably. jacobi works, but subdsyevv is way more eficient.)
    if (diagSub .eq. 0) then
        call subdsyevv(mat,nAtom,evec)
    else if (diagSub .eq. 1) then
        call jacobi(mat, nAtom, evec) 
    else
        write(*,*) 'W A R N I N G !!!'
        write(*,*) 'Diagonalization method choosed not defined'
        stop
    endif

    !NORMALIZE THE EIGENVETORS (coefficients)
    call normvec(evec, nAtom)

    allocate(ener(nAtom)) !vector to store the energies

    !STORE THE EIGENVALUES (energies) IN ONE VECTOR
    do i = 1, nAtom
        ener(i) = mat(i,i)
    enddo

    !SORT THE ENERGIES FROM LOWER TO HIGHER
    call lowToHigh(ener, nAtom)

    !WRITE THE ENERGIES (output)
    call enLvl(ener, nAtom)

    deallocate(ener)

    !WRITE THE EIGENVECTORS (output)
    call eigenvectors(nAtom,mat,evec)

end program huckel
