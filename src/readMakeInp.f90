module start 

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         January, 2022          +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Module that contains the subroutines needed to start the Huckel calculation
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    use types
    implicit none

    contains

    subroutine input(nAtom,nBond,nAlpha, nBeta,conBond,conAlpha,conBeta,debug, diagSub) 

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         January, 2022          +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Reads from the input.dat file and decides if reads the variables directly or updates input.dat with correct variables
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        integer(ik), intent(out)                               :: diagSub
        integer(ik), intent(out)                               :: nAtom, nBond !number of atoms, number of bonds
        integer(ik), intent(out)                               :: nAlpha, nBeta !number of diferent alphas and betas
        integer(ik), allocatable, dimension(:,:), intent(out)  :: conBond  !pair of atoms that form a bond
        real(rk), allocatable, dimension(:), intent(out)       :: conBeta, conAlpha !alphas and beta for every atom and every bond  
        logical, intent(out)                                   :: debug

        integer(ik)                                            :: inpTyp
        character                                              :: junk

        !Main input file
        open(unit=100, file='Global.inp', action='read', status='old')
        read(100,*) junk
        read(100,*) inpTyp !choose if reads the file or creates the new input file
        read(100,*) debug !do you want to debbug?
        read(100,*) diagSub !which diagonalization subrotuine 
        close(100)

        if (inpTyp .eq. 0 ) then !Explicit input

            write(*,*) 'Calculation done for a system defined explicitly'
            write(*,*) 'Results on the directory output'

            !Reads the explcit input file
            call inpread(nAtom,nBond,nAlpha, nBeta,conBond,conAlpha,conBeta)


        elseif(inpTyp .eq. 1 ) then !Nanotube input 

            write(*,*) 'Calculation done for a NANOTUBE'
            write(*,*) 'Results on the directory output'

            !Creates the explicit input file thath correspond to the nanotube described in 1nanotube.inp
            call nanoTubeInp

            !Reads the explcit input file
            call inpread(nAtom,nBond,nAlpha, nBeta,conBond,conAlpha,conBeta)

            !Copy nanotube input into the output directory
            call system('cp 1nanotube.inp output/.')

        elseif(inpTyp .eq. 2) then !Polyene input

            write(*,*) 'Calculation done for a POLYENE'
            write(*,*) 'Results on the directory output'

            !Creates the explicit input file thath correspond to the polyene described in 2polyene.inp
            call polyeneInp

            !Reads the explcit input file
            call inpread(nAtom,nBond,nAlpha, nBeta,conBond,conAlpha,conBeta)

            !Copy polyene input into the output directory
            call system('cp 2polyene.inp output/.')

        elseif(inpTyp .eq. 3) then !Triangulene input

            write(*,*) 'Calculation done for a TRIANGULENE'
            write(*,*) 'Results on the directory output'

            !Creates the explicit input file thath correspond to the triangulene described in 3triangulene.inp
            call trianguleneInp

            !Reads the explcit input file
            call inpread(nAtom,nBond,nAlpha, nBeta,conBond,conAlpha,conBeta)

            !Copy triangulene input into the output directory
            call system('cp 3triangulene.inp output/.')

        else

            write(*,*)''
            write(*,*)'W A R N I N G !!!'
            write(*,*)''
            write(*,*)'The type selected in Global.inp is not defined'
            write(*,*)''
            stop

        endif

    end subroutine input

    subroutine inpread(nAtom,nBond,nAlpha, nBeta,conBond,conAlpha,conBeta) 

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         January, 2022          +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Continue reading from the input.dat file the variables needed
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        integer(ik), intent(out)                               :: nAtom, nBond !number of atoms, number of bonds
        integer(ik), intent(out)                               :: nAlpha, nBeta !number of diferent alphas and betas
        integer(ik), allocatable, dimension(:,:), intent(out)  :: conBond  !pair of atoms that form a bond
        real(rk), allocatable, dimension(:), intent(out)       :: conBeta, conAlpha !alphas and beta for every atom and every bond  

        integer(ik)                               :: i 
        integer(ik), allocatable, dimension(:)    :: conAtm !type of alpha for each atom
        real(rk), allocatable, dimension(:)       :: alphas !all different alphas
        real(rk), allocatable, dimension(:)       :: betas  !all diferent betas
        character                                 :: junk


        !Remove previous input and output files from the output directory
        call system('rm output/*.inp')
        call system('rm output/Eigenvectors/*.dat')

        !Copy global input files int the output directory
        call system('cp Global.inp output/.')
        call system('cp 0explicit.inp output/.')

        !Open the explicit input file
        open(unit=200, file='0explicit.inp', status='old')

        !Read input data
        read(200,*) junk
        read(200,*) nAtom
        read(200,*) nAlpha
        read(200,*) nBond
        read(200,*) nBeta

        !Allocate the needed matrixes 
        allocate (alphas(nAlpha))
        allocate (betas(nBeta))   
        allocate (conBond(nBond,3)) 
        allocate (conAtm(nAtom)) 
        allocate (conBeta(nBond)) 
        allocate (conAlpha(nAtom)) 

        read(200,*) junk

        !Read alphas values
        do i=1,nAlpha
            read(200,*) alphas(i)
        enddo

        read(200,*) junk

        !Read betas values
        do i=1,nBeta
            read(200,*) betas(i)
        enddo

        read(200,*) junk

        !Read all pairs of atoms bonded
        do i = 1, nBond 

            !atom, atom, index of the beta of the bond
            read(200,*) conBond(i,1), conBond(i,2), conBond(i,3)

            !beta of the bond
            conBeta(i) = betas(conBond(i,3))
        enddo

        read(200,*) junk
        do i=1,nAtom

            !index of the alpha of the atom
            read(200,*) conAtm(i)

            !alpha of the atom
            conAlpha(i) = alphas(conAtm(i))

        enddo

        close(200)

    end subroutine inpread

    subroutine nanoTubeInp                                              
                                                                     
      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         January, 2022          +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !creates the input for a nanotube using the parameters in 1nanotube.inp
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        integer                              :: i, j                     
        integer                              :: nRingSide, nRingHeight 
        integer                              :: alpha, beta 
        integer                              :: nAtm, nBond, atmCol, atm 
        integer, allocatable, dimension(:)   :: con                      
        character                            :: junk
      

        !Open intput file for the nanotube
        open(unit=300, file='1nanotube.inp', status='old')

        read(300,*) junk
        read(300,*)nRingSide !number of rings long in the input                      
        read(300,*)nRingHeight !number of rings circular in the input                                                 
        read(300,*)alpha                                                        
        read(300,*)beta                                                         

        close(300)

        !Open the explicit file to create it
        open(unit=200, file='0explicit.inp', action='write', status='unknown')
                                                                         
        !Calculate the number of atoms                                                                  
        nAtm = nRingHeight * 4 * nRingSide                               
        
        write(200,*) 'input title'                                         
        write(200,*) ''
        write(200,*) nAtm, '#Total number of atoms'                                                  
        write(200,*) 1, '#Diferent atoms types = Number of alpha'                                                     
                                                                         
        !Calculate number of bonds
        nBond = nRingSide*(nRingHeight*5)                                
        nBond = nBond + ((nRingSide-1)*nRingHeight)                      
                                                                         
        write(200,*) nBond, '#Total number of conjugated bonds'                                                 
        write(200,*) 1,'#Total number of types of conjugated atoms = number of beta'                                                     
        write(200,*) ''
                                                                         
        write(200,*) '#ALL ALPHA VALUES (in order)'                        
        write(200,*) alpha        !all alphas are equal                                         
        write(200,*) ''
                                                                         
        write(200,*) '#ALL BETA VALUES (in order)'                         
        write(200,*) beta         !all betas are equal 
        write(200,*) ''
        allocate(con(nBond))
        
        write(200,*) '#BONDS'
        
        !Number of atoms in one big ring (circular)
        atmCol = nRingHeight * 2
        
        !Counter of atoms
        atm=1
        
        !For each circular
        do j = 1, nRingSide*2

            !Create the circular structure 
            do i=1, atmCol-1
                write(200,*) atm, atm+1, '1'
                atm = atm + 1
            enddo

            write(200,*) atm, atm-(2*nRingHeight)+1, '1'
            atm = atm+1

        enddo
        
        !Counter of atoms
        atm=1

        !For each circular minus one
        do i = 1, nRingSide*2-1

            !Make needed bonds between circulars
            do j=1, nRingHeight
                write(200,*) atm, atm+(nRingHeight*2), '1'
                atm=atm+2
            enddo

            !Alternate differently to make corrects bonds between circular strucruers
            if (mod(i,2) .eq.  1 ) then
                atm=atm+1
            else
                atm=atm-1
            endif

        enddo

        !All atoms have the same alphas
        write(200,*) ''
        write(200,*) '#ATOM TYPE (alphas)'
        do i = 1, nAtm
            write(200,*) 1
        enddo

        close(200)
    

    end subroutine nanoTubeInp                                   


    subroutine polyeneInp

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         January, 2022          +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Creates the input for a polyene using the parameters in 2polyene.inp
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      integer         :: nAtom, nBond
      integer         :: i, j
      real            :: alpha, betap, betam
      logical         :: pbc
      character       :: junk


      !Open intput file for the polyene
      open(unit=400, file='2polyene.inp', status='old')

      read(400,*) junk
      read(400,*) pbc !periodic boundary conditions?
      read(400,*) junk

      read(400,*) nAtom !number of atoms
      read(400,*) alpha !just one alpha possible
      read(400,*) betap !beta plus
      read(400,*) betam !beta minus

      close(400)


      !Check if given the contidions is possible to create the structure
      if (pbc .and. mod(nAtom,2) .ne. 0) then
          write(*,*) '* * * W A R N I N G * * *'
          write(*,*) 'If you want to use periodic boundary conditions (pbc), the total number of atoms must be even'
          stop
      endif

      !Open the explicit file to create it
      open(unit=100, file='0explicit.inp', status='unknown')

      write(100,*) 'input title'                                         
      write(100,*) ''
      write(100,*) nAtom, '#Total number of atoms' 
      write(100,*) 1, '#Diferent atoms types = Number of alpha'

      !Calculate number of bonds
      if (pbc) then
          nBond = nAtom 
      else
          nBond = nAtom - 1
      endif

      write(100,*) nBond,'#Total number of conjugated bonds'
      write(100,*) 2, '#Total number of types of conjugated atoms = number of beta'

      write(100,*) ''
      write(100,*) '#ALL ALPHA VALUES (in order)'
      write(100,*) alpha

      write(100,*) ''
      write(100,*) ' #ALL BETA VALUES (in order)'
      write(100,*) betap
      write(100,*) betam

      write(100,*) ''
      write(100,*) '#BONDS'
      !write the pait of atoms that makes all the bonds
      do i = 1, nAtom-1
          write(100,*) i, i+1 , mod(i,2)+1
      enddo
      if (pbc) write(100,*) nAtom, 1, mod(i,2)+1 

      !all alphas are the same
      write(100,*) ''
      write(100,*) '#ATOM TYPE (alphas)'
      do i = 1, nAtom
          write(100,*) 1
      enddo

      close(100)


    end subroutine polyeneInp

    subroutine trianguleneInp

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         January, 2022          +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Creates the input for a triangulene using the parameters in 3triangulene.inp
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      integer                               :: nAtom, nBond, L, atomLess, bondLess, nAtomPast
      integer, dimension(3)                 :: lm
      integer                               :: innAtom, outnAtom
      integer                               :: cont, comp, atm, bondCont
      integer                               :: i, j
      integer, allocatable, dimension(:,:)  :: allBonds
      integer, allocatable, dimension(:)    :: allAtomLess
      real                                  :: alpha, beta
      character                             :: junk
 

      !Open intput file for the triangulene
      open(unit=400, file='3triangulene.inp', status='old')

      !( L | lm(1), lm(2), lm(3) )
      read(400,*) junk
      read(400,*) L
      read(400,*) lm(1)
      read(400,*) lm(2)
      read(400,*) lm(3)

      !All alphas are equal. All betas ate equal
      read(400,*) alpha
      read(400,*) beta

      close(400)

      !Check if the combination of lambdas is possible
      call lambdasCorrect(L,lm)

      !Calculates the total number of Atoms and Bonds
      nAtom = (2*L) + 1 
      nBond = 0
      do i = 1, L
          nAtom = nAtom + (2*i) + 1
          nBond = nBond + 3*(i+1)
      enddo

      allocate(allBonds(nBond,3))

      !Creates the complete triangulene (L|0,0,0)
      call completeTriangulene(allBonds, nBond, nAtom,L)      

      !Calculates the total number of bonds that should be eliminated
      atomLess = 0
      bondLess = 0
      do i = 1, 3
          if (lm(i) .gt. 0) then
              cont = 4
          endif
          do j = 1, lm(i)
              atomLess = atomLess + ((2*j) + 1)
              bondLess = bondLess + cont 
              cont = cont + 3
          enddo
      enddo

      allocate(allAtomLess(atomLess))

      !Search for the atoms that should be eliminated given the numeration
      call atomsToEliminate(atomLess, allAtomLess, lm, L, nAtom)

      !Create the new triangulene by eliminating the needed atoms
      call newTriangulene(nBond, bondLess, allBonds, allAtomLess, atomLess)

      !Shift the numeration to have all the atoms numbered in a continium range
      call shiftTrianguleneNumbering(nBond, nAtom, allBonds, atomLess, allAtomLess)

      !Open the explicit input file to create it
      open(unit=100, file='0explicit.inp', status='unknown')

      write(100,*) 'input title'
      write(100,*) ''
      write(100,*) nAtom-atomLess, '#Total number of atoms'
      write(100,*) 1, '#Diferent atoms types = Number of alpha'
      write(100,*) nBond, '#Total number of conjugated bonds'
      write(100,*) 1, '#Total number of types of conjugated atoms = number of beta'
      write(100,*) ''
      write(100,*) '#ALL ALPHA VALUES (in order)'
      write(100,*) alpha

      write(100,*) ''
      write(100,*) ' #ALL BETA VALUES (in order)'
      write(100,*) beta

      write(100,*) ''
      write(100,*) ' #BONDS'
      call writeBonds(nBond, allBonds)
 
      write(100,*) ''
      write(100,*) '#ATOM TYPE (alphas)'
      do i = 1, nAtom
          write(100,*) 1
      enddo

      deallocate(allBonds)
      deallocate(allAtomLess)

      close(100)

    end subroutine trianguleneInp

    subroutine shiftTrianguleneNumbering (nBond,nAtom, allBonds, atomLess, allAtomLess)

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         January, 2022          +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Displaces the indexing to number all the atoms in a continous range from 1 to the total number of atoms
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        integer                                 :: nBond, nAtom, atomLess
        integer                                 :: i, j, k, l, cont
        integer, allocatable, dimension(:,:)    :: allBonds
        integer, allocatable, dimension(:)      :: dummy
        integer, allocatable, dimension(:)      :: allAtomLess


        !Shift all bonds by substracting the lower atom plus 1 (to start the numbering in 1)
        allBonds(:,1:2) = allBonds(:,1:2) - minval(allBonds(:,1:2))+1

        !Dummy vector to store auxiliar indexes of size of the maximum index of all atoms in the actual numbering
        allocate(dummy(maxval(allBonds(:,1:2))))


        !Initialize dummy matrix
        dummy(:) = 0
        do i = 1, size(dummy)
            do j = 1, nBond
                do k = 1, 2
                    !Fill dummy matrix. 0 if atom eliminated, if not, then the number of the atom
                    if (i .eq. allBonds(j,k)) then
                        dummy(i) = i
                    endif
                enddo
            enddo
        enddo

        !initialize auxiliar counter
        cont = sum(dummy(:))

        !while there are zeros in the indexes that should be different of zeros
        do while( cont .gt. sum(dummy(1:nAtom-atomLess)) )

        do i = 2, size(dummy) !first index should not be displaced because it has already been displaced

            if (dummy(i-1) .eq. 0 ) then

                !Displace atom one place
                dummy(i-1) = dummy(i)
                dummy(i) = 0

            endif
        enddo

        enddo

        !Change the bonds with the new numeration
        do i =1, nBond
            do j = 1, 2
                do k = 1, size(dummy)
                    if (allBonds(i,j) .eq. dummy(k)) allBonds(i,j) = k
                enddo
            enddo
        enddo

        deallocate(dummy)

    end subroutine shiftTrianguleneNumbering 


    subroutine newTriangulene(nBond, bondLess, allBonds, allAtomLess, atomLess)

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         February, 2022         +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Createst the triangulene by eliminating the needed atoms 
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      integer                                   :: nBond, bondLess, atomLess
      integer                                   :: i, j, k
      integer, allocatable, dimension(:,:)      :: allBonds, dummyallBonds
      integer, allocatable, dimension(:)        :: allAtomLess

      !Auxiliar dummy matrix with all the bonds and its betas
      allocate(dummyallBonds(nBond-bondLess,3))

      !Look for the bonds thath should be eliminated dependind on the indexes of the atoms
      do i = 1, nBond
          do j = 1 ,2
              do k = 1, atomLess
                  !if the bond contain one atom that should be eliminated
                  if (allBonds(i,j) .eq. allAtomLess(k) ) then
                      !index of its beta is equal to 0; becuas it is index, it should be > 0 if it exist
                      allBonds(i,3) = 0
                      exit
                  endif
              enddo
          enddo
      enddo

      !auxiliar counter for indexing
      k = 1
      do i = 1, nBond
          !If the bond does exist
          if (allBonds(i,3) .ne. 0) then
              !copy thje bond
              dummyallBonds(k,:) = allBonds(i,:)
              k = k + 1
          endif
      enddo

      deallocate(allBonds)

      !new number of bonds
      nBond = nBond - bondLess

      !Reshape the size of the matrix that stores the bonds
      allocate(allBonds(nBond,3))

      !Copy alll the bonds that does exist
      allBonds(:,:) = dummyallBonds(:,:)

      deallocate(dummyallBonds)

    end subroutine newTriangulene

    subroutine atomsToEliminate(atomLess, allAtomLess, lm, L, nAtom)

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         February, 2022         +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Calculates the index of the atoms that should be eliminted 
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        integer                             :: atomLess, L, nAtom
        integer                             :: i, j, k, cont, inf, sup, supCont, infCont
        integer, allocatable, dimension(:)  :: allAtomLess
        integer, dimension(3)               :: lm

        !Initialize index conuter of atoms that shouls be aliminated
        cont = 1

        !For lm(1). Corner that contains the atom numer 1
        if (lm(1) .gt. 0) then !start to "eliminate" the necesary atoms

            !fisrt loop, form where to where should be eliminated
            inf = 1 !because is the corner that dcontains the atom 1
            sup = 2*lm(1) 

            !initialize counters
            supCont = 0
            infCont = 0

            !loop over diagonal chains that has to be eliminated
            do i= 1, lm(1)+1 
                !loop over the atoms that should be eliminated
                do j = inf, sup 
                    allAtomLess(cont) = j !add tot he list the atom that should be eliminated
                    !write(*,*) inf, sup, cont, j
                    cont = cont + 1
                enddo
                inf = inf + (2*L) + 1 - infCont
                sup = inf + (2*lm(1)) - 2 - supCont
                supCont = supCont + 2
                if (i .gt. 1 ) infCont = infCont + 2
            enddo
        endif 

        !for lm(2). Corner that contains the atom numer 2L+1
        if (lm(2) .gt. 0) then

            !fisrt loop, form where to where should be eliminated
            inf = (2*L)+1 !because is the corner that dcontains the atom 2L + 1 
            sup = (2*L) + 2 - (2*lm(2))

            !initialize counters
            supCont = 2
            infCont = 0

            !loop over diagonal chains that has to be eliminated
            do i= 1, lm(2)+1 
                do j = inf, sup, -1 
                    allAtomLess(cont) = j !add tot he list the atom that should be eliminated
                    !write(*,*) inf, sup, cont, j
                    cont = cont + 1
                enddo
                inf = inf + (2*L) + 1 - infCont
                sup = inf - (2*lm(2)) + supCont
                supCont = supCont + 2
                infCont = infCont + 2
            enddo
        endif

        !for lm(3). Corner that contains the atom numer nAtom-1
        if (lm(3) .gt. 0) then

            inf = nAtom !because the corner starts with nAtom-1
            sup = nAtom-2

            supCont=0
            do i = 2, lm(3)
                supCont =  (2*i) + 1
                sup = sup - supCont 
            enddo

            do i = inf, sup, -1
                    allAtomLess(cont) = i !add tot he list the atom that should be eliminated
                    cont = cont + 1
            enddo
        endif

    end subroutine atomsToEliminate

    subroutine completeTriangulene (allBonds, nBond, nAtom, L)                                   

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         February, 2022         +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Creates the bonds for the triangulene of the shape (L|0,0,0)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      integer                                       :: nBond, nAtom, L
      integer                                       :: cont, comp, atm, bondCont
      integer                                       :: i, j, nb
      integer, allocatable, dimension(:,:)          :: allBonds

      !Bonds for one of the sides of the triangulene
      nb=1 !auxiliar counter for indexing

      do i = 1, (L*2)
          allBonds(nb,1) = i
          allBonds(nb,2) = i+1
          allBonds(nb,3) = 1 !beta of each bond 
          nb = nb + 1
      enddo

      !Auxiliar counters to check if the atoms should be conected 
      cont = 0 
      comp = (L*2) 

      do i = (L*2)+2 , nAtom 
          !Bonds for each chain that is in the tringulene 
          if (cont .lt. comp) then
              allBonds(nb,1) = i
              allBonds(nb,2) = i+1
              allBonds(nb,3) = 1 !beta for each bond
              nb = nb + 1
              cont = cont + 1
          else
              cont = 0
              comp = comp - 2
          end if
      enddo 


      !auxiliar bond counter
      bondCont = (L * 2) + 1

      !atom indexing counter
      atm = -1

      !Bonds of the firts chain with the second chain
      do j = 1, L+1
          atm = atm + 2
          allBonds(nb,1) = atm
          allBonds(nb,2) = atm + bondCont
          allBonds(nb,3) = 1 !Beta for each bond
          nb = nb + 1
      enddo

      !Auciliar counters
      bondCont = bondCont - 1
      cont=L+1

      !Bonds of the following chains with the next chain if exist
      do i = 2,L
          cont = cont - 1
          do j = 1, cont 
              atm = atm + 2
              allBonds(nb,1) = atm
              allBonds(nb,2) = atm + bondCont
              allBonds(nb,3) = 1
              nb = nb + 1
          enddo
          atm = atm + 1
          bondCont = bondCont - 2
      enddo

    end subroutine completeTriangulene                                                                                 

    subroutine writeBonds(nBond, allBonds)

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         February, 2022         +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Write the bonds and its beta values
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        integer                                 :: nBond
        integer                                 :: i , j
        integer, allocatable, dimension(:,:)    :: allBonds

        do i = 1, nBond
            write(100,*) allBonds(i,1), allBonds(i,2), allBonds(i,3)
        enddo

    end subroutine writeBonds

    subroutine lambdasCorrect(L, lm)

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! +         February, 2022         +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Check if the system defined is possible
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        integer                     :: L
        integer, dimension(3)       :: lm
        integer                     :: i, j

        do i = 1, 3
            if (lm(i) .gt. L) then
                write(*,*) 'W A R N I N G !!!'
                write(*,*) 'l1, l2 or l3 cannot be bigger than L'
                stop
            endif
        enddo

        do i= 1, 2
            do j = i+1 , 3
                if (lm(i)+lm(j) .gt. L ) then
                    write(*,*) 'W A R N I N G !!!'
                    write(*,*) 'sum of two l cannot be bigger than L'
                    stop
                endif
            enddo
        enddo 

    end subroutine lambdasCorrect

    subroutine matrixMake(mat, nAtom, alphas, nBond, con, betas)

      ! ++++++++++++++++++++++++++++++++++
      ! +    Osvaldo Hernandez Cuellar   +
      ! + Universidad Autonoma de Madrid +
      ! +         January, 2022         +
      ! ++++++++++++++++++++++++++++++++++

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Creates the H matrix using the bonded atoms and the alpha and beta values
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        integer(ik)                               :: i, j
        integer(ik)                               :: nAtom, nBond !number of atoms, number of bonds, number of pi electrons
        integer(ik), allocatable, dimension(:,:)  :: con  !bonds, atom conections usign indexes for every atom  
        real(rk), allocatable, dimension(:,:)     :: mat
        real(rk), allocatable, dimension(:)       :: alphas, betas

        allocate (mat(nAtom,nAtom))

        mat = 0.0_rk

        !diagonal elements
        do i = 1, nAtom
            mat(i,i) = alphas(i)
        enddo

        !atoms bonded
        do i = 1, nBond
            mat(con(i,1), con(i,2)) = betas(con(i,3))
            mat(con(i,2), con(i,1)) = betas(con(i,3))
        enddo

    end subroutine matrixMake

end module start
