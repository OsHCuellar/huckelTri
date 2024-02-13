module dsyevDiag 

    use types

    implicit none

    contains
!
        subroutine subdsyev(mat,W,d)
!        external dsyev
!         this subroutine is to do the diagnonization of a symmetric matrix 
!         mat is the matrxi and d is the dimension of the mat while W is to store the output eigen value.
       integer ::   k, l,LW,INFO,d   ,i
       real(rk), dimension(d,d) :: mat
       real(rk), dimension(d) :: w
       real(rk), allocatable :: ap(:,:), u(:,:),VL(:,:),VR(:,:)
       real(rk), allocatable :: WORK(:),WR(:),WI(:)     
      allocate (ap(d,d))   
      allocate  (WR(1))
        ap(:,:)=mat(:,:) 
        LW=-1
        
       

!         give the value in the upper part of the mat to AP matrix
                 do k=1,d
                   AP(K,K)=mat(K,K)
                  do l=K,d
                   ap(k,l)=mat(k,l)
                  enddo
                 enddo
                   
                           
!       diagonalization of the mat matrix
        call DSYEV('V','U',d,ap,d,W,WR,LW,INFO)  
        LW=WR(1)
        deallocate (WR)
        allocate (WR(lw))
        call DSYEV('V','U',d,ap,d,W,WR,LW,INFO)
        
!OUTPUT PART(optational)
!                write (6,*) 'orthonomal eigenvector'              
!                write (6,*)
!                write (6,*) 
!               write (6,300) W
!                write (6,*) ""
!                write (6,*) ""
!        DO k=1,d

                
!                write (6,'(*(f10.1))')  ap(k,:)
 
!        enddo
!                write (6,*) 'eigenvalue'
!                 do k=1,d
!                write (6,"(i5,20f10.1)") k, W(k)
!                enddo

!100             format (i5,20f9.1)
!300             format(f10.1)
        !return
        end subroutine subdsyev
!
        !subroutine subdsyevv(mat,W,d,VR)
        subroutine subdsyevv(mat,d,VR)
!        external dsyev
!         this subroutine is to do the diagnonization of a symmetric matrix 
!         mat is the matrxi and d is the dimension of the mat while W is to store the output eigen value.
       integer ::   k, l,LW,INFO,d   
       real(rk), dimension(d,d) :: mat
       real(rk), dimension(d) :: w
       real(rk), dimension(:,:), allocatable ::VR(:,:)
       real(rk), allocatable :: ap(:,:), u(:,:),VL(:,:)
       real(rk), allocatable :: WORK(:),WR(:),WI(:)     
      allocate (ap(d,d))   
      allocate  (WR(1))
        ap(:,:)=mat(:,:) 
        LW=-1
        
       

!         give the value in the upper part of the mat to AP matrix
                 do k=1,d
                   AP(K,K)=mat(K,K)
                  do l=K,d
                   ap(k,l)=mat(k,l)
                  enddo
                 enddo
                   
                           
!       diagonalization of the mat matrix
      
        call DSYEV('V','U',d,ap,d,W,WR,LW,INFO)  
        LW=WR(1)
        deallocate (WR)
        allocate (WR(lw))
        call DSYEV('V','U',d,ap,d,W,WR,LW,INFO)
        VR(:,:)=ap(:,:)
        
        do k=1,d
            mat(k,k) = W(k)
        enddo

!OUTPUT PART(optational)
!                write (6,*) 'orthonomal eigenvector'              
!                write (6,*)
!                write (6,*) 
!               write (6,300) W
!                write (6,*) ""
!                write (6,*) ""
!        DO k=1,d

                
!                write (6,'(*(f10.1))')  ap(k,:)
 
!        enddo
!                write (6,*) 'eigenvalue'
!                 do k=1,d
!                write (6,"(i5,20f10.1)") k, W(k)
!                enddo

!100             format (i5,20f9.1)
!300             format(f10.1)
        return
        end subroutine subdsyevv

end module 
