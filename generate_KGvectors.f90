!Program to generate 'K+G' vectors which will then be used to calculate the 'gk' vectors for the primitive cell
program generate_KGvectors
  implicit none
  !
  real*4,allocatable :: K(:,:),G(:,:),KG(:)
  !K :: the k vectors should be in the cartesian co-ordinates in units of Rydberg
  integer :: nKs,nGs,iter,iter1
  !nKs,nGs :: number of K and G points respectively
  !iter :: iterative term
  allocate(K(3,nKs),G(3,nGs),KG(3))
  !
  open(204,file='./data/sc_Kpoints.in')
  do iter = 1,nKs
    read(204,*)K(:,iter)
  end do
  close(204)
  open(205,file='./data/sc_Gvectors.in')
  do iter = 1,nGs
    read(205,*)G(:,iter)
  end do
  close(205)
  !
  KG = 0
  open(206,file='./data/KplusG.in')
  do iter = 1,nKs
    do iter1 = 1,nGs
      KG(1) = K(1,iter) + G(1,iter1)
      KG(2) = K(2,iter) + G(2,iter1)
      KG(3) = K(3,iter) + G(3,iter1)
      write(206,*)KG
    end do
  end do
  close(206)
  deallocate(K)
  deallocate(G)
  deallocate(KG)
end program generate_KGvectors 
