!Code to generate a list of 'g' vectors for each 'k' vector in the primitve cell
program generate_g_vec_list
  implicit none
  !
  !
  real*4,allocatable :: k(:,:),g(:,:)
  real*4 :: ecutwfc
  integer :: nks,ngs,nk,ng,c
  real*4,dimension(3) :: kg
  !
  allocate(k(3,nks),g(6,ngs))
  !
  open(207,file='./data/vectors_g_k.in')
  do nk = 1,nks
    c = 0
    do ng = 1,ngs
      kg(1) = k(1,nk) + g(1,ng)
      kg(2) = k(2,nk) + g(2,ng)
      kg(3) = k(3,nk) + g(3,ng)
      if ((kg(1)**2 + kg(2)**2 + kg(3)**2) .le. ecutwfc) then
        c = c + 1
        write(207,*),g(4,ng),g(5,ng),g(5,ng)
      end if
    end do
    write(207,*),c
  end do
  deallocate(k)
  deallocate(g)
  close(207)
end program generate_g_vec_list
