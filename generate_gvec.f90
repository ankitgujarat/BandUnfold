program generatevec
  implicit none
  !
  !
  real :: ecutwfc
  !ecutwfc : cutoff energy that is obtained from the total ebergy convergence of the parent unit cell or the primitive cell
  integer :: m1,m2,m3,m1_max,m2_max,m3_max,m_max
  !m1,m2,m3 :: the coefficients of the the reciprocal lattice vector
  !m_max :: the maximum range of the coeeficients of the reciprocal lattice vector
  real*4,allocatable :: bg(:,:),gv(:,:)
  !bg : reciprocal lattice vector or vector describing the Brillouin zone
  !gv : linear combination of the lattice vectors described above
  integer :: ng,i
  !ng : number of gvectors satisfying the cutoff condition
  !i : an iteration variable
  !
  !
  !g vectors are in cartesian co-ordinate
  open(199,file='./data/input.in')
  read(199,*)
  read(199,*)ecutwfc
  read(199,*)
  !
  !
  allocate(bg(3,3),gv(3,1))
  !
  !
  do i = 1,3
    read(199,*)bg(:,i)
  end do
  close(199)
  !
  open(200,file='./data/gvecs.dat')
  !
  m1_max = nint((4*ecutwfc)/((bg(1,1))**2 + (bg(2,1))**2 + (bg(3,1))**2))
  m2_max = nint((4*ecutwfc)/((bg(1,2))**2 + (bg(2,2))**2 + (bg(3,2))**2))
  m3_max = nint((4*ecutwfc)/((bg(1,3))**2 + (bg(2,3))**2 + (bg(3,3))**2))
  !
  !find minimum of the above three values
  !debugging portion
  !PRINT *,M1_MAX,M2_MAX,M3_MAX
  !STOP
  !
m_max = min(m1_max,m2_max,m3_max)
  !
  ng = 0
  do m1 = -m_max,m_max
    do m2 = -m_max,m_max
      do m3 = -m_max,m_max
        gv(:,1) = m1*bg(:,1) + m2*bg(:,2) + m3*bg(:,3)
        if (gv(1,1)**2 + gv(2,1)**2 + gv(3,1)**2 .le. 4*ecutwfc) then
          ng = ng + 1
          write(200,*)gv(:,1),m1,m2,m3
        end if
      end do
    end do
  end do
  print *,'Total number of gvectors generated:',ng
  close(200)
  deallocate(bg)
  deallocate(gv)
end program generategvec
