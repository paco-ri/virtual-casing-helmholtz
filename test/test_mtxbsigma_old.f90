! Script modeled after test_gradcurllap.f90
! See 'Boundary integral equation analysis on the sphere' by Vico et al 

implicit real *8 (a-h,o-z)
  
real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
real *8, allocatable :: wts(:)
character *100 fname
integer ipars(2)

integer, allocatable :: norders(:),ixyzs(:),iptype(:)
integer, allocatable :: ipatch_id(:)
real *8, allocatable :: uvs_src(:,:)
integer, allocatable :: novers(:),ixyzso(:)
integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
real *8, allocatable :: srcover(:,:),wover(:)

integer ndd,ndz,ndi,nnz_gc,nnz_helm
real *8 dpars
complex *16 zpars(3)
integer ipars

real *8, allocatable :: targs(:,:)
integer, allocatable :: ipatch_id_targ(:)
real *8, allocatable :: uvs_targ(:,:)

real *8, allocatable :: hvec(:,:)
real *8, allocatable :: rhs(:)
real *8, allocatable :: hvecdivreal(:),hvecdivimag(:)

real *8, allocatable :: ynmreal(:),ynmimag(:)
real *8, allocatable :: unmreal(:,:),unmimag(:,:),xnmreal(:,:),xnmimag(:,:)
complex *16, allocatable :: ynm(:),unm(:,:),xnm(:,:)

real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
complex *16, allocatable :: rjvec(:,:),rho(:)
complex *16, allocatable :: curlj(:,:),gradrho(:,:)
complex *16, allocatable :: curlj_targ(:,:),gradrho_targ(:,:)
complex *16 zk
complex *16, allocatable :: hnk(:),hnkp(:)
real *8, allocatable ::jnk(:),jnkp(:)
complex *16, allocatable :: rhnk(:),rhnkp(:),rjnk(:),rjnkp(:)
integer mn

complex *16 runc, rxnc, rund, rxnd, rynd

real *8, allocatable :: wnear(:)
real *8, allocatable :: w(:,:)
complex *16 vtmp1(3),vtmp2(3),wtmp
real *8 xyz_out(3),xyz_in(3)
complex *16 zpars
complex *16 ima

data ima/(0.0d0,1.0d0)/

call prini(6,13)

done = 1
pi = atan(done)*4

! set up geometry
igeomtype = 1
ipars(1) = 3!2
npatches = 12*(4**ipars(1)) 
norder = 5
npols = (norder+1)*(norder+2)/2
npts = npatches*npols
allocate(srcvals(12,npts),srccoefs(9,npts))
ifplot = 0
call setup_geom(igeomtype,norder,npatches,ipars, &
   srcvals,srccoefs,ifplot,fname)
allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
do i=1,npatches
  norders(i) = norder
  ixyzs(i) = 1 +(i-1)*npols
  iptype(i) = 1
enddo
ixyzs(npatches+1) = 1+npols*npatches

! get and test quadrature weights
allocate(wts(npts))
call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)
ra = 0
do i=1,npts
  ra = ra + wts(i)
enddo
call prin2('surface area of sphere=*',ra,1)

! compute spherical harmonics 
nn = 3
mm = 2
nmax = nn
allocate(w(0:nmax,0:nmax))
allocate(ynm(npts),ynmreal(npts),ynmimag(npts))
allocate(unm(3,npts),unmreal(3,npts),unmimag(3,npts))
allocate(xnm(3,npts),xnmreal(3,npts),xnmimag(3,npts))
call l3getsph(nmax,mm,nn,12,srcvals,ynm,npts,w)
! zynm <-- Y_n^m, n=1..nn, m=1..mm (each point)
! w <-- Y_nn^mm (last point) 
do i=1,npts
   ynmreal(i) = real(ynm(i))
   ynmimag(i) = aimag(ynm(i))
enddo

! Unm <-- surface gradient of Ynm, scaled
! Xnm <-- cross product of unit normal and Unm
call surf_grad(npatches,norders,ixyzs,iptype,npts,srccoefs,&
     srcvals,ynmreal,unmreal)
call surf_grad(npatches,norders,ixyzs,iptype,npts,srccoefs,&
     srcvals,ynmimag,unmimag)
unm = dcmplx(unmreal,unmimag)
do i=1,npts
  unm(1:3,i) = unm(1:3,i)/sqrt(nn*(nn+1.0d0))
  call dzcross_prod3d(srcvals(10,i),unm(1,i),xnm(1,i))
enddo
xnmreal = real(xnm)
xnmimag = aimag(xnm)

! surface divergence test
allocate(hvecdivreal(npts),hvecdivimag(npts))
call surf_div(npatches,norders,ixyzs,iptype,npts,srccoefs,&
     srcvals,xnmreal,hvecdivreal)
call surf_div(npatches,norders,ixyzs,iptype,npts,srccoefs,&
     srcvals,xnmimag,hvecdivimag)
ra = 0
do i=1,npts
  ra = ra + (hvecdivreal(i)**2 + hvecdivimag(i)**2)*wts(i)
enddo
ra = sqrt(ra)
call prin2('error in surface divergence of xnm=*',ra,1)

! create sources
allocate(rjvec(3,npts),rho(npts),curlj(3,npts),gradrho(3,npts))
rcu = 0.0d0 ! 1.1d0
rcx = 1.0d0 ! 0.0d0
do i=1,npts
   ! check this!!!
   rjvec(1:3,i) = rcu*unm(1:3,i) + rcx*xnm(1:3,i)
   rho(i) = ynm(i)
enddo

eps = 1.0d-8
allocate(ipatch_id(npts),uvs_src(2,npts))
call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, &
     ipatch_id,uvs_src)

! compute layer potentials
zk = 1.0d0
dk = real(zk)
call prin2('zk = *', zk, 2)
call lpcomp_gradcurlhelm(npatches,norders,ixyzs,&
     iptype,npts,srccoefs,srcvals,12,npts,srcvals,ipatch_id,uvs_src,&
     eps,zk,rjvec,rho,curlj,gradrho)


! special function evaluations
allocate(jnk(0:nn),jnkp(0:nn),hnk(0:nn),hnkp(0:nn))
call SPHJ(nn,dk,mn,jnk,jnkp)
call SPHH(nn,dk,mn,hnk,hnkp)

allocate(rjnk(nn),rjnkp(nn),rhnk(0:nn),rhnkp(0:nn))
call riccati_jp(nn,dk,rjnk,rjnkp)
call riccati_hp(nn,dk,rhnk,rhnkp)

! test layer potential computation
errnc = 0
errnd = 0
rnc = 0
rnd = 0

! curl Sk tests

! coefficients, eq. (54), third line
! see fmm3dbie/src/maxwell/analytic_sphere_pw_pec:riccati_jp(), etc.
runc = ima*rjnkp(nn)*rhnk(nn)*rcu
rxnc = -ima*rjnk(nn)*rhnkp(nn)*rcx

! eqn. (62) 
rund = 0
rxnd = -ima*rjnk(nn)*rhnk(nn)*sqrt(nn*(nn+1.0d0))/zk*rcx
call prin2('rxnd = *', rxnd, 2)
rxnd = -ima*zk*jnk(nn)*hnk(nn)*sqrt(nn*(nn+1.0d0))*rcx
call prin2('jnk(nn) = *', jnk(nn), 2)
call prin2('hnk(nn) = *', hnk(nn), 2)
call prin2('rxnd = *', rxnd, 2)

do i=1,npts
   call dzcross_prod3d(srcvals(10,i),curlj(1,i),vtmp1)
   vtmp1(1:3) = vtmp1(1:3) + rjvec(1:3,i)/2

   vtmp2(1) = runc*unm(1,i) + rxnc*xnm(1,i) 
   vtmp2(2) = runc*unm(2,i) + rxnc*xnm(2,i) 
   vtmp2(3) = runc*unm(3,i) + rxnc*xnm(3,i)

   rnc = rnc + (real(vtmp2(1))**2 + aimag(vtmp2(1))**2 + &
        real(vtmp2(2))**2 + aimag(vtmp2(2))**2 + &
        real(vtmp2(3))**2 + aimag(vtmp2(3))**2)*wts(i)
   errnc = errnc + abs(vtmp1(1)-vtmp2(1))**2*wts(i)
   errnc = errnc + abs(vtmp1(2)-vtmp2(2))**2*wts(i)
   errnc = errnc + abs(vtmp1(3)-vtmp2(3))**2*wts(i)
   
   wtmp = 0
   call dzdot_prod3d(srcvals(10,i),curlj(1,i),wtmp)

   rnd = rnd + (real(ynm(i))**2 + aimag(ynm(i))**2)*wts(i)
   errnd = errnd + ((real(wtmp) - real(rxnd*ynm(i)))**2 + &
       (aimag(wtmp) - aimag(rxnd*ynm(i)))**2)*wts(i)
enddo

errnc = sqrt(errnc/rnc)
errnd = sqrt(errnd/rnd)

call prin2('rnc=*',rnc,1)
call prin2('rnd=*',rnd,1)

call prin2('error in n times curl sk = *',errnc,1)
call prin2('error in n dot curl s0 =*',errnd,1)

! grad Sk tests

! eqn. (61)
runc = 0
rxnc = ima*zk*jnk(nn)*hnk(nn)*sqrt(nn*(nn+1.0d0))

! eqn. (29)
rynd = ima*zk**2*jnk(nn)*hnkp(nn)
call prin2('rynd = *', rynd, 2)

errnc = 0
errnd = 0
rnc = 0
rnd = 0

do i=1,npts
    vtmp1 = 0
    vtmp2 = 0
    call dzcross_prod3d(srcvals(10,i),gradrho(1,i),vtmp1)
    vtmp2(1:3) = rxnc*xnm(1:3,i)
    rnc = rnc + abs(vtmp2(1))**2*wts(i)
    rnc = rnc + abs(vtmp2(2))**2*wts(i)
    rnc = rnc + abs(vtmp2(3))**2*wts(i)
    ! rnc = rnc + (vtmp2(1)**2 + vtmp2(2)**2 + vtmp2(3)**2)*wts(i)
    errnc = errnc + aimag(vtmp1(1)-vtmp2(1))**2*wts(i)
    errnc = errnc + aimag(vtmp1(2)-vtmp2(2))**2*wts(i)
    errnc = errnc + aimag(vtmp1(3)-vtmp2(3))**2*wts(i)

    wtmp = 0
    call dzdot_prod3d(srcvals(10,i),gradrho(1,i),wtmp)
    wtmp = wtmp - rho(i)/2
    errnd = errnd + ((real(wtmp)-real(rynd*ynm(i)))**2 &
         + (aimag(wtmp)-aimag(rynd*ynm(i)))**2)*wts(i)
    rnd = rnd + (real(ynm(i))**2+aimag(ynm(i))**2)*wts(i)
enddo

errnc = sqrt(errnc/rnc)
errnd = sqrt(errnd/rnd)
call prin2('error in n times grad s0 = *',errnc,1)
call prin2('error in n dot grad s0 =*',errnd,1)

return
end


subroutine setup_geom(igeomtype,norder,npatches,ipars, & 
   srcvals,srccoefs,ifplot,fname)
implicit real *8 (a-h,o-z)
integer igeomtype,norder,npatches,ipars(*),ifplot
character (len=*) fname
real *8 srcvals(12,*), srccoefs(9,*)
real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

real *8, pointer :: ptr1,ptr2,ptr3,ptr4
integer, pointer :: iptr1,iptr2,iptr3,iptr4
real *8, target :: p1(10),p2(10),p3(10),p4(10)
real *8, allocatable, target :: triaskel(:,:,:)
real *8, allocatable, target :: deltas(:,:)
integer, allocatable :: isides(:)
integer, target :: nmax,mmax

procedure (), pointer :: xtri_geometry


external xtri_stell_eval,xtri_sphere_eval,xtri_wtorus_eval

npols = (norder+1)*(norder+2)/2
allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
allocate(wts(npols))

call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

if(igeomtype.eq.1) then
  itype = 2
  allocate(triaskel(3,3,npatches))
  allocate(isides(npatches))
  npmax = npatches
  ntri = 0
  call xtri_platonic(itype, ipars(1), npmax, ntri, & 
     triaskel, isides)

  xtri_geometry => xtri_sphere_eval
  ptr1 => triaskel(1,1,1)
  ptr2 => p2(1)
  ptr3 => p3(1)
  ptr4 => p4(1)


  if(ifplot.eq.1) then
     call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, & 
        ptr3,ptr4, norder,'Triangulated surface of the sphere')
  endif


  call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4, &
    npols,uvs,umatr,srcvals,srccoefs)
endif

if(igeomtype.eq.2) then
  done = 1
  pi = atan(done)*4
  umin = 0
  umax = 2*pi
  vmin = 2*pi
  vmax = 0
  allocate(triaskel(3,3,npatches))
  nover = 0
  call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2), &
    nover,npatches,npatches,triaskel)

  mmax = 2
  nmax = 1
  xtri_geometry => xtri_stell_eval

  allocate(deltas(-1:mmax,-1:nmax))
  deltas(-1,-1) = 0.17d0
  deltas(0,-1) = 0
  deltas(1,-1) = 0
  deltas(2,-1) = 0

  deltas(-1,0) = 0.11d0
  deltas(0,0) = 1
  deltas(1,0) = 4.5d0
  deltas(2,0) = -0.25d0

  deltas(-1,1) = 0
  deltas(0,1) = 0.07d0
  deltas(1,1) = 0
  deltas(2,1) = -0.45d0

  ptr1 => triaskel(1,1,1)
  ptr2 => deltas(-1,-1)
  iptr3 => mmax
  iptr4 => nmax

  if(ifplot.eq.1) then
     call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,  &
        iptr3,iptr4, norder, &
        'Triangulated surface of the stellarator')
  endif

  call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4, &
    npols,uvs,umatr,srcvals,srccoefs)
endif

if(igeomtype.eq.3) then
  done = 1
  pi = atan(done)*4
  umin = 0
  umax = 2*pi
  vmin = 0
  vmax = 2*pi
  allocate(triaskel(3,3,npatches))
  nover = 0
  call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2), &
    nover,npatches,npatches,triaskel)
  call prinf('npatches=*',npatches,1)

  p1(1) = 1
  p1(2) = 2
  p1(3) = 0.25d0

  p2(1) = 1.2d0
  p2(2) = 1.0d0
  p2(3) = 1.7d0

!
!         numberof oscillations
!
  p4(1) = 5.0d0


  ptr1 => triaskel(1,1,1)
  ptr2 => p1(1)
  ptr3 => p2(1)
  ptr4 => p4(1)
  xtri_geometry => xtri_wtorus_eval
  if(ifplot.eq.1) then
     call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, &
        ptr3,ptr4, norder,&
        'Triangulated surface of the wtorus')
  endif


  call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4, &
    npols,uvs,umatr,srcvals,srccoefs)
endif

if(igeomtype.eq.4) then
  done = 1
  pi = atan(done)*4
  umin = 0
  umax = 2*pi
  vmin = 0
  vmax = 2*pi
  allocate(triaskel(3,3,npatches))
  nover = 0
  call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2), &
    nover,npatches,npatches,triaskel)
  call prinf('npatches=*',npatches,1)

  p1(1) = 1.0d0
  p1(2) = 3.0d0
  p1(3) = 0.25d0

  p2(1) = 1.0d0
  p2(2) = 1.0d0
  p2(3) = 1.0d0

!
!         number of oscillations
!
  p4(1) = 0.0d0


  ptr1 => triaskel(1,1,1)
  ptr2 => p1(1)
  ptr3 => p2(1)
  ptr4 => p4(1)
  xtri_geometry => xtri_wtorus_eval
  if(ifplot.eq.1) then
     call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, &
        ptr3,ptr4, norder, &
        'Triangulated surface of the torus')
  endif


  call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4, &
    npols,uvs,umatr,srcvals,srccoefs)
endif

return  
end





subroutine test_exterior_pt(npatches,norder,npts,srcvals, &
  srccoefs,wts,xyzout,isout)
!
!
!  this subroutine tests whether the pt xyzin, is
!  in the exterior of a surface, and also estimates the error
!  in representing e^{ir/2}/r and \grad e^{ir/2}/r \cdot n
!  centered at the interior point. Whether a point 
!  is in the interior or not is tested using Gauss' 
!  identity for the flux due to a point charge
!
!
!  input:
!    npatches - integer
!       number of patches
!    norder - integer
!       order of discretization
!    npts - integer
!       total number of discretization points on the surface
!    srccoefs - real *8 (9,npts)
!       koornwinder expansion coefficients of geometry info
!    xyzout -  real *8 (3)
!       point to be tested
!
!  output: 
!    isout - boolean
!      whether the target is in the interior or not
!

implicit none
integer npatches,norder,npts,npols
real *8 srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
real *8 tmp(3)
real *8 dpars,done,pi
real *8, allocatable :: rsurf(:),err_p(:,:) 
integer ipars,norderhead,nd
complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
complex *16 zk,val

integer ipatch,j,i
real *8 ra,ds
logical isout

done = 1
pi = atan(done)*4

npols = (norder+1)*(norder+2)/2


zk = 0

ra = 0



do ipatch=1,npatches
  do j=1,npols
    i = (ipatch-1)*npols + j
    call h3d_sprime(xyzout,12,srcvals(1,i),0,dpars,1,zk,0,ipars, &
      val)

    call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
    ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
    ra = ra + real(val)*wts(i)
  enddo
enddo

if(abs(ra+4*pi).le.1.0d-1) isout = .false.
if(abs(ra).le.1.0d-1) isout = .true.

return
end








subroutine l3getsph(nmax,mm,nn,ndx,xyzs,ynms,npts,ynm)
implicit real *8 (a-h,o-z)
real *8 :: xyzs(ndx,npts)
complex *16 ynms(npts),ima
real *8 rat1(10000),rat2(10000)
real *8 ynm(0:nmax,0:nmax)
data ima/(0.0d0,1.0d0)/

call ylgndrini(nmax, rat1, rat2)

do i=1,npts
  x=xyzs(1,i)
  y=xyzs(2,i)
  z=xyzs(3,i)
  r=sqrt(x**2+y**2+z**2)
  call cart2polar(xyzs(1,i),r,theta,phi)
  ctheta = cos(theta)
  call ylgndrf(nmax, ctheta, ynm, rat1, rat2)
  ynms(i) = ynm(nn,abs(mm))*exp(ima*mm*phi)        
enddo

return
end
