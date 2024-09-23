subroutine helm_comb_dir_eval_addsub_vec(npatches, norders, ixyzs, &
     iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, &
     dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
     nquad, nker, wnear, novers, nptso, ixyzso, srcover, whtsover, &
     lwork, work, idensflag, ndim, sigma, ipotflag, ndim_p, pot)

  ! Vectorized version of lpcomp_helm_comb_dir_addsub
  !
  ! Changed input arguments:
  !   ndim: integer
  !     number of densities for which to evaluate layer potential
  !   sigma: complex *16 (ndim,npts)
  !
  ! Changed output argument:
  !   pot: complex *16 (ndim,npts)

  implicit none
  integer npatches,npts,ndd,ndz,ndi,nnz,ndtarg
  integer norders(npatches),ixyzs(npatches+1),iptype(npatches),ipars(ndi)
  real *8 srccoefs(9,npts),srcvals(12,npts),eps,dpars(ndd)
  complex *16 zpars(ndz)
  integer ntarg,nquad,nker,nptso,lwork,ndim
  integer row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1),novers(npatches+1)
  complex *16 wnear(nquad),sigma(ndim,npts)
  integer ixyzso(npatches+1)
  real *8 targs(ndtarg,ntarg),srcover(12,nptso),whtsover(nptso),work(lwork)
  integer idensflag, ipotflag, ndim_p
  complex *16 pot(ndim,ntarg)

  integer i

  do i = 1,ndim
     call helm_comb_dir_eval_addsub(npatches, norders, ixyzs, &
     iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, &
     dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
     nquad, nker, wnear, novers, nptso, ixyzso, srcover, whtsover, &
     lwork, work, idensflag, 1, sigma(i,:), ipotflag, ndim_p, pot(i,:))
  enddo

  return
end subroutine helm_comb_dir_eval_addsub_vec

subroutine mtxbsigma_eval_addsub(npatches,norders,ixyzs, &
     iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs, &
     eps,ndd,dpars,ndz,zpars,ndi,ipars,nnz_gc,nnz_helm, &
     row_ptr_gc,row_ptr_helm,col_ind_gc,col_ind_helm, &
     iquad_gc,iquad_helm,nquad_gc,nquad_helm,nker,&
     wnear_gc,wnear_helm,novers_gc,novers_helm,nptso_gc,&
     nptso_helm,ixyzso_gc,ixyzso_helm,srcover_gc, &
     srcover_helm,whtsover_gc,whtsover_helm,lwork, &
     work,idensflag,ndim,sigma,ipotflag,ndim_p,rjvec,rho, &
     pot,curlj,gradrho)

  implicit none

  integer npatches,npts,ndd,ndz,ndi,nnz_gc,ndtarg
  integer norders(npatches),ixyzs(npatches+1)
  integer iptype(npatches),ipars(ndi)
  real *8 srccoefs(9,npts),srcvals(12,npts),eps(2),dpars(ndd)
  complex *16 zpars(ndz)
  integer ntarg,nquad_gc,nker,nptso_gc,lwork,ndim
  integer row_ptr_gc(ntarg+1),col_ind_gc(nnz_gc)
  integer iquad_gc(nnz_gc+1),novers_gc(npatches+1)
  complex *16 wnear_gc(nquad_gc),sigma(ndim,npts)
  integer ixyzso_gc(npatches+1)
  real *8 targs(ndtarg,ntarg),srcover_gc(12,nptso_gc)
  real *8 whtsover_gc(nptso_gc),work(lwork)
  integer idensflag, ipotflag, ndim_p
  complex *16 pot(ndim,ntarg)

  integer ixyzso_helm(npatches+1)
  integer nnz_helm,row_ptr_helm(ntarg+1)
  integer col_ind_helm(nnz_helm),nquad_helm
  integer iquad_helm(nnz_helm+1)
  complex *16 rjvec(3,npts),rho(npts)
  complex *16 wnear_helm(nquad_helm,3)

  integer novers_helm(npatches),nptso_helm
  real *8 srcover_helm(12,nptso_helm),whtsover_helm(nptso_helm)
  complex *16 curlj(3,ntarg),gradrho(3,ntarg)
  
!   call lpcomp_gradcurlhelm_addsub(npatches,norders,ixyzs,&
!      iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs, &
!      eps(1),zpars(1),nnz_gc,row_ptr_gc,col_ind_gc,iquad_gc,&
!      nquad_gc,wnear_gc,rjvec,rho,novers_gc,nptso_gc,ixyzso_gc,&
!      srcover_gc,whtsover_gc,curlj,gradrho)
! 
!   call helm_comb_dir_eval_addsub_vec(npatches, norders, ixyzs, &
!      iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps(2), ndd, &
!      dpars, ndz, zpars, ndi, ipars, nnz_helm, row_ptr_helm, col_ind_helm, iquad_helm, &
!      nquad_helm, nker, wnear_helm, novers_helm, nptso_helm, ixyzso_helm, srcover_helm, whtsover_helm, &
!      lwork, work, idensflag, ndim, sigma, ipotflag, ndim_p, pot)
! 
  return
end subroutine mtxbsigma_eval_addsub

subroutine testfun(npts, ndim, ntarg, sigma, pot)
  implicit none
  integer npts, ndim, ntarg
  complex *16 sigma(ndim,npts), pot(3,ntarg)

  pot = 0

  return 
end subroutine testfun

! subroutine mtxbsigma_eval_addsub(npatches, norders, ixyzs, &
!      iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, &
!      dpars, ndz, zpars, ndi, ipars, nnz_gc, nnz_helm, row_ptr_gc, &
!      row_ptr_helm, col_ind_gc, col_ind_helm, iquad_gc, iquad_helm, nquad_gc, &
!      nquad_helm, nker, wnear_gc, wnear_helm, novers_gc, novers_helm, &
!      nptso_gc, nptso_helm, ixyzso_gc, ixyzso_helm, srcover_gc, srcover_helm, &
!      whtsover_gc, whtsover_helm, lwork, work, idensflag, ndim, sigma, &
!      ipotflag, ndim_p, rjvec, rho, potsigma, curlj, gradrho)
! 
!   implicit none
!   integer npatches,npts,ndd,ndz,ndi,nnz_gc,nnz_helm,ndtarg
!   integer norders(npatches),ixyzs(npatches+1),iptype(npatches),ipars(ndi)
!   real *8 srccoefs(9,npts),srcvals(12,npts),eps(2),dpars(ndd)
!   complex *16 zpars(ndz)
!   integer ntarg,nquad_gc,nquad_helm,nker,nptso_gc,nptso_helm,lwork,ndim
!   integer row_ptr_gc(ntarg+1),col_ind_gc(nnz_gc)
!   integer iquad_gc(nnz_gc+1),novers_gc(npatches+1)
!   integer row_ptr_helm(ntarg+1),col_ind_helm(nnz_helm)
!   integer iquad_helm(nnz_helm+1),novers_helm(npatches+1)
!   complex *16 wnear_gc(nquad_gc,3),wnear_helm(nquad_helm)
!   complex *16 sigma(ndim,npts),rjvec(3,npts),rho(npts)
!   integer ixyzso_gc(npatches+1),ixyzso_helm(npatches+1)
!   real *8 targs(ndtarg,ntarg)
!   real *8 srcover_gc(12,nptso_gc),srcover_helm(12,nptso_helm)
!   real *8 whtsover_gc(nptso_gc),whtsover_helm(nptso_helm)
!   real *8 work(lwork)
!   integer idensflag,ipotflag,ndim_p
!   complex *16 potsigma(ndim,ntarg),curlj(3,ntarg),gradrho(3,ntarg)
! 
!   call helm_comb_dir_eval_addsub_vec(npatches, norders, ixyzs, &
!        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps(2), ndd, &
!        dpars, ndz, zpars, ndi, ipars, nnz_helm, row_ptr_helm, col_ind_helm, &
!        iquad_helm, nquad_helm, nker, wnear_helm, novers_helm, nptso_helm, &
!        ixyzso_helm, srcover_helm, whtsover_helm, lwork, work, idensflag, &
!        ndim, sigma, ipotflag, ndim_p, potsigma)
!   
!   call lpcomp_gradcurlhelm_addsub(npatches, norders, ixyzs, &
!        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps(1), &
!        zpars(1), nnz_gc, row_ptr_gc, col_ind_gc, iquad_gc, nquad_gc, &
!        wnear_gc, rjvec, rho, novers_gc, nptso_gc, ixyzso_gc, srcover_gc, &
!        whtsover_gc, curlj, gradrho)
! 
!   return
! end subroutine mtxbsigma_eval_addsub
! 
! subroutine mtxbsigma_eval_addsub_stupid(npatches,norders,ixyzs,&
!      iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs, &
!      eps,zk,nnz,row_ptr,col_ind,iquad,nquad,wnear,rjvec,rho, &
!      novers,nptso,ixyzso,srcover,whtsover,curlj,gradrho,&
!      npatches2, norders2, ixyzs2, iptype2, npts2, srccoefs2, &
!      srcvals2, ndtarg2, ntarg2, targs2, eps2, ndd2, &
!      dpars2, ndz2, zpars2, ndi2, ipars2, nnz2, row_ptr2, col_ind2, iquad2, &
!      nquad2, nker2, wnear2, novers2, nptso2, ixyzso2, srcover2, whtsover2, &
!      lwork2, work2, idensflag2, ndim2, sigma, ipotflag2, ndim_p2, pot)
! 
!   implicit none
!   integer npatches,norder,npols,npts
!   integer ndtarg,ntarg
!   integer norders(npatches),ixyzs(npatches+1)
!   integer ixyzso(npatches+1),iptype(npatches)
!   real *8 srccoefs(9,npts),srcvals(12,npts),eps
!   complex *16 zk
!   real *8 targs(ndtarg,ntarg) 
!   integer nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
!   integer iquad(nnz+1)
!   complex *16 rjvec(3,npts),rho(npts)
!   complex *16 wnear(nquad,3)
! 
!   integer novers(npatches)
!   integer nover,npolso,nptso
!   real *8 srcover(12,nptso),whtsover(nptso)
!   complex *16 curlj(3,ntarg),gradrho(3,ntarg)
! 
!   integer npatches2, npts2, ndd2, ndz2, ndi2, nnz2, ndtarg2
!   integer norders2(npatches2), ixyzs2(npatches2+1), iptype2(npatches2), ipars2(ndi2)
!   real *8 srccoefs2(9, npts2), srcvals2(12, npts2), eps2, dpars2(ndd2)
!   complex *16 zpars2(ndz2)
!   integer ntarg2, nquad2, nker2, nptso2, lwork2, ndim2
!   integer row_ptr2(ntarg2+1), col_ind2(nnz2), iquad2(nnz2+1), novers2(npatches2+1)
!   complex *16 wnear2(nquad2), sigma(ndim2, npts2)
!   integer ixyzso2(npatches2+1)
!   real *8 targs2(ndtarg2, ntarg2), srcover2(12, nptso2), whtsover2(nptso2), work2(lwork2)
!   integer idensflag2, ipotflag2, ndim_p2
!   complex *16 pot(ndim2, ntarg2)
! 
!   call lpcomp_gradcurlhelm_addsub(npatches,norders,ixyzs,&
!      iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs, &
!      eps,zk,nnz,row_ptr,col_ind,iquad,nquad,wnear,rjvec,rho, &
!      novers,nptso,ixyzso,srcover,whtsover,curlj,gradrho)
! 
!   call helm_comb_dir_eval_addsub_vec(npatches2, norders2, ixyzs2, &
!      iptype2, npts2, srccoefs2, srcvals2, ndtarg2, ntarg2, targs2, eps2, ndd2, &
!      dpars2, ndz2, zpars2, ndi2, ipars2, nnz2, row_ptr2, col_ind2, iquad2, &
!      nquad2, nker2, wnear2, novers2, nptso2, ixyzso2, srcover2, whtsover2, &
!      lwork2, work2, idensflag2, ndim2, sigma, ipotflag2, ndim_p2, pot)
! 
! end subroutine mtxbsigma_eval_addsub_stupid
! 






