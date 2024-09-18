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

subroutine mtxbsigma_eval_addsub(npatches, norders, ixyzs, &
     iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, &
     dpars, ndz, zpars, ndi, ipars, nnz_gc, nnz_helm, row_ptr_gc, &
     row_ptr_helm, col_ind_gc, col_ind_helm, iquad_gc, iquad_helm, nquad_gc, &
     nquad_helm, nker, wnear_gc, wnear_helm, novers_gc, novers_helm, &
     nptso_gc, nptso_helm, ixyzso_gc, ixyzso_helm, srcover_gc, srcover_helm, &
     whtsover_gc, whtsover_helm, lwork, work, idensflag, ndim, sigma, &
     ipotflag, ndim_p, rjvec, rho, potsigma, curlj, gradrho)

  implicit none
  integer npatches,npts,ndd,ndz,ndi,nnz_gc,nnz_helm,ndtarg
  integer norders(npatches),ixyzs(npatches+1),iptype(npatches),ipars(ndi)
  real *8 srccoefs(9,npts),srcvals(12,npts),eps(2),dpars(ndd)
  complex *16 zpars(ndz)
  integer ntarg,nquad_gc,nquad_helm,nker,nptso_gc,nptso_helm,lwork,ndim
  integer row_ptr_gc(ntarg+1),col_ind_gc(nnz_gc)
  integer iquad_gc(nnz_gc+1),novers_gc(npatches+1)
  integer row_ptr_helm(ntarg+1),col_ind_helm(nnz_helm)
  integer iquad_helm(nnz_helm+1),novers_helm(npatches+1)
  complex *16 wnear_gc(nquad_gc,3),wnear_helm(nquad_helm)
  complex *16 sigma(ndim,npts),rjvec(3,npts),rho(npts)
  integer ixyzso_gc(npatches+1),ixyzso_helm(npatches+1)
  real *8 targs(ndtarg,ntarg)
  real *8 srcover_gc(12,nptso_gc),srcover_helm(12,nptso_helm)
  real *8 whtsover_gc(nptso_gc),whtsover_helm(nptso_helm)
  real *8 work(lwork)
  integer idensflag,ipotflag,ndim_p
  complex *16 potsigma(ndim,ntarg),curlj(3,ntarg),gradrho(3,ntarg)

  call helm_comb_dir_eval_addsub_vec(npatches, norders, ixyzs, &
       iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps(2), ndd, &
       dpars, ndz, zpars, ndi, ipars, nnz_helm, row_ptr_helm, col_ind_helm, &
       iquad_helm, nquad_helm, nker, wnear_helm, novers_helm, nptso_helm, &
       ixyzso_helm, srcover_helm, whtsover_helm, lwork, work, idensflag, &
       ndim, sigma, ipotflag, ndim_p, potsigma)
  
  call lpcomp_gradcurlhelm_addsub(npatches, norders, ixyzs, &
       iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps(1), &
       zpars(1), nnz_gc, row_ptr_gc, col_ind_gc, iquad_gc, nquad_gc, &
       wnear_gc, rjvec, rho, novers_gc, nptso_gc, ixyzso_gc, srcover_gc, &
       whtsover_gc, curlj, gradrho)

  return
end subroutine mtxBsigma_eval_addsub
