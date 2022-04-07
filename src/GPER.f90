
        SUBROUTINE sqrdrv1(tau,nobs,nvars,x,r,vl)
        IMPLICIT NONE
        INTEGER:: nobs
        INTEGER:: nvars
        INTEGER:: i
        DOUBLE PRECISION:: tau
        DOUBLE PRECISION:: dl(nobs)
        DOUBLE PRECISION:: r(nobs)
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION:: vl(nvars)
        vl = 0.0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0-tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
			ENDDO
        vl = matmul(dl, x) / nobs
        END SUBROUTINE sqrdrv1

! --------------------------------------------------! -------------------------
! -------------------------------------------------! --------------------------
! --------quantile regression  with group lasso penalite--------------------
! ---------------------------------------------f1 and f2-----------------------
! --------------------------------------------------! ------------------------- jd , isd, ibeta

SUBROUTINE sqr1lasso (tau,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                                eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::tau
        DOUBLE PRECISION::delta
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION:: dd1
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
        ! - - - some initial setup - - -
        delta=2*max(tau,1-tau)
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
        IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv1(tau,nobs,nvars,x,r,vl)

        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO
!----------------------------------------------------------------------------------------

        DO l=1,nlam

        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0

        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO

        al = al0 * alf
        ENDIF
        ENDIF

        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
        ! --------- outer loop ----------------------------
        DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    u = u + dl(i)*x(i,start:end)/nobs
			ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb

        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))

        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    d = d + dl(i)
			ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    u = u + dl(i)*x(i,start:end)/nobs
			ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
		dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    d = d + dl(i)
			ENDDO
        d = 0.50 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
        jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
  
        IF (jx /= 0) CYCLE
        CALL sqrdrv1(tau,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE sqr1lasso



SUBROUTINE sqr1mcp (gamm,tau,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                                eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::tau
        DOUBLE PRECISION::delta
        DOUBLE PRECISION::gamm
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION:: dd1
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::tt
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
        ! - - - some initial setup - - -
        delta=2*max(tau,1-tau)
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
        IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv1(tau,nobs,nvars,x,r,vl)

        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO
!----------------------------------------------------------------------------------------

        DO l=1,nlam

        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0

        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO

        al = al0 * alf
        ENDIF
        ENDIF

        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
        ! --------- outer loop ----------------------------
        DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    u = u + dl(i)*x(i,start:end)/nobs
			ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))

        t=unorm-pf(g)*al

        IF (t>0.0D0) THEN
				tt=unorm-gam(g)*gamm*al
				IF (tt>0.0D0) THEN
						b(start:end) = u/gam(g)
				ELSE
						b(start:end)=(u*t/((gam(g)-pf(g)/gamm)*unorm))
				END IF
        ELSE
				b(start:end)=0.0D0

	END IF

        dd=b(start:end)-oldb

        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))

        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    d = d + dl(i)
			ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    u = u + dl(i)*x(i,start:end)/nobs
			ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))

        t=unorm-pf(g)*al

        IF (t>0.0D0) THEN
				tt=unorm-gam(g)*gamm*al
				IF (tt>0.0D0) THEN
						b(start:end) = u/gam(g)
				ELSE
						b(start:end)=(u*t/((gam(g)-pf(g)/gamm)*unorm))
				END IF
        ELSE
				b(start:end)=0.0D0

	END IF

        dd=b(start:end)-oldb

        IF(any(dd/=0.0D0)) THEN
		dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    d = d + dl(i)
			ENDDO
        d = 0.50 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
        jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
  
        IF (jx /= 0) CYCLE
        CALL sqrdrv1(tau,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE sqr1mcp


SUBROUTINE sqr1scad (gamm,tau,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                                eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::tau
        DOUBLE PRECISION::delta
        DOUBLE PRECISION::gamm
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION:: dd1
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::tt
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
        ! - - - some initial setup - - -
        delta=2*max(tau,1-tau)
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
        IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv1(tau,nobs,nvars,x,r,vl)

        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO
!----------------------------------------------------------------------------------------

        DO l=1,nlam

        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0

        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO

        al = al0 * alf
        ENDIF
        ENDIF

        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
        ! --------- outer loop ----------------------------
        DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    u = u + dl(i)*x(i,start:end)/nobs
			ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))

	IF (unorm <= ((pf (g)+gam(g))*al)) THEN
        	t=unorm-pf(g)*al
        	IF(t>0.0D0) THEN
        		b(start:end)=u*t/(gam(g)*unorm)
        	ELSE
        		b(start:end)=0.0D0
        	ENDIF
        ELSE IF (unorm <= (gam(g)*al*gamm )) THEN
        	t=unorm-pf(g)*al*(gamm/(gamm-1))
        	IF(t>0.0D0) THEN
        		b(start:end)=u*t/((unorm)*(gam(g)-pf(g)/(gamm-1.0D0)))
        	ELSE
        		b(start:end)=0.0D0
        	ENDIF
        ELSE
        	b(start:end) = u/gam(g)
        END IF

        dd=b(start:end)-oldb

        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))

        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    d = d + dl(i)
			ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    u = u + dl(i)*x(i,start:end)/nobs
			ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))

	IF (unorm <= ((pf (g)+gam(g))*al)) THEN
        	t=unorm-pf(g)*al
        	IF(t>0.0D0) THEN
        		b(start:end)=u*t/(gam(g)*unorm)
        	ELSE
        		b(start:end)=0.0D0
        	ENDIF
        ELSE IF (unorm <= (gam(g)*al*gamm )) THEN
        	t=unorm-pf(g)*al*(gamm/(gamm-1))
        	IF(t>0.0D0) THEN
        		b(start:end)=u*t/((unorm)*(gam(g)-pf(g)/(gamm-1.0D0)))
        	ELSE
        		b(start:end)=0.0D0
        	ENDIF
        ELSE
        	b(start:end) = u/gam(g)
        END IF

        dd=b(start:end)-oldb

        IF(any(dd/=0.0D0)) THEN
		dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    d = d + dl(i)
			ENDDO
        d = 0.50 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
        jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
  
        IF (jx /= 0) CYCLE
        CALL sqrdrv1(tau,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE sqr1scad

!----------------------------------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE asqr1mcp (gamm,tau,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                                eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::tau
        DOUBLE PRECISION::delta
        DOUBLE PRECISION::gamm
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::pfHis(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION:: dd1
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::tt
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::bnorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
		pfHis = pf
        ! - - - some initial setup - - -
        delta=2*max(tau,1-tau)
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
        IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv1(tau,nobs,nvars,x,r,vl)

        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO
!----------------------------------------------------------------------------------------

        DO l=1,nlam

		IF(l>1) THEN
					DO j=1,bn
						g=idx(j)
						bnorm = sqrt(dot_product(b(ix(g):iy(g)),b(ix(g):iy(g))))
						IF  (pf(g)*al*gamm>bnorm)   THEN
								pf(g)=pf(g)-bnorm/(al*gamm)
						ELSE
								pf(g)=0.0D0
						END IF
					ENDDO
		END IF

        ! -------------------------------------


        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0

        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO

        al = al0 * alf
        ENDIF
        ENDIF

        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO

        ! --------- outer loop ----------------------------
        DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    u = u + dl(i)*x(i,start:end)/nobs
			ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF

        dd=b(start:end)-oldb

        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))

        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    d = d + dl(i)
			ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    u = u + dl(i)*x(i,start:end)/nobs
			ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF

        dd=b(start:end)-oldb

        IF(any(dd/=0.0D0)) THEN
		dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    d = d + dl(i)
			ENDDO
        d = 0.50 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
        jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
  
        IF (jx /= 0) CYCLE
        CALL sqrdrv1(tau,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE asqr1mcp

SUBROUTINE asqr1scad (gamm,tau,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                                eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::tau
        DOUBLE PRECISION::delta
        DOUBLE PRECISION::gamm
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::pfHis(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION:: dd1
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::tt
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::bnorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
		pfHis = pf
        ! - - - some initial setup - - -
        delta=2*max(tau,1-tau)
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
        IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv1(tau,nobs,nvars,x,r,vl)

        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO
!----------------------------------------------------------------------------------------

        DO l=1,nlam
        pf = pfHis
        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0

        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO

        al = al0 * alf
        ENDIF
        ENDIF

        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
        ! --------- ----------------------------
        ! --------- ----------------------------
        !pf = pfHis
        IF(l>1) THEN
				DO j=1,bn
					g=idx(j)
					bnorm = sqrt(dot_product(b(ix(g):iy(g)),b(ix(g):iy(g))))
					IF  (pf(g)*al>bnorm)   THEN
							pf(g)=pf(g)
					ELSE IF ((pf(g)*al*gamm>bnorm)) THEN
							pf(g)=(gamm/(gamm-1))*pf(g)-bnorm/(al*(gamm-1))
					ELSE
							pf(g)=0.0D0
					END IF
				ENDDO
		END IF
        ! -------------------------------------
        ! -------------------------------------
        ! --------- outer loop ----------------------------
        DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    u = u + dl(i)*x(i,start:end)/nobs
			ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF

        dd=b(start:end)-oldb

        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))

        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    d = d + dl(i)
			ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    u = u + dl(i)*x(i,start:end)/nobs
			ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF

        dd=b(start:end)-oldb

        IF(any(dd/=0.0D0)) THEN
		dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
			DO i = 1,nobs
				    IF (r(i) < 0) THEN
				        dl (i) = 2.0D0*(1.0D0 - tau)*r(i)
				    ELSE
				        dl (i) = 2.0D0*(tau)*r(i)
				    END IF
				    d = d + dl(i)
			ENDDO
        d = 0.50 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
        jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
  
        IF (jx /= 0) CYCLE
        CALL sqrdrv1(tau,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE asqr1scad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------
!------------------------------Cosales for group lasso-----------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------

SUBROUTINE serdrvBeta(w,tau,nobs,nvars,x,r,rf,vl)

                IMPLICIT NONE
                INTEGER:: nobs
                INTEGER:: nvars
                INTEGER:: i
                DOUBLE PRECISION:: tau
                DOUBLE PRECISION:: w
                DOUBLE PRECISION:: dl(nobs)
                DOUBLE PRECISION:: r(nobs)
                DOUBLE PRECISION:: rf(nobs)
                DOUBLE PRECISION:: x(nobs,nvars)
                DOUBLE PRECISION:: vl(nvars)

                vl = 0.0D0

                DO i = 1, nobs
                    IF (rf(i) < 0) THEN
                        dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)+w*r(i)
                    ELSE
                        dl (i) = 2.0D0*(tau)*rf(i)+w*r(i)
                    END IF
                ENDDO
                vl = matmul(dl, x) / nobs 

END SUBROUTINE serdrvBeta


SUBROUTINE serdrvTh(tau,nobs,nvars,x,rf,vl)

                IMPLICIT NONE
                INTEGER:: nobs
                INTEGER:: nvars
                INTEGER:: i
                DOUBLE PRECISION:: tau
                DOUBLE PRECISION:: dl(nobs)
                DOUBLE PRECISION:: rf(nobs)
                DOUBLE PRECISION:: x(nobs,nvars)
                DOUBLE PRECISION:: vl(nvars)

                vl = 0.0D0

                DO i = 1, nobs
                    IF (rf(i) < 0) THEN
                        dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)
                    ELSE
                        dl (i) = 2.0D0*(tau)*rf(i)
                    END IF
                ENDDO
                vl = matmul(dl, x) / nobs

END SUBROUTINE serdrvTh



SUBROUTINE CserLassoUnifier (w,tau,bn,bs,ix,iy,gam,nobs,nvars,x,y,pfmean,pfscale,dfmax,pmax,nlam,flmin,ulam,&
                                eps,maxit,intr,nalam,b0,beta,t0,theta,idx,idxf,nbeta,ntheta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idxf(pmax)
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        INTEGER::ntheta(nlam)
        DOUBLE PRECISION:: w
        DOUBLE PRECISION:: delta
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::tau
        DOUBLE PRECISION::lambda
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pfmean(bn)
        DOUBLE PRECISION::pfscale(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
	    DOUBLE PRECISION:: t0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::theta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION:: dd1
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al  !al represente mu
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: th
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldtheta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::rf
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldth
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidxf
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::nib
        INTEGER::nif
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        DOUBLE PRECISION:: ga1(bn)
        DOUBLE PRECISION:: ga2(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
        ALLOCATE(b(0:nvars), STAT = jerr)
	jerr=jerr + ierr
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(th(0:nvars), STAT = jerr)
	    jerr=jerr + ierr
        ALLOCATE(oldtheta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(rf(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidxf(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pfmean pfscale- - -
        IF(maxval(pfmean) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pfmean = max(0.0D0, pfmean)

        IF(maxval(pfscale) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pfscale = max(0.0D0, pfscale)

        ! - - - some initial setup - - -
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
	    rf = y
        b = 0.0D0
        delta=2*max(tau,1.0D0-tau)
        oldbeta = 0.0D0
        th = 0.0D0
        oldtheta = 0.0D0
        idx = 0
        idxf = 0
        oidx = 0
        oidxf = 0
        npass = 0
        nib = npass
	    nif = npass
        ! --------- lambda loop ----------------------------
        IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF

!----------------------------------------------------------------------------------------

        DO l=1,nlam
            al0 = al
            IF(flmin>=1.0D0) THEN
                al=ulam(l)
                ELSE
                IF(l > 2) THEN
                al=al*alf
                ELSE IF(l==1) THEN
                       al=big
                ELSE IF(l==2) THEN
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        vl=0.0D0
			CALL serdrvBeta(w,tau,nobs,nvars,x,r,rf,vl)
			DO g = 1,bn
				ALLOCATE(u(bs(g)),STAT=ierr)
				IF(ierr/=0) RETURN
				u = vl(ix(g):iy(g))
				ga2(g) = sqrt(dot_product(u,u))
				DEALLOCATE(u)
			END DO
                        vl=0.0D0
			CALL serdrvTh(tau,nobs,nvars,x,rf,vl)
			DO g = 1,bn
				ALLOCATE(u(bs(g)),STAT=ierr)
				IF(ierr/=0) RETURN
				u = vl(ix(g):iy(g))
				ga1(g) = sqrt(dot_product(u,u))
				DEALLOCATE(u)
			END DO
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    al0 = 0.0D0
                    DO g = 1,bn
                                IF(pfscale(g)>0.0D0) THEN
                                al0 = max(al0, ga1(g) / pfscale(g))
                                ENDIF
                    END DO
                    DO g = 1,bn
                                IF(pfmean(g)>0.0D0) THEN
                                al0 = max(al0, ga2(g) / pfmean(g))
                                ENDIF
                    END DO
                    al = al0 * alf
                ENDIF
            ENDIF
        ! --------- outer loop ----------------------------
        DO  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                oldbeta(0)=b(0)
                IF(nib>0) THEN
		        DO j=1,nib
		                g=idx(j)
		                oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
		        ENDDO
                ENDIF

                oldtheta(0)=th(0)
                IF(nif>0) THEN
		        DO j=1,nif
		                g=idxf(j)
		                oldtheta(ix(g):iy(g))=th(ix(g):iy(g))
		        ENDDO
                ENDIF
        
		DO  ! --middle loop-------------------------------------

		npass=npass+1
		dif=0.0D0
			DO g=1,bn  ! theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta

		                start=ix(g)
		                end=iy(g)
		                ALLOCATE(u(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(dd(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(oldth(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                IF(jerr/=0) RETURN
		                oldth=th(start:end)
		                u = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)
					    END IF
					    u = u + dl(i)*x(i,start:end)/nobs
				ENDDO

				u=delta*gam(g)*th(start:end) + u
				unorm=sqrt(dot_product(u,u))
				t=unorm-pfscale(g)*al
				IF(t>0.0D0)     THEN
							th(start:end)=u*t/(delta*gam(g)*unorm)
						ELSE
							th(start:end)=0.0D0
				ENDIF
				dd=th(start:end)-oldth
		                IF(any(dd/=0.0D0)) THEN

		                    dif=max(dif,delta*gam(g)**2*dot_product(dd,dd))
		                    rf=rf-matmul(x(:,start:end),dd)
		                    IF(oidxf(g)==0) THEN
		                        nif=nif+1
		                        IF(nif>pmax) EXIT
		                                oidxf(g)=nif
		                                idxf(nif)=g
		                        ENDIF
		                  ENDIF
				DEALLOCATE(u,dd,oldth)
			ENDDO  ! theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta

			DO g=1,bn    ! beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta
		                start=ix(g)
		                end=iy(g)
		                ALLOCATE(u(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(dd(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(oldb(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                IF(jerr/=0) RETURN
		                oldb=b(start:end)
		                u = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)+ w*r(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)+ w*r(i)
					    END IF
					    u = u + dl(i)*x(i,start:end)/nobs
				ENDDO
				u=(w+delta)*gam(g)*b(start:end) + u
				unorm=sqrt(dot_product(u,u))
				t=unorm-pfmean(g)*al
		                IF(t>0.0D0) THEN
		                        b(start:end)=u*t/((w+delta)*gam(g)*unorm)
		                        ELSE
		                        b(start:end)=0.0D0
		                ENDIF

		                dd=b(start:end)-oldb

		                IF(any(dd/=0.0D0)) THEN
		                    dif=max(dif,delta*gam(g)**2*dot_product(dd,dd))

		                    rf=rf-matmul(x(:,start:end),dd)
				    r=r-matmul(x(:,start:end),dd)
		                    IF(oidx(g)==0) THEN
		                        nib=nib+1
		                        IF(nib>pmax) EXIT
		                                oidx(g)=nib
		                                idx(nib)=g
		                    ENDIF
		                 ENDIF
		                DEALLOCATE(u,dd,oldb)
			ENDDO  ! ! ! beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta


			IF(intr /= 0) THEN !  intercept theta theta theta intercept theta theta theta intercept theta theta theta intercept theta theta th
				d = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)
					    END IF
					    d = d + dl(i) / nobs
				ENDDO
				d = d / delta
				IF(d/=0.0D0) THEN
					th(0)=th(0)+d
					rf=rf-d
					dif=max(dif, d**2)
				ENDIF
			ENDIF   !  !  intercept theta theta theta intercept theta theta theta intercept theta theta theta intercept theta theta th

			IF(intr /= 0) THEN !  intercept beta beta beta intercept beta beta beta intercept beta beta betaintercept
				d = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)+ w*r(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)+ w*r(i)
					    END IF
					    d = d + dl(i) / nobs
				ENDDO
				d = d /(w+delta)
				IF(d/=0.0D0) THEN
					b(0)=b(0)+d
					r=r-d
					rf=rf-d
					dif=max(dif, d**2)
				ENDIF
			ENDIF   !  intercept beta beta beta intercept beta beta beta intercept beta beta beta intercept beta beta beta                   
			IF (nib > pmax) EXIT
			IF (nif > pmax) EXIT
			IF (dif < eps) EXIT

			IF(npass > maxit) THEN
				jerr=-l
				RETURN
			ENDIF
		

		ENDDO   ! --middle loop-------------------------------------
		IF(nib>pmax) EXIT
		IF(nif>pmax) EXIT

		!--- final check ------------------------
		jx = 0
		max_gam = maxval(gam)
		IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
		IF(any((max_gam*(th-oldtheta)/(1+abs(th)))**2 >= eps)) jx = 1
		IF (jx /= 0) cycle
                EXIT
        ENDDO
        !---------- final update variable and save results------------
            IF(nif>pmax) THEN
                jerr=-10000-l
                EXIT
            ENDIF


            IF(nib>0) THEN  !0000000000000000000000
		    DO j=1,nib
		        g=idx(j)
		        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
		    ENDDO
            ENDIF
            nbeta(l)=nib
            b0(l)=b(0)
            alam(l)=al
            nalam=l      !0000000000000000000000

            IF(nif>0) THEN  !111111111111111111111
		    DO j=1,nif
		        g=idxf(j)
		        theta(ix(g):iy(g),l)=th(ix(g):iy(g))
		    ENDDO
            ENDIF
            ntheta(l)=nif
	    t0(l)=th(0)     !11111111111111111111111

            IF (l < mnl) CYCLE

            me=0
            DO j=1,nib
                g=idx(j)
                IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
            ENDDO
            IF(me>dfmax) EXIT

            me=0
            DO j=1,nif
                g=idxf(j)
                IF(any(theta(ix(g):iy(g),l)/=0.0D0)) me=me+1
            ENDDO
            IF(me>dfmax) EXIT
         

        ENDDO  ! fin de la grnade boucle for l in 1:nlamth 

        DEALLOCATE(b,oldbeta,th,oldtheta,r,rf,oidx,oidxf)
        RETURN
END SUBROUTINE CserLassoUnifier
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CserMcpUnifier (gamm,w,tau,bn,bs,ix,iy,gam,nobs,nvars,x,y,pfmean,pfscale,dfmax,pmax,nlam,flmin,ulam,&
                                eps,maxit,intr,nalam,b0,beta,t0,theta,idx,idxf,nbeta,ntheta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idxf(pmax)
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        INTEGER::ntheta(nlam)
        DOUBLE PRECISION:: w
        DOUBLE PRECISION:: gamm
        DOUBLE PRECISION:: delta
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::tau
        DOUBLE PRECISION::lambda
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pfmean(bn)
        DOUBLE PRECISION::pfscale(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
	DOUBLE PRECISION:: t0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::theta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION:: dd1
        DOUBLE PRECISION::d
        DOUBLE PRECISION::tt
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al  !al represente mu
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: th
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldtheta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::rf
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldth
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidxf
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::nib
        INTEGER::nif
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        DOUBLE PRECISION:: ga1(bn)
        DOUBLE PRECISION:: ga2(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
        ALLOCATE(b(0:nvars), STAT = jerr)
	jerr=jerr + ierr
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(th(0:nvars), STAT = jerr)
	jerr=jerr + ierr
        ALLOCATE(oldtheta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(rf(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidxf(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pfmean pfscale- - -
        IF(maxval(pfmean) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pfmean = max(0.0D0, pfmean)

        IF(maxval(pfscale) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pfscale = max(0.0D0, pfscale)

        ! - - - some initial setup - - -
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
	rf = y
        b = 0.0D0
        delta=2*max(tau,1.0D0-tau)
        oldbeta = 0.0D0
        th = 0.0D0
        oldtheta = 0.0D0
        idx = 0
        idxf = 0
        oidx = 0
        oidxf = 0
        npass = 0
        nib = npass
	nif = npass
        ! --------- lambda loop ----------------------------
        IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF

!----------------------------------------------------------------------------------------

        DO l=1,nlam
            al0 = al
            IF(flmin>=1.0D0) THEN
                al=ulam(l)
                ELSE
                IF(l > 2) THEN
                al=al*alf
                ELSE IF(l==1) THEN
                       al=big

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			CALL serdrvBeta(w,tau,nobs,nvars,x,r,rf,vl)
			DO g = 1,bn
				ALLOCATE(u(bs(g)),STAT=ierr)
				IF(ierr/=0) RETURN
				u = vl(ix(g):iy(g))
				ga2(g) = sqrt(dot_product(u,u))
				DEALLOCATE(u)
			END DO

			CALL serdrvTh(tau,nobs,nvars,x,rf,vl)
			DO g = 1,bn
				ALLOCATE(u(bs(g)),STAT=ierr)
				IF(ierr/=0) RETURN
				u = vl(ix(g):iy(g))
				ga1(g) = sqrt(dot_product(u,u))
				DEALLOCATE(u)
			END DO
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                ELSE IF(l==2) THEN
                    al0 = 0.0D0
                    DO g = 1,bn
                                IF(pfscale(g)>0.0D0) THEN
                                al0 = max(al0, ga1(g) / pfscale(g))
                                ENDIF
                    END DO
                    DO g = 1,bn
                                IF(pfmean(g)>0.0D0) THEN
                                al0 = max(al0, ga2(g) / pfmean(g))
                                ENDIF
                    END DO
                    al = al0 * alf
                ENDIF
            ENDIF
        ! --------- outer loop ----------------------------
        DO  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                oldbeta(0)=b(0)
                IF(nib>0) THEN
		        DO j=1,nib
		                g=idx(j)
		                oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
		        ENDDO
                ENDIF

                oldtheta(0)=th(0)
                IF(nif>0) THEN
		        DO j=1,nif
		                g=idxf(j)
		                oldtheta(ix(g):iy(g))=th(ix(g):iy(g))
		        ENDDO
                ENDIF
        
		DO  ! --middle loop-------------------------------------

		npass=npass+1
		dif=0.0D0
			DO g=1,bn  ! theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta

		                start=ix(g)
		                end=iy(g)
		                ALLOCATE(u(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(dd(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(oldth(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                IF(jerr/=0) RETURN
		                oldth=th(start:end)
		                u = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)
					    END IF
					    u = u + dl(i)*x(i,start:end)/nobs
				ENDDO

				u=delta*gam(g)*th(start:end) + u
				unorm=sqrt(dot_product(u,u))
				t=unorm-delta*gam(g)*gamm*al

				IF(t<0.0D0)     THEN
             				        tt=unorm-pfscale(g)*al
                                                IF (tt>0.0D0) THEN
                                                        th(start:end)= u*tt/((delta*gam(g)-pfscale(g)/gamm)*unorm)
                                                ELSE
							th(start:end)=0.0D0
						END IF
                                ELSE
						th(start:end)= u/(delta*gam(g))
				ENDIF

				dd=th(start:end)-oldth
		                IF(any(dd/=0.0D0)) THEN
		                    dif=max(dif,delta*gam(g)**2*dot_product(dd,dd))
		                    rf=rf-matmul(x(:,start:end),dd)
		                    IF(oidxf(g)==0) THEN
		                        nif=nif+1
		                        IF(nif>pmax) EXIT
		                                oidxf(g)=nif
		                                idxf(nif)=g
		                        ENDIF
		                  ENDIF
				DEALLOCATE(u,dd,oldth)
			ENDDO  ! theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta

			DO g=1,bn    ! beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta
		                start=ix(g)
		                end=iy(g)
		                ALLOCATE(u(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(dd(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(oldb(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                IF(jerr/=0) RETURN
		                oldb=b(start:end)
		                u = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)+ w*r(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)+ w*r(i)
					    END IF
					    u = u + dl(i)*x(i,start:end)/nobs
				ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!
				u=(w+delta)*gam(g)*b(start:end) + u
				unorm=sqrt(dot_product(u,u))
				t=unorm-(w+delta)*gam(g)*gamm*al

				IF(t<0.0D0)     THEN
             				        tt=unorm-al*pfmean(g)
                                                IF (tt>0.0D0) THEN
                                                        b(start:end)= u*tt/(((w+delta)*gam(g)-pfmean(g)/gamm)*unorm)
                                                ELSE

							b(start:end)=0.0D0
						END IF
				ELSE
						        b(start:end)=u/((w+delta)*gam(g))
				ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!



		                dd=b(start:end)-oldb

		                IF(any(dd/=0.0D0)) THEN
		                    dif=max(dif,delta*gam(g)**2*dot_product(dd,dd))

		                    rf=rf-matmul(x(:,start:end),dd)
				    r=r-matmul(x(:,start:end),dd)
		                    IF(oidx(g)==0) THEN
		                        nib=nib+1
		                        IF(nib>pmax) EXIT
		                                oidx(g)=nib
		                                idx(nib)=g
		                    ENDIF
		                 ENDIF
		                DEALLOCATE(u,dd,oldb)
			ENDDO  ! ! ! beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta


			IF(intr /= 0) THEN !  intercept theta theta theta intercept theta theta theta intercept theta theta theta intercept theta theta th
				d = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)
					    END IF
					    d = d + dl(i) / nobs
				ENDDO
				d = d / delta
				IF(d/=0.0D0) THEN
					th(0)=th(0)+d
					rf=rf-d
					dif=max(dif, d**2)
				ENDIF
			ENDIF   !  !  intercept theta theta theta intercept theta theta theta intercept theta theta theta intercept theta theta th

			IF(intr /= 0) THEN !  intercept beta beta beta intercept beta beta beta intercept beta beta betaintercept
				d = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)+ w*r(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)+ w*r(i)
					    END IF
					    d = d + dl(i) / nobs
				ENDDO
				d = d /(w+delta)
				IF(d/=0.0D0) THEN
					b(0)=b(0)+d
					r=r-d
					rf=rf-d
					dif=max(dif, d**2)
				ENDIF
			ENDIF   !  intercept beta beta beta intercept beta beta beta intercept beta beta beta intercept beta beta beta                  
			IF (nib > pmax) EXIT
			IF (nif > pmax) EXIT
			IF (dif < eps) EXIT

			IF(npass > maxit) THEN
				jerr=-l
				RETURN
			ENDIF
		

		ENDDO   ! --middle loop-------------------------------------
		IF(nib>pmax) EXIT
		IF(nif>pmax) EXIT

		!--- final check ------------------------
		jx = 0
		max_gam = maxval(gam)
		IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
		IF(any((max_gam*(th-oldtheta)/(1+abs(th)))**2 >= eps)) jx = 1
		IF (jx /= 0) cycle
                EXIT
        ENDDO
        !---------- final update variable and save results------------
            IF(nif>pmax) THEN
                jerr=-10000-l
                EXIT
            ENDIF


            IF(nib>0) THEN  !0000000000000000000000
		    DO j=1,nib
		        g=idx(j)
		        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
		    ENDDO
            ENDIF
            nbeta(l)=nib
            b0(l)=b(0)
            alam(l)=al
            nalam=l      !0000000000000000000000

            IF(nif>0) THEN  !111111111111111111111
		    DO j=1,nif
		        g=idxf(j)
		        theta(ix(g):iy(g),l)=th(ix(g):iy(g))
		    ENDDO
            ENDIF
            ntheta(l)=nif
	    t0(l)=th(0)     !11111111111111111111111

            IF (l < mnl) CYCLE

            me=0
            DO j=1,nib
                g=idx(j)
                IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
            ENDDO
            IF(me>dfmax) EXIT

            me=0
            DO j=1,nif
                g=idxf(j)
                IF(any(theta(ix(g):iy(g),l)/=0.0D0)) me=me+1
            ENDDO
            IF(me>dfmax) EXIT
         

        ENDDO  ! fin de la grnade boucle for l in 1:nlamth 

        DEALLOCATE(b,oldbeta,th,oldtheta,r,rf,oidx,oidxf)
        RETURN
END SUBROUTINE CserMcpUnifier


SUBROUTINE CserScadUnifier (gamm,w,tau,bn,bs,ix,iy,gam,nobs,nvars,x,y,pfmean,pfscale,dfmax,pmax,nlam,flmin,ulam,&
                                eps,maxit,intr,nalam,b0,beta,t0,theta,idx,idxf,nbeta,ntheta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idxf(pmax)
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        INTEGER::ntheta(nlam)
        DOUBLE PRECISION:: w
        DOUBLE PRECISION:: gamm
        DOUBLE PRECISION:: delta
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::tau
        DOUBLE PRECISION::lambda
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pfmean(bn)
        DOUBLE PRECISION::pfscale(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
	DOUBLE PRECISION:: t0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::theta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION:: dd1
        DOUBLE PRECISION::d
        DOUBLE PRECISION::tt
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al  !al represente mu
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: th
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldtheta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::rf
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldth
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidxf
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::nib
        INTEGER::nif
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        DOUBLE PRECISION:: ga1(bn)
        DOUBLE PRECISION:: ga2(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
        ALLOCATE(b(0:nvars), STAT = jerr)
	jerr=jerr + ierr
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(th(0:nvars), STAT = jerr)
	jerr=jerr + ierr
        ALLOCATE(oldtheta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(rf(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidxf(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pfmean pfscale- - -
        IF(maxval(pfmean) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pfmean = max(0.0D0, pfmean)

        IF(maxval(pfscale) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pfscale = max(0.0D0, pfscale)

        ! - - - some initial setup - - -
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
	rf = y
        b = 0.0D0
        delta=2*max(tau,1.0D0-tau)
        oldbeta = 0.0D0
        th = 0.0D0
        oldtheta = 0.0D0
        idx = 0
        idxf = 0
        oidx = 0
        oidxf = 0
        npass = 0
        nib = npass
	nif = npass
        ! --------- lambda loop ----------------------------
        IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF

!----------------------------------------------------------------------------------------

        DO l=1,nlam
            al0 = al
            IF(flmin>=1.0D0) THEN
                al=ulam(l)
                ELSE
                IF(l > 2) THEN
                al=al*alf
                ELSE IF(l==1) THEN
                       al=big
                ELSE IF(l==2) THEN
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        vl=0.0D0
			CALL serdrvBeta(w,tau,nobs,nvars,x,r,rf,vl)
			DO g = 1,bn
				ALLOCATE(u(bs(g)),STAT=ierr)
				IF(ierr/=0) RETURN
				u = vl(ix(g):iy(g))
				ga2(g) = sqrt(dot_product(u,u))
				DEALLOCATE(u)
			END DO
                        vl=0.0D0
			CALL serdrvTh(tau,nobs,nvars,x,rf,vl)
			DO g = 1,bn
				ALLOCATE(u(bs(g)),STAT=ierr)
				IF(ierr/=0) RETURN
				u = vl(ix(g):iy(g))
				ga1(g) = sqrt(dot_product(u,u))
				DEALLOCATE(u)
			END DO
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    al0 = 0.0D0
                    DO g = 1,bn
                                IF(pfscale(g)>0.0D0) THEN
                                al0 = max(al0, ga1(g) / pfscale(g))
                                ENDIF
                    END DO
                    DO g = 1,bn
                                IF(pfmean(g)>0.0D0) THEN
                                al0 = max(al0, ga2(g) / pfmean(g))
                                ENDIF
                    END DO
                    al = al0 * alf
                ENDIF
            ENDIF
        ! --------- outer loop ----------------------------
        DO  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                oldbeta(0)=b(0)
                IF(nib>0) THEN
		        DO j=1,nib
		                g=idx(j)
		                oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
		        ENDDO
                ENDIF

                oldtheta(0)=th(0)
                IF(nif>0) THEN
		        DO j=1,nif
		                g=idxf(j)
		                oldtheta(ix(g):iy(g))=th(ix(g):iy(g))
		        ENDDO
                ENDIF
        
		DO  ! --middle loop-------------------------------------

		npass=npass+1
		dif=0.0D0
			DO g=1,bn  ! theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta

		                start=ix(g)
		                end=iy(g)
		                ALLOCATE(u(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(dd(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(oldth(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                IF(jerr/=0) RETURN
		                oldth=th(start:end)
		                u = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)
					    END IF
					    u = u + dl(i)*x(i,start:end)/nobs
				ENDDO

!!!!!!!!!!!!!1
				u=delta*gam(g)*th(start:end) + u
				unorm=sqrt(dot_product(u,u))

				IF (unorm <= (delta*gam(g)+pfscale(g))*al) THEN
						t=unorm-pfscale(g)*al
						IF(t>0.0D0) THEN
							th(start:end)=u*t/(delta*gam(g)*unorm)
						ELSE
							th(start:end)=0.0D0
						ENDIF
				ELSE IF (unorm <= (delta*gam(g)*al*gamm)) THEN
						tt=unorm-pfscale(g)*al*(gamm/(gamm-1))
						IF(tt>0.0D0) THEN
							th(start:end)=u*tt/((unorm)*(delta*gam(g)-pfscale(g)/(gamm-1.0D0)))
						ELSE
							th(start:end)=0.0D0
						ENDIF
				ELSE
					th(start:end) = u/(delta*gam(g))
				END IF
!!!!!!!!!!!!!!!!


				dd=th(start:end)-oldth
		                IF(any(dd/=0.0D0)) THEN
		                    dif=max(dif,delta*gam(g)**2*dot_product(dd,dd))
		                    rf=rf-matmul(x(:,start:end),dd)
		                    IF(oidxf(g)==0) THEN
		                        nif=nif+1
		                        IF(nif>pmax) EXIT
		                                oidxf(g)=nif
		                                idxf(nif)=g
		                        ENDIF
		                  ENDIF
				DEALLOCATE(u,dd,oldth)
			ENDDO  ! theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta

			DO g=1,bn    ! beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta
		                start=ix(g)
		                end=iy(g)
		                ALLOCATE(u(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(dd(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(oldb(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                IF(jerr/=0) RETURN
		                oldb=b(start:end)
		                u = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)+ w*r(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)+ w*r(i)
					    END IF
					    u = u + dl(i)*x(i,start:end)/nobs
				ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!

				u=(w+delta)*gam(g)*b(start:end) + u
				unorm=sqrt(dot_product(u,u))

                                IF (unorm <= ((w+delta)*gam(g)+pfmean(g))*al) THEN
					t=unorm-pfmean(g)*al
					IF(t>0.0D0) THEN
						b(start:end)=u*t/((w+delta)*gam(g)*unorm)
					ELSE
						b(start:end)=0.0D0
					ENDIF
				ELSE IF (unorm <= ((w+delta)*gam(g)*al*gamm)) THEN
					tt=unorm-pfmean(g)*al*(gamm/(gamm-1))
					IF(tt>0.0D0) THEN
						b(start:end)=u*tt/((unorm)*((w+delta)*gam(g)-pfmean(g)/(gamm-1.0D0)))
					ELSE
						b(start:end)=0.0D0
					ENDIF
				ELSE
					b(start:end) = u/((w+delta)*gam(g))
				END IF

!!!!!!!!!!!!!!!!!!!!!!!!


		                dd=b(start:end)-oldb

		                IF(any(dd/=0.0D0)) THEN
		                    dif=max(dif,delta*gam(g)**2*dot_product(dd,dd))

		                    rf=rf-matmul(x(:,start:end),dd)
				    r=r-matmul(x(:,start:end),dd)
		                    IF(oidx(g)==0) THEN
		                        nib=nib+1
		                        IF(nib>pmax) EXIT
		                                oidx(g)=nib
		                                idx(nib)=g
		                    ENDIF
		                 ENDIF
		                DEALLOCATE(u,dd,oldb)
			ENDDO  ! ! ! beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta


			IF(intr /= 0) THEN !  intercept theta theta theta intercept theta theta theta intercept theta theta theta intercept theta theta th
				d = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)
					    END IF
					    d = d + dl(i) / nobs
				ENDDO
				d = d / delta
				IF(d/=0.0D0) THEN
					th(0)=th(0)+d
					rf=rf-d
					dif=max(dif, d**2)
				ENDIF
			ENDIF   !  !  intercept theta theta theta intercept theta theta theta intercept theta theta theta intercept theta theta th

			IF(intr /= 0) THEN !  intercept beta beta beta intercept beta beta beta intercept beta beta betaintercept
				d = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)+ w*r(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)+ w*r(i)
					    END IF
					    d = d + dl(i) / nobs
				ENDDO
				d = d /(w+delta)
				IF(d/=0.0D0) THEN
					b(0)=b(0)+d
					r=r-d
					rf=rf-d
					dif=max(dif, d**2)
				ENDIF
			ENDIF   !  intercept beta beta beta intercept beta beta beta intercept beta beta beta intercept beta beta beta                     
			IF (nib > pmax) EXIT
			IF (nif > pmax) EXIT
			IF (dif < eps) EXIT

			IF(npass > maxit) THEN
				jerr=-l
				RETURN
			ENDIF
		

		ENDDO   ! --middle loop-------------------------------------
		IF(nib>pmax) EXIT
		IF(nif>pmax) EXIT

		!--- final check ------------------------
		jx = 0
		max_gam = maxval(gam)
		IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
		IF(any((max_gam*(th-oldtheta)/(1+abs(th)))**2 >= eps)) jx = 1
		IF (jx /= 0) cycle
                EXIT
        ENDDO
        !---------- final update variable and save results------------
            IF(nif>pmax) THEN
                jerr=-10000-l
                EXIT
            ENDIF


            IF(nib>0) THEN  !0000000000000000000000
		    DO j=1,nib
		        g=idx(j)
		        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
		    ENDDO
            ENDIF
            nbeta(l)=nib
            b0(l)=b(0)
            alam(l)=al
            nalam=l      !0000000000000000000000

            IF(nif>0) THEN  !111111111111111111111
		    DO j=1,nif
		        g=idxf(j)
		        theta(ix(g):iy(g),l)=th(ix(g):iy(g))
		    ENDDO
            ENDIF
            ntheta(l)=nif
	    t0(l)=th(0)     !11111111111111111111111

            IF (l < mnl) CYCLE

            me=0
            DO j=1,nib
                g=idx(j)
                IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
            ENDDO
            IF(me>dfmax) EXIT

            me=0
            DO j=1,nif
                g=idxf(j)
                IF(any(theta(ix(g):iy(g),l)/=0.0D0)) me=me+1
            ENDDO
            IF(me>dfmax) EXIT
         

        ENDDO  ! fin de la grnade boucle for l in 1:nlamth 

        DEALLOCATE(b,oldbeta,th,oldtheta,r,rf,oidx,oidxf)
        RETURN
END SUBROUTINE CserScadUnifier

!---------------------------------------------------------------------
!---------------------------------------------------------------------

SUBROUTINE aCserMcpUnifier (gamm,w,tau,bn,bs,ix,iy,gam,nobs,nvars,x,y,pfmean,pfscale,dfmax,pmax,nlam,flmin,ulam,&
                                eps,maxit,intr,nalam,b0,beta,t0,theta,idx,idxf,nbeta,ntheta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idxf(pmax)
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        INTEGER::ntheta(nlam)
        DOUBLE PRECISION:: w
        DOUBLE PRECISION:: gamm
        DOUBLE PRECISION:: delta
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::tau
        DOUBLE PRECISION::lambda
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pfmean(bn)
        DOUBLE PRECISION::pfscale(bn)
	DOUBLE PRECISION::bnorm
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
	DOUBLE PRECISION:: t0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::theta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION:: dd1
        DOUBLE PRECISION::d
        DOUBLE PRECISION::tt
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al  !al represente mu
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: th
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldtheta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::rf
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldth
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidxf
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::nib
        INTEGER::nif
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        DOUBLE PRECISION:: ga1(bn)
        DOUBLE PRECISION:: ga2(bn)
        DOUBLE PRECISION:: pfmeanHis(bn)
	DOUBLE PRECISION:: pfscaleHis(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
        ALLOCATE(b(0:nvars), STAT = jerr)
	jerr=jerr + ierr
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(th(0:nvars), STAT = jerr)
	jerr=jerr + ierr
        ALLOCATE(oldtheta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(rf(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidxf(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pfmean pfscale- - -
        IF(maxval(pfmean) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pfmean = max(0.0D0, pfmean)
	pfmeanHis  = pfmean
        IF(maxval(pfscale) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pfscale = max(0.0D0, pfscale)
	pfscaleHis  = pfscale
        ! - - - some initial setup - - -
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
	rf = y
        b = 0.0D0
        delta=2*max(tau,1.0D0-tau)
        oldbeta = 0.0D0
        th = 0.0D0
        oldtheta = 0.0D0
        idx = 0
        idxf = 0
        oidx = 0
        oidxf = 0
        npass = 0
        nib = npass
	nif = npass
        ! --------- lambda loop ----------------------------
        IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF

!----------------------------------------------------------------------------------------

        DO l=1,nlam
            al0 = al
            IF(flmin>=1.0D0) THEN
                al=ulam(l)
                ELSE
                IF(l > 2) THEN
                al=al*alf
                ELSE IF(l==1) THEN
                       al=big

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			CALL serdrvBeta(w,tau,nobs,nvars,x,r,rf,vl)
			DO g = 1,bn
				ALLOCATE(u(bs(g)),STAT=ierr)
				IF(ierr/=0) RETURN
				u = vl(ix(g):iy(g))
				ga2(g) = sqrt(dot_product(u,u))
				DEALLOCATE(u)
			END DO

			CALL serdrvTh(tau,nobs,nvars,x,rf,vl)
			DO g = 1,bn
				ALLOCATE(u(bs(g)),STAT=ierr)
				IF(ierr/=0) RETURN
				u = vl(ix(g):iy(g))
				ga1(g) = sqrt(dot_product(u,u))
				DEALLOCATE(u)
			END DO
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                ELSE IF(l==2) THEN
                    al0 = 0.0D0
                    DO g = 1,bn
                                IF(pfscale(g)>0.0D0) THEN
                                al0 = max(al0, ga1(g) / pfscale(g))
                                ENDIF
                    END DO
                    DO g = 1,bn
                                IF(pfmean(g)>0.0D0) THEN
                                al0 = max(al0, ga2(g) / pfmean(g))
                                ENDIF
                    END DO
                    al = al0 * alf
                ENDIF
            ENDIF
        ! --------- ----------------------------
        ! --------- ----------------------------
        !pf = pfHis
                IF(l>1) THEN
				DO j=1,bn

					g=idx(j)
					bnorm = sqrt(dot_product(b(ix(g):iy(g)),b(ix(g):iy(g))))
					IF  (pfmean(g)*al*gamm>bnorm)   THEN
							pfmean(g) = pfmean(g) - bnorm/(al*gamm)
					ELSE
							pfmean(g) = 0.0D0
					END IF

					bnorm = sqrt(dot_product(th(ix(g):iy(g)),th(ix(g):iy(g))))
					IF  (pfscale(g)*al*gamm>bnorm)   THEN
							pfscale(g)=pfscale(g) - bnorm/(al*gamm)
					ELSE 
							pfscale(g)=0.0D0
					END IF

				ENDDO
		END IF
        ! -------------------------------------
        ! --------- outer loop ----------------------------
        DO  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                oldbeta(0)=b(0)
                IF(nib>0) THEN
		        DO j=1,nib
		                g=idx(j)
		                oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
		        ENDDO
                ENDIF

                oldtheta(0)=th(0)
                IF(nif>0) THEN
		        DO j=1,nif
		                g=idxf(j)
		                oldtheta(ix(g):iy(g))=th(ix(g):iy(g))
		        ENDDO
                ENDIF
        
		DO  ! --middle loop-------------------------------------

		npass=npass+1
		dif=0.0D0
			DO g=1,bn  ! theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta

		                start=ix(g)
		                end=iy(g)
		                ALLOCATE(u(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(dd(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(oldth(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                IF(jerr/=0) RETURN
		                oldth=th(start:end)
		                u = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)
					    END IF
					    u = u + dl(i)*x(i,start:end)/nobs
				ENDDO

				u=delta*gam(g)*th(start:end) + u
				unorm=sqrt(dot_product(u,u))
				t=unorm-pfscale(g)*al
				IF(t>0.0D0)     THEN
							th(start:end)=u*t/(delta*gam(g)*unorm)
						ELSE
							th(start:end)=0.0D0
				ENDIF

				dd=th(start:end)-oldth
		                IF(any(dd/=0.0D0)) THEN
		                    dif=max(dif,delta*gam(g)**2*dot_product(dd,dd))
		                    rf=rf-matmul(x(:,start:end),dd)
		                    IF(oidxf(g)==0) THEN
		                        nif=nif+1
		                        IF(nif>pmax) EXIT
		                                oidxf(g)=nif
		                                idxf(nif)=g
		                        ENDIF
		                  ENDIF
				DEALLOCATE(u,dd,oldth)
			ENDDO  ! theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta

			DO g=1,bn    ! beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta
		                start=ix(g)
		                end=iy(g)
		                ALLOCATE(u(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(dd(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(oldb(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                IF(jerr/=0) RETURN
		                oldb=b(start:end)
		                u = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)+ w*r(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)+ w*r(i)
					    END IF
					    u = u + dl(i)*x(i,start:end)/nobs
				ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!

				u=(w+delta)*gam(g)*b(start:end) + u
				unorm=sqrt(dot_product(u,u))
				t=unorm-pfmean(g)*al
		                IF(t>0.0D0) THEN
		                        b(start:end)=u*t/((w+delta)*gam(g)*unorm)
		                        ELSE
		                        b(start:end)=0.0D0
		                ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!

		                dd=b(start:end)-oldb

		                IF(any(dd/=0.0D0)) THEN
		                    dif=max(dif,delta*gam(g)**2*dot_product(dd,dd))

		                    rf=rf-matmul(x(:,start:end),dd)
				    r=r-matmul(x(:,start:end),dd)
		                    IF(oidx(g)==0) THEN
		                        nib=nib+1
		                        IF(nib>pmax) EXIT
		                                oidx(g)=nib
		                                idx(nib)=g
		                    ENDIF
		                 ENDIF
		                DEALLOCATE(u,dd,oldb)
			ENDDO  ! ! ! beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta


			IF(intr /= 0) THEN !  intercept theta theta theta intercept theta theta theta intercept theta theta theta intercept theta theta th
				d = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)
					    END IF
					    d = d + dl(i) / nobs
				ENDDO
				d = d / delta
				IF(d/=0.0D0) THEN
					th(0)=th(0)+d
					rf=rf-d
					dif=max(dif, d**2)
				ENDIF
			ENDIF   !  !  intercept theta theta theta intercept theta theta theta intercept theta theta theta intercept theta theta th

			IF(intr /= 0) THEN !  intercept beta beta beta intercept beta beta beta intercept beta beta betaintercept
				d = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)+ w*r(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)+ w*r(i)
					    END IF
					    d = d + dl(i) / nobs
				ENDDO
				d = d /(w+delta)
				IF(d/=0.0D0) THEN
					b(0)=b(0)+d
					r=r-d
					rf=rf-d
					dif=max(dif, d**2)
				ENDIF
			ENDIF   !  intercept beta beta beta intercept beta beta beta intercept beta beta beta intercept beta beta beta                  
			IF (nib > pmax) EXIT
			IF (nif > pmax) EXIT
			IF (dif < eps) EXIT

			IF(npass > maxit) THEN
				jerr=-l
				RETURN
			ENDIF
		

		ENDDO   ! --middle loop-------------------------------------
		IF(nib>pmax) EXIT
		IF(nif>pmax) EXIT

		!--- final check ------------------------
		jx = 0
		max_gam = maxval(gam)
		IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
		IF(any((max_gam*(th-oldtheta)/(1+abs(th)))**2 >= eps)) jx = 1
		IF (jx /= 0) cycle
                EXIT
        ENDDO
        !---------- final update variable and save results------------
            IF(nif>pmax) THEN
                jerr=-10000-l
                EXIT
            ENDIF


            IF(nib>0) THEN  !0000000000000000000000
		    DO j=1,nib
		        g=idx(j)
		        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
		    ENDDO
            ENDIF
            nbeta(l)=nib
            b0(l)=b(0)
            alam(l)=al
            nalam=l      !0000000000000000000000

            IF(nif>0) THEN  !111111111111111111111
		    DO j=1,nif
		        g=idxf(j)
		        theta(ix(g):iy(g),l)=th(ix(g):iy(g))
		    ENDDO
            ENDIF
            ntheta(l)=nif
	    t0(l)=th(0)     !11111111111111111111111

            IF (l < mnl) CYCLE

            me=0
            DO j=1,nib
                g=idx(j)
                IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
            ENDDO
            IF(me>dfmax) EXIT

            me=0
            DO j=1,nif
                g=idxf(j)
                IF(any(theta(ix(g):iy(g),l)/=0.0D0)) me=me+1
            ENDDO
            IF(me>dfmax) EXIT
         

        ENDDO  ! fin de la grnade boucle for l in 1:nlamth 

        DEALLOCATE(b,oldbeta,th,oldtheta,r,rf,oidx,oidxf)
        RETURN
END SUBROUTINE aCserMcpUnifier

!-------------------------------------------------------------------------------

SUBROUTINE aCserScadUnifier (gamm,w,tau,bn,bs,ix,iy,gam,nobs,nvars,x,y,pfmean,pfscale,dfmax,pmax,nlam,flmin,ulam,&
                                eps,maxit,intr,nalam,b0,beta,t0,theta,idx,idxf,nbeta,ntheta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idxf(pmax)
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        INTEGER::ntheta(nlam)
        DOUBLE PRECISION:: w
        DOUBLE PRECISION:: gamm
        DOUBLE PRECISION:: delta
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::tau
        DOUBLE PRECISION::lambda
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pfmean(bn)
        DOUBLE PRECISION::pfscale(bn)
	DOUBLE PRECISION::bnorm
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
	DOUBLE PRECISION:: t0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::theta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION:: dd1
        DOUBLE PRECISION::d
        DOUBLE PRECISION::tt
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al  !al represente mu
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: th
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldtheta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::rf
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldth
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidxf
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::nib
        INTEGER::nif
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        DOUBLE PRECISION:: ga1(bn)
        DOUBLE PRECISION:: ga2(bn)
        DOUBLE PRECISION:: pfmeanHis(bn)
	DOUBLE PRECISION:: pfscaleHis(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
        ALLOCATE(b(0:nvars), STAT = jerr)
	jerr=jerr + ierr
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(th(0:nvars), STAT = jerr)
	jerr=jerr + ierr
        ALLOCATE(oldtheta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(rf(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidxf(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pfmean pfscale- - -
        IF(maxval(pfmean) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pfmean = max(0.0D0, pfmean)
	pfmeanHis  = pfmean
        IF(maxval(pfscale) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pfscale = max(0.0D0, pfscale)
	pfscaleHis  = pfscale
        ! - - - some initial setup - - -
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
	rf = y
        b = 0.0D0
        delta=2*max(tau,1.0D0-tau)
        oldbeta = 0.0D0
        th = 0.0D0
        oldtheta = 0.0D0
        idx = 0
        idxf = 0
        oidx = 0
        oidxf = 0
        npass = 0
        nib = npass
	nif = npass
        ! --------- lambda loop ----------------------------
        IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF

!----------------------------------------------------------------------------------------

        DO l=1,nlam
            al0 = al
            IF(flmin>=1.0D0) THEN
                al=ulam(l)
                ELSE
                IF(l > 2) THEN
                al=al*alf
                ELSE IF(l==1) THEN
                       al=big

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			CALL serdrvBeta(w,tau,nobs,nvars,x,r,rf,vl)
			DO g = 1,bn
				ALLOCATE(u(bs(g)),STAT=ierr)
				IF(ierr/=0) RETURN
				u = vl(ix(g):iy(g))
				ga2(g) = sqrt(dot_product(u,u))
				DEALLOCATE(u)
			END DO

			CALL serdrvTh(tau,nobs,nvars,x,rf,vl)
			DO g = 1,bn
				ALLOCATE(u(bs(g)),STAT=ierr)
				IF(ierr/=0) RETURN
				u = vl(ix(g):iy(g))
				ga1(g) = sqrt(dot_product(u,u))
				DEALLOCATE(u)
			END DO
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                ELSE IF(l==2) THEN
                    al0 = 0.0D0
                    DO g = 1,bn
                                IF(pfscale(g)>0.0D0) THEN
                                al0 = max(al0, ga1(g) / pfscale(g))
                                ENDIF
                    END DO
                    DO g = 1,bn
                                IF(pfmean(g)>0.0D0) THEN
                                al0 = max(al0, ga2(g) / pfmean(g))
                                ENDIF
                    END DO
                    al = al0 * alf
                ENDIF
            ENDIF
        ! --------- ----------------------------
        ! --------- ----------------------------
        !pf = pfHis
        IF(l>1) THEN
				DO j=1,bn

					g=idx(j)
					bnorm = sqrt(dot_product(b(ix(g):iy(g)),b(ix(g):iy(g))))
					IF  (pfmean(g)*al>bnorm)   THEN
							pfmean(g)=pfmean(g)
					ELSE IF ((pfmean(g)*al*gamm>bnorm)) THEN
							pfmean(g)=(gamm/(gamm-1))*pfmean(g)-bnorm/(al*(gamm-1))
					ELSE
							pfmean(g)=0.0D0
					END IF

					bnorm = sqrt(dot_product(th(ix(g):iy(g)),th(ix(g):iy(g))))
					IF  (pfscale(g)*al>bnorm)   THEN
							pfscale(g)=pfscale(g)
					ELSE IF ((pfscale(g)*al*gamm>bnorm)) THEN
							pfscale(g)=(gamm/(gamm-1))*pfscale(g)-bnorm/(al*(gamm-1))
					ELSE
							pfscale(g)=0.0D0
					END IF

				ENDDO
		END IF

        ! -------------------------------------
        ! -------------------------------------
        ! --------- outer loop ----------------------------
        DO  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                oldbeta(0)=b(0)
                IF(nib>0) THEN
		        DO j=1,nib
		                g=idx(j)
		                oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
		        ENDDO
                ENDIF

                oldtheta(0)=th(0)
                IF(nif>0) THEN
		        DO j=1,nif
		                g=idxf(j)
		                oldtheta(ix(g):iy(g))=th(ix(g):iy(g))
		        ENDDO
                ENDIF
        
		DO  ! --middle loop-------------------------------------

		npass=npass+1
		dif=0.0D0
			DO g=1,bn  ! theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta

		                start=ix(g)
		                end=iy(g)
		                ALLOCATE(u(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(dd(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(oldth(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                IF(jerr/=0) RETURN
		                oldth=th(start:end)
		                u = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)
					    END IF
					    u = u + dl(i)*x(i,start:end)/nobs
				ENDDO

				u=delta*gam(g)*th(start:end) + u
				unorm=sqrt(dot_product(u,u))
				t=unorm-pfscale(g)*al
				IF(t>0.0D0)     THEN
							th(start:end)=u*t/(delta*gam(g)*unorm)
						ELSE
							th(start:end)=0.0D0
				ENDIF

				dd=th(start:end)-oldth
		                IF(any(dd/=0.0D0)) THEN
		                    dif=max(dif,delta*gam(g)**2*dot_product(dd,dd))
		                    rf=rf-matmul(x(:,start:end),dd)
		                    IF(oidxf(g)==0) THEN
		                        nif=nif+1
		                        IF(nif>pmax) EXIT
		                                oidxf(g)=nif
		                                idxf(nif)=g
		                        ENDIF
		                  ENDIF
				DEALLOCATE(u,dd,oldth)
			ENDDO  ! theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta theta

			DO g=1,bn    ! beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta
		                start=ix(g)
		                end=iy(g)
		                ALLOCATE(u(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(dd(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                ALLOCATE(oldb(bs(g)),STAT=ierr)
		                jerr=jerr+ierr
		                IF(jerr/=0) RETURN
		                oldb=b(start:end)
		                u = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)+ w*r(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)+ w*r(i)
					    END IF
					    u = u + dl(i)*x(i,start:end)/nobs
				ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!

				u=(w+delta)*gam(g)*b(start:end) + u
				unorm=sqrt(dot_product(u,u))
				t=unorm-pfmean(g)*al
		                IF(t>0.0D0) THEN
		                        b(start:end)=u*t/((w+delta)*gam(g)*unorm)
		                        ELSE
		                        b(start:end)=0.0D0
		                ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!

		                dd=b(start:end)-oldb

		                IF(any(dd/=0.0D0)) THEN
		                    dif=max(dif,delta*gam(g)**2*dot_product(dd,dd))

		                    rf=rf-matmul(x(:,start:end),dd)
				    r=r-matmul(x(:,start:end),dd)
		                    IF(oidx(g)==0) THEN
		                        nib=nib+1
		                        IF(nib>pmax) EXIT
		                                oidx(g)=nib
		                                idx(nib)=g
		                    ENDIF
		                 ENDIF
		                DEALLOCATE(u,dd,oldb)
			ENDDO  ! ! ! beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta beta


			IF(intr /= 0) THEN !  intercept theta theta theta intercept theta theta theta intercept theta theta theta intercept theta theta th
				d = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)
					    END IF
					    d = d + dl(i) / nobs
				ENDDO
				d = d / delta
				IF(d/=0.0D0) THEN
					th(0)=th(0)+d
					rf=rf-d
					dif=max(dif, d**2)
				ENDIF
			ENDIF   !  !  intercept theta theta theta intercept theta theta theta intercept theta theta theta intercept theta theta th

			IF(intr /= 0) THEN !  intercept beta beta beta intercept beta beta beta intercept beta beta betaintercept
				d = 0.0D0
				DO i = 1,nobs
					    IF (rf(i) < 0) THEN
						dl (i) = 2.0D0*(1.0D0 - tau)*rf(i)+ w*r(i)
					    ELSE
						dl (i) = 2.0D0*(tau)*rf(i)+ w*r(i)
					    END IF
					    d = d + dl(i) / nobs
				ENDDO
				d = d /(w+delta)
				IF(d/=0.0D0) THEN
					b(0)=b(0)+d
					r=r-d
					rf=rf-d
					dif=max(dif, d**2)
				ENDIF
			ENDIF   !  intercept beta beta beta intercept beta beta beta intercept beta beta beta intercept beta beta beta                  
			IF (nib > pmax) EXIT
			IF (nif > pmax) EXIT
			IF (dif < eps) EXIT

			IF(npass > maxit) THEN
				jerr=-l
				RETURN
			ENDIF
		

		ENDDO   ! --middle loop-------------------------------------
		IF(nib>pmax) EXIT
		IF(nif>pmax) EXIT

		!--- final check ------------------------
		jx = 0
		max_gam = maxval(gam)
		IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
		IF(any((max_gam*(th-oldtheta)/(1+abs(th)))**2 >= eps)) jx = 1
		IF (jx /= 0) cycle
                EXIT
        ENDDO
        !---------- final update variable and save results------------
            IF(nif>pmax) THEN
                jerr=-10000-l
                EXIT
            ENDIF


            IF(nib>0) THEN  !0000000000000000000000
		    DO j=1,nib
		        g=idx(j)
		        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
		    ENDDO
            ENDIF
            nbeta(l)=nib
            b0(l)=b(0)
            alam(l)=al
            nalam=l      !0000000000000000000000

            IF(nif>0) THEN  !111111111111111111111
		    DO j=1,nif
		        g=idxf(j)
		        theta(ix(g):iy(g),l)=th(ix(g):iy(g))
		    ENDDO
            ENDIF
            ntheta(l)=nif
	    t0(l)=th(0)     !11111111111111111111111

            IF (l < mnl) CYCLE

            me=0
            DO j=1,nib
                g=idx(j)
                IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
            ENDDO
            IF(me>dfmax) EXIT

            me=0
            DO j=1,nif
                g=idxf(j)
                IF(any(theta(ix(g):iy(g),l)/=0.0D0)) me=me+1
            ENDDO
            IF(me>dfmax) EXIT
         

        ENDDO  ! fin de la grnade boucle for l in 1:nlamth 

        DEALLOCATE(b,oldbeta,th,oldtheta,r,rf,oidx,oidxf)
        RETURN
END SUBROUTINE aCserScadUnifier

