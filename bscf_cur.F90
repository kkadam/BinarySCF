       subroutine binary_scf(iam,ra,rb,rc,initial,model_num,rhom1,rhom2) 
       implicit none 
       include 'startup.h' 
 
       real, dimension(numr,numz,numphi) :: rho, pot, pot_it, pot_old, h
       real, dimension(numr,numz,numphi) :: temp, rchpot
       real, dimension(numr,numphi) :: psi
       real, dimension(maxit) :: c1, c2, mass1, mass2, omsq, hm1, hm2
       real, dimension(numr) :: r, rhf
       real, dimension(numz) :: z, zhf
       real, dimension(numphi) :: phi, sine, cosine
       real :: com, pi, rhom1, rhom2, ret1, ret2, dr, dz, dphi, dpot, dpsi
       real :: factor, gamma, n1, cnvgom, cnvgc1, cnvgc2, cnvgh1, cnvgh2
       real :: s1, s2, stot, t1, t2, ttot, w1, w2, wtot, j1, j2, jtot
       real :: e1, e2, etot, en1, en2, entot, pm1, pm2, virialerr, virialerr1
       real :: virialerr2, eps, position, separation, xavg1, xavg2
       real :: yavg1, yavg2, kappa1, kappa2, vol1, vol2, reff1, reff2
       real :: rchmax, rchmin, xcrit, rpotcrit, volr1, volr2, reffr1, reffr2
       real :: rchtest, curvature, pottmp1, pottmp2, psitmp1, psitmp2, omega
       real :: kepler, period, epsilon, rho_max_temp
       real :: timef, stime, ftime
       integer :: i, j, k, q, iam
       integer :: ra, phia, za, rb, phib, zb, rc, phic, zc, rm1, phim1, zm1
       integer :: rm2, phim2, zm2, qfinal, rmax, star2maxr
       integer :: flag, isave, rochemax1, rochemax2
       integer, dimension(3) :: rminloc, rmaxloc
       integer :: primary, initial, start_center_1, start_center_2
       integer :: louter1, louter2, model_num
       real :: radius_start_1, radius_start_2, d
       character(len=50) :: dens_file, dens_template

!  Variables for poisson solver
       real dummy2(numphi), dummy3(numphi), grav
       real dd3(numr), dd4(numz), dd5, dd6, dd7, dd8
       common /pois/ pot, rho
       common /blok6/ dphi, dummy2, dummy3, pi, grav
       common /grid/ r, z, rhf, zhf, dd3, dd4, dd5, dd6, dd7, dd8

       call cpu_time(stime)
       
!  Initialize the grid
       grav = 1.0
       pi = acos(-1.0)
       dr = 1.0/(numr-1.0)
       dz = dr
       dphi = 2.0*pi/numphi
       r(1) = - dr
       rhf(1) = - dr/2.0
       z(1) = - dz
       zhf(1) = - dz/2.0 
       phi(1) = 0.0
       do i = 2,numr
          r(i) = r(i-1) + dr 
          rhf(i) = rhf(i-1) + dr
       enddo
       do i = 2,numz
          z(i) = z(i-1) + dz
          zhf(i) = zhf(i-1) + dz
       enddo
       do i = 2,numphi
          phi(i) = phi(i-1) + dphi
       enddo
       do i = 1,numphi
          cosine(i) = cos(phi(i))
          sine(i) = sin(phi(i))
       enddo
       factor = 2.0*dr*dz*dphi

       phia = 1
       phib = 1
       phic = numphi/2 + 1
       za = 2
       zb = 2
       zc = 2
       eps = 2.0e-4
       n1 = 1.5

       rmax = numr - 8
       gamma = 1.0 + 1.0/n1
       epsilon = 1.0e-5

       dens_template = 'density'

!  Re-Initialize variables
       omsq = 0.0
       hm1 = 0.0
       hm2 = 0.0
       c1 = 0.0
       c2 = 0.0


! generate the initial model
       select case(initial) 
          case (1)		!  roll your own initial Gaussian rho
          start_center_1 = ( ra - rb )  / 2 + rb
          start_center_2 = ( rmax - rc ) / 2 + rc
          do i = 1, phi1+1
             do j = 2, numz
                do k = 2, rmax
                   d = (rhf(k)*cosine(i) -                           &
                              rhf(start_center_1)*cosine(phia))*     &
                      (rhf(k)*cosine(i) -                            &
                              rhf(start_center_1)*cosine(phia)) +    &
                      (rhf(k)*sine(i) -                              &
                              rhf(start_center_1)*sine(phia))*       &
                      (rhf(k)*sine(i) -                              &
                              rhf(start_center_1)*sine(phia)) +      &
                      zhf(j)*zhf(j) 
                   rho(k,j,i) = rhom1*exp(-20.0*d)
                enddo 
             enddo 
          enddo 
          do i = phi2, phi3+1
             do j = 2, numz
                do k = 2, rmax
                   d = (rhf(k)*cosine(i) -                          &
                              rhf(start_center_2)*cosine(phic))*    &
                      (rhf(k)*cosine(i) -                           &
                              rhf(start_center_2)*cosine(phic)) +   &
                      (rhf(k)*sine(i) -                             &
                              rhf(start_center_2)*sine(phic))*      &
                      (rhf(k)*sine(i) -                             &
                              rhf(start_center_2)*sine(phic)) +     &
                      zhf(j)*zhf(j)
                   rho(k,j,i) = rhom2*exp(-20.0*d)
                enddo 
             enddo 
          enddo 
          do i = phi4, numphi
             do j = 2, numz
                do k = 2, rmax
                   d = (rhf(k)*cosine(i) -                          &
                              rhf(start_center_1)*cosine(phia))*    &
                      (rhf(k)*cosine(i) -                           &
                              rhf(start_center_1)*cosine(phia)) +   &
                      (rhf(k)*sine(i) -                             &
                              rhf(start_center_1)*sine(phia))*      &
                      (rhf(k)*sine(i) -                             &
                              rhf(start_center_1)*sine(phia)) +     &
                      zhf(j)*zhf(j)
                   rho(k,j,i) = rhom1*exp(-20.0*d)
                enddo 
             enddo 
          enddo 
          rho(:,1,:) = rho(:,2,:)	!  equatoiral bounary condition
          rho(1,:,:) = cshift(rho(2,:,:),dim=2,shift=numphi/2) ! axial b.c.

          case (2)		!  roll your own uniform rho
          start_center_1 = ( ra - rb ) / 2 + rb
          start_center_2 = (rmax - rc ) /2 + rc
          radius_start_1 = (rhf(ra) - rhf(start_center_1))*        &
                          (rhf(ra) - rhf(start_center_1))
          radius_start_2 = (rhf(start_center_2) - rhf(rc))*        &
                          (rhf(start_center_2) - rhf(rc))
          do i = 1, phi1+1
             do j = 2, numz
                do k = 2, rmax
                   d = (rhf(k)*cosine(i) -                         &
                              rhf(start_center_1)*cosine(phia))*   &
                      (rhf(k)*cosine(i) -                          &
                              rhf(start_center_1)*cosine(phia)) +  &
                      (rhf(k)*sine(i) -                            &
                              rhf(start_center_1)*sine(phia))*     &
                      (rhf(k)*sine(i) -                            &
                              rhf(start_center_1)*sine(phia)) +    &
                      zhf(j)*zhf(j)
                   if( d <= radius_start_1 ) then
                      rho(k,j,i) = rhom1
                   else
                      rho(k,j,i) = 0.0
                   endif
                enddo 
             enddo 
          enddo 
          do i = phi2, phi3+1
             do j = 2, numz
                do k = 2, rmax
                   d = (rhf(k)*cosine(i) -                         &
                              rhf(start_center_2)*cosine(phic))*   &
                      (rhf(k)*cosine(i) -                          &
                              rhf(start_center_2)*cosine(phic)) +  &
                      (rhf(k)*sine(i) -                            & 
                              rhf(start_center_2)*sine(phic))*     &
                      (rhf(k)*sine(i) -                            &
                              rhf(start_center_2)*sine(phic)) +    & 
                      zhf(j)*zhf(j)
                   if( d <= radius_start_2 ) then
                      rho(k,j,i) = rhom2
                   else
                      rho(k,j,i) = 0.0
                   endif 
                enddo 
             enddo 
          enddo 
          do i = phi4, numphi
             do j = 2, numz
                do k = 2, rmax
                   d = (rhf(k)*cosine(i) -                         &
                              rhf(start_center_1)*cosine(phia))*   &
                      (rhf(k)*cosine(i) -                          &
                              rhf(start_center_1)*cosine(phia)) +  &
                      (rhf(k)*sine(i) -                            &
                              rhf(start_center_1)*sine(phia))*     &
                      (rhf(k)*sine(i) -                            &
                              rhf(start_center_1)*sine(phia)) +    &
                      zhf(j)*zhf(j)
                   if( d <= radius_start_1 ) then
                      rho(k,j,i) = rhom1
                   else
                      rho(k,j,i) = 0.0
                   endif 
                enddo 
             enddo 
          enddo 
          rho(:,1,:) = rho(:,2,:)	!  equatoiral bounary condition
          rho(1,:,:) = cshift(rho(2,:,:),dim=2,shift=numphi/2) ! axial b.c.

       end select

!  Calculate the initial total mass for each star
       do i = 1,numphi
          do j = 2, numz
             do k = 2, numr
                temp(k,j,i) = rhf(k)*rho(k,j,i)
             enddo
          enddo
       enddo
       call bin_sum(temp, ret1, ret2)
       mass1(1) = factor*ret1
       mass2(1) = factor*ret2

!  Calculate the initial center of mass
       do i = 1, numphi
          do j = 2,numz
             do k = 2,numr
                temp(k,j,i) = rhf(k)*rhf(k)*cosine(i)*rho(k,j,i)
             enddo
          enddo
       enddo
       call bin_sum(temp,ret1, ret2)
       xavg1 = factor*ret1/mass1(1)
       xavg2 = factor*ret2/mass2(1)
       separation = xavg1 - xavg2
       
       do i = 1, numphi
          do j = 2, numz
             do k = 2, numr
                com = com + rhf(k)*rhf(k)*cosine(i)*rho(k,j,i)
             enddo
          enddo
       enddo
       com = factor * com / (mass1(1) + mass2(1))

       
!  START OF THE ITERATION
       do q = 2,maxit		
          print*, "Iteration number = ",q
!  Solve Poisson's equation for gravitational potential, calculate
!  the new centrifugal potential, angular frequency and integration
!  constants
          if( q == 2 ) then
             call setbdy(icall,isym)
             call bdygen(maxterm,isym,redge)
             call pot3(npoint,iprint,isym)
             pot_it = pot
             pot_old = pot
          else
             call bdygen(maxterm,isym,redge)
             call pot3(npoint,iprint,isym)
             pot_it = 0.5*(pot + pot_old)
             pot_old = pot
          endif 
          
          pottmp1 = 0.5*(pot_it(ra,za,phia) + pot_it(ra-1,za,phia))
          pottmp2 = 0.5*(pot_it(rb,zb,phib) + pot_it(rb+1,zb,phib))
          dpot = pottmp1 - pottmp2
          
          do i = 1,numphi
             do j = 2,numr
                psi(j,i) = -0.5*( (rhf(j)*cosine(i)-com)**2 +        &
                           rhf(j)*rhf(j)*sine(i)*sine(i))
             enddo
          enddo
          
          psitmp1 = 0.5*(psi(rb,phib) + psi(rb+1,phib))
          psitmp2 = 0.5*(psi(ra,phia) + psi(ra-1,phia))
          dpsi = psitmp1 - psitmp2 
          omsq(q) = dpot/dpsi
          
          pottmp1 = 0.5*(pot_it(rb,zb,phib) + pot_it(rb+1,zb,phib))
          psitmp1 = 0.5*(psi(rb,phib) + psi(rb+1,phib))
          pottmp2 = 0.5*(pot_it(rc,zc,phic) + pot_it(rc+1,zc,phic))
          psitmp2 = 0.5*(psi(rc,phic) + psi(rc+1,phic))
          c1(q) = pottmp1 + omsq(q)*psitmp1
          c2(q) = pottmp2 + omsq(q)*psitmp2

!  Calculate new enthalpy field and max enthalpies and their positions
          do i = 1,phi1+1
             do j = 2,numz
                do k = 2,rmax
                   h(k,j,i) = c1(q) - pot_it(k,j,i) - omsq(q)*psi(k,i)
                   if(h(k,j,i).gt.hm1(q)) then
                      hm1(q) = h(k,j,i)
                      rm1 = k
                      zm1 = j
                      phim1 = i
                   endif
                enddo
             enddo
          enddo
          do i = phi2,phi3+1
             do j = 2,numz
                do k = 2,rmax
                   h(k,j,i) = c2(q) - pot_it(k,j,i) - omsq(q)*psi(k,i)
                   if(h(k,j,i).gt.hm2(q)) then
                      hm2(q) = h(k,j,i)
                      rm2 = k
                      zm2 = j
                      phim2 = i
                   endif
                enddo
             enddo
          enddo
          do i = phi4,numphi
             do j = 2,numz
                do k = 2,rmax
                   h(k,j,i) = c1(q) - pot_it(k,j,i) - omsq(q)*psi(k,i)
                   if(h(k,j,i).gt.hm1(q)) then 
                      hm1(q) = h(k,j,i)
                      rm1 = k
                      zm1 = j
                      phim1 = i
                   endif
                enddo
             enddo
          enddo 

!  Calculate new density field from enthalpy field
          do i = 1,phi1+1
             do j = 2,numz
                do k = 2,numr
                   if(h(k,j,i).gt.0.0) then
                      rho(k,j,i) = rhom1*(h(k,j,i)/hm1(q))**n1
                   else
                      rho(k,j,i) = 0.0
                   endif
                enddo
             enddo
          enddo
          do i = phi2,phi3+1
             do j = 2,numz
                do k = 2,numr
                   if(h(k,j,i).gt.0.0) then
                      rho(k,j,i) = rhom2*(h(k,j,i)/hm2(q))**n1
                   else
                      rho(k,j,i) = 0.0
                   endif
                enddo
             enddo
          enddo
          do i = phi4,numphi
             do j = 2,numz
                do k = 2,numr
                   if(h(k,j,i).gt.0.0) then
                      rho(k,j,i) = rhom1*(h(k,j,i)/hm1(q))**n1
                   else
                      rho(k,j,i) = 0.0
                   endif
                enddo
             enddo
          enddo
          
!  zero out grid between inner boundary points for both stars and
!  the z axis 
          do i = phi4, numphi
             do j = 2, numz
                do k = 2, rb
                   rho(k,j,i) = 0.0
                enddo
             enddo
          enddo
          do i = 1, phi1+1
             do j = 2, numz
                do k = 2, rb
                   rho(k,j,i) = 0.0
                enddo
             enddo
          enddo 
          do i = phi2, phi3+1
             do j = 2, numz
                do k = 2, rc
                   rho(k,j,i) = 0.0
                enddo
             enddo
          enddo
          rho(:,1,:) = rho(:,2,:)	! equatorial boundary condition
          rho(1,:,:) = cshift(rho(2,:,:),dim=2,shift=numphi/2)	! axial b.c. 

!  Calculate total mass in each star
          do i = 1, numphi
             do j = 2, numz
                do k = 2, numr
                   temp(k,j,i) = rhf(k)*rho(k,j,i)
                enddo
             enddo
          enddo
          call bin_sum(temp, ret1, ret2)
          mass1(q) = factor*ret1
          mass2(q) = factor*ret2

!  Calculate center of mass for current density field
          xavg1 = 0.0
          xavg2 = 0.0
          separation = 0.0
          com = 0.0
          do i = 1, numphi
             do j = 2, numz
                do k = 2, numr
                   temp(k,j,i) = rhf(k)*rhf(k)*cosine(i)*rho(k,j,i)
                enddo
             enddo
          enddo
          call bin_sum(temp, ret1, ret2)
          xavg1 = factor*ret1/mass1(q)
          xavg2 = factor*ret2/mass2(q)
          separation = xavg1 - xavg2
          do i = 1, numphi
             do j = 2, numz
                do k = 2, numr
                   com = com + rhf(k)*rhf(k)*cosine(i)*rho(k,j,i)
                 enddo
             enddo
          enddo
          com = factor * com / (mass1(q)+mass2(q))

          !  Has solution converged sufficiently? 
          cnvgom = abs(omsq(q)-omsq(q-1))/abs(omsq(q))
          cnvgc1 = abs(c1(q)-c1(q-1))/abs(c1(q))
          cnvgc2 = abs(c2(q)-c2(q-1))/abs(c2(q))
          cnvgh1 = abs(hm1(q)-hm1(q-1))/abs(hm1(q))
          cnvgh2 = abs(hm2(q)-hm2(q-1))/abs(hm2(q))
          if(cnvgom.lt.eps.and.cnvgc1.lt.eps.and.cnvgc2.lt.eps) then
             if(cnvgh1.lt.eps.and.cnvgh2.lt.eps) then
                qfinal = q
                goto 10
             else
                if(q.eq.maxit) then 
                   goto 10
                endif
             endif
          else
             if(q.eq.maxit) then 
                goto 10
             endif
          endif
       enddo			!  END OF ITERATION CYCLE
 10    continue
       qfinal = q

        !   Set an integer flag to indicate which star, 1 or 2, is the more massive
        primary = 1
        if (mass2(qfinal) .gt. mass1(qfinal)) then
           primary = 2
        endif

        omega = sqrt(omsq(qfinal))

       !  Compute virial preasure for each star
       do i = 1, numphi
          do j = 2, numz
             do k = 2, numr
                temp(k,j,i) = rhf(k)*h(k,j,i)*rho(k,j,i)
             enddo
          enddo
       enddo
       call bin_sum(temp, ret1, ret2)
       s1 = factor*ret1/(n1+1.0)
       s2 = factor*ret2/(n1+1.0)
       stot = s1 + s2
       
       !  Compute total potential energy of each star
       do i = 1, numphi
          do j = 2, numz
             do k = 2, numr
                temp(k,j,i) = rhf(k)*pot_it(k,j,i)*rho(k,j,i)
             enddo
          enddo
       enddo
       call bin_sum(temp, ret1, ret2)
       w1 = 0.5*factor*ret1
       w2 = 0.5*factor*ret2
       wtot = w1 + w2
       
       !  Compute total kinetic energy of rotation for each star
       do i = 1, numphi
          do j = 2, numz
             do k = 2, numr
                temp(k,j,i) = rhf(k)*psi(k,i)*rho(k,j,i)
             enddo
          enddo
       enddo
       call bin_sum(temp, ret1, ret2)
       t1 = -omega*omega*factor*ret1
       t2 = -omega*omega*factor*ret2
       ttot = t1 + t2
      
       virialerr = abs(2.0*ttot + 3.0*stot + wtot) / abs(wtot)
       virialerr1 = abs(2.0*t1 + 3.0*s1 + w1) / abs(w1)
       virialerr2 = abs(2.0*t2 + 3.0*s2 + w2) / abs(w2)
       pm1 = rho(rm1,zm1,phim1)*hm1(qfinal)/(n1+1.0)
       pm2 = rho(rm2,zm2,phim2)*hm2(qfinal)/(n1+1.0)
       kappa1 = pm1/rhom1**gamma
       kappa2 = pm2/rhom2**gamma
       period = 2.0*pi/omega
       kepler = (separation**3)*omega*omega /                          &
               (mass1(qfinal)+mass2(qfinal))

       !  Calculate angular momentum
       do i = 1, numphi
          do j = 2, numz
             do k = 2, numr
                temp(k,j,i) = rhf(k)*psi(k,i)*rho(k,j,i)
             enddo
          enddo
       enddo
       call bin_sum(temp, ret1, ret2)
       j1 = -2.0*omega*factor*ret1
       j2 = -2.0*omega*factor*ret2
       jtot = j1 + j2

       !  Calculate internal energies
       do i = 1, numphi
          do j = 2, numz
             do k = 2, numr
                temp(k,j,i) = rhf(k)*rho(k,j,i)**gamma
             enddo
          enddo
       enddo
       call bin_sum(temp, ret1, ret2)
       e1 = n1*factor*kappa1*ret1
       e2 = n1*factor*kappa2*ret2
       etot = e1 + e2

       !  Calculate total energoies
       en1 = t1 + e1 + w1
       en2 = t2 + e2 + w2
       entot  = ttot + etot + wtot

       !  Compute total volume for each star
       do i = 1, numphi
          do j = 2, numz
             do k = 2, numr
                if(rho(k,j,i) > 0.0) then
                   temp(k,j,i) = rhf(k)
                endif
             enddo
          enddo
       enddo
       call bin_sum(temp, ret1, ret2)
       vol1 = factor*ret1
       vol2 = factor*ret2
       reff1 = (0.75*vol1/pi)**0.3333333
       reff2 = (0.75*vol2/pi)**0.3333333
 
       !  Compute y moment of density distribution 
       do i = 1, numphi
          do j = 2, numz
             do k = 2, numr
                temp(k,j,i) = rhf(k)*rhf(k)*sine(i)*rho(k,j,i)
             enddo
          enddo
       enddo
       call bin_sum(temp, ret1, ret2)
       yavg1 = factor*ret1/mass1(qfinal)
       yavg2 = factor*ret2/mass2(qfinal)

       !  Calculate the Roche Potential
       do i = 1,numphi
          do j = 2,numz
             do k = 2,numr
                rchpot(k,j,i) = pot_it(k,j,i) + omega*omega*psi(k,i)
             enddo
          enddo
       enddo
       rchpot(1,:,:) = cshift(rchpot(2,:,:),dim=2,shift=numphi/2)

       ! Find min and max values of the Roche potential
       rminloc = minloc(rchpot)
       rmaxloc = maxloc(rchpot)
       rchmin = rchpot(rminloc(1),rminloc(2),rminloc(3))
       rchmax = rchpot(rmaxloc(1),rmaxloc(2),rmaxloc(3))

       ! Find Inner Lagrange Point and Roche Potential at L1
       ! flag = 0 for L1 in *2 (-ve x), flag = 1 for L1 in *1 (+ve x)
       flag = 0
       isave = 0
       j = 2 
       i = phic
       do k = rm2,2,-1
          rchtest = (rchpot(k,j,i)-rchpot(k+1,j,i))*                     &
                   (rchpot(k-1,j,i)-rchpot(k,j,i))
          if(rchtest.lt.0.0) then
             curvature = rchpot(k+1,j,i) + rchpot(k-1,j,i) -             &
                        2.0*rchpot(k,j,i)
             if(curvature.lt.0.0) then 
                  isave = k
             endif 
          endif
       enddo
       !  L1 not on -ve x axis, try on +ve x axis
       if(isave.eq.0) then
          i = 1
          do k = 2,rm1
             rchtest = (rchpot(k+1,j,i)-rchpot(k,j,i))*                  &
             (rchpot(k,j,i)-rchpot(k-1,j,i))
             if(rchtest.lt.0.0) then
                curvature = rchpot(k+1,j,i) + rchpot(k-1,j,i)            &
                -2.0*rchpot(k,j,i)
                if(curvature.lt.0.0) then
                   isave = k
                   flag = 1
                endif
             endif
          enddo
       endif
       !  If isave is still equal to zero then L1 point at x = 0.  If
       !  isave is not equal to zero flag tells which side of x=0 the
       !  L1 point resides at, flag = 0 means -ve x, flag = 1 means +ve x
       if(isave.eq.0) then
          rpotcrit = 0.5*(rchpot(1,2,1) + rchpot(2,2,1))
          xcrit = 0.0
       else
          if(flag.eq.0) then
             rpotcrit = rchpot(isave,2,phic)
             xcrit = - rhf(isave)
          else
             rpotcrit = rchpot(isave,2,1)
             xcrit = rhf(isave)
          endif
       endif 

       ! find the L2 and L3 points along the line of centers if
       ! they are on the computational grid
       louter1 = numr
       louter2 = numr
       j = 2
       i = phia
       do k = rm1, numr - 2
          rchtest = (rchpot(k+1,j,i) - rchpot(k,j,i))*                  &
                   (rchpot(k,j,i) - rchpot(k-1,j,i))
          if( rchtest < 0.0 ) then
             curvature = rchpot(k+1,j,i) + rchpot(k-1,j,i)              &
                        - 2.0 * rchpot(k,j,i)
             if( curvature < 0.0 ) then
                louter1 = k
                exit
             endif 
          endif 
       enddo 
       i = phic
       do k = rm2, numr-2
          rchtest = (rchpot(k+1,j,i) - rchpot(k,j,i))*                  &
                   (rchpot(k,j,i) - rchpot(k-1,j,i))
          if( rchtest < 0.0 ) then
             curvature = rchpot(k+1,j,i) + rchpot(k-1,j,i)              &
                        -2.0*rchpot(k,j,i)
             if( curvature < 0.0 ) then
                louter2 = k
                exit
             endif	 
          endif	
       enddo 

       ! find the outer edge of each Roche lobe along the line of centers
       ! if it is on the computationsl grid.
       rochemax1 = numr
       rochemax2 = numr
       i = phia
       j = 2
       do k = rm1, numr-2
          if( (rchpot(k,j,i) <= rpotcrit) .and.                        &
             (rpotcrit <= rchpot(k+1,j,i)) ) then
             rochemax1 = k+1
             exit
          endif 
       enddo 
       i = phic
       do k = rm2, numr-2
          if( (rchpot(k,j,i) <= rpotcrit) .and.                        &
             (rpotcrit <= rchpot(k+1,j,i)) ) then
             rochemax2 = k+1
             exit
          endif 
       enddo 

       ! now add the volume of all the cells where the local potantial
       ! is less than the potential at the L1 point.
       volr1 = 0.0
       volr2 = 0.0
       if( xcrit >= 0.0 ) then
          do i = 1, phi1+1
             do j = 2, numz-1
                do k = 2, numr-1 
                   if( rchpot(k,j,i) <= rpotcrit ) then
                      if( rhf(k)*cosine(i) < xcrit ) then
                         volr2 = volr2 + rhf(k)
                      else
                         volr1 = volr1 + rhf(k)
                      endif 
                   endif 
                enddo
             enddo 
          enddo 
          do i = phi4, numphi
             do j = 2, numz-1
                do k = 2, numr-1
                   if( rchpot(k,j,i) <= rpotcrit ) then
                      if( rhf(k)*cosine(i) < xcrit ) then
                         volr2 = volr2 + rhf(k)
                      else
                         volr1 = volr1 + rhf(k)
                      endif 
                   endif 
                enddo 
             enddo 
          enddo 
          do i = phi2, phi3+1
             do j = 2, numz-1
                do k = 2, numr-1
                   if( rchpot(k,j,i) <= rpotcrit ) then
                      volr2 = volr2 + rhf(k)
                   endif 
                enddo
             enddo 
          enddo 
       else
          do i = 1, phi1+1
             do j = 2, numz-1
                do k = 2, numr-1
                   if( rchpot(k,j,i) <= rpotcrit ) then
                      volr1 = volr1 + rhf(k)
                   endif 
                enddo 
             enddo 
          enddo
          do i = phi4, numphi
             do j = 2, numz-1
                do k = 2, numr-1
                   if( rchpot(k,j,i) <= rpotcrit  ) then
                      volr1 = volr1 + rhf(k)
                   endif 
                enddo
             enddo
          enddo
          do i = phi2, phi3+1
             do j = 2, numz-1
                do k = 2, numr-1
                   if( rchpot(k,j,i) <= rpotcrit ) then
                      if( rhf(k)*cosine(i) > xcrit ) then
                         volr1 = volr1 + rhf(k)
                      else
                         volr2 = volr2 + rhf(k)
                      endif 
                   endif
                enddo
             enddo
          enddo
       endif 
       volr1 = factor*volr1
       volr2 = factor*volr2
       reffr1 = (0.75*volr1/pi)**0.3333333
       reffr2 = (0.75*volr2/pi)**0.3333333

       call cpu_time(ftime)

       i = phic
       j = zc
       do k = rm2, numr
          if( (rho(k,j,i) > epsilon) .and.                            &
             (rho(k+1,j,i) < epsilon) ) then
             star2maxr = k
          endif
       enddo 

       ! write the verbose output file 
       write(11,*) 'Model Number: ',model_num
       write(11,*)
       write(11,*) 'For Star 1:'
       write(11,*) 'Mass 1: ',mass1(qfinal)
       write(11,*) '(x,y): ',xavg1, yavg1
       write(11,*) 'Maximum density: ',rhom1
       write(11,*) 'Polytropic Constant: ',kappa1
       write(11,*) 'Virial Preasure:',s1
       write(11,*) 'Potential Energy:',w1
       write(11,*) 'Kinetic Energy:',t1
       write(11,*) 'Virial Error:',virialerr1
       write(11,*) 'Preasure and Enthalpy Maximums:',pm1, hm1(qfinal)
       write(11,*) 'Maximum at (r,z,phi):',rhf(rm1),zhf(zm1),phi(phim1)
       write(11,*) 'Maximum as integers: ',rm1,zm1,phim1
       write(11,*) 'Inner boundary point:',rhf(rb), zhf(zb), phi(phib)
       write(11,*) 'As integers: ',rb, zb, phib
       write(11,*) 'Outer boundary point:',rhf(ra), zhf(za), phi(phia)
       write(11,*) 'As integers: ', ra, za, phia
       write(11,*) 'Volume: ',vol1
       write(11,*) 'Effective radius: ',reff1
       write(11,*) 'Roche Volume: ',volr1
       write(11,*) 'Effective Radius of Roche lobe: ',reffr1
       write(11,*) 'Angular Momentum: ',j1
       write(11,*) 'Internal Energy: ',e1
       write(11,*) 'Total Energy: ',en1
       write(11,*) 'Outer Lagrange Point1 ',louter1, rhf(louter1)
       write(11,*) 'Outer edge of Roche lobe1: ',rochemax1,rhf(rochemax1)
       if( rho(ra+1,za,phia) > 0.0 ) then
          write(11,*) '*1 extends beyond outer boundary point'
       endif
       write(11,*)
       write(11,*) 'For Star 2:'
       write(11,*) 'Mass:',mass2(qfinal)
       write(11,*) '(x,y): ',xavg2, yavg2
       write(11,*) 'Maximum density: ',rhom2
       write(11,*) 'Polytropic Constant: ',kappa2
       write(11,*) 'Virial Preasure:',s2
       write(11,*) 'Potantial Energy:',w2
       write(11,*) 'Kinetic Energy:',t2
       write(11,*) 'Virial Error:',virialerr2
       write(11,*) 'Preasure and Enthalpy Maximums:',pm2,hm2(qfinal)
       write(11,*) 'Maximum at (r,z,phi):',rhf(rm2),zhf(zm2),phi(phim2)
       write(11,*) 'Maximum as integers: ',rm2,zm2,phim2
       write(11,*) 'Inner boundary point:',rhf(rc),zhf(zc),phi(phic)
       write(11,*) 'As integers: ', rc, zc, phic
       write(11,*) 'Volume: ',vol2
       write(11,*) 'Effective radius: ',reff2
       write(11,*) 'Roche Volume: ',volr2
       write(11,*) 'Effective Radius of Roche lobe: ',reffr2
       write(11,*) 'Density at outer edge',rho(rmax,2,phic)
       write(11,*) 'Star 2 outer extent: ',star2maxr
       write(11,*) 'Angular Momentum: ',j2
       write(11,*) 'Internal Energy: ',e2
       write(11,*) 'Total Energy: ',en2
       write(11,*) 'Outer Lagrange Point2: ',louter2, rhf(louter2)
       write(11,*) 'Outer edge of Roche lobe2: ',rochemax2,rhf(rochemax2)
       write(11,*)
       write(11,*) 'Mass Ratio', mass1(qfinal)/mass2(qfinal)
       write(11,*) 'Primary is star: ', primary
       write(11,*) 'Virial Error:', virialerr
       write(11,*) 'Center of mass:',com
       write(11,*) 'Seperation: ', separation
       write(11,*) 'Angular frequncy: ', sqrt(omsq(qfinal))
       write(11,*) 'Period: ',period
       write(11,*) 'Keplers 3rd constant: ',kepler
       write(11,*) 'Integration constant 1: ', c1(qfinal)
       write(11,*) 'Integration constant 2: ', c2(qfinal)
       write(11,*) 'Total Angular Momentum: ',jtot
       write(11,*) 'Total Energy: ',entot
       write(11,*) 
       write(11,*) 'Roche potential at L1: ',rpotcrit
       write(11,*) 'x coord of L1: ', xcrit
       write(11,*) 'Max value of Roche potential: ',rchmax
       write(11,*) 'Min value of Roche potential: ',rchmin
       write(11,*) 'Convergence criterion: ',eps
       write(11,*) 'Number of iterations: ', qfinal
       write(11,*) 'Total Potential @ a: ',rchpot(ra,za,phia)
       write(11,*) 'Total Potential @ b: ',rchpot(rb,zb,phib)
       write(11,*) 'Total Potential @ c: ',rchpot(rc,zc,phic)
       if (rochemax1.eq.numr) then
           write(11,*) 'Roche Lobe of *1 off grid'
       endif
       if (rochemax2.eq.numr) then
           write(11,*) 'Roche Lobe of *2 off grid'
       endif
       write(11,*)
       write(11,*) 'Execution Time = ',(ftime-stime)/60.0,' minutes'
       write(11,*)
       write(11,*) '==========================================================================================='
       write(11,*)


       write(12,*) model_num,ra,rb,rc,rhom1,rhom2,kappa1,kappa2,           &
                  mass1(qfinal),mass2(qfinal),mass1(qfinal)/mass2(qfinal), &
                  vol1/volr1,vol2/volr2,virialerr,com,xcrit,qfinal

       write(dens_file,'(a,i5)') trim(dens_template), model_num
       open(unit=13,file=trim(dens_file),form='unformatted',status='unknown')
       write(13) rho(:,2,:)
       close(13)

!FOR WRITING OUT A MODEL TO BE RUN BY THE HYDRODYNAMICS CODE
       open(unit=14,file='roche_pot',form='unformatted',status='unknown')
       write(14) rchpot
       close(14)

       open(unit=15,file='density',form='unformatted',status='unknown')
       write(15) rho
       close(15) 

       
       open(unit=10,file="star1")
         do j=1,numz
           do i=1,numr  
             write(10,*) i,j,rho(i,j,1) 
           enddo
           write(10,*)
         enddo
       close(10)       
       
       open(unit=10,file="pot1")
         do j=1,numz
           do i=1,numr  
             write(10,*) i,j,pot(i,j,1) 
           enddo
           write(10,*)
         enddo
       close(10)  
       
       
       
       open(unit=10,file="star2")
         do j=1,numz
           do i=1,numr  
             write(10,*) i,j,rho(i,j,256) 
           enddo
           write(10,*)
         enddo
       close(10) 
       
       open(unit=10,file="pot2")
         do j=1,numz
           do i=1,numr  
             write(10,*) i,j,pot(i,j,256) 
           enddo
           write(10,*)
         enddo
       close(10)
       
       end subroutine binary_scf

       subroutine bin_sum(temp, ret1, ret2)
!!Finds summation of density for each component of the binary
       implicit none
       include 'startup.h'
       real, intent(in), dimension(numr,numz,numphi) :: temp
       real, intent(out) :: ret1, ret2
       integer :: i, j, k
       ret1 = 0.0
       ret2 = 0.0
       do i = phi2, phi3
          do j = 2, numz
             do k = 2, numr
                ret2 = ret2 + temp(k,j,i)
             enddo
          enddo
       enddo
       do i = phi4, numphi
          do j = 2, numz
             do k = 2, numr
                ret1 = ret1 + temp(k,j,i)
             enddo
          enddo
       enddo
       do i = 1, phi1
          do j = 2, numz
             do k = 2, numr
                ret1 = ret1 + temp(k,j,i)
             enddo
          enddo
       enddo
       return 
       end subroutine bin_sum
