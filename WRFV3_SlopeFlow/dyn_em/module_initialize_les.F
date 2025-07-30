! sed -e "s/grid%mu/gridmu/g" -e "s/grid%Mu/gridMu/g" module_initialize_hill2d_x.F | cpp -DHYBRID_COORD | sed -e "s/gridmu/grid%mu/g" -e "s/gridMu/grid%Mu/g" >> module_initialize_hill2d_x.next
#if ( HYBRID_COORD==1 )
#  define gridmu_1(...) (grid%c1h(k)*XXPC1HXX(__VA_ARGS__))
#  define XXPC1HXX(...) grid%mu_1(__VA_ARGS__)

#  define gridMu_1(...) (grid%c1f(k)*XXPC1FXX(__VA_ARGS__))
#  define XXPC1FXX(...) grid%Mu_1(__VA_ARGS__)

#  define gridmub(...) (grid%c1h(k)*XXPCBHXX(__VA_ARGS__)+grid%c2h(k))
#  define XXPCBHXX(...) grid%mub(__VA_ARGS__)

#  define gridMub(...) (grid%c1f(k)*XXPCBFXX(__VA_ARGS__)+grid%c2f(k))
#  define XXPCBFXX(...) grid%Mub(__VA_ARGS__)
#endif
!IDEAL:MODEL_LAYER:INITIALIZATION
!

!  This MODULE holds the routines which are used to perform various initializations
!  for the individual domains.  

!  This MODULE CONTAINS the following routines:

!  initialize_field_test - 1. Set different fields to different constant
!                             values.  This is only a test.  If the correct
!                             domain is not found (based upon the "id")
!                             then a fatal error is issued.               

!-----------------------------------------------------------------------

MODULE module_initialize_ideal

   USE module_domain
   USE module_io_domain
   USE module_state_description
   USE module_model_constants
   USE module_bc
   USE module_timing
   USE module_configure
   USE module_init_utilities
   !DJW 140408 added USE module_soil_pre along with subroutines 
!           monthly_min_max and monthly_interp_to_date. module_soil_pre
!           is necessary since several variables are defined in it that 
!           are needed for monthly_interp_to_date.
   USE module_soil_pre
#ifdef DM_PARALLEL
   USE module_dm
#endif


CONTAINS


!-------------------------------------------------------------------
! this is a wrapper for the solver-specific init_domain routines.
! Also dereferences the grid variables and passes them down as arguments.
! This is crucial, since the lower level routines may do message passing
! and this will get fouled up on machines that insist on passing down
! copies of assumed-shape arrays (by passing down as arguments, the 
! data are treated as assumed-size -- ie. f77 -- arrays and the copying
! business is avoided).  Fie on the F90 designers.  Fie and a pox.

   SUBROUTINE init_domain ( grid )

   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid 
   !  Local data.
   INTEGER :: idum1, idum2

   CALL set_scalar_indices_from_config ( head_grid%id , idum1, idum2 )

     CALL init_domain_rk( grid &
!
#include "actual_new_args.inc"
!
                        )

   END SUBROUTINE init_domain

!-------------------------------------------------------------------

   SUBROUTINE init_domain_rk ( grid &
!
# include "dummy_new_args.inc"
!
)
   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid

# include "dummy_new_decl.inc"

   TYPE (grid_config_rec_type)              :: config_flags

   !  Local data
   INTEGER                             ::                       &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte, &
                                  i, j, k, kk

   ! Local data

   INTEGER, PARAMETER :: nl_max = 1000
   REAL, DIMENSION(nl_max) :: zk, p_in, theta, rho, u, v, qv, pd_in
   INTEGER :: nl_in


   INTEGER :: icm,jcm, ii, im1, jj, jm1, loop, error, fid, nxc, nyc
   REAL    :: u_mean,v_mean, f0, p_surf, p_level, qvf, z_at_v, z_at_u
   REAL    :: z_scale, xrad, yrad, zrad, rad, delt, cof1, cof2
!   REAL, EXTERNAL :: interp_0
   REAL    :: hm
   REAL    :: pi

!  stuff from original initialization that has been dropped from the Registry 
   REAL    :: vnu, xnu, xnus, dinit0, cbh, p0_temp, t0_temp, zd, zt
   REAL    :: qvf1, qvf2, pd_surf, theta_surf
   INTEGER :: it
   real :: thtmp, ptmp, temp(3)

   LOGICAL :: moisture_init
   LOGICAL :: stretch_grid, dry_sounding

  INTEGER :: xs , xe , ys , ye
  REAL :: mtn_ht,x_dist,y_dist,L1,L2,H1,H,tantheta1,tantheta2
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
!  For LES, add randx
   real :: randx
   !DJW 140723 added for ideal_soil code
   INTEGER :: ns


!DJW added for specifying different eta levels for each domain
  INTEGER :: ks, ke, id
  LOGICAL :: vnest !DJW T if using vertical nesting, otherwise F

! for the hybrid coordinate
   REAL    :: B1, B2, B3, B4, B5, sin_arg

#ifdef DM_PARALLEL
#    include "data_calls.inc"
#endif


   SELECT CASE ( model_data_order )
         CASE ( DATA_ORDER_ZXY )
   kds = grid%sd31 ; kde = grid%ed31 ;
   ids = grid%sd32 ; ide = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   kms = grid%sm31 ; kme = grid%em31 ;
   ims = grid%sm32 ; ime = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   kts = grid%sp31 ; kte = grid%ep31 ;   ! note that tile is entire patch
   its = grid%sp32 ; ite = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XYZ )
   ids = grid%sd31 ; ide = grid%ed31 ;
   jds = grid%sd32 ; jde = grid%ed32 ;
   kds = grid%sd33 ; kde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   jms = grid%sm32 ; jme = grid%em32 ;
   kms = grid%sm33 ; kme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   jts = grid%sp32 ; jte = grid%ep32 ;   ! note that tile is entire patch
   kts = grid%sp33 ; kte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XZY )
   ids = grid%sd31 ; ide = grid%ed31 ;
   kds = grid%sd32 ; kde = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   kms = grid%sm32 ; kme = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   kts = grid%sp32 ; kte = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch

   END SELECT


!  stretch_grid = .true.
!  FOR LES, set stretch to false
   stretch_grid = .false.
   delt = 3.
!   z_scale = .50
   z_scale = .40
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = (jde-jds)/2

   CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )

! here we check to see if the boundary conditions are set properly

   CALL boundary_condition_check( config_flags, bdyzone, error, grid%id )

   moisture_init = .true.

    grid%itimestep=0

#ifdef DM_PARALLEL
   CALL wrf_dm_bcast_bytes( icm , IWORDSIZE )
   CALL wrf_dm_bcast_bytes( jcm , IWORDSIZE )
#endif

    CALL nl_set_mminlu(1, '    ')
    CALL nl_set_iswater(1,0)
    CALL nl_set_cen_lat(1,40.)
    CALL nl_set_cen_lon(1,-105.)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)


!  here we initialize data we currently is not initialized 
!  in the input data

    DO j = jts, jte
      DO i = its, ite
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%msfvy(i,j)    = 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.
         grid%e(i,j)        = 0.
!  for LES, include Coriolis force
         grid%f(i,j)        = 1.e-4 

      END DO
   END DO

    DO j = jts, jte
    DO k = kts, kte
      DO i = its, ite
         grid%ww(i,k,j)     = 0.
      END DO
   END DO
   END DO

   grid%step_number = 0

! set up the grid

    !DJW Added code for specifying multiple domains' eta_levels.
    !First check to make sure that we've not specified more
    !eta_levels than the dimensionality of eta_levels can handle! This
    !issue will most likely cause a break sometime before we real this
    !check, however it doesn't hurt to include it. To increase max_eta,
    !go to frame/module_driver_constants.F.
   vnest = .FALSE.
   DO id=1,model_config_rec%max_dom
     IF (model_config_rec%vert_refine_method(id) .NE. 0) THEN
       vnest = .TRUE.
     ENDIF
   ENDDO

       IF (model_config_rec%eta_levels(1) .EQ. -1) THEN !we do not have eta_levels from namelist      
         !DJW start of original code to set eta levels
   IF (stretch_grid) THEN ! exponential stretch for eta (nearly constant dz)
	    CALL wrf_debug(0, "module_initialize_les: eta_levels is not specified in the namelist, setting levels with stretched spacing in eta.")
     DO k=1, kde
      grid%znw(k) = (exp(-(k-1)/float(kde-1)/z_scale) - exp(-1./z_scale))/ &
                                (1.-exp(-1./z_scale))
     ENDDO
   ELSE
	    CALL wrf_debug(0,"module_initialize_les: eta_levels is not specified in the namelist, setting levels with constant spacing in eta.")
     DO k=1, kde
      grid%znw(k) = 1. - float(k-1)/float(kde-1)
     ENDDO
   ENDIF
       ELSE !we have specified eta levels from the namelist
          CALL wrf_debug(0,"module_initialize_les: vertical nesting is enabled, using eta_levels specified in namelist.input")
          ks = 0
          DO id=1,grid%id
             ks = ks+model_config_rec%e_vert(id)
          ENDDO
          IF (ks .GT. max_eta) THEN
             CALL wrf_error_fatal("too many vertical levels, increase max_eta in frame/module_driver_constants.F")
          ENDIF
	  !Now set the eta_levels to what we specified in the namelist. We've
	  !packed all the domains' eta_levels into a 'vector' and now we need
	  !to pull only the section of the vector associated with our domain
	  !of interest, which is between indicies ks and ke.
	  IF (grid%id .EQ. 1) THEN
            ks = 1
            ke = model_config_rec%e_vert(1)
	  ELSE
            id = 1
            ks = 1
            ke = 0
            DO WHILE (grid%id .GT. id)
              id = id+1
              ks = ks+model_config_rec%e_vert(id-1)
              ke = ks+model_config_rec%e_vert(id)
            ENDDO
	  ENDIF
	  DO k=1,kde
            grid%znw(k) = model_config_rec%eta_levels(ks+k-1)
	  ENDDO
	  !Check the value of the first and last eta level for our domain,
	  !then check that the vector of eta levels is only decreasing
	  IF (grid%znw(1) .NE. 1.0) THEN
            CALL wrf_error_fatal("error with specified eta_levels, first level is not 1.0")
	  ENDIF
	  IF (grid%znw(kde) .NE. 0.0) THEN
            CALL wrf_error_fatal("error with specified eta_levels, last level is not 0.0")
	  ENDIF
	  DO k=2,kde
            IF (grid%znw(k) .GT. grid%znw(k-1)) THEN
              CALL wrf_error_fatal("eta_levels are not uniformly decreasing from 1.0 to 0.0")
            ENDIF
	  ENDDO
       ENDIF

   DO k=1, kde-1
    grid%dnw(k) = grid%znw(k+1) - grid%znw(k)
    grid%rdnw(k) = 1./grid%dnw(k)
    grid%znu(k) = 0.5*(grid%znw(k+1)+grid%znw(k))
   ENDDO
   DO k=2, kde-1
    grid%dn(k) = 0.5*(grid%dnw(k)+grid%dnw(k-1))
    grid%rdn(k) = 1./grid%dn(k)
    grid%fnp(k) = .5* grid%dnw(k  )/grid%dn(k)
    grid%fnm(k) = .5* grid%dnw(k-1)/grid%dn(k)
   ENDDO

   cof1 = (2.*grid%dn(2)+grid%dn(3))/(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(2) 
   cof2 =     grid%dn(2)        /(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(3) 
   grid%cf1  = grid%fnp(2) + cof1
   grid%cf2  = grid%fnm(2) - cof1 - cof2
   grid%cf3  = cof2       

   grid%cfn  = (.5*grid%dnw(kde-1)+grid%dn(kde-1))/grid%dn(kde-1)
   grid%cfn1 = -.5*grid%dnw(kde-1)/grid%dn(kde-1)
   grid%rdx = 1./config_flags%dx
   grid%rdy = 1./config_flags%dy

!  get the sounding from the ascii sounding file, first get dry sounding and 
!  calculate base state

  dry_sounding = .true.
  IF ( wrf_dm_on_monitor() ) THEN
  write(6,*) ' getting dry sounding for base state '

  CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, nl_max, nl_in, theta_surf )
  ENDIF
  CALL wrf_dm_bcast_real( zk , nl_max )
  CALL wrf_dm_bcast_real( p_in , nl_max )
  CALL wrf_dm_bcast_real( pd_in , nl_max )
  CALL wrf_dm_bcast_real( theta , nl_max )
  CALL wrf_dm_bcast_real( rho , nl_max )
  CALL wrf_dm_bcast_real( u , nl_max )
  CALL wrf_dm_bcast_real( v , nl_max )
  CALL wrf_dm_bcast_real( qv , nl_max )
  CALL wrf_dm_bcast_integer ( nl_in , 1 ) 

  write(6,*) ' returned from reading sounding, nl_in is ',nl_in

!  find ptop for the desired ztop (ztop is input from the namelist),
!  and find surface pressure

  grid%p_top = interp_0( p_in, zk, config_flags%ztop, nl_in )

!  Fill in the hybrid coordinate coefficients

   DO k=1, kde
      IF      ( config_flags%hybrid_opt .EQ. 0 ) THEN
         grid%c3f(k) = grid%znw(k)
      ELSE IF ( config_flags%hybrid_opt .EQ. 1 ) THEN
         grid%c3f(k) = grid%znw(k)
      ELSE IF ( config_flags%hybrid_opt .EQ. 2 ) THEN
         B1 = 2. * grid%etac**2 * ( 1. - grid%etac )
         B2 = -grid%etac * ( 4. - 3. * grid%etac - grid%etac**3 )
         B3 = 2. * ( 1. - grid%etac**3 )
         B4 = - ( 1. - grid%etac**2 )
         B5 = (1.-grid%etac)**4
         grid%c3f(k) = ( B1 + B2*grid%znw(k) + B3*grid%znw(k)**2 + B4*grid%znw(k)**3 ) / B5
         IF ( grid%znw(k) .LT. grid%etac ) THEN
            grid%c3f(k) = 0.
         END IF
         IF ( k .EQ. kds ) THEN
            grid%c3f(k) = 1.
         ELSE IF ( k .EQ. kde ) THEN
            grid%c3f(k) = 0.
         END IF
      ELSE IF ( config_flags%hybrid_opt .EQ. 3 ) THEN
         IF ( grid%znw(k) .GE. grid%etac ) THEN
            sin_arg = (1./(1.-grid%etac))*(grid%znw(k)-1.)+1
            grid%c3f(k) = (sin(sin_arg*3.14159265358/2.))**2
         ELSE
            grid%c3f(k) = 0.
         END IF
         IF ( k .EQ. kds ) THEN
            grid%c3f(k) = 1.
         ELSE IF ( k .EQ. kds ) THEN
            grid%c3f(kde) = 0.
         END IF
      ELSE
         CALL wrf_error_fatal ( 'ERROR: --- hybrid_opt=0 ===> Standard WRF Coordinate; hybrid_opt>=1 ===> Hybrid Vertical Coordinate' )
      END IF
   END DO

   DO k=1, kde
      grid%c4f(k) = ( grid%znw(k) - grid%c3f(k) ) * ( p1000mb - grid%p_top )
   ENDDO

   !  Now on half levels, just add up and divide by 2 (for c3h).  Use (eta-c3)*(p00-pt) for c4 on half levels.

   DO k=1, kde-1
      grid%c3h(k) = ( grid%c3f(k+1) + grid%c3f(k) ) * 0.5
      grid%c4h(k) = ( grid%znu(k) - grid%c3h(k) ) * ( p1000mb - grid%p_top )
   ENDDO

   !  c1 = d(B)/d(eta).  We define c1f as c1 on FULL levels.  For a vertical difference,
   !  we need to use B and eta on half levels.  The k-loop ends up referring to the 
   !  full levels, neglecting the top and bottom.

   DO k=kds+1, kde-1
      grid%c1f(k) = ( grid%c3h(k) - grid%c3h(k-1) ) / ( grid%znu(k) - grid%znu(k-1) )
   ENDDO

   !  The boundary conditions to get the coefficients:
   !  1) At k=kts: define d(B)/d(eta) = 1.  This gives us the same value of B and d(B)/d(eta)
   !     when doing the sigma-only B=eta.
   !  2) At k=kte: with the new vertical coordinate, define d(B)/d(eta) = 0.  The curve B SMOOTHLY 
   !     goes to zero, and at the very top, B continues to SMOOTHLY go to zero.  Note that for 
   !     almost all cases of non B=eta, B is ALREADY=ZERO at the top, so this is a reasonable BC to 
   !     assume.
   !  3) At k=kte: when trying to mimic the original vertical coordinate, since B = eta, then
   !     d(B)/d(eta) = 1.

   grid%c1f(kds) = 1.
   IF      ( ( config_flags%hybrid_opt .EQ. 0 ) .OR. ( config_flags%hybrid_opt .EQ. 1 ) ) THEN
      grid%c1f(kde) = 1.
   ELSE
      grid%c1f(kde) = 0.
   END IF

   !  c2 = ( 1. - c1(k) ) * (p00 - pt).  There is no vertical differencing, so we can do the
   !  full kds to kde looping.

   DO k=kds, kde
      grid%c2f(k) = ( 1. - grid%c1f(k) ) * ( p1000mb - grid%p_top )
   END DO

   !  Now on half levels for c1 and c2.  The c1h will result from the full level c3 and full
   !  level eta differences.  The c2 value use the half level c1(k).

   DO k=1, kde-1
      grid%c1h(k) = ( grid%c3f(k+1) - grid%c3f(k) ) / ( grid%znw(k+1) - grid%znw(k) )
      grid%c2h(k) = ( 1. - grid%c1h(k) ) * ( p1000mb - grid%p_top )
   END DO

#if 0
   DO k=1, kde
      grid%c3f(k) = grid%znw(k)
      grid%c4f(k) = 0.
      grid%c3h(k) = grid%znu(k)
      grid%c4h(k) = 0.
      grid%c1f(k) = 1.
      grid%c2f(k) = 0.
      grid%c1h(k) = 1.
      grid%c2h(k) = 0.
   END DO
#endif

!KAL  ------------change ht for trex-------------------------------------
  if (config_flags%ideal_terrain == 1) then
      grid%ht = 0.
  elseif (config_flags%ideal_terrain == 2) then
      DO j=jts,jte-1
      DO i=its,ite-1
         ! calculate x with x=0 as the middle of the domain
         x_dist = (-ide/2.)*config_flags%dx + i*config_flags%dx
         if (ABS(x_dist)<=500.) then
            grid%ht(i,j) = 0.
         elseif ((ABS(x_dist)>500.).AND.(ABS(x_dist)<9500.)) then
            grid%ht(i,j) = 1500.*(0.5 - 0.5*COS((pi*(ABS(x_dist)-500.))/9000.))
         elseif ((ABS(x_dist)>=9500.).AND.(ABS(x_dist)<=10500.)) then
            grid%ht(i,j) = 1500.
         elseif ((ABS(x_dist)>10500.).AND.(ABS(x_dist)<19500.)) then
            grid%ht(i,j) = 1500.*(0.5 +0.5*COS((pi*(ABS(x_dist)-10500.))/9000.))
         elseif (ABS(x_dist)>=19500.) then
            grid%ht(i,j) = 0.
         endif
      ENDDO
      ENDDO
  elseif (config_flags%ideal_terrain == 3) then
      DO j=jts,jte-1
      DO i=its,ite-1
         ! calculate x with x=0 as the middle of the domain
         x_dist = (-ide/2.)*config_flags%dx + i*config_flags%dx
         if (ABS(x_dist)<=500.) then
            grid%ht(i,j) = 0.
         elseif ((ABS(x_dist)>500.).AND.(ABS(x_dist)<9500.)) then
            grid%ht(i,j) = 1500.*(0.5 - 0.5*COS((pi*(ABS(x_dist)-500.))/9000.))
         elseif ((ABS(x_dist)>=9500.).AND.(ABS(x_dist)<=10500.)) then
            grid%ht(i,j) = 1500.
         elseif ((ABS(x_dist)>10500.).AND.(ABS(x_dist)<19500.)) then
            grid%ht(i,j) = 1500.*(0.5 +0.5*COS((pi*(ABS(x_dist)-10500.))/9000.))
         elseif (ABS(x_dist)>=19500.) then
            grid%ht(i,j) = 0.
         endif
         ! calculate the y distance with y=0 as the middle of the domain
         y_dist = (-jde/2.)*config_flags%dy + j*config_flags%dy
         if (ABS(y_dist)<=100000.) then
            grid%ht(i,j) = grid%ht(i,j)
         elseif ((ABS(y_dist)>100000.).AND.(ABS(y_dist)<109000.)) then
            grid%ht(i,j) = grid%ht(i,j)*(0.5 +0.5*COS((pi*(ABS(y_dist)-100000.))/9000.))
         elseif (ABS(y_dist)>=109000.) then
            grid%ht(i,j) = 0.
         endif
      ENDDO
      ENDDO
  elseif (config_flags%ideal_terrain == 4) then
      DO j=jts,MIN(jte,jde-1)
      DO i=its,MIN(ite,ide-1)
         ! calculate x with x=0 as the middle of the domain
         x_dist = ABS((-ide/2.)*config_flags%dx + i*config_flags%dx)
         tantheta1=TAN(config_flags%slope1*pi/180.0)
         tantheta2=TAN(config_flags%slope2*pi/180.0)
         L1=config_flags%fetch1
         L2=config_flags%fetch2
         H1=tantheta1*L1
         H=H1+tantheta2*L2
         if (x_dist>=L1+L2) then
            grid%ht(i,j) = 0.
         elseif (x_dist>L2) then
            grid%ht(i,j) = H1-tantheta1*(x_dist-L2)
         else
            grid%ht(i,j) = H-tantheta2*x_dist
         endif
      ENDDO
      ENDDO

  endif
 





  xs=ide/2 -3
  xs=ids   -3
  xe=xs + 6
  ys=jde/2 -3
  ye=ys + 6
  mtn_ht = 500
#ifdef MTN
  DO j=max(ys,jds),min(ye,jde-1)
  DO i=max(xs,ids),min(xe,ide-1)
     grid%ht(i,j) = mtn_ht * 0.25 * &
               ( 1. + COS ( 2*pi/(xe-xs) * ( i-xs ) + pi ) ) * &
               ( 1. + COS ( 2*pi/(ye-ys) * ( j-ys ) + pi ) )
  ENDDO
  ENDDO
#endif
#ifdef EW_RIDGE
  DO j=max(ys,jds),min(ye,jde-1)
  DO i=ids,ide
     grid%ht(i,j) = mtn_ht * 0.50 * &
               ( 1. + COS ( 2*pi/(ye-ys) * ( j-ys ) + pi ) )
  ENDDO
  ENDDO
#endif
#ifdef NS_RIDGE
  DO j=jds,jde
  DO i=max(xs,ids),min(xe,ide-1)
     grid%ht(i,j) = mtn_ht * 0.50 * &
               ( 1. + COS ( 2*pi/(xe-xs) * ( i-xs ) + pi ) )
  ENDDO
  ENDDO
#endif
  DO j=jts,jte
  DO i=its,ite
    grid%phb(i,1,j) = g * grid%ht(i,j)
    grid%ph0(i,1,j) = g * grid%ht(i,j)
  ENDDO
  ENDDO

  DO J = jts, jte
  DO I = its, ite

    p_surf = interp_0( p_in, zk, grid%phb(i,1,j)/g, nl_in )
    grid%MUB(i,j) = p_surf-grid%p_top

!  this is dry hydrostatic sounding (base state), so given grid%p (coordinate),
!  interp theta (from interp) and compute 1/rho from eqn. of state

    DO K = 1, kte-1
#if  !( HYBRID_COORD==1 )
      p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
#elif ( HYBRID_COORD==1 )
      p_level = grid%c3h(k)*(p_surf - grid%p_top) + grid%c4h(k) + grid%p_top
#endif
      grid%pb(i,k,j) = p_level
      grid%t_init(i,k,j) = interp_0( theta, p_in, p_level, nl_in ) - t0
      grid%alb(i,k,j) = (r_d/p1000mb)*(grid%t_init(i,k,j)+t0)*(grid%pb(i,k,j)/p1000mb)**cvpm
    ENDDO

!  calc hydrostatic balance (alternatively we could interp the geopotential from the
!  sounding, but this assures that the base state is in exact hydrostatic balance with
!  respect to the model eqns.

    DO kk  = 2,kte
      k=kk - 1
      grid%phb(i,kk,j) = grid%phb(i,kk-1,j) - grid%dnw(kk-1)*grid%mub(i,j)*grid%alb(i,kk-1,j)
    ENDDO

  ENDDO
  ENDDO

  IF ( wrf_dm_on_monitor() ) THEN
    write(6,*) ' ptop is ',grid%p_top
#if  !( HYBRID_COORD==1 )
    write(6,*) ' base state grid%mub(1,1), p_surf is ',grid%mub(1,1),grid%mub(1,1)+grid%p_top
#elif ( HYBRID_COORD==1 )
  write(6,*) ' base state grid%MUB(1,1), p_surf is ',grid%MUB(1,1),grid%c3f(kts)*grid%MUB(1,1)+grid%c4f(kts)+grid%p_top
#endif
  ENDIF

!  calculate full state for each column - this includes moisture.

  write(6,*) ' getting moist sounding for full state '
  dry_sounding = .false.
  CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, nl_max, nl_in, theta_surf )

  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)

!  At this point grid%p_top is already set. find the DRY mass in the column 
!  by interpolating the DRY pressure.  

   pd_surf = interp_0( pd_in, zk, grid%phb(i,1,j)/g, nl_in )

!  compute the perturbation mass and the full mass

    grid%MU_1(i,j) = pd_surf-grid%p_top - grid%MUB(i,j)
    grid%MU_2(i,j) = grid%MU_1(i,j)
    grid%MU0(i,j) = grid%MU_1(i,j) + grid%MUB(i,j)

! given the dry pressure and coordinate system, interp the potential
! temperature and qv

    do k=1,kde-1

#if  !( HYBRID_COORD==1 )
      p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top
#elif ( HYBRID_COORD==1 )
      p_level = grid%c3h(k)*(pd_surf - grid%p_top) + grid%c4h(k) + grid%p_top
#endif

      moist(i,k,j,P_QV) = interp_0( qv, pd_in, p_level, nl_in )
      grid%t_1(i,k,j)          = interp_0( theta, pd_in, p_level, nl_in ) - t0
      grid%t_2(i,k,j)          = grid%t_1(i,k,j)
      

    enddo

!  integrate the hydrostatic equation (from the RHS of the bigstep
!  vertical momentum equation) down from the top to get grid%p.
!  first from the top of the model to the top pressure

    kk = kte-1  ! top level
    k=kk+1

    qvf1 = 0.5*(moist(i,kk,j,P_QV)+moist(i,kk,j,P_QV))
    qvf2 = 1./(1.+qvf1)
    qvf1 = qvf1*qvf2

    grid%p(i,kk,j) = - 0.5*(grid%Mu_1(i,j)+qvf1*grid%Mub(i,j))/grid%rdnw(kk)/qvf2
    qvf = 1. + rvovrd*moist(i,kk,j,P_QV)
    grid%alt(i,kk,j) = (r_d/p1000mb)*(grid%t_1(i,kk,j)+t0)*qvf* &
                (((grid%p(i,kk,j)+grid%pb(i,kk,j))/p1000mb)**cvpm)
    grid%al(i,kk,j) = grid%alt(i,kk,j) - grid%alb(i,kk,j)

!  down the column

    do kk=kte-2,1,-1
      k = kk + 1
      qvf1 = 0.5*(moist(i,kk,j,P_QV)+moist(i,kk+1,j,P_QV))
      qvf2 = 1./(1.+qvf1)
      qvf1 = qvf1*qvf2
      grid%p(i,kk,j) = grid%p(i,kk+1,j) - (grid%Mu_1(i,j) + qvf1*grid%Mub(i,j))/qvf2/grid%rdn(kk+1)
      qvf = 1. + rvovrd*moist(i,kk,j,P_QV)
      grid%alt(i,kk,j) = (r_d/p1000mb)*(grid%t_1(i,kk,j)+t0)*qvf* &
                  (((grid%p(i,kk,j)+grid%pb(i,kk,j))/p1000mb)**cvpm)
      grid%al(i,kk,j) = grid%alt(i,kk,j) - grid%alb(i,kk,j)
    enddo

!  this is the hydrostatic equation used in the model after the
!  small timesteps.  In the model, grid%al (inverse density)
!  is computed from the geopotential.


    grid%ph_1(i,1,j) = 0.
    DO k  = 2,kte
      grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - &
                    grid%dnw(k-1)*((grid%Mub(i,j)+grid%Mu_1(i,j))*grid%al(i,k-1,j)+ &
                    grid%Mu_1(i,j)*grid%alb(i,k-1,j)  )
                                                   
      grid%ph_2(i,k,j) = grid%ph_1(i,k,j) 
      grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
    ENDDO

    IF ( wrf_dm_on_monitor() ) THEN
    if((i==2) .and. (j==2)) then
     write(6,*) ' grid%ph_1 calc ',grid%ph_1(2,1,2),grid%ph_1(2,2,2),&
                              grid%mu_1(2,2)+grid%mub(2,2),grid%mu_1(2,2), &
                              grid%alb(2,1,2),grid%al(1,2,1),grid%rdnw(1)
    endif
    ENDIF

  ENDDO
  ENDDO

!#if 0

!  thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',delt

! For LES, change the initial random perturbations
! For 2D test, call randx outside I-loop
! For 3D runs, call randx inside both I-J loops

  DO J = jts, min(jde-1,jte)
!   yrad = config_flags%dy*float(j-nyc)/10000.
    yrad = 0.
    DO I = its, min(ide-1,ite)
!     xrad = config_flags%dx*float(i-nxc)/10000.
      xrad = 0.
      call random_number (randx)
      randx = randx - 0.5
!     DO K = 1, kte-1
      DO K = 1, 4 

!  No bubbles for LES!
!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

!       zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
!                  +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
!       zrad = (zrad-1500.)/1500.
        zrad = 0.
        RAD=SQRT(xrad*xrad+yrad*yrad+zrad*zrad)
        IF(RAD <= 1.) THEN
!          grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*COS(.5*PI*RAD)**2
           grid%t_1(i,k,j)=grid%t_1(i,k,j)+ 0.1 *randx
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF
      ENDDO

!  rebalance hydrostatically

      DO kk  = 2,kte
        k = kk - 1
        grid%ph_1(i,kk,j) = grid%ph_1(i,kk-1,j) - (grid%dnw(kk-1))*(       &
                     (grid%mub(i,j)+grid%mu_1(i,j))*grid%al(i,kk-1,j)+ &
                      grid%mu_1(i,j)*grid%alb(i,kk-1,j)  )
                                                   
        grid%ph_2(i,kk,j) = grid%ph_1(i,kk,j) 
        grid%ph0(i,kk,j) = grid%ph_1(i,kk,j) + grid%phb(i,kk,j)
      ENDDO

    ENDDO
  ENDDO

!#endif

   IF ( wrf_dm_on_monitor() ) THEN
   write(6,*) ' grid%mu_1 from comp ', grid%mu_1(1,1)
   write(6,*) ' full state sounding from comp, ph, grid%p, grid%al, grid%t_1, qv '
   do k=1,kde-1
     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,1)+grid%phb(1,k,1), &
                                      grid%p(1,k,1)+grid%pb(1,k,1), grid%alt(1,k,1), &
                                      grid%t_1(1,k,1)+t0, moist(1,k,1,P_QV)
   enddo

   write(6,*) ' pert state sounding from comp, grid%ph_1, pp, alp, grid%t_1, qv '
   do k=1,kde-1
     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,1), &
                                      grid%p(1,k,1), grid%al(1,k,1), &
                                      grid%t_1(1,k,1), moist(1,k,1,P_QV)
   enddo
   ENDIF

! interp v

  DO J = jts, jte
  DO I = its, min(ide-1,ite)

    IF (j == jds) THEN
      z_at_v = grid%phb(i,1,j)/g
    ELSE IF (j == jde) THEN
      z_at_v = grid%phb(i,1,j-1)/g
    ELSE
      z_at_v = 0.5*(grid%phb(i,1,j)+grid%phb(i,1,j-1))/g
    END IF
    p_surf = interp_0( p_in, zk, z_at_v, nl_in )

    DO K = 1, kte-1
#if  !( HYBRID_COORD==1 )
      p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
#elif ( HYBRID_COORD==1 )
      p_level = grid%c3h(k)*(p_surf - grid%p_top) + grid%c4h(k) + grid%p_top
#endif
      grid%v_1(i,k,j) = interp_0( v, p_in, p_level, nl_in )
      grid%v_2(i,k,j) = grid%v_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO

! interp u

  DO J = jts, min(jde-1,jte)
  DO I = its, ite

    IF (i == ids) THEN
      z_at_u = grid%phb(i,1,j)/g
    ELSE IF (i == ide) THEN
      z_at_u = grid%phb(i-1,1,j)/g
    ELSE
      z_at_u = 0.5*(grid%phb(i,1,j)+grid%phb(i-1,1,j))/g
    END IF

    p_surf = interp_0( p_in, zk, z_at_u, nl_in )

    DO K = 1, kte-1
#if  !( HYBRID_COORD==1 )
      p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
#elif ( HYBRID_COORD==1 )
      p_level = grid%c3h(k)*(p_surf - grid%p_top) + grid%c4h(k) + grid%p_top
#endif
      grid%u_1(i,k,j) = interp_0( u, p_in, p_level, nl_in )
      grid%u_2(i,k,j) = grid%u_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO

!  set w

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte
  DO I = its, min(ide-1,ite)
    grid%w_1(i,k,j) = 0.
    grid%w_2(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO

!  set a few more things

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte-1
  DO I = its, min(ide-1,ite)
    grid%h_diabatic(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO

  IF ( wrf_dm_on_monitor() ) THEN
  DO k=1,kte-1
    grid%t_base(k) = grid%t_1(1,k,1)
    grid%qv_base(k) = moist(1,k,1,P_QV)
    grid%u_base(k) = grid%u_1(1,k,1)
    grid%v_base(k) = grid%v_1(1,k,1)
    grid%z_base(k) = 0.5*(grid%phb(1,k,1)+grid%phb(1,k+1,1)+grid%ph_1(1,k,1)+grid%ph_1(1,k+1,1))/g
  ENDDO
  ENDIF
  CALL wrf_dm_bcast_real( grid%t_base , kte )
  CALL wrf_dm_bcast_real( grid%qv_base , kte )
  CALL wrf_dm_bcast_real( grid%u_base , kte )
  CALL wrf_dm_bcast_real( grid%v_base , kte )
  CALL wrf_dm_bcast_real( grid%z_base , kte )

  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)
     thtmp   = grid%t_2(i,1,j)+t0
     ptmp    = grid%p(i,1,j)+grid%pb(i,1,j)
     temp(1) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,2,j)+t0
     ptmp    = grid%p(i,2,j)+grid%pb(i,2,j)
     temp(2) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,3,j)+t0
     ptmp    = grid%p(i,3,j)+grid%pb(i,3,j)
     temp(3) = thtmp * (ptmp/p1000mb)**rcp

!    For LES-CBL, add 5 degrees to the surface temperature!
!      - Removed in 3.6
!
!    grid%tsk(I,J)=grid%cf1*temp(1)+grid%cf2*temp(2)+grid%cf3*temp(3)
!    grid%tsk(I,J)=grid%cf1*temp(1)+grid%cf2*temp(2)+grid%cf3*temp(3)+5.
     grid%tsk(I,J)=theta_surf * (p_surf/p1000mb)**rcp
     grid%tmn(I,J)=grid%tsk(I,J)-0.5
  ENDDO
  ENDDO

!DJW 140407 ---set up variables used for LSM---
  IF (config_flags%ideal_soil) THEN
    write(*,'(A)') "DJW: using ideal_soil initialization"
    CALL nl_set_mminlu(1,'USGS')
    DO j=jts,min(jde-1,jte)
    DO i=its,min(ide-1,ite)
      grid%ivgtyp(i,j) = 10 !corresponds to savanna
      grid%lu_index(i,j) = grid%ivgtyp(i,j)
      grid%isltyp(i,j) = 3 !corresponds to sandy loam
      grid%landmask(i,j) = 1.0 !1 for land, 2 for water
      grid%xland(i,j) = 1.0
      grid%greenfrac = 0.10000001 !monthly greenness fraction between 0 and 1
!---2m values---
      grid%q2(i,j)  = interp_0(qv,   zk,grid%ht(i,j),nl_in)
      grid%th2(i,j) = interp_0(theta,zk,grid%ht(i,j),nl_in)
      grid%t2(i,j) = grid%th2(i,j)*(((grid%p(i,1,j)+grid%pb(i,1,j))/p1000mb)**rcp)
      grid%psfc(i,j) = grid%p(i,1,j)+grid%pb(i,1,j)
!---10m values---
      grid%u10(i,j) = interp_0(u,zk,grid%ht(i,j)+10.0,nl_in)
      grid%v10(i,j) = interp_0(v,zk,grid%ht(i,j)+10.0,nl_in)
    ENDDO
    ENDDO

    CALL monthly_interp_to_date ( grid%greenfrac , current_date , grid%vegfra ,&
                                  ids , ide , jds , jde , kds , kde , &
                                  ims , ime , jms , jme , kms , kme , &
                                  its , ite , jts , jte , kts , kte )

    CALL monthly_min_max ( grid%greenfrac , grid%shdmin , grid%shdmax , &
                           ids , ide , jds , jde , kds , kde , &
                           ims , ime , jms , jme , kms , kme , &
                           its , ite , jts , jte , kts , kte )

    DO j=jts,MIN(jte,jde-1)
    DO i=its,MIN(ite,ide-1)
      grid%vegfra(i,j) = grid%vegfra(i,j)*100.0 !vegetation fraction
      grid%shdmax(i,j) = grid%shdmax(i,j)*100.0 !annual max veg fraction
      grid%shdmin(i,j) = grid%shdmin(i,j)*100.0 !annual min veg fraction
      grid%snoalb(i,j) = 60.0 !average from owens valley data
      DO ns=1,config_flags%num_soil_layers
        grid%tslb(i,ns,j) = grid%t2(i,j) !initialize soil temp as equal to atm
        grid%smois(i,ns,j) = 0.0868 !twenty percent saturation rate*porosity for sandy loam
!        grid%smois(i,ns,j) = 0.0 !zero percent saturation rate*porosity for
!        sandy loam
      ENDDO
      grid%sst(i,j) = grid%tsk(i,j) !sst is sea surface temperature, tsk is surface skin temperature
      grid%albbck(i,j) = 0.27 !background albedo
      grid%xlat(i,j) = 0.0
      grid%xlong(i,j) = 0.0
    END DO
    END DO

    DO j=jts,MIN(jte,jde-1)
    DO i=its,ite
      grid%xlat_u(i,j) = 0.0
      grid%xlong_u(i,j) = 0.0
    ENDDO
    ENDDO

    DO j=jts,jte
    DO i=its,MIN(ite,ide-1)
      grid%xlat_v(i,j) = 0.0
      grid%xlong_v(i,j) = 0.0
    ENDDO
    ENDDO

    grid%zs = [0.05, 0.25, 0.70, 1.5] !from real initialization, depths ofcenters of soil layers
    grid%dzs = [0.1, 0.3, 0.6, 1.0] !from real initialization, thicknesses ofsoil layers

  ENDIF !ideal_soil 

 END SUBROUTINE init_domain_rk

   SUBROUTINE init_module_initialize
   END SUBROUTINE init_module_initialize

!DJW 140408 copied subroutine monthly_min_max from module_initialize_real
   SUBROUTINE monthly_min_max ( field_in , field_min , field_max , &
                      ids , ide , jds , jde , kds , kde , &
                      ims , ime , jms , jme , kms , kme , &
                      its , ite , jts , jte , kts , kte )

   !  Plow through each month, find the max, min values for each i,j.

      IMPLICIT NONE

      INTEGER , INTENT(IN)        :: ids , ide , jds , jde , kds , kde , &
                                     ims , ime , jms , jme , kms , kme , &
                                     its , ite , jts , jte , kts , kte

      REAL , DIMENSION(ims:ime,12,jms:jme) , INTENT(IN)  :: field_in
      REAL , DIMENSION(ims:ime,   jms:jme) , INTENT(OUT) :: field_min ,field_max

      !  Local vars

      INTEGER :: i , j , l
      REAL :: minner , maxxer

      DO j = jts , MIN(jde-1,jte)
         DO i = its , MIN(ide-1,ite)
            minner = field_in(i,1,j)
            maxxer = field_in(i,1,j)
            DO l = 2 , 12
               IF ( field_in(i,l,j) .LT. minner ) THEN
                  minner = field_in(i,l,j)
               END IF
               IF ( field_in(i,l,j) .GT. maxxer ) THEN
                  maxxer = field_in(i,l,j)
               END IF
            END DO
            field_min(i,j) = minner
            field_max(i,j) = maxxer
         END DO
      END DO

   END SUBROUTINE monthly_min_max
!DJW end of copy of monthly_min_max

!DJW 140408 copied subroutine monthly_interp_to_date from module_initialize_real
   SUBROUTINE monthly_interp_to_date ( field_in , date_str , field_out , &
                      ids , ide , jds , jde , kds , kde , &
                      ims , ime , jms , jme , kms , kme , &
                      its , ite , jts , jte , kts , kte )

   !  Linrarly in time interpolate data to a current valid time.  The data is
   !  assumed to come in "monthly", valid at the 15th of every month.

      IMPLICIT NONE

      INTEGER , INTENT(IN)        :: ids , ide , jds , jde , kds , kde , &
                                     ims , ime , jms , jme , kms , kme , &
                                     its , ite , jts , jte , kts , kte

      CHARACTER (LEN=24) , INTENT(IN) :: date_str
      REAL , DIMENSION(ims:ime,12,jms:jme) , INTENT(IN)  :: field_in
      REAL , DIMENSION(ims:ime,   jms:jme) , INTENT(OUT) :: field_out

      !  Local vars

      INTEGER :: i , j , l
      INTEGER , DIMENSION(0:13) :: middle
      INTEGER :: target_julyr , target_julday , target_date
      INTEGER :: julyr , julday , int_month , month1 , month2
      REAL :: gmt
      CHARACTER (LEN=4) :: yr
      CHARACTER (LEN=2) :: mon , day15


      WRITE(day15,FMT='(I2.2)') 15
      DO l = 1 , 12
         WRITE(mon,FMT='(I2.2)') l
         CALL get_julgmt(date_str(1:4)//'-'//mon//'-'//day15//'_'//'00:00:00.0000' , julyr , julday ,gmt)
         middle(l) = julyr*1000 + julday
      END DO

      l = 0
      middle(l) = middle( 1) - 31

      l = 13
      middle(l) = middle(12) + 31

      CALL get_julgmt ( date_str , target_julyr , target_julday , gmt )
      target_date = target_julyr * 1000 + target_julday
      find_month : DO l = 0 , 12
         IF ( ( middle(l) .LT. target_date ) .AND. ( middle(l+1) .GE.target_date ) ) THEN
            DO j = jts , MIN ( jde-1 , jte )
               DO i = its , MIN (ide-1 , ite )
                  IF ( skip_middle_points_t ( ids , ide , jds , jde , i , j , em_width , hold_ups ) ) CYCLE
                  int_month = l
                  IF ( ( int_month .EQ. 0 ) .OR. ( int_month .EQ. 12 ) ) THEN
                     month1 = 12
                     month2 =  1
                  ELSE
                     month1 = int_month
                     month2 = month1 + 1
                  END IF
                  field_out(i,j) =  ( field_in(i,month2,j) * ( target_date -middle(l)   ) + &
                                      field_in(i,month1,j) * ( middle(l+1) -target_date ) ) / &
                                    ( middle(l+1) - middle(l) )
               END DO
            END DO
            EXIT find_month
         END IF
      END DO find_month

   END SUBROUTINE monthly_interp_to_date
!DJW end of copy of monthly_interp_to_date


!---------------------------------------------------------------------

!  test driver for get_sounding
!
!      implicit none
!      integer n
!      parameter(n = 1000)
!      real zk(n),p(n),theta(n),rho(n),u(n),v(n),qv(n),pd(n)
!      logical dry
!      integer nl,k
!
!      dry = .false.
!      dry = .true.
!      call get_sounding( zk, p, pd, theta, rho, u, v, qv, dry, n, nl )
!      write(6,*) ' input levels ',nl
!      write(6,*) ' sounding '
!      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
!      do k=1,nl
!        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), pd(k), theta(k), rho(k), u(k), v(k), qv(k)
!      enddo
!      end
!
!---------------------------------------------------------------------------

      subroutine get_sounding( zk, p, p_dry, theta, rho, &
                               u, v, qv, dry, nl_max, nl_in, th_surf )
      implicit none

      integer nl_max, nl_in
      real zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
           u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)
      logical dry

      integer n
      parameter(n=1000)
      logical debug
      parameter( debug = .true.)

! input sounding data

      real p_surf, th_surf, qv_surf
      real pi_surf, pi(n)
      real h_input(n), th_input(n), qv_input(n), u_input(n), v_input(n)

! diagnostics

      real rho_surf, p_input(n), rho_input(n)
      real pm_input(n)  !  this are for full moist sounding

! local data

      real r
      parameter (r = r_d)
      integer k, it, nl
      real qvf, qvf1, dz

!  first, read the sounding

      call read_sounding( p_surf, th_surf, qv_surf, &
                          h_input, th_input, qv_input, u_input, v_input,n, nl, debug )

      if(dry) then
       do k=1,nl
         qv_input(k) = 0.
       enddo
      endif

      if(debug) write(6,*) ' number of input levels = ',nl

        nl_in = nl
        if(nl_in .gt. nl_max ) then
          write(6,*) ' too many levels for input arrays ',nl_in,nl_max
          call wrf_error_fatal ( ' too many levels for input arrays ' )
        end if

!  compute diagnostics,
!  first, convert qv(g/kg) to qv(g/g)

      do k=1,nl
        qv_input(k) = 0.001*qv_input(k)
      enddo

      p_surf = 100.*p_surf  ! convert to pascals
      qvf = 1. + rvovrd*qv_input(1) 
      rho_surf = 1./((r/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
      pi_surf = (p_surf/p1000mb)**(r/cp)

      if(debug) then
        write(6,*) ' surface density is ',rho_surf
        write(6,*) ' surface pi is      ',pi_surf
      end if


!  integrate moist sounding hydrostatically, starting from the
!  specified surface pressure
!  -> first, integrate from surface to lowest level

          qvf = 1. + rvovrd*qv_input(1) 
          qvf1 = 1. + qv_input(1)
          rho_input(1) = rho_surf
          dz = h_input(1)
          do it=1,10
            pm_input(1) = p_surf &
                    - 0.5*dz*(rho_surf+rho_input(1))*g*qvf1
            rho_input(1) = 1./((r/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm))
          enddo

! integrate up the column

          do k=2,nl
            rho_input(k) = rho_input(k-1)
            dz = h_input(k)-h_input(k-1)
            qvf1 = 0.5*(2.+(qv_input(k-1)+qv_input(k)))
            qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here
 
            do it=1,10
              pm_input(k) = pm_input(k-1) &
                      - 0.5*dz*(rho_input(k)+rho_input(k-1))*g*qvf1
              rho_input(k) = 1./((r/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
            enddo
          enddo

!  we have the moist sounding

!  next, compute the dry sounding using p at the highest level from the
!  moist sounding and integrating down.

        p_input(nl) = pm_input(nl)

          do k=nl-1,1,-1
            dz = h_input(k+1)-h_input(k)
            p_input(k) = p_input(k+1) + 0.5*dz*(rho_input(k)+rho_input(k+1))*g
          enddo


        do k=1,nl

          zk(k) = h_input(k)
          p(k) = pm_input(k)
          p_dry(k) = p_input(k)
          theta(k) = th_input(k)
          rho(k) = rho_input(k)
          u(k) = u_input(k)
          v(k) = v_input(k)
          qv(k) = qv_input(k)

        enddo

     if(debug) then
      write(6,*) ' sounding '
      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
      do k=1,nl
        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k)
      enddo

     end if

      end subroutine get_sounding

!-------------------------------------------------------

      subroutine read_sounding( ps,ts,qvs,h,th,qv,u,v,n,nl,debug )
      implicit none
      integer n,nl
      real ps,ts,qvs,h(n),th(n),qv(n),u(n),v(n)
      logical end_of_file
      logical debug

      integer k

      open(unit=10,file='input_sounding',form='formatted',status='old')
      rewind(10)
      read(10,*) ps, ts, qvs
      if(debug) then
        write(6,*) ' input sounding surface parameters '
        write(6,*) ' surface pressure (mb) ',ps
        write(6,*) ' surface pot. temp (K) ',ts
        write(6,*) ' surface mixing ratio (g/kg) ',qvs
      end if

      end_of_file = .false.
      k = 0

      do while (.not. end_of_file)

        read(10,*,end=100) h(k+1), th(k+1), qv(k+1), u(k+1), v(k+1)
        k = k+1
        if(debug) write(6,'(1x,i3,5(1x,e10.3))') k, h(k), th(k), qv(k), u(k), v(k)
        go to 110
 100    end_of_file = .true.
 110    continue
      enddo

      nl = k

      close(unit=10,status = 'keep')

      end subroutine read_sounding

END MODULE module_initialize_ideal
