!=============================================================================
! geom_module.f90
! Geometry intersection tests:
!   - Ray vs Axis-Aligned Bounding Box (AABB)
!   - Ray vs Finite Cylinder (with caps, arbitrary orientation)
!=============================================================================
module geom_module
  implicit none
  private
  public :: aabb_t, cyl_t
  public :: ray_hits_aabb, ray_hits_cylinder

  real(8), parameter :: GEOM_EPS = 1.0d-14

  !---------------------------------------------------------------------------
  ! Rectangular box (axis-aligned)
  !---------------------------------------------------------------------------
  type :: aabb_t
    real(8) :: xmin, xmax
    real(8) :: ymin, ymax
    real(8) :: zmin, zmax
    real(8) :: margin = 0.0d0   ! safety inflation [cm]
  end type aabb_t

  !---------------------------------------------------------------------------
  ! Finite capped cylinder (arbitrary axis A->B)
  !---------------------------------------------------------------------------
  type :: cyl_t
    real(8) :: ax, ay, az       ! axis bottom point A [cm]
    real(8) :: bx, by, bz       ! axis top point    B [cm]
    real(8) :: r                ! radius [cm]
    real(8) :: margin = 0.0d0   ! safety inflation [cm]
    logical :: caps = .true.    ! include end caps
  end type cyl_t

contains

  !==========================================================================
  ! Ray-AABB intersection (slab method)
  ! Returns hit=T if ray r(t)=o+t*d (t>=0) intersects the inflated box.
  !==========================================================================
  subroutine ray_hits_aabb(ox,oy,oz, dx,dy,dz, box, hit, t_enter, t_exit)
    real(8), intent(in)  :: ox,oy,oz, dx,dy,dz
    type(aabb_t), intent(in)  :: box
    logical, intent(out) :: hit
    real(8), intent(out) :: t_enter, t_exit

    real(8) :: xlo,xhi, ylo,yhi, zlo,zhi
    real(8) :: tmin, tmax, t1, t2
    real(8) :: mg

    mg  = box%margin
    xlo = box%xmin - mg;  xhi = box%xmax + mg
    ylo = box%ymin - mg;  yhi = box%ymax + mg
    zlo = box%zmin - mg;  zhi = box%zmax + mg

    tmin = -huge(1.0d0)
    tmax =  huge(1.0d0)

    ! X slab
    if (abs(dx) < GEOM_EPS) then
      if (ox < xlo .or. ox > xhi) then
        hit=.false.; t_enter=0d0; t_exit=0d0; return
      end if
    else
      t1 = (xlo-ox)/dx;  t2 = (xhi-ox)/dx
      tmin = max(tmin, min(t1,t2))
      tmax = min(tmax, max(t1,t2))
    end if

    ! Y slab
    if (abs(dy) < GEOM_EPS) then
      if (oy < ylo .or. oy > yhi) then
        hit=.false.; t_enter=0d0; t_exit=0d0; return
      end if
    else
      t1 = (ylo-oy)/dy;  t2 = (yhi-oy)/dy
      tmin = max(tmin, min(t1,t2))
      tmax = min(tmax, max(t1,t2))
    end if

    ! Z slab
    if (abs(dz) < GEOM_EPS) then
      if (oz < zlo .or. oz > zhi) then
        hit=.false.; t_enter=0d0; t_exit=0d0; return
      end if
    else
      t1 = (zlo-oz)/dz;  t2 = (zhi-oz)/dz
      tmin = max(tmin, min(t1,t2))
      tmax = min(tmax, max(t1,t2))
    end if

    if (tmax >= tmin .and. tmax >= 0d0) then
      hit    = .true.
      t_enter = max(tmin, 0d0)
      t_exit  = tmax
    else
      hit    = .false.
      t_enter = 0d0
      t_exit  = 0d0
    end if
  end subroutine ray_hits_aabb


  !==========================================================================
  ! Ray-finite-cylinder intersection
  ! Tests wall (quadratic) + both caps (ray-plane + circle check).
  ! Works for any axis orientation A->B.
  !==========================================================================
  subroutine ray_hits_cylinder(ox,oy,oz, dx,dy,dz, cyl, hit, t_hit)
    real(8), intent(in)  :: ox,oy,oz, dx,dy,dz
    type(cyl_t), intent(in)  :: cyl
    logical, intent(out) :: hit
    real(8), intent(out) :: t_hit

    real(8) :: vx,vy,vz, h, invh
    real(8) :: ocx,ocy,ocz
    real(8) :: dpar, ocpar
    real(8) :: dpx,dpy,dpz, ocpx,ocpy,ocpz
    real(8) :: a,b,c, disc, s1, s2, t1, t2
    real(8) :: tbest, tcap
    real(8) :: sx,sy,sz, sproj, rr2, reff
    logical :: ok

    reff = cyl%r + cyl%margin           ! inflated radius

    ! Axis unit vector
    vx = cyl%bx-cyl%ax;  vy = cyl%by-cyl%ay;  vz = cyl%bz-cyl%az
    h  = sqrt(vx*vx + vy*vy + vz*vz)
    if (h < GEOM_EPS) then
      hit=.false.; t_hit=0d0; return
    end if
    invh = 1d0/h
    vx = vx*invh;  vy = vy*invh;  vz = vz*invh

    ! Vector from A to ray origin
    ocx = ox-cyl%ax;  ocy = oy-cyl%ay;  ocz = oz-cyl%az

    ! Projections onto axis
    dpar  = dx*vx  + dy*vy  + dz*vz
    ocpar = ocx*vx + ocy*vy + ocz*vz

    ! Perpendicular components
    dpx  = dx  - dpar*vx;   dpy  = dy  - dpar*vy;   dpz  = dz  - dpar*vz
    ocpx = ocx - ocpar*vx;  ocpy = ocy - ocpar*vy;  ocpz = ocz - ocpar*vz

    tbest = huge(1.0d0)
    ok    = .false.

    !--- Wall: |oc_perp + t*d_perp|^2 = reff^2 ---
    a = dpx*dpx + dpy*dpy + dpz*dpz
    b = 2d0*(ocpx*dpx + ocpy*dpy + ocpz*dpz)
    c = (ocpx*ocpx + ocpy*ocpy + ocpz*ocpz) - reff*reff

    if (a > GEOM_EPS) then
      disc = b*b - 4d0*a*c
      if (disc >= 0d0) then
        disc = sqrt(disc)
        t1 = (-b - disc)/(2d0*a)
        t2 = (-b + disc)/(2d0*a)
        if (t1 >= 0d0) then
          s1 = ocpar + t1*dpar
          if (s1 >= -cyl%margin .and. s1 <= h+cyl%margin) then
            tbest=min(tbest,t1); ok=.true.
          end if
        end if
        if (t2 >= 0d0) then
          s2 = ocpar + t2*dpar
          if (s2 >= -cyl%margin .and. s2 <= h+cyl%margin) then
            tbest=min(tbest,t2); ok=.true.
          end if
        end if
      end if
    end if

    !--- Caps ---
    if (cyl%caps .and. abs(dpar) > GEOM_EPS) then

      ! Bottom cap (s=0): ocpar + t*dpar = 0
      tcap = -ocpar / dpar
      if (tcap >= 0d0) then
        sx = ocx + tcap*dx;  sy = ocy + tcap*dy;  sz = ocz + tcap*dz
        sproj = sx*vx + sy*vy + sz*vz
        rr2 = (sx-sproj*vx)**2 + (sy-sproj*vy)**2 + (sz-sproj*vz)**2
        if (rr2 <= reff*reff) then
          tbest=min(tbest,tcap); ok=.true.
        end if
      end if

      ! Top cap (s=h): ocpar + t*dpar = h
      tcap = (h - ocpar) / dpar
      if (tcap >= 0d0) then
        sx = ocx + tcap*dx;  sy = ocy + tcap*dy;  sz = ocz + tcap*dz
        sproj = sx*vx + sy*vy + sz*vz
        rr2 = (sx-sproj*vx)**2 + (sy-sproj*vy)**2 + (sz-sproj*vz)**2
        if (rr2 <= reff*reff) then
          tbest=min(tbest,tcap); ok=.true.
        end if
      end if

    end if

    hit   = ok
    t_hit = merge(tbest, 0d0, ok)

  end subroutine ray_hits_cylinder

end module geom_module
