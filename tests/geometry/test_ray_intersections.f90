!=============================================================================
! test_ray_intersections.f90 -- ray_hits_cylinder and ray_hits_aabb
! (src/generator/geom_module.f90) against brute force.
!
! The detector filter keeps a muon when its straight line from the source
! crosses the detector inflated by the safety margin. For a cylinder that
! inflated volume is radius r+m over the axial range [-m, h+m] (the wall test
! extends axially by m, the caps sit at s = 0 and s = h with radius r+m, and
! together they catch every ray through that volume); for a box, every face
! moved out by m. Random rays (t >= 0) and random detectors, including tilted
! axes, rays starting inside, and near-axial rays that only a cap can catch,
! are tested against dense sampling of the ray: a point inside the inflated
! volume means a hit. Rays whose closest approach lies within a thin band of
! the surface are ambiguous for the sampler and are skipped (and counted).
!
!   gfortran -O2 src/generator/geom_module.f90 tests/geometry/test_ray_intersections.f90
!=============================================================================
program test_ray_intersections
  use geom_module
  implicit none
  integer, parameter :: NRAY = 200000, NSTEP = 20000
  type(cyl_t)  :: cyl
  type(aabb_t) :: box
  integer :: i, k, n_bad_c, n_bad_b, n_amb_c, n_amb_b, n_hit_c, n_hit_b, n_cap_only
  real(8) :: o(3), d(3), a(3), b(3), u(3), r, m, h, tmax, t, p(3), w(3), s, rr
  real(8) :: dmin, sd, lo(3), hi(3), q(3), t_hit, t1, t2, band
  logical :: hit, bf
  call random_seed(put=[(12345 + k, k = 1, 64)])
  n_bad_c = 0; n_bad_b = 0; n_amb_c = 0; n_amb_b = 0; n_hit_c = 0; n_hit_b = 0
  n_cap_only = 0
  do i = 1, NRAY
    ! --- random cylinder: axis anywhere in a 4 m cube, any orientation -------
    a = rnd3(-200d0, 200d0)
    u = unit(rnd3(-1d0, 1d0))
    if (mod(i, 4) == 0) u = [0d0, 0d0, 1d0]            ! many vertical boreholes
    h = urand(10d0, 300d0);  r = urand(2d0, 80d0);  m = urand(0d0, 40d0)
    if (mod(i, 5) == 0) m = 0d0
    b = a + h * u
    cyl%ax = a(1); cyl%ay = a(2); cyl%az = a(3)
    cyl%bx = b(1); cyl%by = b(2); cyl%bz = b(3)
    cyl%r = r;  cyl%margin = m;  cyl%caps = .true.
    ! --- random ray aimed near the cylinder ----------------------------------
    o = rnd3(-600d0, 600d0)
    if (mod(i, 7) == 0) o = a + urand(0d0, h) * u       ! some start inside
    d = unit(a + urand(-0.2d0, 1.2d0) * h * u + rnd3(-2d0*r, 2d0*r) - o)
    if (mod(i, 9) == 0) d = u * merge(1d0, -1d0, urand(0d0,1d0) < 0.5d0)  ! axial
    call ray_hits_cylinder(o(1),o(2),o(3), d(1),d(2),d(3), cyl, hit, t_hit)
    ! brute force: signed distance to the inflated cylinder along the ray
    tmax = 2000d0;  dmin = huge(1d0)
    do k = 0, NSTEP
      t = tmax * dble(k) / dble(NSTEP)
      p = o + t * d;  w = p - a;  s = dot_product(w, u)
      rr = sqrt(max(0d0, dot_product(w, w) - s*s))
      sd = max(rr - (r + m), -m - s, s - (h + m))       ! < 0 inside
      dmin = min(dmin, sd)
    end do
    band = 0.15d0                                          ! step 0.1 cm
    if (abs(dmin) < band) then
      n_amb_c = n_amb_c + 1
    else
      bf = dmin < 0d0
      if (bf) n_hit_c = n_hit_c + 1
      if (bf .neqv. hit) n_bad_c = n_bad_c + 1
      if (bf .and. mod(i, 9) == 0) n_cap_only = n_cap_only + 1
    end if

    ! --- random box ----------------------------------------------------------
    lo = rnd3(-200d0, 100d0);  hi = lo + rnd3(5d0, 300d0);  m = urand(0d0, 30d0)
    box%xmin = lo(1); box%xmax = hi(1); box%ymin = lo(2); box%ymax = hi(2)
    box%zmin = lo(3); box%zmax = hi(3); box%margin = m
    o = rnd3(-600d0, 600d0)
    q = lo + rnd3(-0.2d0, 1.2d0) * (hi - lo)
    d = unit(q - o)
    if (mod(i, 11) == 0) d = [0d0, 0d0, -1d0]
    call ray_hits_aabb(o(1),o(2),o(3), d(1),d(2),d(3), box, hit, t1, t2)
    dmin = huge(1d0)
    do k = 0, NSTEP
      t = tmax * dble(k) / dble(NSTEP)
      p = o + t * d
      sd = maxval(max(lo - m - p, p - hi - m))
      dmin = min(dmin, sd)
    end do
    if (abs(dmin) < band) then
      n_amb_b = n_amb_b + 1
    else
      bf = dmin < 0d0
      if (bf) n_hit_b = n_hit_b + 1
      if (bf .neqv. hit) n_bad_b = n_bad_b + 1
    end if
  end do
  print '(A,I7,A,I7,A,I6,A,I5,A)', ' cylinder: ', NRAY - n_amb_c, ' rays checked (', n_hit_c, &
        ' hits, ', n_cap_only, ' axial), ', n_bad_c, ' mismatches'
  print '(A,I7,A,I7,A,I5,A)', ' box:      ', NRAY - n_amb_b, ' rays checked (', n_hit_b, &
        ' hits), ', n_bad_b, ' mismatches'
  print '(A,I6,A,I6,A)', ' skipped as tangent (|distance| < 1.5 mm): ', n_amb_c, ' cyl, ', n_amb_b, ' box'
  if (n_bad_c + n_bad_b == 0) then
    print '(A)', ' ALL AGREE'
  else
    print '(A)', ' FAILED'
    error stop 1
  end if
contains
  real(8) function urand(x0, x1)
    real(8), intent(in) :: x0, x1
    real(8) :: v
    call random_number(v);  urand = x0 + (x1 - x0) * v
  end function urand
  function rnd3(x0, x1) result(v)
    real(8), intent(in) :: x0, x1
    real(8) :: v(3)
    call random_number(v);  v = x0 + (x1 - x0) * v
  end function rnd3
  function unit(v) result(uu)
    real(8), intent(in) :: v(3)
    real(8) :: uu(3)
    uu = v / sqrt(dot_product(v, v))
  end function unit
end program test_ray_intersections
