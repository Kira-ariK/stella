!###############################################################################
!################## IMPLICIT TERMS OF THE GYROKINETIC EQUATION #################
!###############################################################################
! 
! This module adds the implicit gyrokinetic terms to the right-hand-side of the 
! gyrokinetic equation, in order to advance the distribution function in time.
! 
! First advance the distribution function <g> in time using the gyrokinetic equation.
! Next, adance the fields (electrostatic potential <phi>, as well as the electromagnetic
! fields <apar> and <bpar>) in time using the quasi-neutrality condition.
! 
!###############################################################################
module gyrokinetic_equation_implicit

   ! Load debug flags
   use debug_flags, only: debug => time_advance_debug

   implicit none
   
   ! Make routines available to stella.f90
   public :: advance_distribution_function_using_implicit_gyrokinetic_terms

   ! Number of remaps written by dump_remap_comparison (remap_em_ndump)
   integer :: ndump_done = 0

   private

contains

   !****************************************************************************
   !                      IMPLICIT TIME ADVANCE SUBROUTINE                      
   !****************************************************************************
   subroutine advance_distribution_function_using_implicit_gyrokinetic_terms(istep, phi, apar, bpar, g)

      ! Parallelisation
      use mp, only: proc0
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo
      use timers, only: time_gke
      
      ! Fields
      use field_equations, only: advance_fields
      use field_equations, only: fields_updated

      ! Grids
      use grids_z, only: nzgrid

      ! Physics flags
      use parameters_physics, only: include_parallel_streaming
      use parameters_physics, only: include_mirror
      use gk_flow_shear, only: prp_shear_enabled
      
      ! Numerical flags
      use parameters_numerical, only: stream_implicit
      use parameters_numerical, only: mirror_implicit
      use parameters_numerical, only: flip_flop
      
      ! Radial variation
      use parameters_multibox, only: rk_step
      use gk_radial_variation, only: mb_communicate

      ! Collisions
      use dissipation_and_collisions, only: hyper_dissipation
      use dissipation_and_collisions, only: include_collisions
      use dissipation_and_collisions, only: collisions_implicit
      use dissipation_and_collisions, only: advance_collisions_implicit
      
      ! Advance implicit terms of the gyrokinetic equation
      use dissipation_hyper, only: advance_hyper_dissipation
      use gk_implicit_terms, only: advance_implicit_terms
      use gk_mirror, only: advance_mirror_implicit
      use gk_flow_shear, only: advance_perp_flow_shear

      implicit none

      ! Arguments
      integer, intent(in) :: istep
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      
      !-------------------------------------------------------------------------

      ! Start the timer for the implicit part of the solve
      if (proc0) call time_message(.false., time_gke(:, 9), ' implicit')

      ! Reverse the order of operations every time step
      ! as part of alternating direction operator splitting
      ! this is needed to ensure 2nd order accuracy in time
      ! if (mod(istep,2)==0) then
      ! g^{*} (coming from explicit solve) is input
      ! get g^{**}, with g^{**}-g^{*} due to mirror term

      if (rk_step) call mb_communicate(g)

      if (mod(istep, 2) == 1 .or. .not. flip_flop) then

         if (prp_shear_enabled) then
            call remap_flow_shear(g, phi, apar, bpar)
            fields_updated = .false.
         end if

         if (hyper_dissipation) then
            call advance_hyper_dissipation(g)
            fields_updated = .false.
         end if

         if (collisions_implicit .and. include_collisions) then
            call advance_fields(g, phi, apar, bpar, dist='g')
            call advance_collisions_implicit(mirror_implicit, phi, apar, bpar, g)
            fields_updated = .false.
         end if

         if (mirror_implicit .and. include_mirror) then
            call advance_mirror_implicit(collisions_implicit, g, apar)
            fields_updated = .false.
         end if

         ! If the distribution function has been updated due to the addition
         ! of any implicit terms of the gyrokinetic equation, we need to use the 
         ! quasi-neutrality condition to update the fields <phi>, <apar> and <bpar>
         call advance_fields(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 
         
         ! g^{**} is input
         ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
         if (stream_implicit .and. include_parallel_streaming) then
            call advance_implicit_terms(g, phi, apar, bpar)
            fields_updated = .false.
         end if

         ! Update the fields if not already updated
         call advance_fields(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 
         
      else

         ! Get updated fields corresponding to advanced g
         ! note that hyper-dissipation and mirror advances
         ! depended only on g and so did not need field update
         call advance_fields(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 

         ! g^{**} is input
         ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
         if (stream_implicit .and. include_parallel_streaming) then
            call advance_implicit_terms(g, phi, apar, bpar)
            fields_updated = .false.            
         end if

         if (mirror_implicit .and. include_mirror) then
            call advance_mirror_implicit(collisions_implicit, g, apar)
            fields_updated = .false.
         end if

         if (collisions_implicit .and. include_collisions) then
            call advance_fields(g, phi, apar, bpar, dist='g')
            call advance_collisions_implicit(mirror_implicit, phi, apar, bpar, g)
            fields_updated = .false.
         end if

         if (hyper_dissipation) then
            call advance_hyper_dissipation(g)
            fields_updated = .false.
         end if

         if (prp_shear_enabled) then
            call remap_flow_shear(g, phi, apar, bpar)
            fields_updated = .false.
         end if

      end if

      ! Stop the timer for the implict part of the solve
      if (proc0) call time_message(.false., time_gke(:, 9), ' implicit')

   end subroutine advance_distribution_function_using_implicit_gyrokinetic_terms

   !****************************************************************************
   !                  Hammett remap in electromagnetic runs (remap_em)
   !****************************************************************************
   ! The Hammett remap moves each k_y by one k_x cell once its drifted k_x has moved half a cell.  With c = (Z F0/T)
   ! 2 v_par v_th and gbar = g + c J0 A_par:
   !   remap_em = 0 (original): g -> shift(g); A_par = j/k_perp^2 with the new k_perp, so A_par jumps by
   !     k_perp,old^2/k_perp,new^2 (~2.6 at the smallest k_y and k_x = 0): an impulsive inductive E_par.
   !   remap_em = 1: shift gbar built at the old cell; A_par = j/(k_perp^2 + skin term) stays nearly continuous,
   !     but g keeps c [J0(k_old) A_par,old - J0(k_new) A_par,new], mixing the J0 of two cells.
   !   remap_em = 2: shift g, carry A_par with its mode, build gbar at the new cell from the carried A_par, solve
   !     A_par from gbar, rebuild g: g -> shift(g) + c J0(k_new) [A_par,carried - A_par,new].  A_par nearly
   !     continuous, one J0 only.
   ! advance_fields returns at once if fields_updated is set, so the flag is cleared before each solve.
   ! Only the k_y shifted on this call go through g -> gbar -> g; every other k_y (including the zonal ones and
   ! k_x = k_y = 0, where A_par from g is 0 but A_par from gbar is j/skin) is restored unchanged afterwards, so the
   ! options act on the shifted modes only.  The field solve is local in (k_y, k_x), so this is exact.
   ! Shifted k_y outside [remap_em_kymin, remap_em_kymax] get the plain shift of g (remap_em = 0): the kept copy is
   ! shifted with shift_field_hammett (same shift, shift_state untouched) before the remap and restored afterwards.
   !****************************************************************************
   subroutine remap_flow_shear(g, phi, apar, bpar)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid
      use parameters_physics, only: include_apar
      use field_equations, only: advance_fields, fields_updated
      use calculations_tofrom_ghf, only: gbar_to_g
      use gk_flow_shear, only: advance_perp_flow_shear, perp_shift_pending, remap_em, shift_field_hammett
      use gk_flow_shear, only: perp_shift_mask, remap_em_kymin, remap_em_kymax, remap_em_ndump
      use grids_kxky, only: naky, aky

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
      complex, dimension(:, :, :, :), allocatable :: apar_s
      complex, dimension(:, :, :, :, :), allocatable :: gkeep
      logical, dimension(:), allocatable :: shifted
      integer :: iky, ivmu

      if (remap_em == 0 .or. .not. include_apar .or. .not. perp_shift_pending()) then
         call advance_perp_flow_shear(g)
         return
      end if

      ! Diagnostic: compare A_par from the original and the gbar remap of this same state (leaves g unchanged)
      if (ndump_done < remap_em_ndump) call dump_remap_comparison(g, phi, apar, bpar)

      ! Remember which k_y move this call, and keep the others untouched
      allocate (shifted(naky))
      call perp_shift_mask(shifted)
      allocate (gkeep, source=g)
      ! Shifted k_y outside the remap_em range: the kept copy gets the plain shift of g (before shift_state resets)
      if (any(shifted .and. (aky < remap_em_kymin .or. aky > remap_em_kymax))) then
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            call shift_field_hammett(gkeep(:, :, :, :, ivmu))
         end do
         where (aky < remap_em_kymin .or. aky > remap_em_kymax) shifted = .false.
      end if

      fields_updated = .false.
      call advance_fields(g, phi, apar, bpar, dist='g')

      select case (remap_em)
      case (1)
         call gbar_to_g(g, apar, -1.0)
         call advance_perp_flow_shear(g)
      case (2)
         allocate (apar_s, source=apar)
         call shift_field_hammett(apar_s)
         call advance_perp_flow_shear(g)
         call gbar_to_g(g, apar_s, -1.0)
         deallocate (apar_s)
      case default
         call advance_perp_flow_shear(g)
         deallocate (gkeep, shifted)
         return
      end select

      fields_updated = .false.
      call advance_fields(g, phi, apar, bpar, dist='gbar')
      call gbar_to_g(g, apar, 1.0)

      ! Restore every k_y that was not shifted (the remap leaves them unchanged; the round trip would not), and the
      ! shifted k_y outside the remap_em range (plain shift of g, done on the kept copy above)
      do iky = 1, naky
         if (.not. shifted(iky)) g(iky, :, :, :, :) = gkeep(iky, :, :, :, :)
      end do
      deallocate (gkeep, shifted)
      fields_updated = .false.

   end subroutine remap_flow_shear

   !****************************************************************************
   !        Diagnostic: A_par after the original and after the gbar remap
   !****************************************************************************
   ! From the same incoming g, on copies (g and shift_state are left unchanged):
   !   A_pre   A_par solved from g before the remap, moved one cell with the shift (the value each cell inherits);
   !   A_g     original remap: shift g, solve A_par from g            (A_par = j / k_perp^2);
   !   A_gbar  remap_em = 1:   build gbar, shift it, solve from gbar  (A_par = j(gbar) / apar_denom);
   !   A_cons  no shift: build gbar from g and solve from gbar; must equal A_par from g (checks that the g <-> gbar
   !           conversion and apar_denom use consistent velocity integrals).
   ! proc0 appends one line per shifted (k_y, k_x) to remap_dump.txt: field-line sums of |A|^2 and the values at z = 0;
   ! for the smallest non-zero k_y it also writes the whole field A_par(k_x, z) to remap_dump_field.txt.
   subroutine dump_remap_comparison(g, phi, apar, bpar)

      use mp, only: proc0
      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid, ntubes, zed
      use grids_kxky, only: naky, nakx, aky, akx, zonal_mode
      use grids_time, only: code_time
      use arrays, only: kperp2
      use field_equations, only: advance_fields, fields_updated
      use calculations_tofrom_ghf, only: gbar_to_g
      use gk_flow_shear, only: shift_field_hammett, perp_shift_mask

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
      complex, dimension(:, :, :, :, :), allocatable :: gw
      complex, dimension(:, :, :, :), allocatable :: a_pre, a_pre0, a_g, a_gbar, a_cons, kp_src
      logical, dimension(:), allocatable :: shifted
      integer :: iky, ikx, ivmu, unit, iz, iky_min
      real :: cons
      logical, save :: first = .true., first_field = .true.

      allocate (shifted(naky))
      call perp_shift_mask(shifted)
      allocate (gw, source=g)
      allocate (a_pre, a_pre0, a_g, a_gbar, a_cons, kp_src, mold=apar)

      ! A_par from g before the remap
      fields_updated = .false.
      call advance_fields(gw, phi, apar, bpar, dist='g')
      a_pre = apar

      ! Consistency: gbar from g, A_par from gbar, no shift
      call gbar_to_g(gw, a_pre, -1.0)
      fields_updated = .false.
      call advance_fields(gw, phi, apar, bpar, dist='gbar')
      a_cons = apar
      a_pre0 = a_pre
      cons = sqrt(sum(abs(a_cons - a_pre)**2) / max(sum(abs(a_pre)**2), tiny(1.0)))

      ! remap_em = 1: shift that gbar, A_par from gbar
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         call shift_field_hammett(gw(:, :, :, :, ivmu))
      end do
      fields_updated = .false.
      call advance_fields(gw, phi, apar, bpar, dist='gbar')
      a_gbar = apar

      ! Original: shift g, A_par from g
      gw = g
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         call shift_field_hammett(gw(:, :, :, :, ivmu))
      end do
      fields_updated = .false.
      call advance_fields(gw, phi, apar, bpar, dist='g')
      a_g = apar
      fields_updated = .false.

      ! k_perp^2 of the cell each mode comes from, and A_pre moved with its mode
      kp_src = spread(cmplx(kperp2(:, :, 1, :)), 4, ntubes)
      call shift_field_hammett(kp_src)
      call shift_field_hammett(a_pre)

      if (proc0) then
         open (newunit=unit, file='remap_dump.txt', position='append', action='write')
         if (first) then
            write (unit, '(A)') '# remap_em_ndump: A_par from the same state after the original remap (shift g) and the'
            write (unit, '(A)') '# gbar remap; one line per shifted (ky,kx); sums are over z (tube 1); z0 = values at z = 0.'
            write (unit, '(A)') '# event t ky kx_new kp2_src_z0 kp2_new_z0 sumA2_pre sumA2_g sumA2_gbar' // &
               ' ReA_pre_z0 ImA_pre_z0 ReA_g_z0 ImA_g_z0 ReA_gbar_z0 ImA_gbar_z0 cons_rel'
            first = .false.
         end if
         do iky = 1, naky
            if (.not. shifted(iky) .or. zonal_mode(iky)) cycle
            do ikx = 1, nakx
               write (unit, '(I4, 1X, ES14.7, 2(1X, F9.5), 2(1X, ES12.5), 3(1X, ES12.5), 6(1X, ES13.5), 1X, ES10.3)') &
                  ndump_done + 1, code_time, aky(iky), akx(ikx), real(kp_src(iky, ikx, 0, 1)), kperp2(iky, ikx, 1, 0), &
                  sum(abs(a_pre(iky, ikx, :, 1))**2), sum(abs(a_g(iky, ikx, :, 1))**2), &
                  sum(abs(a_gbar(iky, ikx, :, 1))**2), &
                  real(a_pre(iky, ikx, 0, 1)), aimag(a_pre(iky, ikx, 0, 1)), real(a_g(iky, ikx, 0, 1)), &
                  aimag(a_g(iky, ikx, 0, 1)), real(a_gbar(iky, ikx, 0, 1)), aimag(a_gbar(iky, ikx, 0, 1)), cons
            end do
         end do
         close (unit)

         ! The whole field for the smallest non-zero k_y (its remap is the one that matters most)
         iky_min = minloc(aky, dim=1, mask=aky > 0.)
         if (shifted(iky_min)) then
            open (newunit=unit, file='remap_dump_field.txt', position='append', action='write')
            if (first_field) then
               write (unit, '(A)') '# A_par(k_x, z) of the smallest non-zero k_y at its remaps (tube 1).  pre0, cons: before'
               write (unit, '(A)') '# the remap at the cell itself (from g; from gbar, no shift; must agree).  pre: the value'
               write (unit, '(A)') '# the cell inherits (pre0 moved with the shift).  g: original remap.  gbar: gbar remap.'
               write (unit, '(A)') '# event t ky ikx kx iz zed Re_pre0 Im_pre0 Re_cons Im_cons Re_pre Im_pre Re_g Im_g' // &
                  ' Re_gbar Im_gbar kp2_src kp2_new'
               first_field = .false.
            end if
            do ikx = 1, nakx
               do iz = -nzgrid, nzgrid
                  write (unit, '(I4, 1X, ES14.7, 1X, F9.5, 1X, I4, 1X, F9.5, 1X, I4, 1X, F10.5, 10(1X, ES13.5), 2(1X, ES12.5))') &
                     ndump_done + 1, code_time, aky(iky_min), ikx, akx(ikx), iz, zed(iz), &
                     real(a_pre0(iky_min, ikx, iz, 1)), aimag(a_pre0(iky_min, ikx, iz, 1)), &
                     real(a_cons(iky_min, ikx, iz, 1)), aimag(a_cons(iky_min, ikx, iz, 1)), &
                     real(a_pre(iky_min, ikx, iz, 1)), aimag(a_pre(iky_min, ikx, iz, 1)), &
                     real(a_g(iky_min, ikx, iz, 1)), aimag(a_g(iky_min, ikx, iz, 1)), &
                     real(a_gbar(iky_min, ikx, iz, 1)), aimag(a_gbar(iky_min, ikx, iz, 1)), &
                     real(kp_src(iky_min, ikx, iz, 1)), kperp2(iky_min, ikx, 1, iz)
               end do
            end do
            close (unit)
         end if
      end if
      ndump_done = ndump_done + 1

      deallocate (gw, a_pre, a_pre0, a_g, a_gbar, a_cons, kp_src, shifted)

   end subroutine dump_remap_comparison

end module gyrokinetic_equation_implicit
