module oas_sendReceiveMod
  use shr_kind_mod     , only: r8 => shr_kind_r8
  use abortutils       , only: endrun
  use clm_time_manager , only: get_nstep, get_step_size
  use decompMod        , only: bounds_type
  use clm_varpar       , only: nlevgrnd
  use clm_varctl       , only: iulog
  use oas_vardefMod
  use mod_oasis
  implicit none
  save
  private

#ifdef COUP_OAS_PFL
  public  :: oas_send
  public  :: oas_receive
#endif

#ifdef COUP_OAS_ICON
  public  :: oas_send_icon
  public  :: oas_receive_icon
#endif 

contains

#ifdef COUP_OAS_PFL
  subroutine oas_receive(bounds, seconds_elapsed, atm2lnd_inst)
    use atm2lndType, only: atm2lnd_type

    type(bounds_type),  intent(in)    :: bounds
    integer          ,  intent(in)    :: seconds_elapsed
    type(atm2lnd_type), intent(inout) :: atm2lnd_inst
    integer                           :: info

    call oasis_get(oas_psi_id, seconds_elapsed, atm2lnd_inst%pfl_psi_grc, info)
    call oasis_get(oas_sat_id, seconds_elapsed, atm2lnd_inst%pfl_h2osoi_liq_grc, info)

  end subroutine oas_receive

  subroutine oas_send(bounds, seconds_elapsed, lnd2atm_inst)
    use lnd2atmType, only : lnd2atm_type
    use spmdMod,     only : mpicom
    use shr_mpi_mod, only: shr_mpi_barrier

    type(bounds_type),  intent(in)    :: bounds
    integer          ,  intent(in)    :: seconds_elapsed
    type(lnd2atm_type), intent(inout) :: lnd2atm_inst

    integer                           :: info
    
    call oasis_put(oas_et_loss_id, seconds_elapsed, lnd2atm_inst%qflx_parflow_grc, info)
    
  end subroutine oas_send
#endif

#ifdef COUP_OAS_ICON
  subroutine oas_receive_icon(bounds, seconds_elapsed, atm2lnd_inst)
    use atm2lndType, only: atm2lnd_type
    use clm_varctl      , only: co2_type, co2_ppmv, iulog, use_c13
    use clm_varcon      , only: rair, o2_molar_const, c13ratio
    use shr_const_mod   , only: SHR_CONST_TKFRZ
    use clm_time_manager   , only : get_nstep, get_step_size
    use m_MCTWorld      , only: ThisMCTWorld
    use clm_cpl_indices , only: index_x2l_Sa_co2prog, index_x2l_Sa_co2diag

    type(bounds_type),  intent(in)    :: bounds
    integer          ,  intent(in)    :: seconds_elapsed
    type(atm2lnd_type), intent(inout) :: atm2lnd_inst
    real(kind=r8),      allocatable   :: buffer(:,:)
    integer                           :: num_grid_points
    integer                           :: info
    integer                           :: g

    real(r8) :: forc_rainc           ! rainxy Atm flux mm/s
    real(r8) :: e                    ! vapor pressure (Pa)
    real(r8) :: qsat                 ! saturation specific humidity (kg/kg)
    real(r8) :: forc_t               ! atmospheric temperature (Kelvin)
    real(r8) :: forc_q               ! atmospheric specific humidity (kg/kg)
    real(r8) :: forc_pbot            ! atmospheric pressure (Pa)
    real(r8) :: forc_rainl           ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc           ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl           ! snowfxl Atm flux  mm/s
    real(r8) :: co2_ppmv_diag        ! temporary
    real(r8) :: co2_ppmv_prog        ! temporary
    real(r8) :: co2_ppmv_val         ! temporary
    integer  :: co2_type_idx         ! integer flag for co2_type options
    real(r8) :: esatw                ! saturation vapor pressure over water (Pa)
    real(r8) :: esati                ! saturation vapor pressure over ice (Pa)
    real(r8) :: a0,a1,a2,a3,a4,a5,a6 ! coefficients for esat over water
    real(r8) :: b0,b1,b2,b3,b4,b5,b6 ! coefficients for esat over ice
    real(r8) :: tdc, t               ! Kelvins to Celcius function and its input
    character(len=32), parameter :: sub = 'lnd_import'

    ! Constants to compute vapor pressure
    parameter (a0=6.107799961_r8    , a1=4.436518521e-01_r8, &
         a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
         a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
         a6=6.136820929e-11_r8)

    parameter (b0=6.109177956_r8    , b1=5.034698970e-01_r8, &
         b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
         b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
         b6=1.838826904e-10_r8)
    !
    ! function declarations
    !
    tdc(t) = min( 50._r8, max(-50._r8,(t-SHR_CONST_TKFRZ)) )
    esatw(t) = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    esati(t) = 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
    !---------------------------------------------------------------------------

    co2_type_idx = 0
    if (co2_type == 'prognostic') then
       co2_type_idx = 1
    else if (co2_type == 'diagnostic') then
       co2_type_idx = 2
    end if
    if (co2_type == 'prognostic' .and. index_x2l_Sa_co2prog == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2prog for co2_type equal to prognostic' )
    else if (co2_type == 'diagnostic' .and. index_x2l_Sa_co2diag == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2diag for co2_type equal to diagnostic' )
    end if

    num_grid_points = (bounds%endg - bounds%begg) + 1
    allocate(buffer(num_grid_points, 1))

    !write(iulog,*) 'BEFORE oasis_get(): ', ThisMCTWorld%mygrank, 'nstep = ', get_nstep()
    !write(iulog,*) 'BEFORE oasis_get(): ', ThisMCTWorld%mygrank, '<Faxa_lwdn> = ', sum( atm2lnd_inst%forc_lwrad_not_downscaled_grc(:) )/size( atm2lnd_inst%forc_lwrad_not_downscaled_grc(:) )
    !atm2lnd_inst%forc_lwrad_not_downscaled_grc(:) = 333.333
    write(iulog,*) 'BEFOR oasis_get():  ', get_nstep(), ThisMCTWorld%mygrank, &
    & 'lwrad =', sum( atm2lnd_inst%forc_lwrad_not_downscaled_grc(:) )/size( atm2lnd_inst%forc_lwrad_not_downscaled_grc(:) )

    !    call oasis_get(oas_id_t, seconds_elapsed, oas_rcv_meta(:,:,oas_id_t), info)
    call oasis_get(oas_id_t,  seconds_elapsed, atm2lnd_inst%forc_t_not_downscaled_grc, info)
    call oasis_get(oas_id_u,  seconds_elapsed, atm2lnd_inst%forc_u_grc, info)
    call oasis_get(oas_id_v,  seconds_elapsed, atm2lnd_inst%forc_v_grc, info)
    call oasis_get(oas_id_qv, seconds_elapsed, atm2lnd_inst%forc_q_not_downscaled_grc, info)
    call oasis_get(oas_id_ht, seconds_elapsed, atm2lnd_inst%forc_hgt_grc, info)
    call oasis_get(oas_id_pr, seconds_elapsed, atm2lnd_inst%forc_pbot_not_downscaled_grc, info)
    call oasis_get(oas_id_rs, seconds_elapsed, atm2lnd_inst%forc_solad_grc(:,1), info)
    call oasis_get(oas_id_fs, seconds_elapsed, atm2lnd_inst%forc_solai_grc(:,1), info)
    call oasis_get(oas_id_lw, seconds_elapsed, atm2lnd_inst%forc_lwrad_not_downscaled_grc, info)
    call oasis_get(oas_id_cr, seconds_elapsed, atm2lnd_inst%forc_rain_not_downscaled_grc, info)
    call oasis_get(oas_id_gr, seconds_elapsed, atm2lnd_inst%forc_snow_not_downscaled_grc, info)

    !SPo: some postprocessing of atm2lnd is missing; may better use x2l 

    do g=bounds%begg,bounds%endg
       forc_t = atm2lnd_inst%forc_t_not_downscaled_grc(g)
       forc_q = atm2lnd_inst%forc_q_not_downscaled_grc(g)
       forc_pbot = atm2lnd_inst%forc_pbot_not_downscaled_grc(g)

       atm2lnd_inst%forc_hgt_u_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of wind [m]
       atm2lnd_inst%forc_hgt_t_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of temperature [m]
       atm2lnd_inst%forc_hgt_q_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of humidity [m]
       atm2lnd_inst%forc_vp_grc(g)    = forc_q * forc_pbot  / (0.622_r8 + 0.378_r8 * forc_q)
       atm2lnd_inst%forc_rho_not_downscaled_grc(g) = &
            (forc_pbot - 0.378_r8 * atm2lnd_inst%forc_vp_grc(g)) / (rair * forc_t)
       atm2lnd_inst%forc_po2_grc(g)   = o2_molar_const * forc_pbot
       atm2lnd_inst%forc_wind_grc(g)  = sqrt(atm2lnd_inst%forc_u_grc(g)**2 + atm2lnd_inst%forc_v_grc(g)**2)
       atm2lnd_inst%forc_solar_grc(g) = atm2lnd_inst%forc_solad_grc(g,1) + atm2lnd_inst%forc_solai_grc(g,1) + &
                                        atm2lnd_inst%forc_solad_grc(g,2) + atm2lnd_inst%forc_solai_grc(g,2)

       forc_rainc = 0
       forc_rainl = atm2lnd_inst%forc_rain_not_downscaled_grc(g)  ! mm/s
       forc_snowc = 0
       forc_snowl = atm2lnd_inst%forc_snow_not_downscaled_grc(g)  ! mm/s

       !atm2lnd_inst%forc_rain_not_downscaled_grc(g)  = forc_rainc + forc_rainl
       !atm2lnd_inst%forc_snow_not_downscaled_grc(g)  = forc_snowc + forc_snowl

       if (forc_t > SHR_CONST_TKFRZ) then
          e = esatw(tdc(forc_t))
       else
          e = esati(tdc(forc_t))
       end if
       qsat           = 0.622_r8*e / (forc_pbot - 0.378_r8*e)

       !modify specific humidity if precip occurs
       if(1==2) then
          if((forc_rainc+forc_rainl) > 0._r8) then
             forc_q = 0.95_r8*qsat
             !           forc_q = qsat
             atm2lnd_inst%forc_q_not_downscaled_grc(g) = forc_q
          endif
       endif

       atm2lnd_inst%forc_rh_grc(g) = 100.0_r8*(forc_q / qsat)

       ! Determine derived quantities for optional fields
       ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
       ! Note that forc_pbot is in Pa

       if (co2_type_idx == 1) then
          co2_ppmv_val = co2_ppmv_prog
       else if (co2_type_idx == 2) then
          co2_ppmv_val = co2_ppmv_diag
       else
          co2_ppmv_val = co2_ppmv
       end if
       if ( (co2_ppmv_val < 10.0_r8) .or. (co2_ppmv_val > 15000.0_r8) )then
          call endrun( sub//' ERROR: CO2 is outside of an expected range' )
       end if
       atm2lnd_inst%forc_pco2_grc(g)   = co2_ppmv_val * 1.e-6_r8 * forc_pbot
       if (use_c13) then
          atm2lnd_inst%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
       end if

       atm2lnd_inst%forc_solad_grc(g,1) = 0.5_r8 * atm2lnd_inst%forc_solad_grc(g,1)
       atm2lnd_inst%forc_solad_grc(g,2) = atm2lnd_inst%forc_solad_grc(g,1)
       atm2lnd_inst%forc_solai_grc(g,1) = 0.5_r8 * atm2lnd_inst%forc_solai_grc(g,1)
       atm2lnd_inst%forc_solai_grc(g,2) = atm2lnd_inst%forc_solai_grc(g,1)
       atm2lnd_inst%forc_solar_grc(g)   =  atm2lnd_inst%forc_solad_grc(g,2) + atm2lnd_inst%forc_solad_grc(g,1) &
                                        +  atm2lnd_inst%forc_solai_grc(g,2) + atm2lnd_inst%forc_solai_grc(g,1)
    enddo

    write(iulog,*) 'AFTER oasis_get():  ', get_nstep(), ThisMCTWorld%mygrank, &
    & 'lwrad =', sum( atm2lnd_inst%forc_lwrad_not_downscaled_grc(:) )/size( atm2lnd_inst%forc_lwrad_not_downscaled_grc(:) )

  end subroutine oas_receive_icon

  subroutine oas_send_icon(bounds, seconds_elapsed, lnd2atm_inst)
    use lnd2atmType, only : lnd2atm_type
    use spmdMod,     only : mpicom
    use shr_mpi_mod, only: shr_mpi_barrier

    type(bounds_type),  intent(in)    :: bounds
    integer          ,  intent(in)    :: seconds_elapsed
    type(lnd2atm_type), intent(inout) :: lnd2atm_inst
    real(kind=r8),      allocatable   :: aux_buffer(:,:)
    integer                           :: num_grid_points
    integer                           :: info

    num_grid_points = (bounds%endg - bounds%begg) + 1
    allocate(aux_buffer(num_grid_points, 1))

     call oasis_put(oas_id_it, seconds_elapsed,lnd2atm_inst%t_rad_grc, info)         ! "CLM_INFRA"
     call oasis_put(oas_id_ad, seconds_elapsed,lnd2atm_inst%albd_grc(:,1), info)     ! "CLM_ALBED"
     call oasis_put(oas_id_ai, seconds_elapsed,lnd2atm_inst%albi_grc(:,1), info)     ! "CLM_ALBEI"
     call oasis_put(oas_id_tx, seconds_elapsed,lnd2atm_inst%taux_grc, info)          ! "CLM_TAUX"
     call oasis_put(oas_id_ty, seconds_elapsed,lnd2atm_inst%tauy_grc, info)          ! "CLM_TAUY"
     call oasis_put(oas_id_sh, seconds_elapsed,lnd2atm_inst%eflx_sh_tot_grc, info)   ! "CLM_SHFLX"
     call oasis_put(oas_id_lh, seconds_elapsed,lnd2atm_inst%eflx_lh_tot_grc, info)   ! "CLM_LHFLX"
     call oasis_put(oas_id_st, seconds_elapsed,lnd2atm_inst%t_sf_grc, info)          ! "CLM_TGRND"

  end subroutine oas_send_icon
#endif

end module oas_sendReceiveMod
