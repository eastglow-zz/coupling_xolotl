# Length unit: nm
# Time unit: s
# Mass unit: ?

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 120
  ny = 120
  xmin = 0
  xmax = 1200
  ymin = 0
  ymax = 1200
  #uniform_refine = 3
[]

[GlobalParams]
  op_num = 5
  grain_num = 5
  var_name_base = etam
  numbub = 1
  bubspac = 150
  radius = 44
  int_width = 50
  # displacements = 'disp_x disp_y'
[]

[Variables]
  [./wv]
  [../]
  [./wg]
  [../]
  [./etab0]
  [../]
  [./PolycrystalVariables]
  [../]
  # Displacement fields in x and y directions
  # [./disp_x]
  #   order = FIRST
  #   family = LAGRANGE
  # [../]
  # [./disp_y]
  #   order = FIRST
  #   family = LAGRANGE
  # [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./XolotlXeRate]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalVoronoiVoidIC]
      invalue = 1.0
      outvalue = 0.0
    [../]
  [../]
  [./bubble_IC]
    variable = etab0
    type = PolycrystalVoronoiVoidIC
    structure_type = voids
    invalue = 1.0
    outvalue = 0.0
  [../]
  [./IC_wv]
    variable = wv
    type = PolycrystalVoronoiVoidIC
    structure_type = voids
    invalue = 0.0
    outvalue = 0.0
  [../]
  [./IC_wg]
    variable = wg
    type = PolycrystalVoronoiVoidIC
    structure_type = voids
    invalue = 0.0
    outvalue = 0.0
  [../]
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Kernels]
  # [./TensorMechanics]
  # [../]

# Order parameter eta_b0 for bubble phase
  [./ACb0_bulk]
    type = ACGrGrMulti
    variable = etab0
    v =           'etam0 etam1 etam2 etam3 etam4'
    gamma_names = 'gmb   gmb   gmb   gmb   gmb'
  [../]
  [./ACb0_sw]
    type = ACSwitching
    variable = etab0
    Fj_names  = 'omega_total_bubble   omega_total_matrix'
    hj_names  = 'hb                   hm'
    args = 'etam0 etam1 etam2 etam3 etam4 wv wg'
  [../]
  [./ACb0_int]
    type = ACInterface
    variable = etab0
    kappa_name = kappa
  [../]
  [./eb0_dot]
    type = TimeDerivative
    variable = etab0
  [../]
# Order parameter eta_m0 for matrix grain 0
  [./ACm0_bulk]
    type = ACGrGrMulti
    variable = etam0
    v =           'etab0 etam1 etam2 etam3 etam4'
    gamma_names = 'gmb   gmm   gmm   gmm   gmm'
  [../]
  [./ACm0_sw]
    type = ACSwitching
    variable = etam0
    Fj_names  = 'omega_total_bubble   omega_total_matrix'
    hj_names  = 'hb                   hm'
    args = 'etab0 etam1 etam2 etam3 etam4 wv wg'
  [../]
  [./ACm0_int]
    type = ACInterface
    variable = etam0
    kappa_name = kappa
  [../]
  [./em0_dot]
    type = TimeDerivative
    variable = etam0
  [../]
# Order parameter eta_m1 for matrix grain 1
  [./ACm1_bulk]
    type = ACGrGrMulti
    variable = etam1
    v =           'etab0 etam0 etam2 etam3 etam4'
    gamma_names = 'gmb   gmm   gmm   gmm   gmm'
  [../]
  [./ACm1_sw]
    type = ACSwitching
    variable = etam1
    Fj_names  = 'omega_total_bubble   omega_total_matrix'
    hj_names  = 'hb                   hm'
    args = 'etab0 etam0 etam2 etam3 etam4 wv wg'
  [../]
  [./ACm1_int]
    type = ACInterface
    variable = etam1
    kappa_name = kappa
  [../]
  [./em1_dot]
    type = TimeDerivative
    variable = etam1
  [../]
# Order parameter eta_m2 for matrix grain 2
  [./ACm2_bulk]
    type = ACGrGrMulti
    variable = etam2
    v =           'etab0 etam0 etam1 etam3 etam4'
    gamma_names = 'gmb   gmm   gmm   gmm   gmm'
  [../]
  [./ACm2_sw]
    type = ACSwitching
    variable = etam2
    Fj_names  = 'omega_total_bubble   omega_total_matrix'
    hj_names  = 'hb                   hm'
    args = 'etab0 etam0 etam1 etam3 etam4 wv wg'
  [../]
  [./ACm2_int]
    type = ACInterface
    variable = etam2
    kappa_name = kappa
  [../]
  [./em2_dot]
    type = TimeDerivative
    variable = etam2
  [../]
# Order parameter eta_m3 for matrix grain 3
  [./ACm3_bulk]
    type = ACGrGrMulti
    variable = etam3
    v =           'etab0 etam0 etam1 etam2 etam4'
    gamma_names = 'gmb   gmm   gmm   gmm   gmm'
  [../]
  [./ACm3_sw]
    type = ACSwitching
    variable = etam3
    Fj_names  = 'omega_total_bubble   omega_total_matrix'
    hj_names  = 'hb                   hm'
    args = 'etab0 etam0 etam1 etam2 etam4 wv wg'
  [../]
  [./ACm3_int]
    type = ACInterface
    variable = etam3
    kappa_name = kappa
  [../]
  [./em3_dot]
    type = TimeDerivative
    variable = etam3
  [../]
# Order parameter eta_m4 for matrix grain 4
  [./ACm4_bulk]
    type = ACGrGrMulti
    variable = etam4
    v =           'etab0 etam0 etam1 etam2 etam3'
    gamma_names = 'gmb   gmm   gmm   gmm   gmm'
  [../]
  [./ACm4_sw]
    type = ACSwitching
    variable = etam4
    Fj_names  = 'omega_total_bubble   omega_total_matrix'
    hj_names  = 'hb                   hm'
    args = 'etab0 etam0 etam1 etam2 etam3 wv wg'
  [../]
  [./ACm4_int]
    type = ACInterface
    variable = etam4
    kappa_name = kappa
  [../]
  [./em4_dot]
    type = TimeDerivative
    variable = etam4
  [../]

#Chemical potential for vacancies
  [./wv_dot]
    type = SusceptibilityTimeDerivative
    variable = wv
    f_name = chiv
    args = '' # in this case chi (the susceptibility) is simply a constant
  [../]
  [./Diffusion_v]
    type = MatDiffusion
    variable = wv
    D_name = Dchiv
    args = ''
  [../]
  # [./Source_v]
  #   type = MaskedBodyForce
  #   variable = wv
  #   value = 1e-10
  #   mask = hm
  # [../]
  [./Source_v]
    type = MaskedBodyForce
    variable = wv
    value = 1
    mask = VacRate
  [../]
  [./coupled_v_etab0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wv
    v = etab0
    Fj_names = 'rhovbub rhovmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1 etam2 etam3 etam4'
  [../]
  [./coupled_v_etam0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wv
    v = etam0
    Fj_names = 'rhovbub rhovmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1 etam2 etam3 etam4'
  [../]
  [./coupled_v_etam1dot]
    type = CoupledSwitchingTimeDerivative
    variable = wv
    v = etam1
    Fj_names = 'rhovbub rhovmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1 etam2 etam3 etam4'
  [../]
  [./coupled_v_etam2dot]
    type = CoupledSwitchingTimeDerivative
    variable = wv
    v = etam2
    Fj_names = 'rhovbub rhovmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1 etam2 etam3 etam4'
  [../]
  [./coupled_v_etam3dot]
    type = CoupledSwitchingTimeDerivative
    variable = wv
    v = etam3
    Fj_names = 'rhovbub rhovmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1 etam2 etam3 etam4'
  [../]
  [./coupled_v_etam4dot]
    type = CoupledSwitchingTimeDerivative
    variable = wv
    v = etam4
    Fj_names = 'rhovbub rhovmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1 etam2 etam3 etam4'
  [../]

#Chemical potential for gas atoms
  [./wg_dot]
    type = SusceptibilityTimeDerivative
    variable = wg
    f_name = chig
    args = '' # in this case chi (the susceptibility) is simply a constant
  [../]
  [./Diffusion_g]
    type = MatDiffusion
    variable = wg
    D_name = Dchig
    args = ''
  [../]
  # [./Source_g]
  #   type = MaskedBodyForce
  #   variable = wg
  #   value = 1e-10
  #   mask = hm
  # [../]
  [./Source_g]
    type = MaskedBodyForce
    variable = wg
    value = 1 # for unit conversion between PF app and Xolotl
    mask = XeRate
  [../]
  [./coupled_g_etab0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wg
    v = etab0
    Fj_names = 'rhogbub rhogmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1 etam2 etam3 etam4'
  [../]
  [./coupled_g_etam0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wg
    v = etam0
    Fj_names = 'rhogbub rhogmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1 etam2 etam3 etam4'
  [../]
  [./coupled_g_etam1dot]
    type = CoupledSwitchingTimeDerivative
    variable = wg
    v = etam1
    Fj_names = 'rhogbub rhogmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1 etam2 etam3 etam4'
  [../]
  [./coupled_g_etam2dot]
    type = CoupledSwitchingTimeDerivative
    variable = wg
    v = etam2
    Fj_names = 'rhogbub rhogmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1 etam2 etam3 etam4'
  [../]
  [./coupled_g_etam3dot]
    type = CoupledSwitchingTimeDerivative
    variable = wg
    v = etam3
    Fj_names = 'rhogbub rhogmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1 etam2 etam3 etam4'
  [../]
  [./coupled_g_etam4dot]
    type = CoupledSwitchingTimeDerivative
    variable = wg
    v = etam4
    Fj_names = 'rhogbub rhogmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1 etam2 etam3 etam4'
  [../]

[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
[]

[Materials]
  [./hb]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = hb
    all_etas = 'etab0 etam0 etam1 etam2 etam3 etam4'
    phase_etas = 'etab0'
    #outputs = exodus
  [../]
  [./hm]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = hm
    all_etas = 'etab0 etam0 etam1 etam2 etam3 etam4'
    phase_etas = 'etam0 etam1 etam2 etam3 etam4'
    #outputs = exodus
  [../]
# Chemical contribution to grand potential of bubble
  [./omegab]
    type = DerivativeParsedMaterial
    args = 'wv wg'
    f_name = omegab
    material_property_names = 'Va kvbub cvbubeq kgbub cgbubeq'
    function = '-0.5*wv^2/Va^2/kvbub-wv/Va*cvbubeq-0.5*wg^2/Va^2/kgbub-wg/Va*cgbubeq'
    derivative_order = 2
    #outputs = exodus
  [../]
  # [./elastic_energy_bubble]
  #   type = ElasticEnergyMaterial
  #   base_name = bubble
  #   f_name = fe_bub
  #   args = ' '
  # [../]
# Total free energy of the bubble
  [./Total_energy_bubble]
    type = DerivativeSumMaterial
    f_name = omega_total_bubble
    # sum_materials = 'omegab fe_bub'
    sum_materials = 'omegab'
    args = 'wv wg'
  [../]

# Chemical contribution to grand potential of matrix
  [./omegam]
    type = DerivativeParsedMaterial
    args = 'wv wg'
    f_name = omegam
    material_property_names = 'Va kvmatrix cvmatrixeq kgmatrix cgmatrixeq'
    function = '-0.5*wv^2/Va^2/kvmatrix-wv/Va*cvmatrixeq-0.5*wg^2/Va^2/kgmatrix-wg/Va*cgmatrixeq'
    derivative_order = 2
    #outputs = exodus
  [../]
  # [./elastic_energy_matrix]
  #   type = ElasticEnergyMaterial
  #   base_name = matrix
  #   f_name = fe_m
  #   args = ' '
  # [../]
# Total free energy of the matrix
  [./Total_energy_matrix]
    type = DerivativeSumMaterial
    block = 0
    f_name = omega_total_matrix
    # sum_materials = 'omegam fe_m'
    sum_materials = 'omegam'
    args = 'wv wg'
  [../]
# Densities
  [./rhovbub]
    type = DerivativeParsedMaterial
    args = 'wv'
    f_name = rhovbub
    material_property_names = 'Va kvbub cvbubeq'
    function = 'wv/Va^2/kvbub + cvbubeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./rhovmatrix]
    type = DerivativeParsedMaterial
    args = 'wv'
    f_name = rhovmatrix
    material_property_names = 'Va kvmatrix cvmatrixeq'
    function = 'wv/Va^2/kvmatrix + cvmatrixeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./rhogbub]
    type = DerivativeParsedMaterial
    args = 'wg'
    f_name = rhogbub
    material_property_names = 'Va kgbub cgbubeq'
    function = 'wg/Va^2/kgbub + cgbubeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./rhogmatrix]
    type = DerivativeParsedMaterial
    args = 'wg'
    f_name = rhogmatrix
    material_property_names = 'Va kgmatrix cgmatrixeq'
    function = 'wg/Va^2/kgmatrix + cgmatrixeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./const]
    type = GenericConstantMaterial
    prop_names =  'kappa   mu       L   D    Va      cvbubeq cgbubeq gmb 	  gmm T    tgrad_corr_mult YXe'
    prop_values = '0.5273  0.004688 1.0 0.01 0.04092 0.61    0.39    0.9218 1.5 1200 0.0             0.2156'
  [../]
  [./cvmatrixeq]    #For values, see Li et al., Nuc. Inst. Methods in Phys. Res. B, 303, 62-27 (2013).
    type = ParsedMaterial
    f_name = cvmatrixeq
    material_property_names = 'T'
    constant_names        = 'kB           Efv'
    constant_expressions  = '8.6173324e-5 3.0'
    function = 'exp(-Efv/(kB*T))'
  [../]
  [./cgmatrixeq]
    type = ParsedMaterial
    f_name = cgmatrixeq
    material_property_names = 'T'
    constant_names        = 'kB           Efg'
    constant_expressions  = '8.6173324e-5 3.0'
    function = 'exp(-Efg/(kB*T))'
  [../]
  [./kvmatrix_parabola]
    type = ParsedMaterial
    f_name = kvmatrix
    material_property_names = 'T  cvmatrixeq'
    constant_names        = 'c0v  c0g  a1                                               a2'
    constant_expressions  = '0.01 0.01 0.178605-0.0030782*log(1-c0v)+0.0030782*log(c0v) 0.178605-0.00923461*log(1-c0v)+0.00923461*log(c0v)'
    function = '((-a2+3*a1)/(4*(c0v-cvmatrixeq))+(a2-a1)/(2400*(c0v-cvmatrixeq))*T)'
    #outputs = exodus
  [../]
  [./kgmatrix_parabola]
    type = ParsedMaterial
    f_name = kgmatrix
    material_property_names = 'kvmatrix'
    function = 'kvmatrix'
  [../]
  [./kgbub_parabola]
    type = ParsedMaterial
    f_name = kgbub
    material_property_names = 'kvmatrix  cvmatrixeq cvbubeq'
    constant_names        = 'fcross'
    constant_expressions  = '0.5'   #Scaled by C44
    function = 'kvmatrix * fcross/(sqrt(kvmatrix)*(cvmatrixeq-cvbubeq) + sqrt(fcross))^2'
    #outputs = exodus
  [../]
  [./kvbub_parabola]
    type = ParsedMaterial
    f_name = kvbub
    material_property_names = 'kgbub'
    function = 'kgbub'
  [../]
  [./Mobility_v]
    type = DerivativeParsedMaterial
    f_name = Dchiv
    material_property_names = 'D chiv'
    function = 'D*chiv'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./Mobility_g]
    type = DerivativeParsedMaterial
    f_name = Dchig
    material_property_names = 'D chig'
    function = 'D*chig'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./chiv]
    type = DerivativeParsedMaterial
    f_name = chiv
    material_property_names = 'Va hb kvbub hm kvmatrix '
    function = '(hm/kvmatrix + hb/kvbub) / Va^2'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./chig]
    type = DerivativeParsedMaterial
    f_name = chig
    material_property_names = 'Va hb kgbub hm kgmatrix '
    function = '(hm/kgmatrix + hb/kgbub) / Va^2'
    derivative_order = 2
    #outputs = exodus
  [../]

  #Mechanical properties
  # [./Stiffness_matrix]
  #   type = ComputeElasticityTensor
  #   C_ijkl = '6.17 1.89 1.89 6.17 1.89 6.17 1 1 1'
  #   # Elastic constants are in eV/atom to match other energies in the problem
  #   # For symmetric9 with cubic anisotropy the constants are in order
  #   # C11 C12 C12 C11 C12 C11 C44 C44 C44
  #   base_name = matrix
  #   fill_method = symmetric9
  # [../]
  # [./strain_matrix]
  #   type = ComputeSmallStrain
  #   base_name = matrix
  # [../]
  # [./stress_matrix]
  #   type = ComputeLinearElasticStress
  #   base_name = matrix
  # [../]
  #
  # [./Stiffness_bubble]
  #   type = ComputeElasticityTensor
  #   C_ijkl = '6.17e-4 1.89e-4 1.89e-4 6.17e-4 1.89e-4 6.17e-4 1e-4 1e-4 1e-4'
  #   base_name = bubble
  #   fill_method = symmetric9
  # [../]
  # [./strain_bubble]
  #   type = ComputeSmallStrain
  #   base_name = bubble
  # [../]
  # [./stress_bubble]
  #   type = ComputeLinearElasticStress
  #   base_name = bubble
  # [../]
  # [./const_stress]
  #   type = ComputeExtraStressConstant
  #   base_name = bubble
  #   extra_stress_tensor = '-6.25e-3 -6.25e-3 -6.25e-3 0 0 0'
  # [../]
  #
  # [./global_stress]
  #   type = MultiPhaseStressMaterial
  #   phase_base = 'bubble matrix'
  #   h          = 'hb     hm'
  # [../]

  [./XeRate]
    type = ParsedMaterial
    f_name = XeRate
    args = 'XolotlXeRate'
    function = 'XolotlXeRate'
    outputs = exodus
  [../]

  [./VacRate]
    type = ParsedMaterial
    f_name = VacRate
    material_property_names = 'YXe'
    args = 'XolotlXeRate'
    function = 'XolotlXeRate / YXe'
    outputs = exodus
  [../]
[]

#[Adaptivity]
#  marker = errorfrac
#  max_h_level = 3
#  [./Indicators]
#    [./error]
#      type = GradientJumpIndicator
#      variable = bnds
#    [../]
#  [../]
#  [./Markers]
#    [./bound_adapt]
#      type = ValueThresholdMarker
#      third_state = DO_NOTHING
#      coarsen = 1.0
#      refine = 0.99
#      variable = bnds
#      invert = true
#    [../]
#    [./errorfrac]
#      type = ErrorFractionMarker
#      coarsen = 0.1
#      indicator = error
#      refine = 0.7
#    [../]
#  [../]
#[]

[Postprocessors]
  [./number_DOFs]
    type = NumDOFs
  [../]
  [./dt]
    type = TimestepSize
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient
  nl_max_its = 15
  scheme = bdf2
  #solve_type = NEWTON
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -sub_pc_type'
  petsc_options_value = 'asm      lu'
  l_max_its = 15
  l_tol = 1.0e-3
  nl_rel_tol = 1.0e-8
  start_time = 0.0
  #num_steps = 1000
  end_time = 1e9
  nl_abs_tol = 1e-10
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 0.5
    adapt_log = true
  [../]
[]

[Outputs]
  [./exodus]
    type = Exodus
    # interval = 10
    interval = 1
  [../]
  checkpoint = true
  csv = true
[]

[MultiApps]
  [./XolotlWrapper]
    type = TransientMultiApp
    app_type = coupling_xolotlApp
    execute_on = TIMESTEP_BEGIN
    positions = '0 0 0'
    input_files = 'xolotl_userobj.i'
  [../]
[]

[Transfers]
  [./Xolotl_to_PF]
    type = MultiAppInterpolationTransfer
    direction = from_multiapp
    multi_app = XolotlWrapper
    source_variable = Auxv
    variable = XolotlXeRate
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./PF_to_Xolotl]
    type = MultiAppInterpolationTransfer
    direction = to_multiapp
    multi_app = XolotlWrapper
    source_variable = bnds
    variable = Auxv_gb
    execute_on = SAME_AS_MULTIAPP
  [../]
[]
