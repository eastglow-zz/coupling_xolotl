# Length unit: nm
# Time unit: s
# Mass unit: ?

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1000
  ny = 10
  xmin = 0
  xmax = 100000
  ymin = 0
  ymax = 100000
[]

[GlobalParams]
  op_num = 2
  grain_num = 2
  var_name_base = etam
  numbub = 1
  bubspac = 150
  radius = 44
  int_width = 20
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
  [./bnds]
    type = ConstantIC
    variable = bnds
    value = 1
  [../]
  [./etam0_IC]
    type = BoundingBoxIC
    variable = etam0
    inside = 1
    outside = 0
    x1 = 0
    x2 = 49998
    y1 = 0
    y2 = 100000
  [../]
  [./etam1_IC]
    type = BoundingBoxIC
    variable = etam1
    inside = 1
    outside = 0
    x1 = 50002
    x2 = 100000
    y1 = 0
    y2 = 100000
  [../]
  [./bubble_IC]
    type = BoundingBoxIC
    variable = etab0
    inside = 1
    outside = 0
    x1 = 49998
    x2 = 50002
    y1 = 0
    y2 = 100000
  [../]
  [./IC_wv]
    type = ConstantIC
    variable = wv
    value = 0
  [../]
  [./IC_wg]
    type = ConstantIC
    variable = wg
    value = 0
  [../]
[]

[BCs]
  [./etam0_adiabatic]
    type = NeumannBC
    boundary = 'left right top bottom'
    variable = etam0
    value = 0
  [../]
  [./etam1_adiabatic]
    type = NeumannBC
    boundary = 'left right top bottom'
    variable = etam1
    value = 0
  [../]
  [./etab0_adiabatic]
    type = NeumannBC
    boundary = 'left right top bottom'
    variable = etab0
    value = 0
  [../]
  [./wg_adiabatic]
    type = NeumannBC
    boundary = 'left right top bottom'
    variable = wg
    value = 0
  [../]
  [./wb_adiabatic]
    type = NeumannBC
    boundary = 'left right top bottom'
    variable = wv
    value = 0
  [../]
[]

[Kernels]
  [./ACb0_bulk]
    type = ACGrGrMulti
    variable = etab0
    v =           'etam0 etam1'
    gamma_names = 'gmb   gmb  '
  [../]
  [./ACb0_sw]
    type = ACSwitching
    variable = etab0
    Fj_names  = 'omega_total_bubble   omega_total_matrix'
    hj_names  = 'hb                   hm'
    args = 'etam0 etam1 wv wg'
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
  [./ACm0_bulk]
    type = ACGrGrMulti
    variable = etam0
    v =           'etab0 etam1'
    gamma_names = 'gmb   gmm  '
  [../]
  [./ACm0_sw]
    type = ACSwitching
    variable = etam0
    Fj_names  = 'omega_total_bubble   omega_total_matrix'
    hj_names  = 'hb                   hm'
    args = 'etab0 etam1 wv wg'
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
  [./ACm1_bulk]
    type = ACGrGrMulti
    variable = etam1
    v =           'etab0 etam0'
    gamma_names = 'gmb   gmm  '
  [../]
  [./ACm1_sw]
    type = ACSwitching
    variable = etam1
    Fj_names  = 'omega_total_bubble   omega_total_matrix'
    hj_names  = 'hb                   hm'
    args = 'etab0 etam0 wv wg'
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
  [./Source_v]
    type = CoupledForce
    coef = 1
    variable = wv
    v = VacRate
  [../]
  [./coupled_v_etab0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wv
    v = etab0
    Fj_names = 'rhovbub rhovmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
  [./coupled_v_etam0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wv
    v = etam0
    Fj_names = 'rhovbub rhovmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
  [./coupled_v_etam1dot]
    type = CoupledSwitchingTimeDerivative
    variable = wv
    v = etam1
    Fj_names = 'rhovbub rhovmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
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
  [./Source_g]
    type = CoupledForce
    variable = wg
    coef = 1 # for unit conversion between PF app and Xolotl
    v = XolotlXeRate
  [../]
  [./coupled_g_etab0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wg
    v = etab0
    Fj_names = 'rhogbub rhogmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
  [./coupled_g_etam0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wg
    v = etam0
    Fj_names = 'rhogbub rhogmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
  [./coupled_g_etam1dot]
    type = CoupledSwitchingTimeDerivative
    variable = wg
    v = etam1
    Fj_names = 'rhogbub rhogmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
[]

[AuxKernels]
  [./BndsCalcInit]
    type = BndsCalcAux
    variable = bnds
    execute_on = initial
  [../]
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
    all_etas = 'etab0 etam0 etam1'
    phase_etas = 'etab0'
  [../]
  [./hm]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = hm
    all_etas = 'etab0 etam0 etam1'
    phase_etas = 'etam0 etam1'
  [../]
  [./omegab]
    type = DerivativeParsedMaterial
    args = 'wv wg'
    f_name = omegab
    material_property_names = 'Va kvbub cvbubeq kgbub cgbubeq'
    function = '-0.5*wv^2/Va^2/kvbub-wv/Va*cvbubeq-0.5*wg^2/Va^2/kgbub-wg/Va*cgbubeq'
    derivative_order = 2
  [../]
  [./Total_energy_bubble]
    type = DerivativeSumMaterial
    f_name = omega_total_bubble
    sum_materials = 'omegab'
    args = 'wv wg'
  [../]
  [./omegam]
    type = DerivativeParsedMaterial
    args = 'wv wg'
    f_name = omegam
    material_property_names = 'Va kvmatrix cvmatrixeq kgmatrix cgmatrixeq'
    function = '-0.5*wv^2/Va^2/kvmatrix-wv/Va*cvmatrixeq-0.5*wg^2/Va^2/kgmatrix-wg/Va*cgmatrixeq'
    derivative_order = 2
  [../]
  [./Total_energy_matrix]
    type = DerivativeSumMaterial
    block = 0
    f_name = omega_total_matrix
    sum_materials = 'omegam'
    args = 'wv wg'
  [../]
  [./rhovbub]
    type = DerivativeParsedMaterial
    args = 'wv'
    f_name = rhovbub
    material_property_names = 'Va kvbub cvbubeq'
    function = 'wv/Va^2/kvbub + cvbubeq/Va'
    derivative_order = 2
  [../]
  [./rhovmatrix]
    type = DerivativeParsedMaterial
    args = 'wv'
    f_name = rhovmatrix
    material_property_names = 'Va kvmatrix cvmatrixeq'
    function = 'wv/Va^2/kvmatrix + cvmatrixeq/Va'
    derivative_order = 2
  [../]
  [./rhogbub]
    type = DerivativeParsedMaterial
    args = 'wg'
    f_name = rhogbub
    material_property_names = 'Va kgbub cgbubeq'
    function = 'wg/Va^2/kgbub + cgbubeq/Va'
    derivative_order = 2
  [../]
  [./rhogmatrix]
    type = DerivativeParsedMaterial
    args = 'wg'
    f_name = rhogmatrix
    material_property_names = 'Va kgmatrix cgmatrixeq'
    function = 'wg/Va^2/kgmatrix + cgmatrixeq/Va'
    derivative_order = 2
  [../]
  [./const]
    type = GenericConstantMaterial
    prop_names =  'kappa   mu       L    Dm    Db     Va      cvbubeq cgbubeq gmb 	  gmm T    tgrad_corr_mult YXe'
    prop_values = '1.0     0.004688 0.01 2.3398 2.3398 0.04092 0.61    0.39    0.9218 1.5 1800 0.0             0.2156'
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
    outputs = exodus
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
    outputs = exodus
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
    material_property_names = 'Dg chiv'
    function = 'Dg*chiv'
    derivative_order = 2
    outputs = exodus
  [../]
  [./Mobility_g]
    type = DerivativeParsedMaterial
    f_name = Dchig
    material_property_names = 'Dm chig'
    function = 'Dm*chig'
    derivative_order = 2
    outputs = exodus
  [../]
  [./chiv]
    type = DerivativeParsedMaterial
    f_name = chiv
    material_property_names = 'Va hb kvbub hm kvmatrix '
    function = '(hm/kvmatrix + hb/kvbub) / Va^2'
    derivative_order = 2
    outputs = exodus
  [../]
  [./chig]
    type = DerivativeParsedMaterial
    f_name = chig
    material_property_names = 'Va hb kgbub hm kgmatrix '
    function = '(hm/kgmatrix + hb/kgbub) / Va^2'
    derivative_order = 2
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

[Functions]
  [./dts]
    type = PiecewiseLinear
    x = '0 1 10 100 1000 10000 100000000'
    y = '0.1 0.1 1 10 100 1000 5000'
  [../]
[]

[Executioner]
  type = Transient
  nl_max_its = 15
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -sub_pc_type'
  petsc_options_value = 'asm      lu'
  l_max_its = 15
  l_tol = 1.0e-3
  nl_rel_tol = 1.0e-8
  start_time = 0.0
  num_steps = 10000
  end_time = 1e8
  nl_abs_tol = 1e-10
  [./TimeStepper]
    type = FunctionDT
    function = dts
    min_dt = 0.01
  [../]
[]

[Outputs]
  [./exodus]
    type = Exodus
    # interval = 10
    interval = 200
  [../]
  checkpoint = true
  csv = true
[]

[MultiApps]
  [./sub_app]
    type = TransientMultiApp
    input_files = 'xolotl_2D.i'
    app_type = coupling_xolotlApp
    library_path = 'lib'
    sub_cycling = true
  [../]
[]

[Transfers]
  [./fromsub]
    type = MultiAppNearestNodeTransfer
    direction = from_multiapp
    multi_app = sub_app
    source_variable = Auxv
    variable = XolotlXeRate
    fixed_meshes = true
  [../]
  [./tosub]
    type = MultiAppInterpolationTransfer
    direction = to_multiapp
    multi_app = sub_app
    source_variable = bnds
    variable = AuxGB
    fixed_meshes = true
  [../]
[]
