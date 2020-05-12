# Length unit: nm
# Time unit: s
# Energy unit: eV

XolotlWrapperPath = '/home/donguk.kim/gcc_moose/projects/coupling_xolotl/Demo3D/graintracker/small_cube3D/XolotlWrapp.i'

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 59
  ny = 59
  nz = 59
  xmin = 0
  xmax = 10620
  ymin = 0
  ymax = 10620
  zmin = 0
  zmax = 10620
  #uniform_refine = 2
[]

[GlobalParams]
  op_num = 25
  grain_num = 32
  var_name_base = etam
  numbub = 8
  bubspac = 2500
  radius = 733
  int_width = 480
  invalue = 1
  outvalue = 0
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
  [./time]
  [../]
  [./cg]
    order = FIRST
    family = MONOMIAL
  [../]
  [./cv]
    order = FIRST
    family = MONOMIAL
  [../]

  #For use of Grain Tracker
  [./unique_grains]
    order = FIRST
    family = MONOMIAL
  [../]
  [./var_indices]
    order = FIRST
    family = MONOMIAL
  [../]

  [./halos]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo0] #We need the 'op_num' number of halos here; so we need 8 halo variables
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo1]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo2]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo3]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo4]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo5]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo6]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo7]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo8]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo9]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo10]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo11]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo12]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo13]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo14]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo15]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo16]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo17]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo18]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo19]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo20]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo21]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo22]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo23]
    order = FIRST
    family = MONOMIAL
  [../]
  [./halo24]
    order = FIRST
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./time]
    type = FunctionAux
    variable = time
    function = 't'
  [../]
  [./cg]
    type = MaterialRealAux
    variable = cg
    property = cg_from_rhog
  [../]
  [./cv]
    type = MaterialRealAux
    variable = cv
    property = cv_from_rhov
  [../]

  #For use of Grain Tracker
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  [../]
  [./halo0]
    type = FeatureFloodCountAux
    variable = halo0
    map_index = 0
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo1]
    type = FeatureFloodCountAux
    variable = halo1
    map_index = 1
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo2]
    type = FeatureFloodCountAux
    variable = halo2
    map_index = 2
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo3]
    type = FeatureFloodCountAux
    variable = halo3
    map_index = 3
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo4]
    type = FeatureFloodCountAux
    variable = halo4
    map_index = 4
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo5]
    type = FeatureFloodCountAux
    variable = halo5
    map_index = 5
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo6]
    type = FeatureFloodCountAux
    variable = halo6
    map_index = 6
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo7]
    type = FeatureFloodCountAux
    variable = halo7
    map_index = 7
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo8]
    type = FeatureFloodCountAux
    variable = halo8
    map_index = 8
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo9]
    type = FeatureFloodCountAux
    variable = halo9
    map_index = 9
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo10]
    type = FeatureFloodCountAux
    variable = halo10
    map_index = 10
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo11]
    type = FeatureFloodCountAux
    variable = halo11
    map_index = 11
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo12]
    type = FeatureFloodCountAux
    variable = halo12
    map_index = 12
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo13]
    type = FeatureFloodCountAux
    variable = halo13
    map_index = 13
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo14]
    type = FeatureFloodCountAux
    variable = halo14
    map_index = 14
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo15]
    type = FeatureFloodCountAux
    variable = halo15
    map_index = 15
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo16]
    type = FeatureFloodCountAux
    variable = halo16
    map_index = 16
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo17]
    type = FeatureFloodCountAux
    variable = halo17
    map_index = 17
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo18]
    type = FeatureFloodCountAux
    variable = halo18
    map_index = 18
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo19]
    type = FeatureFloodCountAux
    variable = halo19
    map_index = 19
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo20]
    type = FeatureFloodCountAux
    variable = halo20
    map_index = 20
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo21]
    type = FeatureFloodCountAux
    variable = halo21
    map_index = 21
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo22]
    type = FeatureFloodCountAux
    variable = halo22
    map_index = 22
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo23]
    type = FeatureFloodCountAux
    variable = halo23
    map_index = 23
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
  [./halo24]
    type = FeatureFloodCountAux
    variable = halo24
    map_index = 24
    field_display = HALOS
    flood_counter = grain_tracker
  [../]
[]

[ICs]
  # [./bnds]
  #   type = ConstantIC
  #   variable = bnds
  #   value = 1
  # [../]
  [./PolycrystalICs]
    [./PolycrystalVoronoiVoidIC]
      invalue = 1.0
      outvalue = 0.0
      polycrystal_ic_uo = voronoi
    [../]
  [../]
  [./bubble_IC]
    variable = etab0
    type = PolycrystalVoronoiVoidIC
    structure_type = voids
    invalue = 1.0
    outvalue = 0.0
    polycrystal_ic_uo = voronoi
  [../]
  [./IC_wv]
    variable = wv
    type = PolycrystalVoronoiVoidIC
    structure_type = voids
    invalue = 0.0
    outvalue = 0.0
    polycrystal_ic_uo = voronoi
  [../]
  [./IC_wg]
    variable = wg
    type = PolycrystalVoronoiVoidIC
    structure_type = voids
    invalue = 0.0
    outvalue = 0.0
    polycrystal_ic_uo = voronoi
  [../]

  #[./etab0_IC]
  #  type = SmoothCircleIC
  #  variable = etab0
  #  invalue = 1
  #  outvalue = 0
  #  x1 = 4500
  #  y1 = 4500
  #  radius = 733
  #  int_width = 480
  #[../]
  #[./etam0_IC]
  #  type = SmoothCircleIC
  #  variable = etam0
  #  invalue = 0
  #  outvalue = 1
  #  x1 = 4500
  #  y1 = 4500
  #  radius = 733
  #  int_width = 480
  #[../]
  #[./IC_wv]
  #  type = ConstantIC
  #  variable = wv
  #  value = 0
  #[../]
  #[./IC_wg]
  #  type = ConstantIC
  #  variable = wg
  #  value = 0
  #[../]
[]



[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x y'
    [../]
  [../]
  # [./etam0_adiabatic]
  #   type = NeumannBC
  #   boundary = 'left right top bottom'
  #   variable = etam0
  #   value = 0
  # [../]
  # [./etam1_adiabatic]
  #   type = NeumannBC
  #   boundary = 'left right top bottom'
  #   variable = etam1
  #   value = 0
  # [../]
  # [./etam2_adiabatic]
  #   type = NeumannBC
  #   boundary = 'left right top bottom'
  #   variable = etam2
  #   value = 0
  # [../]
  # [./etam3_adiabatic]
  #   type = NeumannBC
  #   boundary = 'left right top bottom'
  #   variable = etam3
  #   value = 0
  # [../]
  # [./etam4_adiabatic]
  #   type = NeumannBC
  #   boundary = 'left right top bottom'
  #   variable = etam4
  #   value = 0
  # [../]
  # [./etab0_adiabatic]
  #   type = NeumannBC
  #   boundary = 'left right top bottom'
  #   variable = etab0
  #   value = 0
  # [../]
  # [./wg_adiabatic]
  #   type = NeumannBC
  #   boundary = 'left right top bottom'
  #   variable = wg
  #   value = 0
  # [../]
  # [./wb_adiabatic]
  #   type = NeumannBC
  #   boundary = 'left right top bottom'
  #   variable = wv
  #   value = 0
  # [../]
[]

[Modules]
  [./PhaseField]
    [./GrandPotential]
      switching_function_names = 'hb hm'
      anisotropic = 'false false'

      chemical_potentials = 'wv wg'
      mobilities = 'Dchiv Dchig'
      susceptibilities = 'chiv chig'
      free_energies_w = 'rhovbub rhovmatrix rhogbub rhogmatrix'

      # chemical_potentials = 'wg'
      # mobilities = 'Dchig'
      # susceptibilities = 'chig'
      # free_energies_w = 'rhogbub rhogmatrix'

      # gamma_gr = gamma
      gamma_gr = gmm
      mobility_name_gr = L
      kappa_gr = kappa
      free_energies_gr = 'omegab omegam'

      additional_ops = 'etab0'
      gamma_grxop = gmb
      mobility_name_op = L
      kappa_op = kappa
      free_energies_op = 'omegab omegam'
    [../]
  [../]
[]

[Kernels]


  #[./Source_v]
  #  type = MaskedBodyForce
  #  variable = wv
  #  value = 1
  #  mask = VacRate0
  #[../]
  [./Source_v]
    type = MaskedBodyForce
    variable = wv
    value = 1
    mask = VacRate
  [../]

  #[./Source_g]
  #  type = MaskedBodyForce
  #  variable = wg
  #  value = 1
  #  mask = XeRate0
  #[../]
  [./Source_g]
    type = MaskedBodyForce
    variable = wg
    value = 1 # for unit conversion between PF app and Xolotl
    mask = XeRate
  [../]
[]

[Materials]
  [./hb]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = hb
    all_etas = 'etab0 etam0 etam1 etam2 etam3 etam4 etam5 etam6 etam7 etam8 etam9 etam10 etam11 etam12 etam13 etam14 etam15 etam16 etam17 etam18 etam19 etam20 etam21 etam22 etam23 etam24'
    #all_etas = 'etab0 etam0 etam1 etam2 etam3 etam4 etam5 etam6 etam7 etam8 etam9'
    phase_etas = 'etab0'
    #outputs = exodus
  [../]
  [./hm]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = hm
    #all_etas = 'etab0 etam0 etam1 etam2 etam3 etam4 etam5 etam6 etam7 etam8 etam9'
    # phase_etas = 'etam0 etam1 etam2 etam3 etam4'
    all_etas = 'etab0 etam0 etam1 etam2 etam3 etam4 etam5 etam6 etam7 etam8 etam9 etam10 etam11 etam12 etam13 etam14 etam15 etam16 etam17 etam18 etam19 etam20 etam21 etam22 etam23 etam24'
    phase_etas = 'etam0 etam1 etam2 etam3 etam4 etam5 etam6 etam7 etam8 etam9 etam10 etam11 etam12 etam13 etam14 etam15 etam16 etam17 etam18 etam19 etam20 etam21 etam22 etam23 etam24'
    #outputs = exodus
  [../]
# Chemical contribution to grand potential of bubble
  [./omegab]
    type = DerivativeParsedMaterial
    args = 'wv wg time'
    f_name = omegab
    material_property_names = 'Va kvbub cvbubeq kgbub cgbubeq'
    function = 'if(time < 0, 0, -0.5*wv^2/Va^2/kvbub-wv/Va*cvbubeq-0.5*wg^2/Va^2/kgbub-wg/Va*cgbubeq)'
    derivative_order = 2
    #outputs = exodus
  [../]

# Chemical contribution to grand potential of matrix
  [./omegam]
    type = DerivativeParsedMaterial
    args = 'wv wg time'
    f_name = omegam
    material_property_names = 'Va kvmatrix cvmatrixeq kgmatrix cgmatrixeq'
    function = 'if(time < 0, 0, -0.5*wv^2/Va^2/kvmatrix-wv/Va*cvmatrixeq-0.5*wg^2/Va^2/kgmatrix-wg/Va*cgmatrixeq)'
    derivative_order = 2
    #outputs = exodus
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
    prop_names =  'kappa     mu      L   Dm      Db     Va      cvbubeq  cgbubeq gmb 	gmm T     YXe'
    #prop_values = '0.84368e2 0.00293 0.1   2.3398  2.3398 0.0409  0.546    0.454   0.922  1.5 1200  0.2156'
    #prop_values = '0.84368e1 0.00293e-1 0.5   4.6796  4.6796  0.0409  0.546    0.454   0.922  1.5 1200  0.2156'
    prop_values = '0.84368e1 0.00293e-1 1    2.3398 2.3398  0.0409  0.546    0.454   0.922  1.5 1200  0.2156'
    #prop_values = '2.21125e2 1.875  0.975e-3 0.1  0.1  0.0409  0.546    0.454   0.922  1.5 1200  0.2156'
  [../]
  [./cvmatrixeq]
    type = ParsedMaterial
    f_name = cvmatrixeq
    material_property_names = 'T'
    constant_names        = 'kB           Efv'  # in eV/atom
    constant_expressions  = '8.6173324e-5 3.0'
    function = 'exp(-Efv/(kB*T))'
  [../]
  [./cgmatrixeq]
    type = ParsedMaterial
    f_name = cgmatrixeq
    material_property_names = 'T'
    constant_names        = 'kB           Efg'  # in eV/atom
    constant_expressions  = '8.6173324e-5 3.0'
    function = 'exp(-Efg/(kB*T))'
  [../]
  [./kvmatrix_parabola]
    type = ParsedMaterial
    f_name = kvmatrix
    args = 'time'
    #function = '3.00625e3' # in eV/nm^3
    #function = '7.516' # in eV/nm^3
    function = '7.516e-2' # in eV/nm^3, decreased
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
    #function = '0.5625e3' # in eV/nm^3
    #function = '1.406' # in eV/nm^3
    function = '1.406e-2' # in eV/nm^3, decreased
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
    material_property_names = 'Db chiv'
    args = 'time'
    function = 'if(time < 0, 100*Db*chiv, Db*chiv)'
    #function = 'Db*chiv'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./Mobility_g]
    type = DerivativeParsedMaterial
    f_name = Dchig
    material_property_names = 'Dm chig'
    args = 'time'
    function = 'if(time < 0, 100*Dm*chig, Dm*chig)'
    #function = 'Dm*chig'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./chiv]
    type = DerivativeParsedMaterial
    f_name = chiv
    material_property_names = 'Va hb kvbub hm kvmatrix '
    function = '(hm/kvmatrix + hb/kvbub) / Va^2'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./chig]
    type = DerivativeParsedMaterial
    f_name = chig
    material_property_names = 'Va hb kgbub hm kgmatrix '
    function = '(hm/kgmatrix + hb/kgbub) / Va^2'
    derivative_order = 2
    # outputs = exodus
  [../]

  [./XeRate]
    type = ParsedMaterial
    f_name = XeRate
    material_property_names = 'hm'
    args = 'time XolotlXeRate'  # XolotlXeRate is in Xe/(nm^3 * s) & Va is in Xe/nm^3
    # function = 'if(time < 0, 0, XolotlXeRate * hm)'
    function = 'if(time < 0, 0, XolotlXeRate)'
    outputs = exodus
  [../]

  [./VacRate]
    type = ParsedMaterial
    f_name = VacRate
    material_property_names = 'XeRate YXe'
    function = 'XeRate / YXe'
    outputs = exodus
  [../]

  [./XeRate_ref]
    type = ParsedMaterial
    f_name = XeRate0
    material_property_names = 'Va hm'
    constant_names = 's0'
    constant_expressions = '2.35e-9'  # in atoms/(nm^3 * s)
    args = 'time'
    function = 'if(time < 0, 0, s0 * hm)'
    outputs = exodus
  [../]
  [./VacRate_ref]
    type = ParsedMaterial
    f_name = VacRate0
    material_property_names = 'YXe XeRate0'
    args = 'time'
    function = 'if(time < 0, 0, XeRate0 / YXe)'
    outputs = exodus
  [../]

  [./cg]
    type = ParsedMaterial
    f_name = cg_from_rhog
    material_property_names = 'Va rhogbub rhogmatrix hm hb'
    function = 'hb*Va*rhogbub + hm*Va*rhogmatrix'
    outputs = exodus
  [../]
  [./cv]
    type = ParsedMaterial
    f_name = cv_from_rhov
    material_property_names = 'Va rhovbub rhovmatrix hm hb'
    function = 'hb*Va*rhovbub + hm*Va*rhovmatrix'
    outputs = exodus
  [../]
[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    rand_seed = 2
    # int_width = 480
    coloring_algorithm = jp
    output_adjacency_matrix = false
  [../]

  [./grain_tracker]
    type = GrainTracker
    threshold = 0.5
    connecting_threshold = 0.5
    compute_halo_maps = true
    remap_grains = true
  [../]
[]

#[Adaptivity]
#  marker = errorfrac
#  max_h_level = 2
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
  # [./number_DOFs]
  #   type = NumDOFs
  # [../]
  # [./dt]
  #   type = TimestepSize
  # [../]
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
  solve_type = NEWTON
  #solve_type = PJFNK
  #petsc_options_iname = '-pc_type -sub_pc_type'
  #petsc_options_value = 'asm      lu'
  # Use with solve_type = NEWTON in Executioner block
  petsc_options_iname = '-pc_type -ksp_type -ksp_gmres_restart'
  petsc_options_value = 'bjacobi  gmres     30'  # default is 30, the higher the higher resolution but the slower

  l_max_its = 15
  l_tol = 1.0e-5
  nl_rel_tol = 1.0e-8
  start_time = 0
  #num_steps = 1000
  #end_time = 1.2e8
  end_time = 1.0e4
  nl_abs_tol = 1e-10
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5
    growth_factor = 1.2
    cutback_factor = 0.8
    # adapt_log = true
  [../]
[]

[Outputs]
  [./exodus]
    type = Exodus
    #interval = 10
    interval = 1
    #sync_times = '0 1.2e8'
    sync_times = '0 1.0e4'
  [../]
  # checkpoint = true
  csv = true
  perf_graph = true
[]

[MultiApps]
  [./XolotlWrapper]
    type = TransientMultiApp
    app_type = coupling_xolotlApp
    execute_on = TIMESTEP_END
    positions = '0 0 0'
    input_files = ${XolotlWrapperPath}
  [../]
[]

[Transfers]
  [./Xolotl_to_PF]
    type = MultiAppInterpolationTransfer
    # type = MultiAppNearestNodeTransfer
    direction = from_multiapp
    multi_app = XolotlWrapper
    source_variable = Auxv
    variable = XolotlXeRate
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./PF_to_Xolotl]
    type = MultiAppInterpolationTransfer
    # type = MultiAppNearestNodeTransfer
    direction = to_multiapp
    multi_app = XolotlWrapper
    source_variable = bnds
    variable = Auxv_gb
    execute_on = SAME_AS_MULTIAPP
  [../]
[]
