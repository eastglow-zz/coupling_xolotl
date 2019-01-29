[Mesh]
  # type = GeneratedMesh
  # dim = 3
  # nx = 144
  # ny = 20
  # nz = 20
  # xmin = 0
  # xmax = 21
  # ymin = 0
  # ymax = 21
  # zmin = 0
  # zmax = 21
  type = XolotlMesh
  dim = 3
  library_path_name ='/Users/donguk.kim/projects/xolotl-build/lib/libxolotlInter.dylib'
  XolotlInput_path_name = './params_NE_3D.txt'
[]
[Variables]
  [./d]
  [../]
[]
[AuxVariables]
  [./Auxv]
  [../]
  [./Auxv_gb]
  [../]
[]
[ICs]
  [./Init_Aux_const]
    type = ConstantIC
    variable = Auxv
    value = 0
  [../]
  [./Init_Aux_gb]
    type = BoundingBoxIC
    variable = Auxv_gb
    inside = 1
    outside = 0
    x1 = 71999
    x2 = 72001
    y1 = 0
    y2 = 100000
    z1 = 0
    z2 = 100000
  [../]
[]

[AuxKernels]
  [./Spatial_Auxv]
    type = SpatialUserObjectAux
    variable = Auxv
    user_object = 'Xolotl_driver'
  [../]
[]
[Executioner]
  type = Transient
  [./TimeStepper]
    type = ConstantDT
    dt = 5000.0
  [../]
  start_time = 0
  end_time = 7e7
  # end_time = 20000.0
[]
[Problem]
  kernel_coverage_check = false
  solve = false
[]
[UserObjects]
  [./Xolotl_driver]
    type = XolotlUserObject
    variable = Auxv
    variable_gb = Auxv_gb
    library_path_name ='/Users/donguk.kim/projects/xolotl-build/lib/libxolotlInter.dylib'
    XolotlInput_path_name = './params_NE_3D.txt'
  [../]
[]
[Outputs]
  exodus = true
[]
