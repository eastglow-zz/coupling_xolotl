[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 21
  xmax = 10
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
    x1 = 4
    x2 = 6
    y1 = -0.1
    y2 = 0.1
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
    dt = 1
  [../]
  start_time = 0
  end_time = 2
[]
[Problem]
  kernel_coverage_check = false
  solve = false
[]
[UserObjects]
  [./Xolotl_driver]
    type = XolotlUserObject
    dim = 1
    nx = 21
    xmax = 10
    variable = Auxv
    variable_gb = Auxv_gb
    library_path_name ='/Users/donguk.kim/projects/xolotl-build/lib/libxolotlInter.dylib'
    XolotlInput_path_name = './params_NE.txt'
  [../]
[]
[Outputs]
  exodus = true
[]
