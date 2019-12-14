# default length unit: nm
# default time unit: s
# default mass unit: ?

XolotlLibPath = '/home/donguk.kim/gcc_moose/projects/xolotl-dlload-gcc/lib/libxolotlInter.so'    #for UF HPG2
XolotlInpPath = '/home/donguk.kim/gcc_moose/projects/coupling_xolotl/Demo3D/graintracker/small_cube3D/XolotlInput_NoCL.txt'   #for UF HPG2

[Mesh]
  type = XolotlMeshSynced
  # type = XolotlMesh
  dim = 2
  library_path_name = ${XolotlLibPath}
  XolotlInput_path_name = ${XolotlInpPath}
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
  # [./Init_Aux_gb]
  #   type = BoundingBoxIC
  #   variable = Auxv_gb
  #   inside = 0
  #   outside = 1
  #   # x1 = 39999.5
  #   # x2 = 60000.5
  #   # y1 = 0
  #   # y2 = 100000
  #   # z1 = 0
  #   # z2 = 100000
  #   x1 = 0
  #   x2 = 1200
  #   y1 = 399.5
  #   y2 = 600.5
  #   # z1 = 0
  #   # z2 = 100000
  # [../]
  [./Init_Aux_gb]
    type = ConstantIC
    variable = Auxv_gb
    value = 1
  [../]
[]

[AuxKernels]
  [./Spatial_Auxv]
    type = SpatialUserObjectAux
    variable = Auxv
    user_object = 'Xolotl_driver'
    execute_on = 'TIMESTEP_END'
  [../]
[]
[Executioner]
  type = Transient
  [./TimeStepper]
    type = ConstantDT
    dt = 5000.0
  [../]
  start_time = 0
  #end_time = 2e8
  end_time = 1e6
  # end_time = 20000.0
[]
[Problem]
  kernel_coverage_check = false
  solve = false
[]
[UserObjects]
  [./Xolotl_driver]
    # type = XolotlUserObject
    type = XolotlUserObjectSynced
    variable = Auxv
    variable_gb = Auxv_gb
    gb_marker_threshold = 0.9
    library_path_name = ${XolotlLibPath}
    XolotlInput_path_name = ${XolotlInpPath}
  [../]
[]
[Outputs]
  exodus = true
[]
